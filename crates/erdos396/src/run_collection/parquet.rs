use crate::false_positive::RunWindowClassification;
use crate::run_collection::schema::{
    maximal_runs_batch, run_window_classifications_batch, run_windows_batch,
};
use crate::run_collection::{RunRecord, RunWindowRecord};
use anyhow::{Context, Result};
use arrow::array::{Array, BooleanArray, StringArray, UInt32Array, UInt64Array};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

const DEFAULT_BATCH_SIZE: usize = 8_192;

#[derive(Debug, Clone, Default)]
pub struct RunWindowFilter {
    pub window_length: Option<u32>,
    pub start_n: Option<u64>,
    pub end_n: Option<u64>,
    pub limit: Option<usize>,
}

#[derive(Debug, Clone, Default)]
pub struct ClassificationFilter {
    pub window_length: Option<u32>,
}

pub fn prepare_build_output(output_root: &Path) -> Result<()> {
    for dataset in ["maximal_runs", "run_windows"] {
        let dataset_root = output_root.join(dataset);
        if dataset_root.exists() {
            fs::remove_dir_all(&dataset_root)
                .with_context(|| format!("remove {}", dataset_root.display()))?;
        }
    }
    fs::create_dir_all(output_root).with_context(|| format!("create {}", output_root.display()))?;
    Ok(())
}

pub fn maximal_run_writer(output_root: &Path) -> PartitionedDatasetWriter<RunRecord> {
    PartitionedDatasetWriter::new(
        output_root.to_path_buf(),
        DEFAULT_BATCH_SIZE,
        maximal_partition_dir,
        maximal_runs_batch,
        false,
        None,
    )
}

pub fn run_window_writer(output_root: &Path) -> PartitionedDatasetWriter<RunWindowRecord> {
    PartitionedDatasetWriter::new(
        output_root.to_path_buf(),
        DEFAULT_BATCH_SIZE,
        run_windows_partition_dir,
        run_windows_batch,
        false,
        None,
    )
}

pub fn classification_writer(
    output_root: &Path,
    invocation_id: &str,
) -> PartitionedDatasetWriter<RunWindowClassification> {
    PartitionedDatasetWriter::new(
        output_root.to_path_buf(),
        DEFAULT_BATCH_SIZE,
        classifications_partition_dir,
        run_window_classifications_batch,
        true,
        Some(invocation_id.to_string()),
    )
}

pub fn read_run_windows(
    output_root: &Path,
    filter: &RunWindowFilter,
) -> Result<Vec<RunWindowRecord>> {
    let root = match filter.window_length {
        Some(length) => output_root
            .join("run_windows")
            .join(format!("window_length={length}")),
        None => output_root.join("run_windows"),
    };

    let mut out = Vec::new();
    for path in parquet_files_under(&root)? {
        let file = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .with_context(|| format!("read parquet {}", path.display()))?;
        let reader = builder.build()?;
        for batch in reader {
            let batch = batch?;
            for record in run_window_records_from_batch(&batch)? {
                if !matches_run_window_filter(&record, filter) {
                    continue;
                }
                out.push(record);
                if let Some(limit) = filter.limit {
                    if out.len() >= limit {
                        return Ok(out);
                    }
                }
            }
        }
    }
    Ok(out)
}

pub fn latest_classifications(
    output_root: &Path,
    filter: &ClassificationFilter,
) -> Result<HashMap<String, RunWindowClassification>> {
    let root = match filter.window_length {
        Some(length) => output_root
            .join("run_window_classifications")
            .join(format!("window_length={length}")),
        None => output_root.join("run_window_classifications"),
    };
    let mut latest = HashMap::new();
    if !root.exists() {
        return Ok(latest);
    }

    for path in parquet_files_under(&root)? {
        let file = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .with_context(|| format!("read parquet {}", path.display()))?;
        let reader = builder.build()?;
        for batch in reader {
            let batch = batch?;
            for record in run_window_classifications_from_batch(&batch)? {
                let replace = latest
                    .get(&record.window_id)
                    .map(|existing: &RunWindowClassification| {
                        record.verified_at_utc > existing.verified_at_utc
                    })
                    .unwrap_or(true);
                if replace {
                    latest.insert(record.window_id.clone(), record);
                }
            }
        }
    }

    Ok(latest)
}

type BatchBuilder<T> = fn(&[T]) -> Result<RecordBatch>;
type PartitionDir<T> = fn(&Path, &T) -> PathBuf;

pub struct PartitionedDatasetWriter<T> {
    output_root: PathBuf,
    batch_size: usize,
    partition_dir: PartitionDir<T>,
    batch_builder: BatchBuilder<T>,
    append_mode: bool,
    invocation_id: Option<String>,
    buffers: HashMap<PathBuf, Vec<T>>,
    part_counters: HashMap<PathBuf, usize>,
}

impl<T> PartitionedDatasetWriter<T>
where
    T: Clone,
{
    fn new(
        output_root: PathBuf,
        batch_size: usize,
        partition_dir: PartitionDir<T>,
        batch_builder: BatchBuilder<T>,
        append_mode: bool,
        invocation_id: Option<String>,
    ) -> Self {
        Self {
            output_root,
            batch_size,
            partition_dir,
            batch_builder,
            append_mode,
            invocation_id,
            buffers: HashMap::new(),
            part_counters: HashMap::new(),
        }
    }

    pub fn push(&mut self, record: T) -> Result<()> {
        let partition = (self.partition_dir)(&self.output_root, &record);
        let buffer = self.buffers.entry(partition.clone()).or_default();
        buffer.push(record);
        if buffer.len() >= self.batch_size {
            self.flush_partition(&partition)?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        let partitions: Vec<PathBuf> = self.buffers.keys().cloned().collect();
        for partition in partitions {
            self.flush_partition(&partition)?;
        }
        Ok(())
    }

    fn flush_partition(&mut self, partition: &PathBuf) -> Result<()> {
        let Some(buffer) = self.buffers.get_mut(partition) else {
            return Ok(());
        };
        if buffer.is_empty() {
            return Ok(());
        }

        fs::create_dir_all(partition).with_context(|| format!("create {}", partition.display()))?;
        let index = self
            .part_counters
            .entry(partition.clone())
            .or_insert(0usize);
        let filename = match (&self.invocation_id, self.append_mode) {
            (Some(invocation_id), true) => format!("part-{invocation_id}-{index:05}.parquet"),
            _ => format!("part-{index:05}.parquet"),
        };
        *index += 1;

        let rows = std::mem::take(buffer);
        let batch = (self.batch_builder)(&rows)?;
        write_batch(&partition.join(filename), &batch)?;
        Ok(())
    }
}

fn write_batch(path: &Path, batch: &RecordBatch) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), None)
        .with_context(|| format!("open writer {}", path.display()))?;
    writer
        .write(batch)
        .with_context(|| format!("write {}", path.display()))?;
    writer
        .close()
        .with_context(|| format!("close {}", path.display()))?;
    Ok(())
}

fn maximal_partition_dir(root: &Path, record: &RunRecord) -> PathBuf {
    root.join("maximal_runs")
        .join(format!("campaign_target_k={}", record.campaign_target_k))
        .join(format!("range_bucket_t={}", record.range_bucket_t))
}

fn run_windows_partition_dir(root: &Path, record: &RunWindowRecord) -> PathBuf {
    root.join("run_windows")
        .join(format!("window_length={}", record.window_length))
        .join(format!("range_bucket_t={}", record.range_bucket_t))
}

fn classifications_partition_dir(root: &Path, record: &RunWindowClassification) -> PathBuf {
    root.join("run_window_classifications")
        .join(format!("window_length={}", record.window_length))
        .join(format!("audit_status={}", record.audit_status))
}

fn parquet_files_under(root: &Path) -> Result<Vec<PathBuf>> {
    if !root.exists() {
        return Ok(Vec::new());
    }
    let mut files = Vec::new();
    for entry in WalkDir::new(root) {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() && path.extension().and_then(|ext| ext.to_str()) == Some("parquet") {
            files.push(path.to_path_buf());
        }
    }
    files.sort();
    Ok(files)
}

fn matches_run_window_filter(record: &RunWindowRecord, filter: &RunWindowFilter) -> bool {
    if let Some(length) = filter.window_length {
        if record.window_length != length {
            return false;
        }
    }
    if let Some(start_n) = filter.start_n {
        if record.window_n < start_n {
            return false;
        }
    }
    if let Some(end_n) = filter.end_n {
        if record.window_n > end_n {
            return false;
        }
    }
    true
}

fn run_window_records_from_batch(batch: &RecordBatch) -> Result<Vec<RunWindowRecord>> {
    let schema_version = string_array(batch, "schema_version")?;
    let window_id = string_array(batch, "window_id")?;
    let parent_run_id = string_array(batch, "parent_run_id")?;
    let checkpoint_id = string_array(batch, "checkpoint_id")?;
    let source_label = string_array(batch, "source_label")?;
    let archive_role = string_array(batch, "archive_role")?;
    let server = string_array(batch, "server")?;
    let source_path = string_array(batch, "source_path")?;
    let campaign_target_k = u32_array(batch, "campaign_target_k")?;
    let worker_id = u32_array(batch, "worker_id")?;
    let max_run_start = u64_array(batch, "max_run_start")?;
    let max_run_end = u64_array(batch, "max_run_end")?;
    let max_run_length = u32_array(batch, "max_run_length")?;
    let window_length = u32_array(batch, "window_length")?;
    let window_k = u32_array(batch, "window_k")?;
    let window_offset = u32_array(batch, "window_offset")?;
    let window_start = u64_array(batch, "window_start")?;
    let window_end = u64_array(batch, "window_end")?;
    let window_n = u64_array(batch, "window_n")?;
    let range_bucket_t = u64_array(batch, "range_bucket_t")?;
    let recorded_status = string_array(batch, "recorded_status")?;
    let recorded_failing_prime = u64_array(batch, "recorded_failing_prime")?;
    let recorded_demand = u64_array(batch, "recorded_demand")?;
    let recorded_supply = u64_array(batch, "recorded_supply")?;
    let is_unique_coverage = bool_array(batch, "is_unique_coverage")?;

    let mut out = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        out.push(RunWindowRecord {
            schema_version: schema_version.value(row).to_string(),
            window_id: window_id.value(row).to_string(),
            parent_run_id: parent_run_id.value(row).to_string(),
            checkpoint_id: checkpoint_id.value(row).to_string(),
            source_label: source_label.value(row).to_string(),
            archive_role: archive_role.value(row).to_string(),
            server: server.value(row).to_string(),
            source_path: source_path.value(row).to_string(),
            campaign_target_k: campaign_target_k.value(row),
            worker_id: worker_id.value(row),
            max_run_start: max_run_start.value(row),
            max_run_end: max_run_end.value(row),
            max_run_length: max_run_length.value(row),
            window_length: window_length.value(row),
            window_k: window_k.value(row),
            window_offset: window_offset.value(row),
            window_start: window_start.value(row),
            window_end: window_end.value(row),
            window_n: window_n.value(row),
            range_bucket_t: range_bucket_t.value(row),
            recorded_status: nullable_string(recorded_status, row),
            recorded_failing_prime: nullable_u64(recorded_failing_prime, row),
            recorded_demand: nullable_u64(recorded_demand, row),
            recorded_supply: nullable_u64(recorded_supply, row),
            is_unique_coverage: is_unique_coverage.value(row),
        });
    }
    Ok(out)
}

fn run_window_classifications_from_batch(
    batch: &RecordBatch,
) -> Result<Vec<RunWindowClassification>> {
    let schema_version = string_array(batch, "schema_version")?;
    let classification_id = string_array(batch, "classification_id")?;
    let window_id = string_array(batch, "window_id")?;
    let window_length = u32_array(batch, "window_length")?;
    let window_k = u32_array(batch, "window_k")?;
    let window_start = u64_array(batch, "window_start")?;
    let window_end = u64_array(batch, "window_end")?;
    let window_n = u64_array(batch, "window_n")?;
    let audit_status = string_array(batch, "audit_status")?;
    let failing_prime = u64_array(batch, "failing_prime")?;
    let demand = u64_array(batch, "demand")?;
    let supply = u64_array(batch, "supply")?;
    let verifier_mode = string_array(batch, "verifier_mode")?;
    let verified_at_utc = string_array(batch, "verified_at_utc")?;
    let build_git_hash = string_array(batch, "build_git_hash")?;

    let mut out = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        out.push(RunWindowClassification {
            schema_version: schema_version.value(row).to_string(),
            classification_id: classification_id.value(row).to_string(),
            window_id: window_id.value(row).to_string(),
            window_length: window_length.value(row),
            window_k: window_k.value(row),
            window_start: window_start.value(row),
            window_end: window_end.value(row),
            window_n: window_n.value(row),
            audit_status: audit_status.value(row).to_string(),
            failing_prime: nullable_u64(failing_prime, row),
            demand: nullable_u64(demand, row),
            supply: nullable_u64(supply, row),
            verifier_mode: verifier_mode.value(row).to_string(),
            verified_at_utc: verified_at_utc.value(row).to_string(),
            build_git_hash: build_git_hash.value(row).to_string(),
        });
    }
    Ok(out)
}

fn string_array<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<StringArray>()
        .with_context(|| format!("column {name} is not StringArray"))
}

fn u32_array<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a UInt32Array> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .with_context(|| format!("column {name} is not UInt32Array"))
}

fn u64_array<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a UInt64Array> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<UInt64Array>()
        .with_context(|| format!("column {name} is not UInt64Array"))
}

fn bool_array<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a BooleanArray> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<BooleanArray>()
        .with_context(|| format!("column {name} is not BooleanArray"))
}

fn nullable_string(array: &StringArray, row: usize) -> Option<String> {
    if array.is_null(row) {
        None
    } else {
        Some(array.value(row).to_string())
    }
}

fn nullable_u64(array: &UInt64Array, row: usize) -> Option<u64> {
    if array.is_null(row) {
        None
    } else {
        Some(array.value(row))
    }
}

#[derive(Debug, Serialize)]
pub struct QuerySummary {
    pub count: usize,
}
