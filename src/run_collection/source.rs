use crate::checkpoint::{Checkpoint, RunInfo};
use crate::run_collection::expand::{expand_run_windows, RecordedEvent};
use crate::run_collection::parquet::{maximal_run_writer, prepare_build_output, run_window_writer};
use crate::run_collection::{
    BuildRunCorpusConfig, BuildRunCorpusStats, RunRecord, RANGE_BUCKET_SIZE, SCHEMA_VERSION,
};
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceKind {
    CheckpointV2,
    RunLogV3,
}

impl SourceKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::CheckpointV2 => "checkpoint_v2",
            Self::RunLogV3 => "run_log_v3",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArchiveRole {
    Canonical,
    Overlap,
    Provenance,
}

impl ArchiveRole {
    fn as_str(self) -> &'static str {
        match self {
            Self::Canonical => "canonical",
            Self::Overlap => "overlap",
            Self::Provenance => "provenance",
        }
    }

    fn is_unique_coverage(self) -> bool {
        matches!(self, Self::Canonical)
    }
}

#[derive(Debug, Clone)]
struct SourceSpec {
    rel_path: &'static str,
    source_label: &'static str,
    archive_role: ArchiveRole,
    source_kind: SourceKind,
    server: &'static str,
    recursive: bool,
    checkpoint_prefix: Option<&'static str>,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct EventKey {
    checkpoint_id: String,
    window_start: u64,
    window_length: u32,
    window_n: u64,
}

pub fn build_run_corpus(config: &BuildRunCorpusConfig) -> Result<BuildRunCorpusStats> {
    if config.source_root.as_os_str().is_empty() {
        return Err(anyhow!("source_root must not be empty"));
    }
    if config.min_length == 0 || config.max_length < config.min_length {
        return Err(anyhow!(
            "invalid length bounds: min={} max={}",
            config.min_length,
            config.max_length
        ));
    }

    prepare_build_output(&config.output_root)?;
    let mut maximal_writer = maximal_run_writer(&config.output_root);
    let mut window_writer = run_window_writer(&config.output_root);
    let mut stats = BuildRunCorpusStats::default();

    for spec in source_catalog() {
        stats.sources_scanned += 1;

        if !config.include_overlaps && spec.archive_role != ArchiveRole::Canonical {
            continue;
        }

        let source_dir = config.source_root.join(spec.rel_path);
        if !source_dir.exists() {
            stats.skipped_sources.push(spec.rel_path.to_string());
            continue;
        }

        stats.sources_used += 1;
        match spec.source_kind {
            SourceKind::CheckpointV2 => ingest_v2_source(
                config,
                &spec,
                &source_dir,
                &mut maximal_writer,
                &mut window_writer,
                &mut stats,
            )?,
            SourceKind::RunLogV3 => ingest_v3_source(
                config,
                &spec,
                &source_dir,
                &mut maximal_writer,
                &mut window_writer,
                &mut stats,
            )?,
        }
    }

    for extra_dir in &config.extra_v3_dirs {
        ingest_extra_v3_source(
            config,
            extra_dir,
            &mut maximal_writer,
            &mut window_writer,
            &mut stats,
        )?;
    }

    maximal_writer.finish()?;
    window_writer.finish()?;
    Ok(stats)
}

fn ingest_extra_v3_source(
    config: &BuildRunCorpusConfig,
    source_dir: &Path,
    maximal_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<RunRecord>,
    window_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<
        crate::run_collection::RunWindowRecord,
    >,
    stats: &mut BuildRunCorpusStats,
) -> Result<()> {
    stats.sources_scanned += 1;
    if !source_dir.exists() {
        stats
            .skipped_sources
            .push(source_dir.to_string_lossy().into_owned());
        return Ok(());
    }
    stats.sources_used += 1;

    let source_label = source_dir
        .file_name()
        .and_then(|name| name.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| "local_v3".to_string());
    let spec = SourceSpec {
        rel_path: "",
        source_label: Box::leak(source_label.into_boxed_str()),
        archive_role: ArchiveRole::Canonical,
        source_kind: SourceKind::RunLogV3,
        server: "local",
        recursive: false,
        checkpoint_prefix: None,
    };

    ingest_v3_source(
        config,
        &spec,
        source_dir,
        maximal_writer,
        window_writer,
        stats,
    )
}

fn source_catalog() -> Vec<SourceSpec> {
    vec![
        SourceSpec {
            rel_path: "server_a/verify_k10",
            source_label: "server_a_k10_verify_1M_18B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_a/k11_search",
            source_label: "server_a_k11_18B_200B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_b_backup_2026-01-30/k11_search_200_400",
            source_label: "server_b_k11_200B_400B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-b",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_a/k11_search_400_600",
            source_label: "server_a_k11_400B_600B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_b_backup_2026-01-30/k11_search_600_800",
            source_label: "server_b_k11_600B_800B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-b",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_b_backup_2026-01-30/k11_search_800_850",
            source_label: "server_b_k11_800B_850B",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-b",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k11_850B-1T",
            source_label: "server_a_k11_850B_1T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_b_k11_1T-2T",
            source_label: "server_b_k11_1T_2T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-b",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k11_2T-2.5T",
            source_label: "server_a_k11_2T_2_5T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k12_2.5T-5T",
            source_label: "server_a_k12_2_5T_5T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_b_k12_5T-10T",
            source_label: "server_b_k12_5T_10T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-b",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k12_10T-15T",
            source_label: "server_a_k12_10T_15T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k13_15T-25T",
            source_label: "server_a_k13_15T_25T",
            archive_role: ArchiveRole::Canonical,
            source_kind: SourceKind::RunLogV3,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k12_gap_5T-6.15T",
            source_label: "server_a_k12_gap_5T_6_15T",
            archive_role: ArchiveRole::Overlap,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: true,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "archives/server_a_k13_15T-25T_v2_final",
            source_label: "server_a_k13_15T_25T_v2_final",
            archive_role: ArchiveRole::Provenance,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: None,
        },
        SourceSpec {
            rel_path: "server_a",
            source_label: "server_a_k9_history",
            archive_role: ArchiveRole::Provenance,
            source_kind: SourceKind::CheckpointV2,
            server: "server-a",
            recursive: false,
            checkpoint_prefix: Some("checkpoint_k9_"),
        },
    ]
}

fn ingest_v2_source(
    config: &BuildRunCorpusConfig,
    spec: &SourceSpec,
    source_dir: &Path,
    maximal_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<RunRecord>,
    window_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<
        crate::run_collection::RunWindowRecord,
    >,
    stats: &mut BuildRunCorpusStats,
) -> Result<()> {
    for checkpoint_path in checkpoint_files(source_dir, spec.recursive, spec.checkpoint_prefix)? {
        let checkpoint = Checkpoint::load(&checkpoint_path)
            .with_context(|| checkpoint_path.display().to_string())?;
        let checkpoint_id = relative_path_string(&config.source_root, &checkpoint_path)?;
        let events = recorded_event_map(&checkpoint, &checkpoint_id);

        for (ordinal, run) in checkpoint.significant_runs.iter().enumerate() {
            let run_record = run_record_from_checkpoint_run(
                spec,
                &checkpoint,
                &checkpoint_path,
                &checkpoint_id,
                ordinal as u64,
                run,
            )?;
            let windows = expand_run_windows(
                &run_record,
                config.min_length,
                config.max_length,
                |window_start, window_length, window_n| {
                    events
                        .get(&EventKey {
                            checkpoint_id: checkpoint_id.clone(),
                            window_start,
                            window_length,
                            window_n,
                        })
                        .cloned()
                },
            );
            maximal_writer.push(run_record)?;
            stats.maximal_runs_written += 1;
            for window in windows {
                window_writer.push(window)?;
                stats.run_windows_written += 1;
            }
        }
    }

    Ok(())
}

fn ingest_v3_source(
    config: &BuildRunCorpusConfig,
    spec: &SourceSpec,
    source_dir: &Path,
    maximal_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<RunRecord>,
    window_writer: &mut crate::run_collection::parquet::PartitionedDatasetWriter<
        crate::run_collection::RunWindowRecord,
    >,
    stats: &mut BuildRunCorpusStats,
) -> Result<()> {
    let checkpoints = checkpoint_files(source_dir, false, spec.checkpoint_prefix)?;
    let mut contexts = HashMap::new();

    for checkpoint_path in checkpoints {
        let checkpoint = Checkpoint::load(&checkpoint_path)
            .with_context(|| checkpoint_path.display().to_string())?;
        let worker_id = checkpoint
            .worker_id
            .map(|id| id as u32)
            .unwrap_or_else(|| parse_worker_id(&checkpoint_path).unwrap_or(0));
        let checkpoint_id = relative_path_string(&config.source_root, &checkpoint_path)?;
        let events = recorded_event_map(&checkpoint, &checkpoint_id);
        contexts.insert(
            worker_id,
            WorkerContext {
                checkpoint,
                checkpoint_id,
                events,
            },
        );
    }

    let nested_run_log_dir = source_dir.join("run_logs");
    let run_log_dir = if nested_run_log_dir.exists() {
        nested_run_log_dir
    } else {
        source_dir.to_path_buf()
    };
    for run_log_path in run_log_files(&run_log_dir)? {
        let worker_id = parse_worker_id(&run_log_path)
            .ok_or_else(|| anyhow!("could not parse worker id from {}", run_log_path.display()))?;
        let context = contexts
            .get(&worker_id)
            .ok_or_else(|| anyhow!("missing checkpoint for worker {}", worker_id))?;
        let source_path = relative_path_string(&config.source_root, &run_log_path)?;
        let file = File::open(&run_log_path)
            .with_context(|| format!("open {}", run_log_path.display()))?;
        let reader = BufReader::new(file);
        for (ordinal, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let run: RunInfo = serde_json::from_str(&line).with_context(|| {
                format!("parse {} line {}", run_log_path.display(), ordinal + 1)
            })?;
            let run_record = run_record_from_v3_line(
                spec,
                &context.checkpoint,
                &context.checkpoint_id,
                &source_path,
                ordinal as u64,
                &run,
            )?;
            let windows = expand_run_windows(
                &run_record,
                config.min_length,
                config.max_length,
                |window_start, window_length, window_n| {
                    context
                        .events
                        .get(&EventKey {
                            checkpoint_id: context.checkpoint_id.clone(),
                            window_start,
                            window_length,
                            window_n,
                        })
                        .cloned()
                },
            );
            maximal_writer.push(run_record)?;
            stats.maximal_runs_written += 1;
            for window in windows {
                window_writer.push(window)?;
                stats.run_windows_written += 1;
            }
        }
    }

    Ok(())
}

fn checkpoint_files(dir: &Path, recursive: bool, prefix: Option<&str>) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    if recursive {
        for entry in WalkDir::new(dir) {
            let entry = entry?;
            let path = entry.path();
            if is_checkpoint_file(path, prefix) {
                files.push(path.to_path_buf());
            }
        }
    } else {
        for entry in
            std::fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))?
        {
            let entry = entry?;
            let path = entry.path();
            if is_checkpoint_file(&path, prefix) {
                files.push(path);
            }
        }
    }
    files.sort();
    Ok(files)
}

fn run_log_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    if !dir.exists() {
        return Ok(files);
    }
    for entry in std::fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file()
            && path
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name| name.starts_with("runs_k") && name.ends_with(".jsonl"))
                .unwrap_or(false)
        {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn is_checkpoint_file(path: &Path, prefix: Option<&str>) -> bool {
    path.is_file()
        && path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| {
                let expected = prefix.unwrap_or("checkpoint_k");
                name.starts_with(expected) && name.ends_with(".json")
            })
            .unwrap_or(false)
}

fn run_record_from_checkpoint_run(
    spec: &SourceSpec,
    checkpoint: &Checkpoint,
    checkpoint_path: &Path,
    checkpoint_id: &str,
    source_record_ordinal: u64,
    run: &RunInfo,
) -> Result<RunRecord> {
    let run_end = run
        .start
        .checked_add(run.length as u64)
        .and_then(|value| value.checked_sub(1))
        .ok_or_else(|| anyhow!("overflow computing run end for {}", checkpoint_id))?;
    let worker_id = checkpoint
        .worker_id
        .map(|value| value as u32)
        .unwrap_or_else(|| parse_worker_id(checkpoint_path).unwrap_or(0));

    Ok(RunRecord {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: format!(
            "{checkpoint_id}|{source_record_ordinal}|{}|{}",
            run.start, run.length
        ),
        checkpoint_id: checkpoint_id.to_string(),
        source_label: spec.source_label.to_string(),
        archive_role: spec.archive_role.as_str().to_string(),
        source_kind: spec.source_kind.as_str().to_string(),
        server: spec.server.to_string(),
        source_path: checkpoint_id.to_string(),
        source_record_ordinal,
        campaign_target_k: checkpoint.target_k,
        worker_id,
        run_start: run.start,
        run_end,
        run_length: run.length as u32,
        range_bucket_t: run_end / RANGE_BUCKET_SIZE,
        is_unique_coverage: spec.archive_role.is_unique_coverage(),
    })
}

fn run_record_from_v3_line(
    spec: &SourceSpec,
    checkpoint: &Checkpoint,
    checkpoint_id: &str,
    source_path: &str,
    source_record_ordinal: u64,
    run: &RunInfo,
) -> Result<RunRecord> {
    let run_end = run
        .start
        .checked_add(run.length as u64)
        .and_then(|value| value.checked_sub(1))
        .ok_or_else(|| anyhow!("overflow computing run end for {}", source_path))?;
    let worker_id = checkpoint.worker_id.map(|value| value as u32).unwrap_or(0);

    Ok(RunRecord {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: format!(
            "{source_path}|{source_record_ordinal}|{}|{}",
            run.start, run.length
        ),
        checkpoint_id: checkpoint_id.to_string(),
        source_label: spec.source_label.to_string(),
        archive_role: spec.archive_role.as_str().to_string(),
        source_kind: spec.source_kind.as_str().to_string(),
        server: spec.server.to_string(),
        source_path: source_path.to_string(),
        source_record_ordinal,
        campaign_target_k: checkpoint.target_k,
        worker_id,
        run_start: run.start,
        run_end,
        run_length: run.length as u32,
        range_bucket_t: run_end / RANGE_BUCKET_SIZE,
        is_unique_coverage: spec.archive_role.is_unique_coverage(),
    })
}

fn recorded_event_map(
    checkpoint: &Checkpoint,
    checkpoint_id: &str,
) -> HashMap<EventKey, RecordedEvent> {
    let mut events = HashMap::new();
    let target_run_length = checkpoint.target_k + 1;

    for &candidate in &checkpoint.candidates {
        if let Some(key) = event_key(checkpoint_id, target_run_length, candidate) {
            events.insert(
                key,
                RecordedEvent {
                    status: "candidate".to_string(),
                    failing_prime: None,
                    demand: None,
                    supply: None,
                },
            );
        }
    }

    for &witness in &checkpoint.witnesses {
        if let Some(key) = event_key(checkpoint_id, target_run_length, witness) {
            events.insert(
                key,
                RecordedEvent {
                    status: "witness".to_string(),
                    failing_prime: None,
                    demand: None,
                    supply: None,
                },
            );
        }
    }

    for detail in &checkpoint.false_positive_details {
        if let Some(key) = event_key(checkpoint_id, detail.run_length as u32, detail.position) {
            events.insert(
                key,
                RecordedEvent {
                    status: "false_positive".to_string(),
                    failing_prime: Some(detail.failing_prime),
                    demand: Some(detail.demand),
                    supply: Some(detail.supply),
                },
            );
        }
    }

    for &false_positive in &checkpoint.false_positives {
        if let Some(key) = event_key(checkpoint_id, target_run_length, false_positive) {
            events.entry(key).or_insert_with(|| RecordedEvent {
                status: "false_positive".to_string(),
                failing_prime: None,
                demand: None,
                supply: None,
            });
        }
    }

    events
}

fn event_key(checkpoint_id: &str, run_length: u32, n: u64) -> Option<EventKey> {
    n.checked_add(1)
        .and_then(|value| value.checked_sub(run_length as u64))
        .map(|window_start| EventKey {
            checkpoint_id: checkpoint_id.to_string(),
            window_start,
            window_length: run_length,
            window_n: n,
        })
}

fn relative_path_string(root: &Path, path: &Path) -> Result<String> {
    Ok(path
        .strip_prefix(root)
        .with_context(|| format!("strip prefix {} from {}", root.display(), path.display()))?
        .to_string_lossy()
        .replace('\\', "/"))
}

fn parse_worker_id(path: &Path) -> Option<u32> {
    let name = path.file_name()?.to_str()?;
    let marker = "_w";
    let start = name.rfind(marker)? + marker.len();
    let digits: String = name[start..]
        .chars()
        .take_while(|ch| ch.is_ascii_digit())
        .collect();
    digits.parse::<u32>().ok()
}

struct WorkerContext {
    checkpoint: Checkpoint,
    checkpoint_id: String,
    events: HashMap<EventKey, RecordedEvent>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::checkpoint::RunInfo;
    use crate::run_collection::parquet::{read_run_windows, RunWindowFilter};
    use tempfile::tempdir;

    #[test]
    fn canonical_build_excludes_overlap_by_default() {
        let fixture = tempdir().unwrap();
        write_v2_checkpoint(
            fixture
                .path()
                .join("archives/server_b_k12_5T-10T/checkpoint_k12_w00.json"),
            12,
            0,
            vec![RunInfo {
                start: 5_000_000_000_000,
                length: 6,
                is_witness: None,
                failing_prime: None,
            }],
            vec![],
            vec![],
        );
        write_v2_checkpoint(
            fixture
                .path()
                .join("archives/server_a_k12_gap_5T-6.15T/gap_00/checkpoint_k12_w00.json"),
            12,
            0,
            vec![RunInfo {
                start: 5_100_000_000_000,
                length: 6,
                is_witness: None,
                failing_prime: None,
            }],
            vec![],
            vec![],
        );

        let output = tempdir().unwrap();
        build_run_corpus(&BuildRunCorpusConfig {
            source_root: fixture.path().to_path_buf(),
            output_root: output.path().to_path_buf(),
            min_length: 6,
            max_length: 14,
            include_overlaps: false,
        })
        .unwrap();

        let windows = read_run_windows(
            output.path(),
            &RunWindowFilter {
                window_length: Some(6),
                ..RunWindowFilter::default()
            },
        )
        .unwrap();

        assert_eq!(windows.len(), 1);
        assert!(windows[0].is_unique_coverage);
    }

    #[test]
    fn overlap_build_includes_overlap_when_requested() {
        let fixture = tempdir().unwrap();
        write_v2_checkpoint(
            fixture
                .path()
                .join("archives/server_b_k12_5T-10T/checkpoint_k12_w00.json"),
            12,
            0,
            vec![RunInfo {
                start: 5_000_000_000_000,
                length: 6,
                is_witness: None,
                failing_prime: None,
            }],
            vec![],
            vec![],
        );
        write_v2_checkpoint(
            fixture
                .path()
                .join("archives/server_a_k12_gap_5T-6.15T/gap_00/checkpoint_k12_w00.json"),
            12,
            0,
            vec![RunInfo {
                start: 5_100_000_000_000,
                length: 6,
                is_witness: None,
                failing_prime: None,
            }],
            vec![],
            vec![],
        );

        let output = tempdir().unwrap();
        build_run_corpus(&BuildRunCorpusConfig {
            source_root: fixture.path().to_path_buf(),
            output_root: output.path().to_path_buf(),
            min_length: 6,
            max_length: 14,
            include_overlaps: true,
        })
        .unwrap();

        let windows = read_run_windows(
            output.path(),
            &RunWindowFilter {
                window_length: Some(6),
                ..RunWindowFilter::default()
            },
        )
        .unwrap();

        assert_eq!(windows.len(), 2);
        assert!(windows.iter().any(|window| !window.is_unique_coverage));
    }

    #[test]
    fn v3_run_logs_attach_recorded_witness_only_to_matching_window() {
        let fixture = tempdir().unwrap();
        write_v3_fixture(
            fixture.path().join("archives/server_a_k13_15T-25T"),
            18_253_129_921_829,
            14,
            0,
            vec![18_253_129_921_842],
            vec![18_253_129_921_842],
        );

        let output = tempdir().unwrap();
        build_run_corpus(&BuildRunCorpusConfig {
            source_root: fixture.path().to_path_buf(),
            output_root: output.path().to_path_buf(),
            min_length: 6,
            max_length: 14,
            include_overlaps: false,
        })
        .unwrap();

        let len14 = read_run_windows(
            output.path(),
            &RunWindowFilter {
                window_length: Some(14),
                ..RunWindowFilter::default()
            },
        )
        .unwrap();
        assert_eq!(len14.len(), 1);
        assert_eq!(len14[0].recorded_status.as_deref(), Some("witness"));

        let len13 = read_run_windows(
            output.path(),
            &RunWindowFilter {
                window_length: Some(13),
                ..RunWindowFilter::default()
            },
        )
        .unwrap();
        assert!(len13.iter().all(|window| window.recorded_status.is_none()));
    }

    fn write_v2_checkpoint(
        path: PathBuf,
        target_k: u32,
        worker_id: usize,
        runs: Vec<RunInfo>,
        candidates: Vec<u64>,
        witnesses: Vec<u64>,
    ) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let mut checkpoint = Checkpoint::new_worker(target_k, 0, 1_000, worker_id);
        checkpoint.version = 2;
        checkpoint.significant_runs = runs;
        checkpoint.candidates = candidates;
        checkpoint.witnesses = witnesses;
        checkpoint.save(&path).unwrap();
    }

    fn write_v3_fixture(
        root: PathBuf,
        run_start: u64,
        run_length: usize,
        worker_id: usize,
        candidates: Vec<u64>,
        witnesses: Vec<u64>,
    ) {
        std::fs::create_dir_all(root.join("run_logs")).unwrap();
        let checkpoint_path = root.join(format!("checkpoint_k13_w{worker_id:02}.json"));
        let mut checkpoint = Checkpoint::new_worker(13, run_start, run_start + 1_000, worker_id);
        checkpoint.version = 3;
        checkpoint.candidates = candidates;
        checkpoint.witnesses = witnesses;
        checkpoint.save(&checkpoint_path).unwrap();

        let run_path = root.join(format!("run_logs/runs_k13_w{worker_id:02}.jsonl"));
        std::fs::write(
            run_path,
            serde_json::to_string(&RunInfo {
                start: run_start,
                length: run_length,
                is_witness: None,
                failing_prime: None,
            })
            .unwrap()
                + "\n",
        )
        .unwrap();
    }
}
