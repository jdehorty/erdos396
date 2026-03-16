use crate::false_positive::RunWindowClassification;
use crate::run_collection::{RunRecord, RunWindowRecord};
use anyhow::Result;
use arrow::array::{ArrayRef, BooleanBuilder, StringBuilder, UInt32Builder, UInt64Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use std::sync::Arc;

pub fn maximal_runs_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("schema_version", DataType::Utf8, false),
        Field::new("run_id", DataType::Utf8, false),
        Field::new("checkpoint_id", DataType::Utf8, false),
        Field::new("source_label", DataType::Utf8, false),
        Field::new("archive_role", DataType::Utf8, false),
        Field::new("source_kind", DataType::Utf8, false),
        Field::new("server", DataType::Utf8, false),
        Field::new("source_path", DataType::Utf8, false),
        Field::new("source_record_ordinal", DataType::UInt64, false),
        Field::new("campaign_target_k", DataType::UInt32, false),
        Field::new("worker_id", DataType::UInt32, false),
        Field::new("run_start", DataType::UInt64, false),
        Field::new("run_end", DataType::UInt64, false),
        Field::new("run_length", DataType::UInt32, false),
        Field::new("range_bucket_t", DataType::UInt64, false),
        Field::new("is_unique_coverage", DataType::Boolean, false),
    ]))
}

pub fn run_windows_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("schema_version", DataType::Utf8, false),
        Field::new("window_id", DataType::Utf8, false),
        Field::new("parent_run_id", DataType::Utf8, false),
        Field::new("checkpoint_id", DataType::Utf8, false),
        Field::new("source_label", DataType::Utf8, false),
        Field::new("archive_role", DataType::Utf8, false),
        Field::new("server", DataType::Utf8, false),
        Field::new("source_path", DataType::Utf8, false),
        Field::new("campaign_target_k", DataType::UInt32, false),
        Field::new("worker_id", DataType::UInt32, false),
        Field::new("max_run_start", DataType::UInt64, false),
        Field::new("max_run_end", DataType::UInt64, false),
        Field::new("max_run_length", DataType::UInt32, false),
        Field::new("window_length", DataType::UInt32, false),
        Field::new("window_k", DataType::UInt32, false),
        Field::new("window_offset", DataType::UInt32, false),
        Field::new("window_start", DataType::UInt64, false),
        Field::new("window_end", DataType::UInt64, false),
        Field::new("window_n", DataType::UInt64, false),
        Field::new("range_bucket_t", DataType::UInt64, false),
        Field::new("recorded_status", DataType::Utf8, true),
        Field::new("recorded_failing_prime", DataType::UInt64, true),
        Field::new("recorded_demand", DataType::UInt64, true),
        Field::new("recorded_supply", DataType::UInt64, true),
        Field::new("is_unique_coverage", DataType::Boolean, false),
    ]))
}

pub fn run_window_classifications_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("schema_version", DataType::Utf8, false),
        Field::new("classification_id", DataType::Utf8, false),
        Field::new("window_id", DataType::Utf8, false),
        Field::new("window_length", DataType::UInt32, false),
        Field::new("window_k", DataType::UInt32, false),
        Field::new("window_start", DataType::UInt64, false),
        Field::new("window_end", DataType::UInt64, false),
        Field::new("window_n", DataType::UInt64, false),
        Field::new("audit_status", DataType::Utf8, false),
        Field::new("failing_prime", DataType::UInt64, true),
        Field::new("demand", DataType::UInt64, true),
        Field::new("supply", DataType::UInt64, true),
        Field::new("verifier_mode", DataType::Utf8, false),
        Field::new("verified_at_utc", DataType::Utf8, false),
        Field::new("build_git_hash", DataType::Utf8, false),
    ]))
}

pub fn maximal_runs_batch(records: &[RunRecord]) -> Result<RecordBatch> {
    let schema = maximal_runs_schema();
    RecordBatch::try_new(
        schema,
        vec![
            string_col(records.iter().map(|r| r.schema_version.as_str())),
            string_col(records.iter().map(|r| r.run_id.as_str())),
            string_col(records.iter().map(|r| r.checkpoint_id.as_str())),
            string_col(records.iter().map(|r| r.source_label.as_str())),
            string_col(records.iter().map(|r| r.archive_role.as_str())),
            string_col(records.iter().map(|r| r.source_kind.as_str())),
            string_col(records.iter().map(|r| r.server.as_str())),
            string_col(records.iter().map(|r| r.source_path.as_str())),
            u64_col(records.iter().map(|r| r.source_record_ordinal)),
            u32_col(records.iter().map(|r| r.campaign_target_k)),
            u32_col(records.iter().map(|r| r.worker_id)),
            u64_col(records.iter().map(|r| r.run_start)),
            u64_col(records.iter().map(|r| r.run_end)),
            u32_col(records.iter().map(|r| r.run_length)),
            u64_col(records.iter().map(|r| r.range_bucket_t)),
            bool_col(records.iter().map(|r| r.is_unique_coverage)),
        ],
    )
    .map_err(Into::into)
}

pub fn run_windows_batch(records: &[RunWindowRecord]) -> Result<RecordBatch> {
    let schema = run_windows_schema();
    RecordBatch::try_new(
        schema,
        vec![
            string_col(records.iter().map(|r| r.schema_version.as_str())),
            string_col(records.iter().map(|r| r.window_id.as_str())),
            string_col(records.iter().map(|r| r.parent_run_id.as_str())),
            string_col(records.iter().map(|r| r.checkpoint_id.as_str())),
            string_col(records.iter().map(|r| r.source_label.as_str())),
            string_col(records.iter().map(|r| r.archive_role.as_str())),
            string_col(records.iter().map(|r| r.server.as_str())),
            string_col(records.iter().map(|r| r.source_path.as_str())),
            u32_col(records.iter().map(|r| r.campaign_target_k)),
            u32_col(records.iter().map(|r| r.worker_id)),
            u64_col(records.iter().map(|r| r.max_run_start)),
            u64_col(records.iter().map(|r| r.max_run_end)),
            u32_col(records.iter().map(|r| r.max_run_length)),
            u32_col(records.iter().map(|r| r.window_length)),
            u32_col(records.iter().map(|r| r.window_k)),
            u32_col(records.iter().map(|r| r.window_offset)),
            u64_col(records.iter().map(|r| r.window_start)),
            u64_col(records.iter().map(|r| r.window_end)),
            u64_col(records.iter().map(|r| r.window_n)),
            u64_col(records.iter().map(|r| r.range_bucket_t)),
            opt_string_col(records.iter().map(|r| r.recorded_status.as_deref())),
            opt_u64_col(records.iter().map(|r| r.recorded_failing_prime)),
            opt_u64_col(records.iter().map(|r| r.recorded_demand)),
            opt_u64_col(records.iter().map(|r| r.recorded_supply)),
            bool_col(records.iter().map(|r| r.is_unique_coverage)),
        ],
    )
    .map_err(Into::into)
}

pub fn run_window_classifications_batch(
    records: &[RunWindowClassification],
) -> Result<RecordBatch> {
    let schema = run_window_classifications_schema();
    RecordBatch::try_new(
        schema,
        vec![
            string_col(records.iter().map(|r| r.schema_version.as_str())),
            string_col(records.iter().map(|r| r.classification_id.as_str())),
            string_col(records.iter().map(|r| r.window_id.as_str())),
            u32_col(records.iter().map(|r| r.window_length)),
            u32_col(records.iter().map(|r| r.window_k)),
            u64_col(records.iter().map(|r| r.window_start)),
            u64_col(records.iter().map(|r| r.window_end)),
            u64_col(records.iter().map(|r| r.window_n)),
            string_col(records.iter().map(|r| r.audit_status.as_str())),
            opt_u64_col(records.iter().map(|r| r.failing_prime)),
            opt_u64_col(records.iter().map(|r| r.demand)),
            opt_u64_col(records.iter().map(|r| r.supply)),
            string_col(records.iter().map(|r| r.verifier_mode.as_str())),
            string_col(records.iter().map(|r| r.verified_at_utc.as_str())),
            string_col(records.iter().map(|r| r.build_git_hash.as_str())),
        ],
    )
    .map_err(Into::into)
}

fn string_col<'a, I>(values: I) -> ArrayRef
where
    I: Iterator<Item = &'a str>,
{
    let items: Vec<&str> = values.collect();
    let mut builder = StringBuilder::new();
    for value in items {
        builder.append_value(value);
    }
    Arc::new(builder.finish())
}

fn opt_string_col<'a, I>(values: I) -> ArrayRef
where
    I: Iterator<Item = Option<&'a str>>,
{
    let items: Vec<Option<&str>> = values.collect();
    let mut builder = StringBuilder::new();
    for value in items {
        match value {
            Some(value) => builder.append_value(value),
            None => builder.append_null(),
        }
    }
    Arc::new(builder.finish())
}

fn u64_col<I>(values: I) -> ArrayRef
where
    I: Iterator<Item = u64>,
{
    let items: Vec<u64> = values.collect();
    let mut builder = UInt64Builder::new();
    for value in items {
        builder.append_value(value);
    }
    Arc::new(builder.finish())
}

fn opt_u64_col<I>(values: I) -> ArrayRef
where
    I: Iterator<Item = Option<u64>>,
{
    let items: Vec<Option<u64>> = values.collect();
    let mut builder = UInt64Builder::new();
    for value in items {
        match value {
            Some(value) => builder.append_value(value),
            None => builder.append_null(),
        }
    }
    Arc::new(builder.finish())
}

fn u32_col<I>(values: I) -> ArrayRef
where
    I: Iterator<Item = u32>,
{
    let items: Vec<u32> = values.collect();
    let mut builder = UInt32Builder::new();
    for value in items {
        builder.append_value(value);
    }
    Arc::new(builder.finish())
}

fn bool_col<I>(values: I) -> ArrayRef
where
    I: Iterator<Item = bool>,
{
    let items: Vec<bool> = values.collect();
    let mut builder = BooleanBuilder::new();
    for value in items {
        builder.append_value(value);
    }
    Arc::new(builder.finish())
}
