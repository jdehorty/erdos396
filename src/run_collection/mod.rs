use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub mod expand;
pub mod parquet;
pub mod schema;
pub mod source;

pub const SCHEMA_VERSION: &str = "run_corpus_v1";
pub const DEFAULT_OUTPUT_ROOT: &str = "data/run_corpus/v1";
pub const RANGE_BUCKET_SIZE: u64 = 1_000_000_000_000;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RunRecord {
    pub schema_version: String,
    pub run_id: String,
    pub checkpoint_id: String,
    pub source_label: String,
    pub archive_role: String,
    pub source_kind: String,
    pub server: String,
    pub source_path: String,
    pub source_record_ordinal: u64,
    pub campaign_target_k: u32,
    pub worker_id: u32,
    pub run_start: u64,
    pub run_end: u64,
    pub run_length: u32,
    pub range_bucket_t: u64,
    pub is_unique_coverage: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RunWindowRecord {
    pub schema_version: String,
    pub window_id: String,
    pub parent_run_id: String,
    pub checkpoint_id: String,
    pub source_label: String,
    pub archive_role: String,
    pub server: String,
    pub source_path: String,
    pub campaign_target_k: u32,
    pub worker_id: u32,
    pub max_run_start: u64,
    pub max_run_end: u64,
    pub max_run_length: u32,
    pub window_length: u32,
    pub window_k: u32,
    pub window_offset: u32,
    pub window_start: u64,
    pub window_end: u64,
    pub window_n: u64,
    pub range_bucket_t: u64,
    pub recorded_status: Option<String>,
    pub recorded_failing_prime: Option<u64>,
    pub recorded_demand: Option<u64>,
    pub recorded_supply: Option<u64>,
    pub is_unique_coverage: bool,
}

#[derive(Debug, Clone)]
pub struct BuildRunCorpusConfig {
    pub source_root: PathBuf,
    pub output_root: PathBuf,
    pub min_length: usize,
    pub max_length: usize,
    pub include_overlaps: bool,
    pub extra_v3_dirs: Vec<PathBuf>,
}

impl Default for BuildRunCorpusConfig {
    fn default() -> Self {
        Self {
            source_root: PathBuf::new(),
            output_root: PathBuf::from(DEFAULT_OUTPUT_ROOT),
            min_length: 6,
            max_length: 14,
            include_overlaps: false,
            extra_v3_dirs: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct BuildRunCorpusStats {
    pub sources_scanned: usize,
    pub sources_used: usize,
    pub maximal_runs_written: u64,
    pub run_windows_written: u64,
    pub skipped_sources: Vec<String>,
}

pub fn build_run_corpus(config: &BuildRunCorpusConfig) -> Result<BuildRunCorpusStats> {
    source::build_run_corpus(config)
}
