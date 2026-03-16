use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub mod classify;

pub use classify::{classify_run_windows, verify_window, ClassificationStats};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RunWindowClassification {
    pub schema_version: String,
    pub classification_id: String,
    pub window_id: String,
    pub window_length: u32,
    pub window_k: u32,
    pub window_start: u64,
    pub window_end: u64,
    pub window_n: u64,
    pub audit_status: String,
    pub failing_prime: Option<u64>,
    pub demand: Option<u64>,
    pub supply: Option<u64>,
    pub verifier_mode: String,
    pub verified_at_utc: String,
    pub build_git_hash: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VerificationAudit {
    pub is_witness: bool,
    pub failing_prime: Option<u64>,
    pub demand: Option<u64>,
    pub supply: Option<u64>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ClassificationStatusSource {
    Recorded,
    Unclassified,
    All,
}

impl ClassificationStatusSource {
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "recorded" => Ok(Self::Recorded),
            "unclassified" => Ok(Self::Unclassified),
            "all" => Ok(Self::All),
            _ => Err(anyhow::anyhow!(
                "invalid status-source {value}; expected recorded|unclassified|all"
            )),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ClassificationConfig {
    pub output_root: PathBuf,
    pub window_length: Option<u32>,
    pub start_n: Option<u64>,
    pub end_n: Option<u64>,
    pub status_source: ClassificationStatusSource,
    pub workers: usize,
    pub limit: Option<usize>,
}
