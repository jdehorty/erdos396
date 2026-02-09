//! Checkpoint management for long-running searches.
//!
//! Provides automatic checkpointing at configurable intervals,
//! allowing searches to be paused and resumed.

use crate::Result;
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

/// Checkpoint data for a search worker.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Checkpoint {
    /// Target k value (looking for run of k+1)
    pub target_k: u32,

    /// Start of search range
    pub start: u64,

    /// End of search range (exclusive)
    pub end: u64,

    /// Current position in search
    pub current_pos: u64,

    /// Number of integers checked so far
    pub checked: u64,

    /// Number of Governor Set members found
    pub governor_count: u64,

    /// Longest run of consecutive Governor members seen
    pub longest_run: usize,

    /// Starting position of the longest run
    pub longest_run_start: u64,

    /// Positions where runs of target length were found (candidates)
    pub candidates: Vec<u64>,

    /// Verified witnesses (candidates that passed full verification)
    pub witnesses: Vec<u64>,

    /// False positives (candidates that failed verification) - legacy field
    pub false_positives: Vec<u64>,

    /// Detailed false positive information (new in v2)
    #[serde(default)]
    pub false_positive_details: Vec<FalsePositiveInfo>,

    /// Run length distribution: maps run_length -> count
    /// Tracks how many runs of each length were found (for runs >= 2)
    #[serde(default)]
    pub run_distribution: HashMap<usize, u64>,

    /// All significant runs found (runs of length >= significant_run_threshold)
    #[serde(default)]
    pub significant_runs: Vec<RunInfo>,

    /// Threshold for recording individual runs (default: 6)
    #[serde(default = "default_significant_threshold")]
    pub significant_run_threshold: usize,

    /// Counter-examples found by safety-net mode
    #[serde(default)]
    pub counter_examples: Vec<CounterExampleInfo>,

    /// Number of sliding windows checked by safety-net mode
    #[serde(default)]
    pub safety_net_windows_checked: u64,

    /// Number of safety-net alerts (windows passing small-prime test but not governor runs)
    #[serde(default)]
    pub safety_net_alerts: u64,

    /// Timestamp of last checkpoint
    pub timestamp: DateTime<Utc>,

    /// Worker ID (for parallel searches)
    #[serde(default)]
    pub worker_id: Option<usize>,

    /// Version for compatibility checking
    #[serde(default = "default_version")]
    pub version: u32,
}

fn default_version() -> u32 {
    2  // Bumped for new statistics fields
}

fn default_significant_threshold() -> usize {
    6  // Record all runs of 6+ consecutive governors
}

/// Detailed information about a false positive (Governor run that fails verification).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FalsePositiveInfo {
    /// Ending position n of the run
    pub position: u64,
    /// Length of the Governor run
    pub run_length: usize,
    /// The prime that caused failure
    pub failing_prime: u64,
    /// p-adic demand (exponent needed in product)
    pub demand: u64,
    /// p-adic supply (exponent available in C(2n,n))
    pub supply: u64,
}

/// Information about a potential counter-example found by safety-net mode.
///
/// A counter-example is a window of k+1 consecutive numbers that passes
/// the small-prime divisibility test but is NOT a full governor run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CounterExampleInfo {
    /// Starting n of the window (highest element, matching witness convention)
    pub n: u64,
    /// How many of the k+1 positions are NOT governors
    pub non_governor_count: usize,
    /// Which positions in the window are not governors (offsets from n)
    pub non_governor_positions: Vec<u64>,
    /// Result of full verification (None if not yet run)
    pub full_verify_result: Option<bool>,
    /// If full verification failed, which prime
    pub full_verify_failing_prime: Option<u64>,
}

/// Information about a significant run (for scientific analysis).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunInfo {
    /// Starting position of the run
    pub start: u64,
    /// Length of the run (number of consecutive governors)
    pub length: usize,
    /// Whether this run was verified as a witness
    pub is_witness: Option<bool>,
    /// If verified and failed, which prime
    pub failing_prime: Option<u64>,
}

impl Checkpoint {
    /// Create a new checkpoint for starting a fresh search.
    pub fn new(target_k: u32, start: u64, end: u64) -> Self {
        Self {
            target_k,
            start,
            end,
            current_pos: start,
            checked: 0,
            governor_count: 0,
            longest_run: 0,
            longest_run_start: start,
            candidates: Vec::new(),
            witnesses: Vec::new(),
            false_positives: Vec::new(),
            false_positive_details: Vec::new(),
            run_distribution: HashMap::new(),
            significant_runs: Vec::new(),
            significant_run_threshold: default_significant_threshold(),
            counter_examples: Vec::new(),
            safety_net_windows_checked: 0,
            safety_net_alerts: 0,
            timestamp: Utc::now(),
            worker_id: None,
            version: 2,
        }
    }

    /// Record a run in the distribution and optionally as a significant run.
    pub fn record_run(&mut self, start: u64, length: usize) {
        // Always update distribution for runs of 2+
        if length >= 2 {
            *self.run_distribution.entry(length).or_insert(0) += 1;
        }

        // Record significant runs individually
        if length >= self.significant_run_threshold {
            self.significant_runs.push(RunInfo {
                start,
                length,
                is_witness: None,
                failing_prime: None,
            });
        }
    }

    /// Record detailed false positive information.
    pub fn record_false_positive_detail(
        &mut self,
        position: u64,
        run_length: usize,
        failing_prime: u64,
        demand: u64,
        supply: u64,
    ) {
        self.false_positive_details.push(FalsePositiveInfo {
            position,
            run_length,
            failing_prime,
            demand,
            supply,
        });
    }

    /// Create a checkpoint with worker ID for parallel searches.
    pub fn new_worker(target_k: u32, start: u64, end: u64, worker_id: usize) -> Self {
        let mut cp = Self::new(target_k, start, end);
        cp.worker_id = Some(worker_id);
        cp
    }

    /// Update the checkpoint timestamp.
    pub fn touch(&mut self) {
        self.timestamp = Utc::now();
    }

    /// Calculate progress as a percentage.
    pub fn progress_percent(&self) -> f64 {
        if self.end <= self.start {
            return 100.0;
        }
        let total = (self.end - self.start) as f64;
        let done = (self.current_pos - self.start) as f64;
        (done / total) * 100.0
    }

    /// Calculate Governor Set density observed so far.
    pub fn observed_density(&self) -> f64 {
        if self.checked == 0 {
            return 0.0;
        }
        (self.governor_count as f64) / (self.checked as f64)
    }

    /// Check if search is complete.
    pub fn is_complete(&self) -> bool {
        self.current_pos >= self.end
    }

    /// Load checkpoint from a JSON file.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let reader = BufReader::new(file);
        let checkpoint: Self = serde_json::from_reader(reader)?;
        Ok(checkpoint)
    }

    /// Save checkpoint to a JSON file.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path.as_ref())?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    /// Save to a temporary file then rename (atomic on most filesystems).
    pub fn save_atomic<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let tmp_path = path.with_extension("tmp");

        // Write to temp file
        self.save(&tmp_path)?;

        // Atomic rename
        fs::rename(&tmp_path, path)?;

        Ok(())
    }
}

/// Manager for automatic checkpointing during searches.
pub struct CheckpointManager {
    /// Path to checkpoint file
    path: PathBuf,

    /// Interval (in numbers checked) between checkpoints
    interval: u64,

    /// Counter for checkpoint interval
    counter: u64,

    /// Current checkpoint state
    checkpoint: Checkpoint,
}

impl CheckpointManager {
    /// Create a new checkpoint manager.
    pub fn new<P: AsRef<Path>>(path: P, interval: u64, checkpoint: Checkpoint) -> Self {
        Self {
            path: path.as_ref().to_path_buf(),
            interval,
            counter: 0,
            checkpoint,
        }
    }

    /// Try to load existing checkpoint, or create new one.
    pub fn load_or_new<P: AsRef<Path>>(
        path: P,
        interval: u64,
        target_k: u32,
        start: u64,
        end: u64,
    ) -> Result<Self> {
        let path = path.as_ref();

        let checkpoint = if path.exists() {
            log::info!("Loading checkpoint from {:?}", path);
            Checkpoint::load(path)?
        } else {
            log::info!("Creating new checkpoint");
            Checkpoint::new(target_k, start, end)
        };

        Ok(Self {
            path: path.to_path_buf(),
            interval,
            counter: 0,
            checkpoint,
        })
    }

    /// Get current checkpoint state.
    pub fn checkpoint(&self) -> &Checkpoint {
        &self.checkpoint
    }

    /// Get mutable reference to checkpoint.
    pub fn checkpoint_mut(&mut self) -> &mut Checkpoint {
        &mut self.checkpoint
    }

    /// Update progress and save if interval reached.
    ///
    /// Returns true if checkpoint was saved.
    pub fn update(&mut self, current_pos: u64, governor_found: bool) -> Result<bool> {
        self.checkpoint.current_pos = current_pos;
        self.checkpoint.checked += 1;
        if governor_found {
            self.checkpoint.governor_count += 1;
        }

        self.counter += 1;

        if self.counter >= self.interval {
            self.save()?;
            self.counter = 0;
            return Ok(true);
        }

        Ok(false)
    }

    /// Force save checkpoint.
    pub fn save(&mut self) -> Result<()> {
        self.checkpoint.touch();
        self.checkpoint.save_atomic(&self.path)?;
        Ok(())
    }

    /// Get the checkpoint file path.
    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Aggregate checkpoints from multiple workers into a summary.
pub fn aggregate_checkpoints(checkpoints: &[Checkpoint]) -> AggregatedStats {
    let mut stats = AggregatedStats::default();

    for cp in checkpoints {
        stats.total_checked += cp.checked;
        stats.total_governor_count += cp.governor_count;
        stats.candidates.extend(cp.candidates.iter().copied());
        stats.witnesses.extend(cp.witnesses.iter().copied());
        stats.false_positives.extend(cp.false_positives.iter().copied());
        stats.false_positive_details.extend(cp.false_positive_details.iter().cloned());
        stats.significant_runs.extend(cp.significant_runs.iter().cloned());
        stats.counter_examples.extend(cp.counter_examples.iter().cloned());
        stats.safety_net_windows_checked += cp.safety_net_windows_checked;
        stats.safety_net_alerts += cp.safety_net_alerts;

        // Merge run distributions
        for (&len, &count) in &cp.run_distribution {
            *stats.run_distribution.entry(len).or_insert(0) += count;
        }

        if cp.longest_run > stats.longest_run {
            stats.longest_run = cp.longest_run;
            stats.longest_run_start = cp.longest_run_start;
        }
    }

    // Sort and dedupe
    stats.candidates.sort();
    stats.candidates.dedup();
    stats.witnesses.sort();
    stats.witnesses.dedup();
    stats.false_positives.sort();
    stats.false_positives.dedup();

    // Sort significant runs by position
    stats.significant_runs.sort_by_key(|r| r.start);

    stats
}

/// Aggregated statistics from multiple checkpoints.
#[derive(Debug, Default)]
pub struct AggregatedStats {
    pub total_checked: u64,
    pub total_governor_count: u64,
    pub longest_run: usize,
    pub longest_run_start: u64,
    pub candidates: Vec<u64>,
    pub witnesses: Vec<u64>,
    pub false_positives: Vec<u64>,
    pub false_positive_details: Vec<FalsePositiveInfo>,
    pub run_distribution: HashMap<usize, u64>,
    pub significant_runs: Vec<RunInfo>,
    pub counter_examples: Vec<CounterExampleInfo>,
    pub safety_net_windows_checked: u64,
    pub safety_net_alerts: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_checkpoint_save_load() {
        let mut cp = Checkpoint::new(9, 1_000_000, 2_000_000);
        cp.current_pos = 1_500_000;
        cp.checked = 500_000;
        cp.governor_count = 57_000;
        cp.longest_run = 7;
        cp.longest_run_start = 1_234_567;

        // Save to temp file
        let tmp = NamedTempFile::new().unwrap();
        cp.save(tmp.path()).unwrap();

        // Load and verify
        let loaded = Checkpoint::load(tmp.path()).unwrap();
        assert_eq!(loaded.target_k, 9);
        assert_eq!(loaded.current_pos, 1_500_000);
        assert_eq!(loaded.longest_run, 7);
    }

    #[test]
    fn test_progress() {
        let mut cp = Checkpoint::new(9, 0, 1000);
        assert_eq!(cp.progress_percent(), 0.0);

        cp.current_pos = 500;
        assert_eq!(cp.progress_percent(), 50.0);

        cp.current_pos = 1000;
        assert_eq!(cp.progress_percent(), 100.0);
    }
}
