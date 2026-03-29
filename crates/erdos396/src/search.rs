//! Parallel search for Governor Set runs.
//!
//! This module implements the core search algorithm that finds runs of
//! consecutive Governor Set members, which are candidates for Erdős 396 witnesses.
//!
//! ## Optimizations
//!
//! - **Fused batch computation**: A single segmented sieve pass divides out all
//!   small primes AND checks v_p(C(2n,n)) inline. This eliminates the redundant
//!   trial division that the two-step (prefilter → governor) approach required.
//!
//! - The sieve visits each (number, prime) pair only once, amortized across the
//!   batch. Trial division would redundantly test non-dividing primes against
//!   every number. For 10M numbers at n=10^10, this saves ~100M+ operations.

use crate::checkpoint::{Checkpoint, RunLogger};
use crate::governor::GovernorChecker;
use crate::int_math::isqrt_u128;
use crate::prefilter::FusedBatchResult;
use crate::sieve::PrimeSieve;
use crate::verify::WitnessVerifier;
use crate::Result;

use rayon::prelude::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Default batch size for fused sieve+governor processing.
/// 1M numbers: remaining array = 8MB, fits in L3 cache.
const DEFAULT_BATCH_SIZE: usize = 1_000_000;

/// Sum of integers in `[a, b)` computed in `u128`.
///
/// This is exact for all `u64` inputs because `Σ_{n < 2^64} n < 2^127`.
#[inline]
fn sum_range_u128(a: u64, b: u64) -> u128 {
    if b <= a {
        return 0;
    }
    let count = (b - a) as u128;
    let first_plus_last = a as u128 + (b - 1) as u128;
    if count % 2 == 0 {
        (count / 2) * first_plus_last
    } else {
        // If count is odd, `first_plus_last` must be even.
        count * (first_plus_last / 2)
    }
}

/// XOR of integers in `[0, n]` (inclusive).
#[inline]
fn xor_upto(n: u64) -> u64 {
    match n & 3 {
        0 => n,
        1 => 1,
        2 => n + 1,
        _ => 0,
    }
}

/// XOR of integers in `[a, b)`.
#[inline]
fn xor_range(a: u64, b: u64) -> u64 {
    if b <= a {
        return 0;
    }
    let hi = xor_upto(b - 1);
    let lo = if a == 0 { 0 } else { xor_upto(a - 1) };
    hi ^ lo
}

/// Configuration for a search operation.
#[derive(Debug, Clone)]
pub struct SearchConfig {
    /// Target k value (searching for run of k+1)
    pub target_k: u32,

    /// Start of search range
    pub start: u64,

    /// End of search range (exclusive)
    pub end: u64,

    /// Number of parallel workers
    pub num_workers: usize,

    /// Checkpoint interval (numbers checked between saves)
    pub checkpoint_interval: u64,

    /// Output directory for checkpoints
    pub output_dir: PathBuf,

    /// Minimum run length to record in the append-only JSONL run logs.
    ///
    /// Target-length candidates are still verified regardless of this threshold.
    pub significant_run_threshold: usize,

    /// Whether to verify candidates immediately
    pub verify_candidates: bool,

    /// Progress report interval (seconds)
    pub report_interval: Duration,

    /// Whether to disable the batch prefilter (for benchmarking)
    pub no_prefilter: bool,

    /// Whether to use full (all-prime) verification instead of fast (small-prime-only)
    pub full_verify: bool,

    /// Whether to enable safety-net mode (detect potential counter-examples)
    pub safety_net: bool,

    /// Startup cross-check samples for the fused sieve path (0 disables).
    ///
    /// When enabled, each worker samples random `n` values in its assigned range
    /// and checks that the fused sieve result matches the direct governor test.
    pub fused_self_check_samples: u32,

    /// Periodically audit fused sieve results during the scan (0 disables).
    ///
    /// Every `fused_audit_interval` checked values per worker, validate that the
    /// fused sieve membership agrees with the direct governor test at a sampled `n`.
    pub fused_audit_interval: u64,

    /// Benchmark mode duration in seconds (0.0 = disabled).
    ///
    /// When > 0, the search processes the entire [start, end) range without
    /// file I/O (no checkpoints, run logs, or search report). A background
    /// timer sets `should_stop` after this many seconds as a safety valve,
    /// but the `--end` bound normally terminates the search first.
    pub bench_secs: f64,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            target_k: 9,
            start: 8_000_000_000,
            end: 16_000_000_000,
            num_workers: num_cpus::get(),
            checkpoint_interval: 1_000_000,
            output_dir: PathBuf::from("checkpoints"),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter: false,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        }
    }
}

/// Result from a completed search.
#[derive(Debug, Clone)]
pub struct SearchResult {
    /// Total numbers checked across all workers
    pub total_checked: u64,

    /// Total Governor Set members found
    pub total_governors: u64,

    /// Longest run of consecutive governors found
    pub longest_run: usize,

    /// Starting position of the longest run
    pub longest_run_start: u64,

    /// All candidate positions (runs of target length)
    pub candidates: Vec<u64>,

    /// Verified witnesses
    pub witnesses: Vec<u64>,

    /// False positives (candidates that failed verification)
    pub false_positives: Vec<u64>,

    /// Detailed false positive information
    pub false_positive_details: Vec<crate::checkpoint::FalsePositiveInfo>,

    /// Run length distribution (length -> count)
    pub run_distribution: std::collections::HashMap<usize, u64>,

    /// All significant runs found
    pub significant_runs: Vec<crate::checkpoint::RunInfo>,

    /// Search duration (includes sieve building)
    pub duration: Duration,

    /// Search-only duration (excludes sieve building).
    /// Used for R-line output where prime-table generation is excluded.
    pub search_duration: Duration,

    /// Processing rate (numbers/second)
    pub rate: f64,

    /// Counter-examples found by safety-net mode
    pub counter_examples: Vec<crate::checkpoint::CounterExampleInfo>,

    /// Safety-net statistics
    pub safety_net_windows_checked: u64,
    pub safety_net_alerts: u64,

    /// Prefilter statistics
    pub prefilter_rejected: u64,
    pub prefilter_rejected_odd_prime: u64,
    pub prefilter_rejected_large_pf: u64,
    pub prefilter_rejected_vp_fail: u64,

    /// Final worker checkpoints (one per worker).
    ///
    /// These include coverage invariants (`sum_checked`, `xor_checked`) and are useful
    /// for producing reviewer-facing reports after long runs.
    pub worker_checkpoints: Vec<Checkpoint>,
}

/// Global progress tracking for parallel search.
struct GlobalProgress {
    best_run: AtomicUsize,
    best_run_pos: AtomicU64,
    total_checked: AtomicU64,
    total_governors: AtomicU64,
    total_rejected_odd_prime: AtomicU64,
    total_rejected_large_pf: AtomicU64,
    total_rejected_vp_fail: AtomicU64,
    should_stop: AtomicBool,
}

impl GlobalProgress {
    fn new() -> Self {
        Self {
            best_run: AtomicUsize::new(0),
            best_run_pos: AtomicU64::new(0),
            total_checked: AtomicU64::new(0),
            total_governors: AtomicU64::new(0),
            total_rejected_odd_prime: AtomicU64::new(0),
            total_rejected_large_pf: AtomicU64::new(0),
            total_rejected_vp_fail: AtomicU64::new(0),
            should_stop: AtomicBool::new(false),
        }
    }

    fn update_best_run(&self, run_len: usize, pos: u64) {
        let mut current = self.best_run.load(Ordering::Relaxed);
        while run_len > current {
            match self.best_run.compare_exchange_weak(
                current,
                run_len,
                Ordering::SeqCst,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    self.best_run_pos.store(pos, Ordering::Release);
                    log::info!("New longest run: {} at n={}", run_len, pos);
                    break;
                }
                Err(x) => current = x,
            }
        }
    }
}

/// A single search worker.
pub struct SearchWorker {
    /// Worker ID
    id: usize,

    /// Governor Set checker (for --no-prefilter path)
    checker: GovernorChecker,

    /// Witness verifier (for candidate verification)
    verifier: WitnessVerifier,

    /// Target run length (k+1)
    target_run_length: usize,

    /// Whether to verify candidates
    verify_candidates: bool,

    /// Whether to use the fused sieve+governor computation
    use_fused: bool,

    /// Whether to use full (all-prime) verification instead of fast
    full_verify: bool,

    /// Whether safety-net mode is enabled
    safety_net: bool,

    /// Primes for the fused sieve (up to √(2·end))
    prefilter_primes: Arc<Vec<u64>>,

    /// Startup cross-check samples for the fused sieve path (0 disables).
    fused_self_check_samples: u32,

    /// Periodic audit interval for fused sieve cross-checks (0 disables).
    fused_audit_interval: u64,

    /// Whether bench mode is active (skip all file I/O).
    bench_mode: bool,
}

impl SearchWorker {
    /// Create a new search worker.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: usize,
        sieve: PrimeSieve,
        target_k: u32,
        verify_candidates: bool,
        use_fused: bool,
        full_verify: bool,
        safety_net: bool,
        prefilter_primes: Arc<Vec<u64>>,
        fused_self_check_samples: u32,
        fused_audit_interval: u64,
        bench_mode: bool,
    ) -> Self {
        Self {
            id,
            checker: GovernorChecker::with_sieve(sieve.clone()),
            verifier: WitnessVerifier::with_sieve(sieve),
            target_run_length: (target_k + 1) as usize,
            verify_candidates,
            use_fused,
            full_verify,
            safety_net,
            prefilter_primes,
            fused_self_check_samples,
            fused_audit_interval,
            bench_mode,
        }
    }

    #[inline]
    fn lcg_next(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state
    }

    fn run_fused_self_check(&self, start: u64, end: u64) -> Result<()> {
        let samples = self.fused_self_check_samples;
        if samples == 0 {
            return Ok(());
        }
        let scan_size = end.saturating_sub(start);
        if scan_size == 0 {
            return Ok(());
        }

        let mut state =
            0x9c9b_a3b2_4a3f_1d77u64 ^ (self.id as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        for _ in 0..samples {
            let n = start + (Self::lcg_next(&mut state) % scan_size);
            let fused = FusedBatchResult::compute(n, 1, &self.prefilter_primes).is_governor[0];
            let direct = self.checker.is_governor(n);
            if fused != direct {
                return Err(crate::Error::Audit(format!(
                    "Fused self-check mismatch (worker {}): n={}, fused={}, direct={}",
                    self.id, n, fused, direct
                )));
            }
        }

        Ok(())
    }

    /// Search a range for Governor Set runs.
    fn search_range(
        &self,
        checkpoint: &mut Checkpoint,
        progress: &GlobalProgress,
        checkpoint_path: &PathBuf,
        checkpoint_interval: u64,
        run_logger: &mut RunLogger,
    ) -> Result<()> {
        let start = checkpoint.start;
        let end = checkpoint.end;
        let search_start = checkpoint.current_pos;

        log::info!(
            "Worker {} starting search [{}, {}), resuming from {}, mode={}",
            self.id,
            start,
            end,
            search_start,
            if self.use_fused { "fused" } else { "linear" }
        );

        if self.use_fused && self.fused_self_check_samples > 0 {
            log::info!(
                "Worker {} running fused self-check ({} samples)...",
                self.id,
                self.fused_self_check_samples
            );
            if let Err(e) = self.run_fused_self_check(start, end) {
                progress.should_stop.store(true, Ordering::Relaxed);
                return Err(e);
            }
            log::info!("Worker {} fused self-check passed.", self.id);
        }

        let result = if self.use_fused {
            self.search_fused(
                search_start,
                checkpoint,
                progress,
                checkpoint_path,
                checkpoint_interval,
                run_logger,
            )
        } else {
            self.search_linear(
                search_start,
                checkpoint,
                progress,
                checkpoint_path,
                checkpoint_interval,
                run_logger,
            )
        };

        // If we weren't asked to stop early, assert full coverage of the worker range.
        // In bench mode, skip assertions (I/O is disabled and timer may fire early).
        if result.is_ok() && !self.bench_mode && !progress.should_stop.load(Ordering::Relaxed) {
            let expected_checked = end.saturating_sub(start);
            let expected_sum = sum_range_u128(start, end);
            let expected_xor = xor_range(start, end);
            if checkpoint.current_pos != end || checkpoint.checked != expected_checked {
                return Err(crate::Error::Audit(format!(
                    "Coverage mismatch (worker {}): checked={} (expected {}), current_pos={} (expected {})",
                    self.id, checkpoint.checked, expected_checked, checkpoint.current_pos, end
                )));
            }
            if checkpoint.sum_checked != expected_sum || checkpoint.xor_checked != expected_xor {
                return Err(crate::Error::Audit(format!(
                    "Coverage invariant mismatch (worker {}): sum_checked={} (expected {}), xor_checked={} (expected {})",
                    self.id, checkpoint.sum_checked, expected_sum, checkpoint.xor_checked, expected_xor
                )));
            }
        }

        result
    }

    /// Search using fused sieve+governor computation (fastest path).
    ///
    /// For each batch of DEFAULT_BATCH_SIZE numbers, the fused sieve computes
    /// exact governor membership in a single pass. The search loop just reads
    /// the precomputed results — no per-number factorization at all.
    fn search_fused(
        &self,
        search_start: u64,
        checkpoint: &mut Checkpoint,
        progress: &GlobalProgress,
        checkpoint_path: &PathBuf,
        checkpoint_interval: u64,
        run_logger: &mut RunLogger,
    ) -> Result<()> {
        let end = checkpoint.end;
        let mut current_run = checkpoint.current_run;
        let mut run_start = if current_run == 0 {
            search_start
        } else {
            checkpoint.current_run_start
        };
        let mut checkpoint_counter = 0u64;
        let mut pos = search_start;

        // Prefix tracking: count consecutive governors from the start of the range.
        // If we're resuming and have already seen a non-governor, prefix is complete.
        let mut prefix_complete =
            checkpoint.checked > 0 && (checkpoint.prefix_run as u64) < checkpoint.checked;

        // Safety-net: ring buffer tracking governor membership for last (k+1) positions
        let k = (self.target_run_length - 1) as u32;
        let window_size = self.target_run_length;
        let mut ring_buf: Vec<bool> = vec![false; window_size];
        let mut ring_idx: usize = 0;
        let mut ring_filled: usize = 0; // how many slots have been written

        // Fused audit scheduling (in terms of `checkpoint.checked`).
        let audit_interval = self.fused_audit_interval;
        let mut next_audit_at = if audit_interval > 0 {
            checkpoint
                .checked
                .div_euclid(audit_interval)
                .saturating_add(1)
                .saturating_mul(audit_interval)
        } else {
            0
        };

        while pos < end {
            if progress.should_stop.load(Ordering::Relaxed) {
                log::info!("Worker {} received stop signal", self.id);
                break;
            }

            // Determine batch boundaries
            let batch_lo = pos;
            let batch_hi = std::cmp::min(pos + DEFAULT_BATCH_SIZE as u64, end);
            let batch_len = (batch_hi - batch_lo) as usize;

            // Fused sieve+governor: compute exact governor membership for entire batch
            let batch = FusedBatchResult::compute(batch_lo, batch_len, &self.prefilter_primes);

            // Optional audit: cross-check the fused result against the direct governor test.
            if audit_interval > 0 && next_audit_at > checkpoint.checked {
                let checked_before = checkpoint.checked;
                let checked_after = checked_before.saturating_add(batch_len as u64);

                while next_audit_at > checked_before && next_audit_at <= checked_after {
                    let idx = (next_audit_at - checked_before - 1) as usize;
                    let n = batch_lo + idx as u64;
                    let fused = batch.is_governor[idx];
                    let direct = self.checker.is_governor(n);
                    if fused != direct {
                        progress.should_stop.store(true, Ordering::Relaxed);
                        return Err(crate::Error::Audit(format!(
                            "Fused audit mismatch (worker {}): n={}, fused={}, direct={}",
                            self.id, n, fused, direct
                        )));
                    }
                    next_audit_at = next_audit_at.saturating_add(audit_interval);
                }
            }

            // Track statistics
            progress
                .total_rejected_odd_prime
                .fetch_add(batch.rejected_odd_prime as u64, Ordering::Relaxed);
            progress
                .total_rejected_large_pf
                .fetch_add(batch.rejected_large_pf as u64, Ordering::Relaxed);
            progress
                .total_rejected_vp_fail
                .fetch_add(batch.rejected_vp_fail as u64, Ordering::Relaxed);

            // Process the precomputed results — just read the array
            for i in 0..batch_len {
                let n = batch_lo + i as u64;
                let is_gov = batch.is_governor[i];

                checkpoint.checked += 1;
                checkpoint.sum_checked += n as u128;
                checkpoint.xor_checked ^= n;
                checkpoint.current_pos = n + 1;
                progress.total_checked.fetch_add(1, Ordering::Relaxed);

                if is_gov {
                    checkpoint.governor_count += 1;
                    progress.total_governors.fetch_add(1, Ordering::Relaxed);

                    if current_run == 0 {
                        run_start = n;
                    }
                    current_run += 1;

                    // Track prefix: consecutive governors from the start of the range
                    if !prefix_complete {
                        checkpoint.prefix_run += 1;
                    }

                    if current_run > checkpoint.longest_run {
                        checkpoint.longest_run = current_run;
                        checkpoint.longest_run_start = run_start;
                        progress.update_best_run(current_run, run_start);
                    }

                    if current_run == self.target_run_length {
                        self.handle_candidate(n, checkpoint)?;
                    }
                } else {
                    prefix_complete = true;
                    if !self.bench_mode && current_run >= 2 {
                        if let Some(run_info) = checkpoint.record_run(run_start, current_run) {
                            run_logger.log_run(&run_info)?;
                        }
                    }
                    current_run = 0;
                }

                // Safety-net: sliding window check
                if self.safety_net {
                    ring_buf[ring_idx] = is_gov;
                    ring_idx = (ring_idx + 1) % window_size;
                    if ring_filled < window_size {
                        ring_filled += 1;
                    }

                    if ring_filled == window_size {
                        // Count non-governors in the window
                        let non_gov_count = ring_buf.iter().filter(|&&g| !g).count();

                        // Skip if it's a full governor run (already handled above)
                        // or if too many non-governors (heuristic: unlikely to pass)
                        if non_gov_count > 0 && non_gov_count <= 2 {
                            checkpoint.safety_net_windows_checked += 1;

                            // Check small primes for this window
                            if self.verifier.check_small_primes(k, n)?.is_none() {
                                // Passes small-prime test but NOT a full governor run!
                                checkpoint.safety_net_alerts += 1;
                                log::warn!(
                                    "SAFETY-NET ALERT: n={} passes small-prime test with {} non-governors in window",
                                    n, non_gov_count
                                );

                                // Find which positions are non-governors
                                let mut non_gov_positions = Vec::new();
                                for offset in 0..window_size {
                                    let buf_idx =
                                        (ring_idx + window_size - 1 - offset) % window_size;
                                    if !ring_buf[buf_idx] {
                                        non_gov_positions.push(offset as u64);
                                    }
                                }

                                // Run full verification to confirm
                                let full_result = self.verifier.verify(k, n)?;
                                let full_valid = full_result.is_valid;
                                let full_failing = full_result.failing_prime;

                                if full_valid {
                                    log::error!(
                                        "*** COUNTER-EXAMPLE CONFIRMED: k={}, n={} is a witness but NOT a governor run! ***",
                                        k, n
                                    );
                                } else {
                                    log::info!(
                                        "Safety-net: n={} fails full verify at p={:?} (false alarm)",
                                        n, full_failing
                                    );
                                }

                                checkpoint.counter_examples.push(
                                    crate::checkpoint::CounterExampleInfo {
                                        n,
                                        non_governor_count: non_gov_count,
                                        non_governor_positions: non_gov_positions,
                                        full_verify_result: Some(full_valid),
                                        full_verify_failing_prime: full_failing,
                                    },
                                );
                            }
                        }
                    }
                }

                if !self.bench_mode {
                    checkpoint_counter += 1;
                    if checkpoint_counter >= checkpoint_interval {
                        checkpoint.current_run = current_run;
                        checkpoint.current_run_start = run_start;
                        run_logger.flush()?;
                        checkpoint.touch();
                        checkpoint.save_atomic_slim(checkpoint_path)?;
                        checkpoint_counter = 0;
                    }
                }
            }

            pos = batch_hi;
        }

        // Always update in-memory state (needed for boundary stitching).
        if !self.bench_mode && current_run >= 2 {
            if let Some(run_info) = checkpoint.record_run(run_start, current_run) {
                run_logger.log_run(&run_info)?;
            }
        }
        checkpoint.current_run = current_run;
        checkpoint.current_run_start = run_start;

        if !self.bench_mode {
            run_logger.flush()?;
            checkpoint.touch();
            checkpoint.save_atomic_slim(checkpoint_path)?;
        }

        log::info!(
            "Worker {} completed: checked {}, governors {}, longest run {} at {}",
            self.id,
            checkpoint.checked,
            checkpoint.governor_count,
            checkpoint.longest_run,
            checkpoint.longest_run_start
        );

        Ok(())
    }

    /// Search linearly without any sieve optimization (for benchmarking baseline).
    fn search_linear(
        &self,
        search_start: u64,
        checkpoint: &mut Checkpoint,
        progress: &GlobalProgress,
        checkpoint_path: &PathBuf,
        checkpoint_interval: u64,
        run_logger: &mut RunLogger,
    ) -> Result<()> {
        let end = checkpoint.end;
        let mut current_run = checkpoint.current_run;
        let mut run_start = if current_run == 0 {
            search_start
        } else {
            checkpoint.current_run_start
        };
        let mut checkpoint_counter = 0u64;

        // Prefix tracking: count consecutive governors from the start of the range.
        let mut prefix_complete =
            checkpoint.checked > 0 && (checkpoint.prefix_run as u64) < checkpoint.checked;

        for n in search_start..end {
            if progress.should_stop.load(Ordering::Relaxed) {
                log::info!("Worker {} received stop signal", self.id);
                break;
            }

            let is_gov = self.checker.is_governor(n);

            checkpoint.checked += 1;
            checkpoint.sum_checked += n as u128;
            checkpoint.xor_checked ^= n;
            checkpoint.current_pos = n + 1;
            progress.total_checked.fetch_add(1, Ordering::Relaxed);

            if is_gov {
                checkpoint.governor_count += 1;
                progress.total_governors.fetch_add(1, Ordering::Relaxed);

                if current_run == 0 {
                    run_start = n;
                }
                current_run += 1;

                // Track prefix: consecutive governors from the start of the range
                if !prefix_complete {
                    checkpoint.prefix_run += 1;
                }

                if current_run > checkpoint.longest_run {
                    checkpoint.longest_run = current_run;
                    checkpoint.longest_run_start = run_start;
                    progress.update_best_run(current_run, run_start);
                }

                if current_run == self.target_run_length {
                    self.handle_candidate(n, checkpoint)?;
                }
            } else {
                prefix_complete = true;
                if !self.bench_mode && current_run >= 2 {
                    if let Some(run_info) = checkpoint.record_run(run_start, current_run) {
                        run_logger.log_run(&run_info)?;
                    }
                }
                current_run = 0;
            }

            if !self.bench_mode {
                checkpoint_counter += 1;
                if checkpoint_counter >= checkpoint_interval {
                    checkpoint.current_run = current_run;
                    checkpoint.current_run_start = run_start;
                    run_logger.flush()?;
                    checkpoint.touch();
                    checkpoint.save_atomic_slim(checkpoint_path)?;
                    checkpoint_counter = 0;
                }
            }
        }

        // Always update in-memory state (needed for boundary stitching).
        if !self.bench_mode && current_run >= 2 {
            if let Some(run_info) = checkpoint.record_run(run_start, current_run) {
                run_logger.log_run(&run_info)?;
            }
        }
        checkpoint.current_run = current_run;
        checkpoint.current_run_start = run_start;

        if !self.bench_mode {
            run_logger.flush()?;
            checkpoint.touch();
            checkpoint.save_atomic_slim(checkpoint_path)?;
        }

        log::info!(
            "Worker {} completed: checked {}, governors {}, longest run {} at {}",
            self.id,
            checkpoint.checked,
            checkpoint.governor_count,
            checkpoint.longest_run,
            checkpoint.longest_run_start
        );

        Ok(())
    }

    /// Handle a candidate run of the target length.
    fn handle_candidate(&self, candidate: u64, checkpoint: &mut Checkpoint) -> Result<()> {
        log::warn!(
            "Worker {} found candidate run of {} at n={}!",
            self.id,
            self.target_run_length,
            candidate
        );

        checkpoint.candidates.push(candidate);

        if self.verify_candidates {
            let k = (self.target_run_length - 1) as u32;
            let result = if self.full_verify {
                self.verifier.verify(k, candidate)?
            } else {
                self.verifier.verify_fast(k, candidate)?
            };

            if result.is_valid {
                log::error!("*** VERIFIED WITNESS FOUND: k={}, n={} ***", k, candidate);
                checkpoint.witnesses.push(candidate);
            } else if let Some(failing_p) = result.failing_prime {
                let demand = *result.demand.get(&failing_p).unwrap_or(&0);
                let supply = *result.supply.get(&failing_p).unwrap_or(&0);
                log::warn!(
                    "Candidate n={} is FALSE POSITIVE: p={} (demand={}, supply={})",
                    candidate,
                    failing_p,
                    demand,
                    supply
                );
                checkpoint.false_positives.push(candidate);
                checkpoint.record_false_positive_detail(
                    candidate,
                    self.target_run_length,
                    failing_p,
                    demand,
                    supply,
                );
            } else {
                log::warn!(
                    "Candidate n={} is a FALSE POSITIVE (unknown reason)",
                    candidate
                );
                checkpoint.false_positives.push(candidate);
            }
        }
        Ok(())
    }
}

/// Run a parallel search with the given configuration.
pub fn parallel_search(config: &SearchConfig) -> Result<SearchResult> {
    let start_time = Instant::now();

    if config.target_k == 0 {
        return Err(crate::Error::InvalidParameter(
            "target_k must be >= 1".to_string(),
        ));
    }
    if config.end <= config.start {
        return Err(crate::Error::InvalidParameter(format!(
            "end ({}) must be greater than start ({})",
            config.end, config.start
        )));
    }
    if config.num_workers == 0 {
        return Err(crate::Error::InvalidParameter(
            "num_workers must be >= 1".to_string(),
        ));
    }
    if config.checkpoint_interval == 0 {
        return Err(crate::Error::InvalidParameter(
            "checkpoint_interval must be >= 1".to_string(),
        ));
    }
    // `GovernorChecker` (Legendre supply) requires that `2n` fits in `u64`.
    let max_safe_end = u64::MAX / 2 + 1;
    if config.end > max_safe_end {
        return Err(crate::Error::InvalidParameter(format!(
            "end too large: end={} must be <= {} to ensure 2n fits in u64",
            config.end, max_safe_end
        )));
    }

    let bench_mode = config.bench_secs > 0.0;

    if !bench_mode {
        std::fs::create_dir_all(&config.output_dir)?;
    }

    // Build prime sieve for factorization (needed for --no-prefilter and verification)
    log::info!("Building prime sieve for range up to {}...", config.end);
    let sieve_start = Instant::now();
    let sieve = PrimeSieve::for_range(config.end);
    log::info!(
        "Sieve built: {} primes in {:?}",
        sieve.len(),
        sieve_start.elapsed()
    );

    // Build prefilter sieve (primes up to √(2·end) for the fused computation)
    let prefilter_limit = (isqrt_u128(2u128 * (config.end as u128)) as u64).saturating_add(1000);
    let prefilter_sieve = PrimeSieve::new(prefilter_limit);
    let prefilter_primes = Arc::new(prefilter_sieve.primes().to_vec());
    log::info!(
        "Prefilter sieve: {} primes up to {} (√(2·end) barrier)",
        prefilter_primes.len(),
        prefilter_limit
    );

    let search_start_time = Instant::now();

    // Calculate per-worker ranges
    let range_size = config.end - config.start;
    let chunk_size = range_size / config.num_workers as u64;

    let ranges: Vec<(usize, u64, u64)> = (0..config.num_workers)
        .map(|i| {
            let w_start = config.start + (i as u64) * chunk_size;
            let w_end = if i == config.num_workers - 1 {
                config.end
            } else {
                config.start + ((i + 1) as u64) * chunk_size
            };
            (i, w_start, w_end)
        })
        .collect();

    log::info!(
        "Starting {} workers for k={} search in range [{}, {}), mode={}",
        config.num_workers,
        config.target_k,
        config.start,
        config.end,
        if config.no_prefilter {
            "linear"
        } else {
            "fused"
        }
    );

    let progress = Arc::new(GlobalProgress::new());

    // In bench mode, spawn a safety timer that sets should_stop after bench_secs.
    // The --end bound normally terminates workers first.
    if bench_mode {
        let p = Arc::clone(&progress);
        let secs = config.bench_secs;
        std::thread::spawn(move || {
            std::thread::sleep(Duration::from_secs_f64(secs));
            p.should_stop.store(true, Ordering::Relaxed);
        });
    }

    let results: Vec<Result<Checkpoint>> = ranges
        .into_par_iter()
        .map(|(worker_id, w_start, w_end)| -> Result<Checkpoint> {
            // In bench mode: use a per-worker tempdir for RunLogger and checkpoints.
            // No checkpoint loading, no run log migration. The tempdir is cleaned
            // up automatically when _bench_tmpdir is dropped.
            let _bench_tmpdir;
            let (checkpoint_path, run_log_path) = if bench_mode {
                let td = tempfile::tempdir()
                    .map_err(crate::Error::Io)?;
                let cp = td.path().join("checkpoint.json");
                let rl = td.path().join("runs.jsonl");
                _bench_tmpdir = Some(td);
                (cp, rl)
            } else {
                _bench_tmpdir = None;
                (
                    config.output_dir.join(format!(
                        "checkpoint_k{}_w{:02}.json",
                        config.target_k, worker_id
                    )),
                    config.output_dir.join(format!(
                        "runs_k{}_w{:02}.jsonl",
                        config.target_k, worker_id
                    )),
                )
            };

            let mut checkpoint = if !bench_mode && checkpoint_path.exists() {
                match Checkpoint::load(&checkpoint_path) {
                    Ok(mut cp) => {
                        let mut ok = true;
                        if cp.target_k != config.target_k || cp.start != w_start || cp.end != w_end {
                            ok = false;
                            log::warn!(
                                "Worker {} checkpoint mismatch: expected k={}, range=[{}, {}), got k={}, range=[{}, {}); starting fresh",
                                worker_id,
                                config.target_k,
                                w_start,
                                w_end,
                                cp.target_k,
                                cp.start,
                                cp.end
                            );
                        }
                        if cp.current_pos < cp.start || cp.current_pos > cp.end {
                            ok = false;
                            log::warn!(
                                "Worker {} checkpoint current_pos {} out of range [{}, {}); starting fresh",
                                worker_id,
                                cp.current_pos,
                                cp.start,
                                cp.end
                            );
                        }
                        if cp.checked != cp.current_pos.saturating_sub(cp.start) {
                            ok = false;
                            log::warn!(
                                "Worker {} checkpoint checked {} inconsistent with current_pos/start ({} - {}); starting fresh",
                                worker_id,
                                cp.checked,
                                cp.current_pos,
                                cp.start
                            );
                        }
                        if cp.current_run > 0
                            && (cp.current_run_start < cp.start || cp.current_run_start >= cp.current_pos)
                        {
                            ok = false;
                            log::warn!(
                                "Worker {} checkpoint run state inconsistent (current_run={}, current_run_start={}, current_pos={}); starting fresh",
                                worker_id,
                                cp.current_run,
                                cp.current_run_start,
                                cp.current_pos
                            );
                        }

                        if ok {
                            let expected_sum = sum_range_u128(cp.start, cp.current_pos);
                            let expected_xor = xor_range(cp.start, cp.current_pos);
                            if cp.checked > 0
                                && cp.sum_checked == 0
                                && cp.xor_checked == 0
                                && (expected_sum != 0 || expected_xor != 0)
                            {
                                log::info!(
                                    "Worker {} initializing coverage invariants for legacy checkpoint",
                                    worker_id
                                );
                                cp.sum_checked = expected_sum;
                                cp.xor_checked = expected_xor;
                            } else if cp.sum_checked != expected_sum || cp.xor_checked != expected_xor {
                                ok = false;
                                log::warn!(
                                    "Worker {} checkpoint coverage invariants inconsistent; starting fresh",
                                    worker_id
                                );
                            }
                        }

                        if ok {
                            if cp.worker_id != Some(worker_id) {
                                log::warn!(
                                    "Worker {} checkpoint has worker_id={:?}; continuing with corrected value",
                                    worker_id,
                                    cp.worker_id
                                );
                                cp.worker_id = Some(worker_id);
                            }
                            if cp.significant_run_threshold != config.significant_run_threshold {
                                log::info!(
                                    "Worker {} updating significant_run_threshold from {} to {}",
                                    worker_id,
                                    cp.significant_run_threshold,
                                    config.significant_run_threshold
                                );
                                cp.significant_run_threshold = config.significant_run_threshold;
                            }
                            // v4 → v5 migration: recompute prefix_run for legacy checkpoints.
                            // A v4 checkpoint loads with prefix_run=0 via #[serde(default)],
                            // which is indistinguishable from "range starts with a non-governor".
                            // Recompute by scanning from the range start.
                            if cp.version < 5 && cp.checked > 0 {
                                let checker = GovernorChecker::with_sieve(sieve.clone());
                                let mut prefix = 0usize;
                                for n in cp.start..cp.current_pos {
                                    if checker.is_governor(n) {
                                        prefix += 1;
                                    } else {
                                        break;
                                    }
                                }
                                cp.prefix_run = prefix;
                                cp.version = 5;
                                log::info!(
                                    "Worker {} migrated to v5: recomputed prefix_run={}",
                                    worker_id,
                                    prefix
                                );
                            }
                            log::info!("Worker {} resuming from checkpoint", worker_id);
                            cp
                        } else {
                            let mut cp =
                                Checkpoint::new_worker(config.target_k, w_start, w_end, worker_id);
                            cp.significant_run_threshold = config.significant_run_threshold;
                            cp
                        }
                    }
                    Err(e) => {
                        log::warn!("Failed to load checkpoint: {}, starting fresh", e);
                        let mut cp =
                            Checkpoint::new_worker(config.target_k, w_start, w_end, worker_id);
                        cp.significant_run_threshold = config.significant_run_threshold;
                        cp
                    }
                }
            } else {
                let mut cp =
                    Checkpoint::new_worker(config.target_k, w_start, w_end, worker_id);
                cp.significant_run_threshold = config.significant_run_threshold;
                cp
            };

            // Create run logger (append-only JSONL)
            let mut run_logger = RunLogger::new(&run_log_path)?;

            // v2 → v3 migration: if checkpoint has existing significant_runs,
            // migrate them to the run log file
            if !checkpoint.significant_runs.is_empty() {
                let count = checkpoint.significant_runs.len();
                log::info!(
                    "Worker {} migrating {} significant runs from v2 checkpoint to run log",
                    worker_id,
                    count
                );
                let migrated = run_logger.migrate_from_v2(&checkpoint.significant_runs)?;
                log::info!(
                    "Worker {} migrated {} runs to {}",
                    worker_id,
                    migrated,
                    run_log_path.display()
                );
                checkpoint.significant_runs.clear();
            }

            let worker = SearchWorker::new(
                worker_id,
                sieve.clone(),
                config.target_k,
                config.verify_candidates,
                !config.no_prefilter,
                config.full_verify,
                config.safety_net,
                prefilter_primes.clone(),
                config.fused_self_check_samples,
                config.fused_audit_interval,
                bench_mode,
            );

            worker.search_range(
                &mut checkpoint,
                &progress,
                &checkpoint_path,
                config.checkpoint_interval,
                &mut run_logger,
            )?;

            Ok(checkpoint)
        })
        .collect();

    let results: Vec<Checkpoint> = results.into_iter().collect::<Result<Vec<_>>>()?;

    let duration = start_time.elapsed();
    let search_duration = search_start_time.elapsed();

    // Aggregate results
    let mut total_checked = 0u64;
    let mut total_governors = 0u64;
    let mut longest_run = 0usize;
    let mut longest_run_start = 0u64;
    let mut candidates = Vec::new();
    let mut witnesses = Vec::new();
    let mut false_positives = Vec::new();
    let mut false_positive_details = Vec::new();
    let mut run_distribution: std::collections::HashMap<usize, u64> =
        std::collections::HashMap::new();
    let mut significant_runs = Vec::new();
    let mut counter_examples = Vec::new();
    let mut safety_net_windows_checked = 0u64;
    let mut safety_net_alerts = 0u64;

    for cp in &results {
        total_checked += cp.checked;
        total_governors += cp.governor_count;
        candidates.extend(cp.candidates.iter().copied());
        witnesses.extend(cp.witnesses.iter().copied());
        false_positives.extend(cp.false_positives.iter().copied());
        false_positive_details.extend(cp.false_positive_details.iter().cloned());
        significant_runs.extend(cp.significant_runs.iter().cloned());
        counter_examples.extend(cp.counter_examples.iter().cloned());
        safety_net_windows_checked += cp.safety_net_windows_checked;
        safety_net_alerts += cp.safety_net_alerts;

        for (&len, &count) in &cp.run_distribution {
            *run_distribution.entry(len).or_insert(0) += count;
        }

        if cp.longest_run > longest_run {
            longest_run = cp.longest_run;
            longest_run_start = cp.longest_run_start;
        }
    }

    // Stitch runs that cross worker boundaries.
    //
    // Each worker's `current_run` is the suffix (consecutive governors at the tail),
    // and `prefix_run` is the prefix (consecutive governors from the start of the range).
    // Adjacent workers' suffix + prefix may form a longer run that neither worker saw.
    //
    // Fragment accounting: each suffix/prefix/full-range fragment was already counted
    // in `run_distribution` by the individual worker. When we stitch, we remove the
    // fragment counts and add one count for the combined length.
    let target_run_length = (config.target_k + 1) as usize;
    if results.len() > 1 {
        let verifier = WitnessVerifier::with_sieve(sieve);

        let mut stitch_len = results[0].current_run;
        let mut stitch_start = results[0].current_run_start;
        // Track fragment lengths so we can fix run_distribution after stitching.
        let mut fragment_lengths: Vec<usize> = if stitch_len > 0 {
            vec![stitch_len]
        } else {
            vec![]
        };

        for cp in &results[1..] {
            let prefix = cp.prefix_run;

            if stitch_len > 0 && prefix > 0 {
                let stitch_len_before = stitch_len;
                stitch_len += prefix;
                fragment_lengths.push(prefix);

                // Update longest run if the stitched run is longer
                if stitch_len > longest_run {
                    longest_run = stitch_len;
                    longest_run_start = stitch_start;
                }

                // Emit a boundary candidate only if the stitch *newly* crosses the
                // target threshold. If the suffix alone already reached target length,
                // that candidate was already handled by the owning worker.
                if stitch_len >= target_run_length && stitch_len_before < target_run_length {
                    // Match single-worker semantics: the candidate position is the n
                    // where the run first reaches target length.
                    let candidate_n = stitch_start + target_run_length as u64 - 1;
                    log::warn!(
                        "Boundary stitch found candidate run of {} at n={} (across workers)",
                        stitch_len,
                        candidate_n
                    );
                    candidates.push(candidate_n);

                    if config.verify_candidates {
                        let k = config.target_k;
                        let result = if config.full_verify {
                            verifier.verify(k, candidate_n)?
                        } else {
                            verifier.verify_fast(k, candidate_n)?
                        };

                        if result.is_valid {
                            log::error!(
                                "*** VERIFIED WITNESS FOUND (boundary stitch): k={}, n={} ***",
                                k,
                                candidate_n
                            );
                            witnesses.push(candidate_n);
                        } else if let Some(failing_p) = result.failing_prime {
                            let demand = *result.demand.get(&failing_p).unwrap_or(&0);
                            let supply = *result.supply.get(&failing_p).unwrap_or(&0);
                            log::warn!(
                                "Boundary candidate n={} is FALSE POSITIVE: p={} (demand={}, supply={})",
                                candidate_n, failing_p, demand, supply
                            );
                            false_positives.push(candidate_n);
                            false_positive_details.push(crate::checkpoint::FalsePositiveInfo {
                                position: candidate_n,
                                run_length: stitch_len,
                                failing_prime: failing_p,
                                demand,
                                supply,
                            });
                        } else {
                            log::warn!(
                                "Boundary candidate n={} is a FALSE POSITIVE (unknown reason)",
                                candidate_n
                            );
                            false_positives.push(candidate_n);
                        }
                    }
                }

                // If this worker is entirely governors, the run continues to the next worker
                if prefix as u64 == cp.checked {
                    continue;
                }

                // Stitch is complete — fix run_distribution by replacing fragments
                // with the stitched whole.
                if fragment_lengths.len() > 1 {
                    for &frag in &fragment_lengths {
                        if frag >= 2 {
                            if let Some(count) = run_distribution.get_mut(&frag) {
                                *count = count.saturating_sub(1);
                            }
                        }
                    }
                    if stitch_len >= 2 {
                        *run_distribution.entry(stitch_len).or_insert(0) += 1;
                    }
                }
                fragment_lengths.clear();
            } else {
                // Chain broken (previous suffix > 0 but this prefix == 0, or
                // previous suffix == 0). Flush any pending multi-fragment stitch
                // before resetting, so the distribution replaces the individual
                // fragment entries with the stitched whole.
                if fragment_lengths.len() > 1 {
                    for &frag in &fragment_lengths {
                        if frag >= 2 {
                            if let Some(count) = run_distribution.get_mut(&frag) {
                                *count = count.saturating_sub(1);
                            }
                        }
                    }
                    if stitch_len >= 2 {
                        *run_distribution.entry(stitch_len).or_insert(0) += 1;
                    }
                }
                fragment_lengths.clear();
            }

            // Reset to this worker's suffix
            stitch_len = cp.current_run;
            stitch_start = cp.current_run_start;
            if stitch_len > 0 {
                fragment_lengths.push(stitch_len);
            }
        }

        // Flush any pending stitch at the end (e.g. chain ending at the last worker).
        if fragment_lengths.len() > 1 {
            for &frag in &fragment_lengths {
                if frag >= 2 {
                    if let Some(count) = run_distribution.get_mut(&frag) {
                        *count = count.saturating_sub(1);
                    }
                }
            }
            if stitch_len >= 2 {
                *run_distribution.entry(stitch_len).or_insert(0) += 1;
            }
        }

        run_distribution.retain(|_, count| *count > 0);
    }

    candidates.sort();
    candidates.dedup();
    witnesses.sort();
    witnesses.dedup();
    false_positives.sort();
    false_positives.dedup();
    false_positive_details.sort_by_key(|fp| (fp.position, fp.run_length, fp.failing_prime));
    false_positive_details.dedup_by_key(|fp| {
        (
            fp.position,
            fp.run_length,
            fp.failing_prime,
            fp.demand,
            fp.supply,
        )
    });
    significant_runs.sort_by_key(|r| r.start);
    significant_runs.dedup_by_key(|r| r.start);
    counter_examples.sort_by_key(|ce| ce.n);
    counter_examples.dedup_by_key(|ce| ce.n);

    let rate = total_checked as f64 / duration.as_secs_f64();

    Ok(SearchResult {
        total_checked,
        total_governors,
        longest_run,
        longest_run_start,
        candidates,
        witnesses,
        false_positives,
        false_positive_details,
        run_distribution,
        significant_runs,
        counter_examples,
        safety_net_windows_checked,
        safety_net_alerts,
        duration,
        search_duration,
        rate,
        prefilter_rejected: progress.total_rejected_odd_prime.load(Ordering::Relaxed)
            + progress.total_rejected_large_pf.load(Ordering::Relaxed)
            + progress.total_rejected_vp_fail.load(Ordering::Relaxed),
        prefilter_rejected_odd_prime: progress.total_rejected_odd_prime.load(Ordering::Relaxed),
        prefilter_rejected_large_pf: progress.total_rejected_large_pf.load(Ordering::Relaxed),
        prefilter_rejected_vp_fail: progress.total_rejected_vp_fail.load(Ordering::Relaxed),
        worker_checkpoints: results,
    })
}

/// Pretty-print search results.
pub fn print_results(result: &SearchResult, config: &SearchConfig) {
    println!("\n{}", "=".repeat(70));
    println!("Search Complete!");
    println!("{}", "=".repeat(70));
    println!("  Target k: {}", config.target_k);
    println!("  Range: [{}, {})", config.start, config.end);
    println!("  Workers: {}", config.num_workers);
    println!(
        "  Mode: {}",
        if config.no_prefilter {
            "linear (no prefilter)"
        } else {
            "fused sieve+governor"
        }
    );
    println!(
        "  Verify: {}",
        if config.full_verify {
            "full (all primes)"
        } else {
            "fast (small primes only)"
        }
    );
    if config.safety_net {
        println!("  Safety-net: ENABLED");
    }
    println!();
    println!("  Total checked: {:>15}", result.total_checked);
    println!("  Governor members: {:>12}", result.total_governors);
    println!(
        "  Observed density: {:>11.2}%",
        (result.total_governors as f64 / result.total_checked as f64) * 100.0
    );
    println!();

    // Fused sieve statistics
    if result.prefilter_rejected > 0 {
        let rejection_rate = result.prefilter_rejected as f64 / result.total_checked as f64 * 100.0;
        println!("*** FUSED SIEVE STATISTICS ***");
        println!(
            "  Not governor: {:>15} ({:.1}%)",
            result.prefilter_rejected, rejection_rate
        );
        println!(
            "    Odd primes:          {:>10}",
            result.prefilter_rejected_odd_prime
        );
        println!(
            "    Large pf (>sqrt(2n)):{:>10}",
            result.prefilter_rejected_large_pf
        );
        println!(
            "    v_p check failed:    {:>10}",
            result.prefilter_rejected_vp_fail
        );
        println!();
    }

    println!("  Duration: {:?}", result.duration);
    println!("  Rate: {:.0} numbers/sec", result.rate);
    println!();
    println!(
        "  Longest run: {} at n={}",
        result.longest_run, result.longest_run_start
    );
    println!(
        "  Runs of {} found: {}",
        config.target_k + 1,
        result.candidates.len()
    );
    println!("  Verified witnesses: {}", result.witnesses.len());
    println!("  False positives: {}", result.false_positives.len());

    if !result.run_distribution.is_empty() {
        println!("\n*** RUN LENGTH DISTRIBUTION ***");
        let mut lengths: Vec<_> = result.run_distribution.keys().copied().collect();
        lengths.sort();
        for len in lengths {
            if let Some(&count) = result.run_distribution.get(&len) {
                println!("  Run of {:>2}: {:>10} occurrences", len, count);
            }
        }
    }

    if !result.significant_runs.is_empty() {
        println!("\n*** SIGNIFICANT RUNS (length >= 6) ***");
        for run in &result.significant_runs {
            let status = match run.is_witness {
                Some(true) => "WITNESS",
                Some(false) => "false positive",
                None => "not verified",
            };
            println!(
                "  n={} (length {}) [{}]",
                run.start + run.length as u64 - 1,
                run.length,
                status
            );
        }
    }

    if !result.witnesses.is_empty() {
        println!("\n*** WITNESSES ***");
        for w in &result.witnesses {
            println!("  k={}, n={}", config.target_k, w);
        }
    }

    if !result.false_positive_details.is_empty() {
        println!("\n*** FALSE POSITIVE DETAILS ***");
        for fp in &result.false_positive_details {
            println!(
                "  n={}: run of {}, fails at p={} (demand={}, supply={})",
                fp.position, fp.run_length, fp.failing_prime, fp.demand, fp.supply
            );
        }
    }

    if result.safety_net_windows_checked > 0 || !result.counter_examples.is_empty() {
        println!("\n*** SAFETY-NET RESULTS ***");
        println!("  Windows checked: {}", result.safety_net_windows_checked);
        println!("  Alerts: {}", result.safety_net_alerts);
        println!("  Counter-examples: {}", result.counter_examples.len());

        for ce in &result.counter_examples {
            let status = match ce.full_verify_result {
                Some(true) => "CONFIRMED COUNTER-EXAMPLE".to_string(),
                Some(false) => format!(
                    "false alarm (fails at p={:?})",
                    ce.full_verify_failing_prime
                ),
                None => "not fully verified".to_string(),
            };
            println!(
                "  n={}: {} non-governors at offsets {:?} [{}]",
                ce.n, ce.non_governor_count, ce.non_governor_positions, status
            );
        }
    }

    if !result.candidates.is_empty() {
        println!("\nCandidate positions:");
        for (i, c) in result.candidates.iter().enumerate() {
            let status = if result.witnesses.contains(c) {
                "VERIFIED"
            } else if result.false_positives.contains(c) {
                "false positive"
            } else {
                "unverified"
            };
            println!("  {}. n={} [{}]", i + 1, c, status);
        }
    }

    println!("{}", "=".repeat(70));
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_small_search_fused() {
        let dir = tempdir().unwrap();

        let config = SearchConfig {
            target_k: 1,
            start: 1,
            end: 100,
            num_workers: 2,
            checkpoint_interval: 50,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: false,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let result = parallel_search(&config).unwrap();

        assert!(
            result.witnesses.contains(&2),
            "Should find k=1 witness at n=2"
        );
        assert!(result.total_governors > 0, "Should find governor members");

        // v3: verify run log files were created
        let run_log_0 = dir.path().join("runs_k1_w00.jsonl");
        let run_log_1 = dir.path().join("runs_k1_w01.jsonl");
        assert!(run_log_0.exists(), "Run log for worker 0 should exist");
        assert!(run_log_1.exists(), "Run log for worker 1 should exist");
    }

    #[test]
    fn test_small_search_linear() {
        let dir = tempdir().unwrap();

        let config = SearchConfig {
            target_k: 1,
            start: 1,
            end: 100,
            num_workers: 2,
            checkpoint_interval: 50,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: true,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let result = parallel_search(&config).unwrap();

        assert!(
            result.witnesses.contains(&2),
            "Should find k=1 witness at n=2"
        );
    }

    #[test]
    fn test_fused_vs_linear_agree() {
        let dir1 = tempdir().unwrap();
        let dir2 = tempdir().unwrap();

        let config_fused = SearchConfig {
            target_k: 2,
            start: 2400,
            end: 2600,
            num_workers: 1,
            checkpoint_interval: 1000,
            output_dir: dir1.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: false,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let config_linear = SearchConfig {
            no_prefilter: true,
            output_dir: dir2.path().to_path_buf(),
            ..config_fused.clone()
        };

        let result_fused = parallel_search(&config_fused).unwrap();
        let result_linear = parallel_search(&config_linear).unwrap();

        assert_eq!(
            result_fused.total_governors, result_linear.total_governors,
            "Governor count should match: fused={}, linear={}",
            result_fused.total_governors, result_linear.total_governors
        );
        assert_eq!(
            result_fused.longest_run, result_linear.longest_run,
            "Longest run should match"
        );
        assert_eq!(
            result_fused.witnesses, result_linear.witnesses,
            "Witnesses should match"
        );
    }

    #[test]
    fn test_fast_vs_full_verify_agree() {
        let dir1 = tempdir().unwrap();
        let dir2 = tempdir().unwrap();

        let config_fast = SearchConfig {
            target_k: 2,
            start: 2400,
            end: 2600,
            num_workers: 1,
            checkpoint_interval: 1000,
            output_dir: dir1.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: false,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let config_full = SearchConfig {
            full_verify: true,
            output_dir: dir2.path().to_path_buf(),
            ..config_fast.clone()
        };

        let result_fast = parallel_search(&config_fast).unwrap();
        let result_full = parallel_search(&config_full).unwrap();

        assert_eq!(
            result_fast.witnesses, result_full.witnesses,
            "Fast and full verify should find same witnesses"
        );
        assert_eq!(
            result_fast.false_positives, result_full.false_positives,
            "Fast and full verify should find same false positives"
        );
    }

    #[test]
    fn test_safety_net_small_range() {
        let dir = tempdir().unwrap();

        let config = SearchConfig {
            target_k: 2,
            start: 1,
            end: 10000,
            num_workers: 1,
            checkpoint_interval: 5000,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: false,
            full_verify: false,
            safety_net: true,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let result = parallel_search(&config).unwrap();

        // Safety-net should run without errors
        // Any counter-examples found should have full_verify_result set
        for ce in &result.counter_examples {
            assert!(
                ce.full_verify_result.is_some(),
                "Counter-example at n={} should have full verify result",
                ce.n
            );
        }
    }

    #[test]
    fn test_v4_checkpoint_is_slim() {
        let dir = tempdir().unwrap();

        let config = SearchConfig {
            target_k: 2,
            start: 1,
            end: 10000,
            num_workers: 1,
            checkpoint_interval: 5000,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(1),
            no_prefilter: false,
            full_verify: false,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let _result = parallel_search(&config).unwrap();

        // Verify checkpoint file is slim (no significant_runs)
        let cp_path = dir.path().join("checkpoint_k2_w00.json");
        let cp = Checkpoint::load(&cp_path).unwrap();
        assert!(
            cp.significant_runs.is_empty(),
            "v4 checkpoint should have empty significant_runs"
        );
        assert_eq!(cp.version, 5, "Checkpoint should be version 5");

        // Coverage invariants should match the assigned worker range.
        assert_eq!(cp.checked, config.end - config.start);
        assert_eq!(cp.sum_checked, sum_range_u128(config.start, config.end));
        assert_eq!(cp.xor_checked, xor_range(config.start, config.end));

        // Verify checkpoint file is small (under 2KB)
        let file_size = std::fs::metadata(&cp_path).unwrap().len();
        assert!(
            file_size < 2048,
            "v4 checkpoint should be under 2KB, got {} bytes",
            file_size
        );

        // Verify run log exists and has content
        let run_log_path = dir.path().join("runs_k2_w00.jsonl");
        assert!(run_log_path.exists(), "Run log should exist");
    }

    #[test]
    fn test_resume_does_not_miss_run_crossing_checkpoint_boundary() {
        // Regression test: if we resume from the *middle* of a governor run,
        // we must not miss candidates that depend on prior run state.
        //
        // Minimal example: k=1 looks for a run of 2 governors. Since 1 and 2
        // are both governors, resuming at n=2 with `current_run=1` must still
        // detect the candidate at n=2.
        let dir = tempdir().unwrap();

        let end = 5u64;
        let sieve = PrimeSieve::for_range(end);
        let prefilter = PrimeSieve::new(100);
        let prefilter_primes = Arc::new(prefilter.primes().to_vec());

        let worker = SearchWorker::new(
            0,
            sieve,
            1,     // k
            true,  // verify candidates
            false, // use fused (irrelevant here, we call search_linear directly)
            true,  // full verify for determinism
            false,
            prefilter_primes,
            0,
            0,
            false, // bench_mode
        );

        let progress = GlobalProgress::new();
        let checkpoint_path = dir.path().join("checkpoint_k1_w00.json");
        let run_log_path = dir.path().join("runs_k1_w00.jsonl");
        let mut run_logger = RunLogger::new(&run_log_path).unwrap();

        // Simulate a saved checkpoint after processing n=1 (a governor).
        let mut checkpoint = Checkpoint::new_worker(1, 1, end, 0);
        checkpoint.current_pos = 2;
        checkpoint.current_run = 1;
        checkpoint.current_run_start = 1;
        checkpoint.checked = 1;
        checkpoint.sum_checked = sum_range_u128(checkpoint.start, checkpoint.current_pos);
        checkpoint.xor_checked = xor_range(checkpoint.start, checkpoint.current_pos);
        checkpoint.governor_count = 1;

        worker
            .search_linear(
                checkpoint.current_pos,
                &mut checkpoint,
                &progress,
                &checkpoint_path,
                10,
                &mut run_logger,
            )
            .unwrap();

        assert!(
            checkpoint.witnesses.contains(&2),
            "Resumed search should still find k=1 witness at n=2"
        );
    }

    /// Helper: run parallel_search with the given worker count on the k=8 witness range.
    fn search_k8_witness_range(num_workers: usize, no_prefilter: bool) -> SearchResult {
        let dir = tempdir().unwrap();
        let config = SearchConfig {
            target_k: 8,
            start: 339_949_244,
            end: 339_949_253,
            num_workers,
            checkpoint_interval: 1_000_000,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter,
            full_verify: true,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };
        parallel_search(&config).unwrap()
    }

    #[test]
    fn test_multiworker_does_not_miss_run_crossing_worker_boundary() {
        // Regression test for GitHub issue #2:
        // The k=8 witness at n=339,949,252 is a run of 9 consecutive governors
        // [339,949,244 .. 339,949,252]. With 2 workers and range [339949244, 339949253),
        // the boundary at 339,949,248 splits the run and neither worker sees target length.

        let result_1w = search_k8_witness_range(1, false);
        assert!(
            result_1w.witnesses.contains(&339_949_252),
            "Single worker should find k=8 witness at n=339,949,252, got witnesses: {:?}",
            result_1w.witnesses
        );

        let result_2w = search_k8_witness_range(2, false);
        assert!(
            result_2w.witnesses.contains(&339_949_252),
            "Two workers should also find k=8 witness at n=339,949,252 via boundary stitching, got witnesses: {:?}",
            result_2w.witnesses
        );
    }

    #[test]
    fn test_multiworker_boundary_stitch_linear() {
        // Same boundary-crossing test but using the linear (--no-prefilter) path.
        let result = search_k8_witness_range(2, true);
        assert!(
            result.witnesses.contains(&339_949_252),
            "Linear mode with 2 workers should find k=8 witness via boundary stitching, got witnesses: {:?}",
            result.witnesses
        );
    }

    #[test]
    fn test_multiworker_witness_parity_across_worker_counts() {
        // w=1, w=2, and w=3 must all produce the same witnesses set on the repro range.
        let w1 = search_k8_witness_range(1, false);
        let w2 = search_k8_witness_range(2, false);
        let w3 = search_k8_witness_range(3, false);

        assert_eq!(
            w1.witnesses, w2.witnesses,
            "w=1 and w=2 should find the same witnesses: {:?} vs {:?}",
            w1.witnesses, w2.witnesses
        );
        assert_eq!(
            w1.witnesses, w3.witnesses,
            "w=1 and w=3 should find the same witnesses: {:?} vs {:?}",
            w1.witnesses, w3.witnesses
        );
    }

    #[test]
    fn test_prefix_run_tracking() {
        // Verify that prefix_run is correctly tracked for each worker.
        // Governors in [1, 10): 1, 2, 4, 6, 8
        // With 2 workers: worker 0 = [1, 5), worker 1 = [5, 10)
        // Worker 0: prefix = 2 (governors 1,2; then 3 is not a governor)
        // Worker 1: prefix = 0 (5 is odd prime, not a governor)
        let dir = tempdir().unwrap();
        let config = SearchConfig {
            target_k: 1,
            start: 1,
            end: 10,
            num_workers: 2,
            checkpoint_interval: 1_000_000,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 2,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter: false,
            full_verify: true,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let result = parallel_search(&config).unwrap();

        assert_eq!(result.worker_checkpoints.len(), 2);
        let w0 = &result.worker_checkpoints[0];
        let w1 = &result.worker_checkpoints[1];

        assert_eq!(
            w0.prefix_run, 2,
            "Worker 0 should have prefix_run=2 (governors 1,2), got {}",
            w0.prefix_run
        );
        assert_eq!(
            w1.prefix_run, 0,
            "Worker 1 should have prefix_run=0 (5 is not a governor), got {}",
            w1.prefix_run
        );
    }

    #[test]
    fn test_resume_preserves_prefix_run() {
        // Verify that prefix_run is correctly preserved across checkpoint resume.
        // Run a search that saves a checkpoint mid-worker, then resume and confirm
        // prefix_run is still correct.
        let dir = tempdir().unwrap();

        let end = 5u64;
        let sieve = PrimeSieve::for_range(end);
        let prefilter = PrimeSieve::new(100);
        let prefilter_primes = Arc::new(prefilter.primes().to_vec());

        let worker = SearchWorker::new(
            0,
            sieve,
            1,     // k
            true,  // verify candidates
            false, // use fused
            true,  // full verify
            false,
            prefilter_primes,
            0,
            0,
            false, // bench_mode
        );

        let progress = GlobalProgress::new();
        let checkpoint_path = dir.path().join("checkpoint_k1_w00.json");
        let run_log_path = dir.path().join("runs_k1_w00.jsonl");
        let mut run_logger = RunLogger::new(&run_log_path).unwrap();

        // Simulate a saved checkpoint after processing n=1 (a governor, prefix still going).
        let mut checkpoint = Checkpoint::new_worker(1, 1, end, 0);
        checkpoint.current_pos = 2;
        checkpoint.current_run = 1;
        checkpoint.current_run_start = 1;
        checkpoint.prefix_run = 1; // one governor seen from start
        checkpoint.checked = 1;
        checkpoint.sum_checked = sum_range_u128(checkpoint.start, checkpoint.current_pos);
        checkpoint.xor_checked = xor_range(checkpoint.start, checkpoint.current_pos);
        checkpoint.governor_count = 1;

        worker
            .search_linear(
                checkpoint.current_pos,
                &mut checkpoint,
                &progress,
                &checkpoint_path,
                10,
                &mut run_logger,
            )
            .unwrap();

        // After resume: 1 and 2 are both governors, 3 is not.
        // prefix_run should be 2 (1 from before + 1 from resume).
        assert_eq!(
            checkpoint.prefix_run, 2,
            "prefix_run should be 2 after resume (1+2 are governors, 3 is not), got {}",
            checkpoint.prefix_run
        );
        assert!(
            checkpoint.witnesses.contains(&2),
            "Should still find k=1 witness at n=2"
        );
    }

    #[test]
    fn test_multiworker_3_worker_chain() {
        // 3-worker search on the k=8 witness range [339949244, 339949253).
        // Range size = 9, so: worker 0 = [339949244, 339949247) (3 numbers),
        // worker 1 = [339949247, 339949250) (3 numbers),
        // worker 2 = [339949250, 339949253) (3 numbers).
        // All 9 numbers are governors, so workers 0 and 1 are entirely governors.
        // The stitching must chain all three to find the run of 9.
        let result = search_k8_witness_range(3, false);
        assert!(
            result.witnesses.contains(&339_949_252),
            "3-worker search should find k=8 witness via chained boundary stitching, got witnesses: {:?}",
            result.witnesses
        );

        // Verify chain: workers 0 and 1 should be entirely governors (prefix == checked).
        assert_eq!(result.worker_checkpoints.len(), 3);
        let w0 = &result.worker_checkpoints[0];
        let w1 = &result.worker_checkpoints[1];
        assert_eq!(
            w0.prefix_run as u64, w0.checked,
            "Worker 0 should be entirely governors (prefix_run={}, checked={})",
            w0.prefix_run, w0.checked
        );
        assert_eq!(
            w1.prefix_run as u64, w1.checked,
            "Worker 1 should be entirely governors (prefix_run={}, checked={})",
            w1.prefix_run, w1.checked
        );
    }

    #[test]
    fn test_boundary_stitch_distribution_correctness() {
        // Regression test: when a chain of entirely-governor workers terminates
        // because the next worker has prefix=0 (the else branch in stitching),
        // the fragment entries must be replaced in run_distribution.
        //
        // Range [339949244, 339949256): 9 consecutive governors followed by numbers
        // that include non-governors.  With 4 workers of size 3:
        //   W0 [339949244, 339949247) - entirely governors
        //   W1 [339949247, 339949250) - entirely governors
        //   W2 [339949250, 339949253) - entirely governors
        //   W3 [339949253, 339949256) - starts with a non-governor (prefix=0)
        //
        // The stitching else branch handles W3, which must flush the pending
        // 3-fragment chain before resetting.  Before the fix, fragment entries
        // remained and the stitched length was never recorded.
        let dir_1w = tempdir().unwrap();
        let dir_4w = tempdir().unwrap();

        let base = SearchConfig {
            target_k: 8,
            start: 339_949_244,
            end: 339_949_256,
            num_workers: 1,
            checkpoint_interval: 1_000_000,
            output_dir: dir_1w.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter: false,
            full_verify: true,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let config_1w = base.clone();
        let config_4w = SearchConfig {
            num_workers: 4,
            output_dir: dir_4w.path().to_path_buf(),
            ..base
        };

        let result_1w = parallel_search(&config_1w).unwrap();
        let result_4w = parallel_search(&config_4w).unwrap();

        assert_eq!(
            result_1w.run_distribution, result_4w.run_distribution,
            "run_distribution must match between 1-worker and 4-worker searches.\n  \
             1-worker: {:?}\n  4-worker: {:?}",
            result_1w.run_distribution, result_4w.run_distribution
        );

        assert_eq!(
            result_1w.witnesses, result_4w.witnesses,
            "Witnesses must match between 1-worker and 4-worker searches"
        );
    }

    #[test]
    fn test_v4_checkpoint_prefix_migration() {
        // Simulate loading a v4 checkpoint (no prefix_run) and verify the v5 migration
        // recomputes prefix_run correctly.
        let dir = tempdir().unwrap();

        // First run: create a v4-style checkpoint by running a small search,
        // then manually downgrade the version and clear prefix_run.
        let config = SearchConfig {
            target_k: 1,
            start: 1,
            end: 10,
            num_workers: 1,
            checkpoint_interval: 3, // checkpoint after 3 numbers
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 2,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter: false,
            full_verify: true,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 0.0,
        };

        let _ = parallel_search(&config).unwrap();

        // Load the checkpoint and downgrade it to v4
        let cp_path = dir.path().join("checkpoint_k1_w00.json");
        let mut cp = Checkpoint::load(&cp_path).unwrap();
        let original_prefix = cp.prefix_run;
        cp.version = 4;
        cp.prefix_run = 0; // simulate v4 (no prefix_run field)
                           // Reset to a mid-scan position so the migration has something to recompute
        cp.current_pos = 5;
        cp.checked = 4;
        cp.sum_checked = sum_range_u128(1, 5);
        cp.xor_checked = xor_range(1, 5);
        cp.save(&cp_path).unwrap();

        // Run search again — it should resume from the v4 checkpoint and migrate prefix_run.
        let result = parallel_search(&config).unwrap();
        let migrated_cp = &result.worker_checkpoints[0];

        // Governors in [1, 5): 1, 2 are governors, 3 is not.
        // So prefix_run should be 2 after migration.
        assert_eq!(
            migrated_cp.prefix_run, original_prefix,
            "v4→v5 migration should recompute prefix_run={}, got {}",
            original_prefix, migrated_cp.prefix_run
        );
        assert_eq!(
            migrated_cp.version, 5,
            "Checkpoint should be upgraded to v5"
        );
    }

    #[test]
    fn test_bench_mode_processes_full_range() {
        let dir = tempdir().unwrap();

        let config = SearchConfig {
            target_k: 8,
            start: 339_949_244,
            end: 339_949_253,
            num_workers: 2,
            checkpoint_interval: 1_000_000,
            output_dir: dir.path().to_path_buf(),
            significant_run_threshold: 6,
            verify_candidates: true,
            report_interval: Duration::from_secs(60),
            no_prefilter: false,
            full_verify: true,
            safety_net: false,
            fused_self_check_samples: 0,
            fused_audit_interval: 0,
            bench_secs: 999.0,
        };

        let result = parallel_search(&config).unwrap();

        // Bench mode should still find the k=8 witness
        assert!(
            result.witnesses.contains(&339_949_252),
            "Bench mode should find k=8 witness at n=339,949,252, got: {:?}",
            result.witnesses
        );

        // Should process entire range
        assert_eq!(
            result.total_checked,
            config.end - config.start,
            "Bench mode should check entire range"
        );

        // No checkpoint files should be written in the real output dir
        let checkpoint_files: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().map_or(false, |ext| ext == "json"))
            .collect();
        assert!(
            checkpoint_files.is_empty(),
            "Bench mode should not write checkpoint files to output dir, found: {:?}",
            checkpoint_files
        );
    }
}
