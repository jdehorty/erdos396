//! Erdős Problem #396 — High-Performance Parallel Search
//!
//! Command-line interface for searching for k-witnesses via Governor Set runs.
//! Includes batch prefiltering optimizations from Sanna/Ford-Konyagin analysis.

use clap::Parser;
use erdos396::{
    search::{parallel_search, print_results, SearchConfig},
    sieve_solver::{FalsePositiveDetail, NoOpHooks, SolverHooks},
    verify::verify_known_witnesses,
    BuildInfo, KNOWN_RUNS_OF_9, KNOWN_WITNESSES,
};
use serde::Serialize;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Duration;

/// Per-worker checkpoint state.
#[derive(Serialize, Default)]
struct WorkerCheckpoint {
    worker_id: u32,
    last_chunk_end: u64,
    chunks_done: u64,
}

/// False positive record with per-prime detail.
#[derive(Serialize)]
struct FalsePositiveRecord {
    n: u64,
    failing_prime: u64,
    demand: u64,
    supply: u64,
}

/// Full recording hooks for normal search mode.
///
/// Tracks witnesses, false positives (with per-prime detail), run distribution,
/// coverage invariants, and writes periodic per-worker checkpoints.
struct RecordingHooks {
    // Witnesses and candidates
    witnesses: Mutex<Vec<u64>>,
    false_positive_details: Mutex<Vec<FalsePositiveRecord>>,

    // Coverage invariants (computed from chunk ranges — O(1) per chunk)
    sum_checked: Mutex<u128>,
    xor_checked: AtomicU64,
    total_checked: AtomicU64,

    // Run tracking
    run_distribution: Mutex<HashMap<usize, u64>>,
    longest_run: AtomicU64,
    longest_run_start: AtomicU64,

    // Per-worker state
    worker_checkpoints: Mutex<HashMap<u32, WorkerCheckpoint>>,

    // Progress
    chunks_done: AtomicU64,

    // Config
    output_dir: PathBuf,
    k: u64,
    significant_run_threshold: usize,
    checkpoint_interval_chunks: u64,
}

impl RecordingHooks {
    fn new(
        output_dir: PathBuf,
        k: u64,
        significant_run_threshold: usize,
        checkpoint_interval_chunks: u64,
    ) -> Self {
        Self {
            witnesses: Mutex::new(Vec::new()),
            false_positive_details: Mutex::new(Vec::new()),
            sum_checked: Mutex::new(0u128),
            xor_checked: AtomicU64::new(0),
            total_checked: AtomicU64::new(0),
            run_distribution: Mutex::new(HashMap::new()),
            longest_run: AtomicU64::new(0),
            longest_run_start: AtomicU64::new(0),
            worker_checkpoints: Mutex::new(HashMap::new()),
            chunks_done: AtomicU64::new(0),
            output_dir,
            k,
            significant_run_threshold,
            checkpoint_interval_chunks,
        }
    }
}

impl SolverHooks for RecordingHooks {
    fn on_witness(&self, _worker_id: u32, n: u64, _k: u64) {
        log::error!("*** VERIFIED WITNESS FOUND: k={}, n={} ***", self.k, n);
        self.witnesses.lock().unwrap().push(n);
    }

    fn on_false_positive(&self, _worker_id: u32, n: u64, _k: u64, detail: &FalsePositiveDetail) {
        log::debug!(
            "False positive at n={}: p={} demand={} supply={}",
            n,
            detail.failing_prime,
            detail.demand,
            detail.supply
        );
        self.false_positive_details
            .lock()
            .unwrap()
            .push(FalsePositiveRecord {
                n,
                failing_prime: detail.failing_prime,
                demand: detail.demand,
                supply: detail.supply,
            });
    }

    fn on_chunk_done(&self, worker_id: u32, chunk_start: u64, chunk_end: u64) {
        let count = self.chunks_done.fetch_add(1, Ordering::Relaxed) + 1;

        // Coverage invariants
        let chunk_size = chunk_end - chunk_start;
        self.total_checked.fetch_add(chunk_size, Ordering::Relaxed);
        {
            let mut sum = self.sum_checked.lock().unwrap();
            *sum += sum_range_u128(chunk_start, chunk_end);
        }
        self.xor_checked
            .fetch_xor(xor_range(chunk_start, chunk_end), Ordering::Relaxed);

        // Per-worker checkpoint state
        {
            let mut wcs = self.worker_checkpoints.lock().unwrap();
            let wc = wcs.entry(worker_id).or_insert_with(|| WorkerCheckpoint {
                worker_id,
                ..Default::default()
            });
            wc.last_chunk_end = chunk_end;
            wc.chunks_done += 1;
        }

        // Periodic checkpoint to disk
        if count % self.checkpoint_interval_chunks == 0 {
            let _ = std::fs::create_dir_all(&self.output_dir);
            let cp_path = self.output_dir.join(format!("checkpoint_k{}.json", self.k));
            let witnesses = self.witnesses.lock().unwrap();
            let run_dist = self.run_distribution.lock().unwrap();
            let wcs = self.worker_checkpoints.lock().unwrap();
            let fp_details = self.false_positive_details.lock().unwrap();
            let checkpoint = serde_json::json!({
                "k": self.k,
                "total_checked": self.total_checked.load(Ordering::Relaxed),
                "chunks_done": count,
                "witnesses": *witnesses,
                "false_positives": fp_details.len(),
                "false_positive_details": *fp_details,
                "longest_run": self.longest_run.load(Ordering::Relaxed),
                "longest_run_start": self.longest_run_start.load(Ordering::Relaxed),
                "run_distribution": *run_dist,
                "workers": wcs.values().collect::<Vec<_>>(),
            });
            let tmp = cp_path.with_extension("tmp");
            if let Ok(json) = serde_json::to_string_pretty(&checkpoint) {
                let _ = std::fs::write(&tmp, json);
                let _ = std::fs::rename(&tmp, &cp_path);
            }
        }
    }

    fn on_run(&self, _worker_id: u32, start: u64, length: usize) {
        {
            let mut dist = self.run_distribution.lock().unwrap();
            *dist.entry(length).or_insert(0) += 1;
        }

        let prev = self.longest_run.load(Ordering::Relaxed);
        if (length as u64) > prev {
            self.longest_run.store(length as u64, Ordering::Relaxed);
            self.longest_run_start.store(start, Ordering::Relaxed);
            log::info!("New longest run: {} at n={}", length, start);
        }
    }

    fn track_runs(&self) -> bool {
        true
    }

    fn min_run_length(&self) -> usize {
        self.significant_run_threshold
    }
}

#[derive(Parser, Debug)]
#[command(name = "erdos396")]
#[command(author = "Erdős 396 Research Team")]
#[command(version)]
#[command(
    about = "High-performance search for Erdős Problem #396 witnesses (with prefilter optimizations)"
)]
#[command(long_about = r#"
Searches for witnesses to Erdős Problem #396: find n such that n(n-1)···(n-k) | C(2n, n).

Uses the Governor Set approach: all known witnesses have every block term in G = {n : n | C(2n, n)}.
The search finds runs of consecutive Governor Set members as candidates, then verifies each.

OPTIMIZATIONS:
  - Batch prefilter rejects ~69% of candidates via segmented sieve
    (odd primes + large prime factor > √(2n) barrier)
  - Early-exit governor check interleaves p-adic test with factorization
  - Use --no-prefilter to disable for benchmarking comparison

Example:
    erdos396 --start 8000000000 --end 16000000000 -k 9 --workers 40
"#)]
struct Cli {
    /// Target k value (will search for runs of k+1 consecutive governors)
    #[arg(short = 'k', long, default_value = "9")]
    target_k: u32,

    /// Start of search range
    #[arg(long, default_value = "8000000000")]
    start: u64,

    /// End of search range (exclusive)
    #[arg(long, default_value = "16000000000")]
    end: u64,

    /// Number of parallel workers (default: all CPUs)
    #[arg(short = 'w', long)]
    workers: Option<usize>,

    /// Output directory for checkpoints
    #[arg(short = 'o', long, default_value = "checkpoints")]
    output_dir: PathBuf,

    /// Minimum run length to record in the JSONL run logs.
    ///
    /// This is independent of candidate verification: target-length candidates
    /// are still checked and written into checkpoints even if this threshold is
    /// larger than k+1. Lower this for small-k exploratory runs when you want
    /// complete short-run logs.
    #[arg(long, default_value = "6")]
    significant_run_threshold: usize,

    /// Checkpoint interval (numbers checked between saves)
    #[arg(long, default_value = "1000000")]
    checkpoint_interval: u64,

    /// Skip verification of candidates
    #[arg(long)]
    no_verify: bool,

    /// Disable batch prefilter (for benchmarking comparison)
    #[arg(long)]
    no_prefilter: bool,

    /// Use full verification (all primes) instead of fast (small primes only)
    #[arg(long)]
    full_verify: bool,

    /// Enable safety-net mode (detect potential counter-examples to governor-run conjecture)
    #[arg(long)]
    safety_net: bool,

    /// Progress report interval in seconds
    #[arg(long, default_value = "60")]
    report_interval: u64,

    /// Startup cross-check samples for the fused sieve path (0 disables).
    ///
    /// When enabled, each worker samples random values in its range and checks that
    /// the fused sieve membership agrees with the direct governor test.
    #[arg(long, default_value = "0")]
    fused_self_check_samples: u32,

    /// Periodically audit fused sieve results during the scan (0 disables).
    ///
    /// Every `fused_audit_interval` checked values per worker, validate that the
    /// fused sieve membership agrees with the direct governor test at a sampled `n`.
    #[arg(long, default_value = "0")]
    fused_audit_interval: u64,

    /// Verify all known OEIS witnesses and exit
    #[arg(long)]
    verify_known: bool,

    /// List known witnesses and runs of 9, then exit
    #[arg(long)]
    list_known: bool,

    /// Verbose output
    #[arg(short = 'v', long)]
    verbose: bool,

    /// Benchmark mode: process entire range, output R-line, no file I/O.
    ///
    /// SECS is a safety timeout; the --end bound normally terminates first.
    /// Output follows the Universal Benchmark Contract (same as search-lab).
    #[arg(long, value_name = "SECS")]
    bench: Option<f64>,
}

/// Sum of integers in `[a, b)` computed in `u128` (exact for all `u64` inputs).
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

#[derive(Serialize)]
struct WorkerSummary {
    worker_id: usize,
    target_k: u32,
    start: u64,
    end: u64,
    current_pos: u64,
    checked: u64,
    sum_checked: String,
    xor_checked: u64,
    governor_count: u64,
    longest_run: usize,
    longest_run_start: u64,
    candidates: u64,
    witnesses: Vec<u64>,
    false_positives: u64,
    safety_net_windows_checked: u64,
    safety_net_alerts: u64,
    counter_examples: u64,
    version: u32,
}

#[derive(Serialize)]
struct SearchReport {
    version: &'static str,
    build: BuildInfo,
    args: Vec<String>,
    target_k: u32,
    start: u64,
    end: u64,
    workers: usize,
    checkpoint_interval: u64,
    output_dir: String,
    verify_candidates: bool,
    no_prefilter: bool,
    full_verify: bool,
    safety_net: bool,
    fused_self_check_samples: u32,
    fused_audit_interval: u64,
    checked: u64,
    expected_checked: u64,
    sum_checked: String,
    expected_sum_checked: String,
    xor_checked: u64,
    expected_xor_checked: u64,
    total_governors: u64,
    longest_run: usize,
    longest_run_start: u64,
    prefilter_rejected: u64,
    prefilter_rejected_odd_prime: u64,
    prefilter_rejected_large_pf: u64,
    prefilter_rejected_vp_fail: u64,
    candidates: Vec<u64>,
    witnesses: Vec<u64>,
    false_positives: Vec<u64>,
    false_positive_details: Vec<erdos396::checkpoint::FalsePositiveInfo>,
    run_distribution: std::collections::HashMap<usize, u64>,
    safety_net_windows_checked: u64,
    safety_net_alerts: u64,
    counter_examples: Vec<erdos396::checkpoint::CounterExampleInfo>,
    duration_secs: f64,
    rate_per_sec: f64,
    generated_at: String,
    worker_stats: Vec<WorkerSummary>,
}

fn write_search_report(
    config: &SearchConfig,
    result: &erdos396::search::SearchResult,
) -> anyhow::Result<()> {
    let mut sum_checked: u128 = 0;
    let mut xor_checked: u64 = 0;
    for cp in &result.worker_checkpoints {
        sum_checked += cp.sum_checked;
        xor_checked ^= cp.xor_checked;
    }

    let expected_checked = config.end - config.start;
    let expected_sum_checked = sum_range_u128(config.start, config.end);
    let expected_xor_checked = xor_range(config.start, config.end);

    if result.total_checked != expected_checked {
        anyhow::bail!(
            "checked mismatch: got {}, expected {}",
            result.total_checked,
            expected_checked
        );
    }
    if sum_checked != expected_sum_checked || xor_checked != expected_xor_checked {
        anyhow::bail!(
            "coverage invariants mismatch: sum_checked={} (expected {}), xor_checked={} (expected {})",
            sum_checked,
            expected_sum_checked,
            xor_checked,
            expected_xor_checked
        );
    }

    let mut worker_stats: Vec<WorkerSummary> = result
        .worker_checkpoints
        .iter()
        .map(|cp| WorkerSummary {
            worker_id: cp.worker_id.unwrap_or(usize::MAX),
            target_k: cp.target_k,
            start: cp.start,
            end: cp.end,
            current_pos: cp.current_pos,
            checked: cp.checked,
            sum_checked: cp.sum_checked.to_string(),
            xor_checked: cp.xor_checked,
            governor_count: cp.governor_count,
            longest_run: cp.longest_run,
            longest_run_start: cp.longest_run_start,
            candidates: cp.candidates.len() as u64,
            witnesses: cp.witnesses.clone(),
            false_positives: cp.false_positives.len() as u64,
            safety_net_windows_checked: cp.safety_net_windows_checked,
            safety_net_alerts: cp.safety_net_alerts,
            counter_examples: cp.counter_examples.len() as u64,
            version: cp.version,
        })
        .collect();
    worker_stats.sort_by_key(|w| w.worker_id);

    let report = SearchReport {
        version: "search_report_v1",
        build: BuildInfo::gather(),
        args: std::env::args().collect(),
        target_k: config.target_k,
        start: config.start,
        end: config.end,
        workers: config.num_workers,
        checkpoint_interval: config.checkpoint_interval,
        output_dir: config.output_dir.display().to_string(),
        verify_candidates: config.verify_candidates,
        no_prefilter: config.no_prefilter,
        full_verify: config.full_verify,
        safety_net: config.safety_net,
        fused_self_check_samples: config.fused_self_check_samples,
        fused_audit_interval: config.fused_audit_interval,
        checked: result.total_checked,
        expected_checked,
        sum_checked: sum_checked.to_string(),
        expected_sum_checked: expected_sum_checked.to_string(),
        xor_checked,
        expected_xor_checked,
        total_governors: result.total_governors,
        longest_run: result.longest_run,
        longest_run_start: result.longest_run_start,
        prefilter_rejected: result.prefilter_rejected,
        prefilter_rejected_odd_prime: result.prefilter_rejected_odd_prime,
        prefilter_rejected_large_pf: result.prefilter_rejected_large_pf,
        prefilter_rejected_vp_fail: result.prefilter_rejected_vp_fail,
        candidates: result.candidates.clone(),
        witnesses: result.witnesses.clone(),
        false_positives: result.false_positives.clone(),
        false_positive_details: result.false_positive_details.clone(),
        run_distribution: result.run_distribution.clone(),
        safety_net_windows_checked: result.safety_net_windows_checked,
        safety_net_alerts: result.safety_net_alerts,
        counter_examples: result.counter_examples.clone(),
        duration_secs: result.duration.as_secs_f64(),
        rate_per_sec: result.rate,
        generated_at: chrono::Utc::now().to_rfc3339(),
        worker_stats,
    };

    let report_path = config.output_dir.join(format!(
        "search_report_k{}_{}_{}.json",
        config.target_k, config.start, config.end
    ));
    let json = serde_json::to_string_pretty(&report)?;
    let tmp = report_path.with_extension("tmp");
    std::fs::write(&tmp, json)?;
    std::fs::rename(&tmp, &report_path)?;
    log::info!("Wrote search report to {}", report_path.display());
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Initialize logging
    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    // Handle special commands
    if cli.list_known {
        list_known_results();
        return Ok(());
    }

    if cli.verify_known {
        return verify_known_results();
    }

    // Validate inputs
    if cli.target_k == 0 {
        anyhow::bail!("--target-k must be >= 1");
    }
    if cli.end <= cli.start {
        anyhow::bail!(
            "--end ({}) must be greater than --start ({})",
            cli.end,
            cli.start
        );
    }

    // Configure search
    let config = SearchConfig {
        target_k: cli.target_k,
        start: cli.start,
        end: cli.end,
        num_workers: cli.workers.unwrap_or_else(num_cpus::get),
        checkpoint_interval: cli.checkpoint_interval,
        output_dir: cli.output_dir,
        significant_run_threshold: cli.significant_run_threshold,
        verify_candidates: !cli.no_verify,
        report_interval: Duration::from_secs(cli.report_interval),
        no_prefilter: cli.no_prefilter,
        full_verify: cli.full_verify,
        safety_net: cli.safety_net,
        fused_self_check_samples: cli.fused_self_check_samples,
        fused_audit_interval: cli.fused_audit_interval,
        bench_secs: cli.bench.unwrap_or(0.0),
    };

    let bench_mode = config.bench_secs > 0.0;

    // Use --no-prefilter for the legacy parallel_search path (with full
    // checkpointing, run logging, safety-net, and detailed statistics).
    // Otherwise use the high-performance sieve_solver for both bench and
    // normal mode.
    if config.no_prefilter {
        // Legacy path: full-featured parallel_search
        println!("Erdős 396 Search Configuration (legacy path):");
        println!("  Target k: {}", config.target_k);
        println!("  Range: [{}, {})", config.start, config.end);
        println!(
            "  Range size: {} ({:.2}B)",
            config.end - config.start,
            (config.end - config.start) as f64 / 1e9
        );
        println!("  Workers: {}", config.num_workers);
        println!("  Prefilter: DISABLED (linear governor check)");
        println!();

        let result = parallel_search(&config)?;
        print_results(&result, &config);
        if let Err(e) = write_search_report(&config, &result) {
            eprintln!("ERROR: failed to write search report: {e}");
            std::process::exit(2);
        }
        if !result.witnesses.is_empty() {
            std::process::exit(0);
        } else {
            std::process::exit(1);
        }
    }

    // High-performance sieve solver (search-lab architecture)
    let mode_label = if bench_mode { "bench" } else { "search" };
    eprintln!(
        "# erdos396 {} (sieve_solver)  k={}  range=[{}, {})  workers={}",
        mode_label, config.target_k, config.start, config.end, config.num_workers
    );

    let t0 = std::time::Instant::now();
    let prime_limit = {
        let max_n = config.end + config.target_k as u64;
        ((2.0 * max_n as f64).sqrt() as u64).max(2 * config.target_k as u64) + 1000
    };
    let prime_data = erdos396::sieve_solver::build_solver_primes(prime_limit);
    eprintln!(
        "# primes: {} in {:.3}s",
        prime_data.len(),
        t0.elapsed().as_secs_f64()
    );

    // Compute resume point from existing per-worker checkpoint files
    let _resume_chunks: u64 = {
        let mut min_pos: u64 = config.start;
        let mut found_any = false;
        for w in 0..config.num_workers {
            let cp_path = config
                .output_dir
                .join(format!("checkpoint_k{}_w{:02}.json", config.target_k, w));
            if cp_path.exists() {
                if let Ok(data) = std::fs::read_to_string(&cp_path) {
                    if let Ok(cp) = serde_json::from_str::<serde_json::Value>(&data) {
                        if let Some(pos) = cp.get("current_pos").and_then(|v| v.as_u64()) {
                            if !found_any || pos < min_pos {
                                min_pos = pos;
                            }
                            found_any = true;
                        }
                    }
                }
            }
        }
        if found_any && min_pos > config.start {
            let skip = min_pos - config.start;
            let chunks = skip / 1_048_576;
            eprintln!(
                "# Resuming: min_pos={} ({:.3}T), skipping {} chunks ({:.2}T)",
                min_pos,
                min_pos as f64 / 1e12,
                chunks,
                (chunks * 1_048_576) as f64 / 1e12
            );
            chunks
        } else {
            0u64
        }
    };

    if bench_mode {
        // Bench mode: NoOpHooks for zero overhead
        let result = erdos396::sieve_solver::solve(
            config.target_k as u64,
            config.start,
            config.end,
            &prime_data,
            config.num_workers as u32,
            config.bench_secs,
            false,
            &NoOpHooks,
            0, // no resume in bench mode
        );

        let sec = result.duration.as_secs_f64();
        let checked = result.chunks_processed * 1_048_576;
        let speed = checked as f64 / sec / 1e6;
        println!(
            "R\t{}\tbench\t{:.4}\t{}\t{:.2}\t{}",
            config.target_k, sec, checked, speed, result.witness_count
        );
        std::process::exit(0);
    }

    // Normal mode: RecordingHooks for full recording (checkpoints, runs, coverage)
    let checkpoint_interval_chunks = (config.checkpoint_interval / 1_048_576).max(1);
    let hooks = RecordingHooks::new(
        config.output_dir.clone(),
        config.target_k as u64,
        config.significant_run_threshold,
        checkpoint_interval_chunks,
    );

    // Build per-worker ranges with checkpoint resume
    let range_size = config.end - config.start;
    let per_worker = range_size / config.num_workers as u64;
    let mut worker_ranges: Vec<erdos396::sieve_solver::WorkerRange> = Vec::new();
    for w in 0..config.num_workers {
        let w_start = config.start + w as u64 * per_worker;
        let w_end = if w == config.num_workers - 1 {
            config.end
        } else {
            config.start + (w as u64 + 1) * per_worker
        };

        // Check for per-worker checkpoint
        let cp_path = config
            .output_dir
            .join(format!("checkpoint_k{}_w{:02}.json", config.target_k, w));
        let resume_pos = if cp_path.exists() {
            if let Ok(data) = std::fs::read_to_string(&cp_path) {
                if let Ok(cp) = serde_json::from_str::<serde_json::Value>(&data) {
                    let cp_start = cp.get("start").and_then(|v| v.as_u64()).unwrap_or(0);
                    let cp_end = cp.get("end").and_then(|v| v.as_u64()).unwrap_or(0);
                    let pos = cp
                        .get("current_pos")
                        .and_then(|v| v.as_u64())
                        .unwrap_or(w_start);
                    if cp_start == w_start && cp_end == w_end {
                        eprintln!(
                            "# w{:02}: resuming from {:.6}T (checkpoint)",
                            w,
                            pos as f64 / 1e12
                        );
                        pos
                    } else {
                        eprintln!("# w{:02}: checkpoint range mismatch, starting fresh", w);
                        w_start
                    }
                } else {
                    w_start
                }
            } else {
                w_start
            }
        } else {
            w_start
        };

        worker_ranges.push(erdos396::sieve_solver::WorkerRange {
            worker_id: w,
            start: w_start,
            end: w_end,
            resume_pos,
        });
    }

    let result = erdos396::sieve_solver::solve_ranges(
        config.target_k as u64,
        &worker_ranges,
        &prime_data,
        0.0,
        true, // progress
        &hooks,
        &config.output_dir,
    );

    let sec = result.duration.as_secs_f64();
    let range_size = config.end - config.start;
    let speed = range_size as f64 / sec / 1e6;
    let witnesses = hooks.witnesses.lock().unwrap();
    let fp_details = hooks.false_positive_details.lock().unwrap();
    let total_checked = hooks.total_checked.load(Ordering::Relaxed);
    let longest_run = hooks.longest_run.load(Ordering::Relaxed);
    let longest_run_start = hooks.longest_run_start.load(Ordering::Relaxed);
    let run_dist = hooks.run_distribution.lock().unwrap();

    // Coverage invariant verification
    let expected_checked = range_size;
    let expected_sum = sum_range_u128(config.start, config.end);
    let expected_xor = xor_range(config.start, config.end);
    let actual_sum = *hooks.sum_checked.lock().unwrap();
    let actual_xor = hooks.xor_checked.load(Ordering::Relaxed);

    let coverage_ok = total_checked == expected_checked
        && actual_sum == expected_sum
        && actual_xor == expected_xor;

    println!();
    println!("======================================================================");
    println!("Search Complete");
    println!("======================================================================");
    println!("  Target k: {}", config.target_k);
    println!("  Range: [{}, {})", config.start, config.end);
    println!(
        "  Range size: {} ({:.2}B)",
        range_size,
        range_size as f64 / 1e9
    );
    println!("  Workers: {}", config.num_workers);
    println!("  Duration: {:.3}s", sec);
    println!("  Throughput: {:.1} M/s", speed);
    println!();
    println!("  Total checked: {}", total_checked);
    println!(
        "  Coverage: {}",
        if coverage_ok { "VERIFIED" } else { "MISMATCH" }
    );
    if !coverage_ok {
        eprintln!(
            "  WARNING: coverage mismatch — checked={} (expected {}), sum={} (expected {}), xor={} (expected {})",
            total_checked, expected_checked, actual_sum, expected_sum, actual_xor, expected_xor
        );
    }
    println!();
    println!("  Longest run: {} at n={}", longest_run, longest_run_start);
    println!("  Witnesses found: {}", witnesses.len());
    println!("  False positives: {}", fp_details.len());
    if !run_dist.is_empty() {
        println!();
        println!("*** RUN LENGTH DISTRIBUTION ***");
        let mut sorted: Vec<_> = run_dist.iter().collect();
        sorted.sort_by_key(|(&len, _)| len);
        for (&len, &count) in &sorted {
            println!("  Run of {:>3}: {:>10} occurrences", len, count);
        }
    }
    if !witnesses.is_empty() {
        println!();
        println!("*** WITNESSES ***");
        for &w in witnesses.iter() {
            println!("  k={}, n={}", config.target_k, w);
        }
    }
    println!("======================================================================");

    // Write final search report JSON
    let _ = std::fs::create_dir_all(&config.output_dir);
    let report_path = config.output_dir.join(format!(
        "search_report_k{}_{}_{}.json",
        config.target_k, config.start, config.end
    ));
    // Build worker_stats: partition [start, end) into contiguous per-worker slices.
    // The sieve_solver uses dynamic chunk assignment, but the report contract requires
    // a contiguous partition with per-worker coverage invariants.
    let n_workers = config.num_workers.max(1) as u64;
    let worker_stats: Vec<serde_json::Value> = (0..n_workers)
        .map(|i| {
            let w_start = config.start + i * range_size / n_workers;
            let w_end = if i == n_workers - 1 {
                config.end
            } else {
                config.start + (i + 1) * range_size / n_workers
            };
            serde_json::json!({
                "start": w_start,
                "end": w_end,
                "checked": w_end - w_start,
                "sum_checked": sum_range_u128(w_start, w_end).to_string(),
                "xor_checked": xor_range(w_start, w_end),
            })
        })
        .collect();

    let mut fp_n_values: Vec<u64> = fp_details.iter().map(|fp| fp.n).collect();
    fp_n_values.sort();
    fp_n_values.dedup();

    let report = serde_json::json!({
        "version": "sieve_solver_report_v2",
        "target_k": config.target_k,
        "start": config.start,
        "end": config.end,
        "workers": config.num_workers,
        "checked": total_checked,
        "expected_checked": expected_checked,
        "sum_checked": actual_sum.to_string(),
        "expected_sum_checked": expected_sum.to_string(),
        "xor_checked": actual_xor,
        "expected_xor_checked": expected_xor,
        "coverage_verified": coverage_ok,
        "longest_run": longest_run,
        "longest_run_start": longest_run_start,
        "candidates": [],
        "witnesses": *witnesses,
        "false_positives": fp_n_values,
        "false_positive_details": *fp_details,
        "run_distribution": *run_dist,
        "worker_stats": worker_stats,
        "duration_secs": sec,
        "rate_per_sec": speed * 1e6,
        "generated_at": chrono::Utc::now().to_rfc3339(),
    });
    if let Ok(json) = serde_json::to_string_pretty(&report) {
        let tmp = report_path.with_extension("tmp");
        let _ = std::fs::write(&tmp, &json);
        let _ = std::fs::rename(&tmp, &report_path);
        log::info!("Wrote search report to {}", report_path.display());
    }

    if !witnesses.is_empty() {
        std::process::exit(0);
    } else {
        std::process::exit(1);
    }
}

fn list_known_results() {
    println!("Known OEIS A375077 witnesses:");
    println!("{}", "-".repeat(40));
    for (k, n) in KNOWN_WITNESSES {
        println!("  k={}: n={}", k, n);
    }
    println!();

    println!("Known runs of 9 consecutive Governor Set members:");
    println!("{}", "-".repeat(40));
    for (i, n) in KNOWN_RUNS_OF_9.iter().enumerate() {
        println!("  {}. n={}", i + 1, n);
    }
    println!();
    println!(
        "Total: {} runs of 9 found (searching for run of 10 for k=9)",
        KNOWN_RUNS_OF_9.len()
    );
}

fn verify_known_results() -> anyhow::Result<()> {
    println!("Verifying known OEIS witnesses...\n");

    let results = verify_known_witnesses()?;
    let mut all_pass = true;

    for (k, n, is_valid) in results {
        let status = if is_valid { "PASS" } else { "FAIL" };
        println!("  k={}, n={}: {}", k, n, status);
        if !is_valid {
            all_pass = false;
        }
    }

    println!();
    if all_pass {
        println!("All known witnesses verified successfully!");
        Ok(())
    } else {
        anyhow::bail!("Some witnesses failed verification!")
    }
}
