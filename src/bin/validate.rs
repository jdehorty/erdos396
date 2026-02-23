//! # Small-Prime Sieve Validator for Erdős Problem 396
//!
//! **Theoretical basis:** Corollary 6 of the Small Prime Barrier Theorem.
//!
//! This binary implements a provably complete method for finding ALL k-witnesses
//! in a given range — including any "non-governor witnesses" that the standard
//! Governor Set search cannot detect.
//!
//! ## How it works
//!
//! For every integer n in [start, end):
//!   1. Check the witness condition at each prime p < 2k+1:
//!        Σ v_p(n-i) for i=0..k  ≤  v_p(C(2n, n))
//!   2. If ALL barrier primes pass → full p-adic verification (all primes)
//!   3. Report confirmed witnesses
//!
//! ## Why this is complete
//!
//! The Small Prime Barrier Theorem (Theorem 1) proves: if a block term n-j
//! fails the governor test at a prime q ≥ 2k+1, then n cannot be a k-witness.
//! Contrapositive: any k-witness satisfies the witness condition at ALL primes
//! q ≥ 2k+1. Therefore checking only the "barrier primes" (p < 2k+1) is
//! sufficient to identify every witness — the screen has zero false negatives.
//!
//! ## Purpose
//!
//! Run this as a second pass over ranges already covered by the Governor Set
//! search. If the validator finds no new witnesses below the known minimum,
//! this constitutes a Certificate of Minimality — independent of the Governor
//! Set Completeness Conjecture.

use clap::Parser;
use erdos396::governor::{
    vp_central_binom_kummer_fast, vp_central_binom_p2, vp_central_binom_p3,
    vp_central_binom_p5, vp_factorial,
};
use erdos396::sieve::PrimeSieve;
use erdos396::verify::{primes_up_to, WitnessVerifier};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

#[derive(Parser, Debug)]
#[command(name = "validate")]
#[command(about = "Small-Prime Sieve Validator — provably complete witness search (Corollary 6)")]
#[command(long_about = r#"
Validates that no Erdős 396 k-witnesses were missed by the Governor Set search.

Uses the Small Prime Barrier Theorem: checks the witness condition at primes
p < 2k+1 for EVERY integer in the range (not just Governor Set members).
Any witness — governor run or not — must pass this screen.

This is a second-pass tool. Run it over ranges already covered by the main
search (erdos396) to produce a Certificate of Minimality for the paper.

Example:
    validate -k 11 --start 0 --end 1070858041585 --workers 40
"#)]
struct Cli {
    /// Target k value
    #[arg(short = 'k', long)]
    target_k: u32,

    /// Start of search range
    #[arg(long)]
    start: u64,

    /// End of search range (exclusive)
    #[arg(long)]
    end: u64,

    /// Number of parallel workers (default: all CPUs)
    #[arg(short = 'w', long)]
    workers: Option<usize>,

    /// Output directory for checkpoints
    #[arg(short = 'o', long, default_value = "validate_checkpoints")]
    output_dir: PathBuf,

    /// Checkpoint interval (numbers between saves)
    #[arg(long, default_value = "10000000")]
    checkpoint_interval: u64,

    /// Progress report interval in seconds
    #[arg(long, default_value = "60")]
    report_interval: u64,

    /// Verbose output
    #[arg(short = 'v', long)]
    verbose: bool,
}

#[derive(Serialize, Deserialize)]
struct ValidateCheckpoint {
    worker_id: usize,
    target_k: u32,
    start: u64,
    end: u64,
    current_pos: u64,
    checked: u64,
    small_prime_candidates: u64,
    confirmed_witnesses: Vec<u64>,
    false_alarms: u64,
    timestamp: String,
}

/// Compute v_p(C(2n, n)) using the best available method for each prime.
#[inline]
fn vp_supply(n: u64, p: u64) -> u64 {
    match p {
        2 => vp_central_binom_p2(n),
        3 => vp_central_binom_p3(n),
        5 => vp_central_binom_p5(n),
        _ => vp_central_binom_kummer_fast(n, p),
    }
}

/// Check the witness condition at all barrier primes (p < 2k+1).
///
/// Returns true if n passes at ALL barrier primes — meaning n is a candidate
/// that requires full verification. Returns false if any barrier prime fails,
/// proving n is not a k-witness.
///
/// Precondition: n > k.
#[inline]
fn passes_small_prime_screen(n: u64, k: u32, barrier_primes: &[u64]) -> bool {
    let block_bottom = n - k as u64;
    for &p in barrier_primes {
        // demand = Σ v_p(n-i) for i=0..k = v_p(n!) - v_p((n-k-1)!)
        let demand = vp_factorial(n, p) - vp_factorial(block_bottom - 1, p);
        let supply = vp_supply(n, p);
        if demand > supply {
            return false;
        }
    }
    true
}

fn save_checkpoint(path: &PathBuf, cp: &ValidateCheckpoint) {
    if let Ok(json) = serde_json::to_string_pretty(cp) {
        let tmp = path.with_extension("tmp");
        let _ = std::fs::write(&tmp, &json);
        let _ = std::fs::rename(&tmp, path);
    }
}

struct WorkerResult {
    checked: u64,
    candidates: u64,
    witnesses: Vec<u64>,
    false_alarms: u64,
}

fn load_or_init(
    worker_id: usize,
    k: u32,
    w_start: u64,
    w_end: u64,
    path: &PathBuf,
) -> ValidateCheckpoint {
    if path.exists() {
        if let Ok(data) = std::fs::read_to_string(path) {
            if let Ok(cp) = serde_json::from_str::<ValidateCheckpoint>(&data) {
                if cp.target_k == k && cp.start == w_start && cp.end == w_end {
                    log::info!(
                        "Worker {} resuming from {} ({:.2}% of worker range)",
                        worker_id,
                        cp.current_pos,
                        (cp.current_pos - w_start) as f64 / (w_end - w_start) as f64 * 100.0
                    );
                    return cp;
                }
            }
        }
    }

    ValidateCheckpoint {
        worker_id,
        target_k: k,
        start: w_start,
        end: w_end,
        current_pos: w_start,
        checked: 0,
        small_prime_candidates: 0,
        confirmed_witnesses: Vec::new(),
        false_alarms: 0,
        timestamp: String::new(),
    }
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    if cli.end <= cli.start {
        anyhow::bail!(
            "--end ({}) must be greater than --start ({})",
            cli.end,
            cli.start
        );
    }

    let k = cli.target_k;
    let barrier_threshold = 2 * k;
    let barrier_primes = primes_up_to(barrier_threshold);
    let num_workers = cli.workers.unwrap_or_else(num_cpus::get);

    println!("═══════════════════════════════════════════════════════════════════");
    println!("  Small-Prime Sieve Validator — Erdős Problem 396");
    println!("  Theoretical basis: Corollary 6, Small Prime Barrier Theorem");
    println!("═══════════════════════════════════════════════════════════════════");
    println!("  Target k:          {}", k);
    println!("  Range:             [{}, {})", cli.start, cli.end);
    println!(
        "  Range size:        {} ({:.3}T)",
        cli.end - cli.start,
        (cli.end - cli.start) as f64 / 1e12
    );
    println!("  Workers:           {}", num_workers);
    println!(
        "  Barrier primes:    {:?}  (all p < 2k+1 = {})",
        barrier_primes,
        2 * k + 1
    );
    println!("  Barrier count:     {} primes", barrier_primes.len());
    println!("  Output dir:        {:?}", cli.output_dir);
    println!();
    println!("  GUARANTEE: This screen has ZERO false negatives. Every k-witness");
    println!("  in the range WILL be detected, whether or not it is a Governor");
    println!("  Set run. (Small Prime Barrier Theorem, Corollary 6)");
    println!("═══════════════════════════════════════════════════════════════════");
    println!();

    std::fs::create_dir_all(&cli.output_dir)?;

    log::info!("Building prime sieve for full verification of candidates...");
    let sieve_start = Instant::now();
    let sieve = PrimeSieve::for_range(cli.end);
    log::info!(
        "Sieve built: {} primes in {:?}",
        sieve.len(),
        sieve_start.elapsed()
    );

    let verifier = WitnessVerifier::with_sieve(sieve);

    let range_size = cli.end - cli.start;
    let chunk_size = range_size / num_workers as u64;

    let ranges: Vec<(usize, u64, u64)> = (0..num_workers)
        .map(|i| {
            let w_start = cli.start + (i as u64) * chunk_size;
            let w_end = if i == num_workers - 1 {
                cli.end
            } else {
                cli.start + ((i + 1) as u64) * chunk_size
            };
            (i, w_start, w_end)
        })
        .collect();

    let total_checked = Arc::new(AtomicU64::new(0));
    let total_candidates = Arc::new(AtomicU64::new(0));
    let total_witnesses = Arc::new(AtomicU64::new(0));
    let stop_reporter = Arc::new(AtomicBool::new(false));

    let reporter = {
        let checked = total_checked.clone();
        let candidates = total_candidates.clone();
        let witnesses = total_witnesses.clone();
        let stop = stop_reporter.clone();
        let interval = cli.report_interval;
        let start_time = Instant::now();

        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                // Sleep in 1-second increments so we can exit promptly
                for _ in 0..interval {
                    if stop.load(Ordering::Relaxed) {
                        return;
                    }
                    std::thread::sleep(Duration::from_secs(1));
                }
                if stop.load(Ordering::Relaxed) {
                    break;
                }

                let c = checked.load(Ordering::Relaxed);
                let cand = candidates.load(Ordering::Relaxed);
                let w = witnesses.load(Ordering::Relaxed);
                let elapsed = start_time.elapsed().as_secs_f64();
                let rate = if elapsed > 0.0 {
                    c as f64 / elapsed
                } else {
                    0.0
                };
                let pct = c as f64 / range_size as f64 * 100.0;
                let eta_hours = if rate > 0.0 {
                    (range_size - c) as f64 / rate / 3600.0
                } else {
                    f64::INFINITY
                };

                log::info!(
                    "Progress: {:.2}% | {:.0}/sec | candidates: {} | witnesses: {} | ETA: {:.1}h",
                    pct,
                    rate,
                    cand,
                    w,
                    eta_hours
                );
            }
        })
    };

    let search_start = Instant::now();

    let results: Vec<WorkerResult> = ranges
        .into_par_iter()
        .map(|(worker_id, w_start, w_end)| {
            let checkpoint_path = cli
                .output_dir
                .join(format!("validate_k{}_w{:02}.json", k, worker_id));

            let mut state = load_or_init(worker_id, k, w_start, w_end, &checkpoint_path);
            let scan_start = std::cmp::max(state.current_pos, k as u64 + 1);
            let mut checkpoint_counter = 0u64;

            for n in scan_start..w_end {
                state.checked += 1;
                state.current_pos = n + 1;
                total_checked.fetch_add(1, Ordering::Relaxed);

                if passes_small_prime_screen(n, k, &barrier_primes) {
                    state.small_prime_candidates += 1;
                    total_candidates.fetch_add(1, Ordering::Relaxed);

                    let result = verifier.verify(k, n);
                    if result.is_valid {
                        log::error!(
                            "*** CONFIRMED WITNESS: k={}, n={} (worker {}) ***",
                            k,
                            n,
                            worker_id
                        );
                        state.confirmed_witnesses.push(n);
                        total_witnesses.fetch_add(1, Ordering::Relaxed);
                    } else {
                        state.false_alarms += 1;
                    }
                }

                checkpoint_counter += 1;
                if checkpoint_counter >= cli.checkpoint_interval {
                    state.timestamp = chrono::Utc::now().to_rfc3339();
                    save_checkpoint(&checkpoint_path, &state);
                    checkpoint_counter = 0;
                }
            }

            state.timestamp = chrono::Utc::now().to_rfc3339();
            save_checkpoint(&checkpoint_path, &state);

            log::info!(
                "Worker {} done: checked={}, screen_pass={}, witnesses={}, false_alarms={}",
                worker_id,
                state.checked,
                state.small_prime_candidates,
                state.confirmed_witnesses.len(),
                state.false_alarms
            );

            WorkerResult {
                checked: state.checked,
                candidates: state.small_prime_candidates,
                witnesses: state.confirmed_witnesses,
                false_alarms: state.false_alarms,
            }
        })
        .collect();

    let duration = search_start.elapsed();

    stop_reporter.store(true, Ordering::Relaxed);
    let _ = reporter.join();

    let agg_checked: u64 = results.iter().map(|r| r.checked).sum();
    let agg_candidates: u64 = results.iter().map(|r| r.candidates).sum();
    let agg_false_alarms: u64 = results.iter().map(|r| r.false_alarms).sum();
    let mut agg_witnesses: Vec<u64> = results
        .iter()
        .flat_map(|r| r.witnesses.iter().copied())
        .collect();
    agg_witnesses.sort();

    let rate = agg_checked as f64 / duration.as_secs_f64();

    println!();
    println!("═══════════════════════════════════════════════════════════════════");
    println!("  VALIDATION COMPLETE");
    println!("═══════════════════════════════════════════════════════════════════");
    println!("  Target k:          {}", k);
    println!("  Range:             [{}, {})", cli.start, cli.end);
    println!("  Total checked:     {}", agg_checked);
    println!("  Duration:          {:?}", duration);
    println!(
        "  Rate:              {:.0} numbers/sec ({:.0}/sec/core)",
        rate,
        rate / num_workers as f64
    );
    println!();
    println!("  Barrier primes:    {:?}", barrier_primes);
    println!(
        "  Screen pass:       {} ({:.6}%)",
        agg_candidates,
        if agg_checked > 0 {
            agg_candidates as f64 / agg_checked as f64 * 100.0
        } else {
            0.0
        }
    );
    println!("  False alarms:      {}", agg_false_alarms);
    println!("  WITNESSES:         {}", agg_witnesses.len());
    println!();

    if agg_witnesses.is_empty() {
        println!("  RESULT: NO WITNESSES FOUND IN RANGE");
        println!();
        println!("  The range [{}, {}) is CLEAR.", cli.start, cli.end);
        println!(
            "  No k={} witnesses exist in this range — proven by the",
            k
        );
        println!("  Small Prime Barrier Theorem (Corollary 6), not conjectured.");
        println!();
        println!("  Combined with the Governor Set search results, this");
        println!("  constitutes a CERTIFICATE OF MINIMALITY for the known");
        println!("  minimum witness.");
    } else {
        println!("  CONFIRMED WITNESSES:");
        for w in &agg_witnesses {
            println!("    k={}, n={}", k, w);
        }
    }

    println!("═══════════════════════════════════════════════════════════════════");

    Ok(())
}
