//! Erdős Problem #396 — High-Performance Parallel Search
//!
//! Command-line interface for searching for k-witnesses via Governor Set runs.
//! Includes batch prefiltering optimizations from Sanna/Ford-Konyagin analysis.

use clap::Parser;
use erdos396::{
    search::{parallel_search, print_results, SearchConfig},
    verify::verify_known_witnesses,
    BuildInfo, KNOWN_RUNS_OF_9, KNOWN_WITNESSES,
};
use serde::Serialize;
use std::path::PathBuf;
use std::time::Duration;

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
        verify_candidates: !cli.no_verify,
        report_interval: Duration::from_secs(cli.report_interval),
        no_prefilter: cli.no_prefilter,
        full_verify: cli.full_verify,
        safety_net: cli.safety_net,
        fused_self_check_samples: cli.fused_self_check_samples,
        fused_audit_interval: cli.fused_audit_interval,
    };

    // Print configuration
    println!("Erdős 396 Search Configuration:");
    println!("  Target k: {}", config.target_k);
    println!("  Range: [{}, {})", config.start, config.end);
    println!(
        "  Range size: {} ({:.2}B)",
        config.end - config.start,
        (config.end - config.start) as f64 / 1e9
    );
    println!("  Workers: {}", config.num_workers);
    println!(
        "  Prefilter: {}",
        if config.no_prefilter {
            "DISABLED"
        } else {
            "ENABLED (sqrt(2n) barrier + odd prime exclusion)"
        }
    );
    if !config.no_prefilter {
        if config.fused_self_check_samples > 0 {
            println!(
                "  Fused self-check: {} samples/worker",
                config.fused_self_check_samples
            );
        }
        if config.fused_audit_interval > 0 {
            println!(
                "  Fused audit interval: {} checked values/worker",
                config.fused_audit_interval
            );
        }
    }
    println!(
        "  Verify mode: {}",
        if config.full_verify {
            "full (all primes)"
        } else {
            "fast (small primes only)"
        }
    );
    println!(
        "  Safety-net: {}",
        if config.safety_net {
            "ENABLED"
        } else {
            "disabled"
        }
    );
    println!("  Checkpoint interval: {}", config.checkpoint_interval);
    println!("  Output directory: {:?}", config.output_dir);
    println!("  Verify candidates: {}", config.verify_candidates);
    println!();

    // Run search
    let result = parallel_search(&config)?;

    // Print results
    print_results(&result, &config);
    if let Err(e) = write_search_report(&config, &result) {
        eprintln!("ERROR: failed to write search report: {e}");
        std::process::exit(2);
    }

    // Exit with appropriate code
    if !result.witnesses.is_empty() {
        std::process::exit(0); // Found witness!
    } else {
        std::process::exit(1); // No witness found
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
