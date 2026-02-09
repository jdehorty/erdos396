//! Erdős Problem #396 — High-Performance Parallel Search
//!
//! Command-line interface for searching for k-witnesses via Governor Set runs.
//! Includes batch prefiltering optimizations from Sanna/Ford-Konyagin analysis.

use clap::Parser;
use erdos396::{
    search::{parallel_search, print_results, SearchConfig},
    verify::verify_known_witnesses,
    KNOWN_RUNS_OF_9, KNOWN_WITNESSES,
};
use std::path::PathBuf;
use std::time::Duration;

#[derive(Parser, Debug)]
#[command(name = "erdos396")]
#[command(author = "Erdős 396 Research Team")]
#[command(version)]
#[command(about = "High-performance search for Erdős Problem #396 witnesses (with prefilter optimizations)")]
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
    if cli.end <= cli.start {
        anyhow::bail!("--end ({}) must be greater than --start ({})", cli.end, cli.start);
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
    println!("  Prefilter: {}", if config.no_prefilter { "DISABLED" } else { "ENABLED (sqrt(2n) barrier + odd prime exclusion)" });
    println!("  Verify mode: {}", if config.full_verify { "full (all primes)" } else { "fast (small primes only)" });
    println!("  Safety-net: {}", if config.safety_net { "ENABLED" } else { "disabled" });
    println!("  Checkpoint interval: {}", config.checkpoint_interval);
    println!("  Output directory: {:?}", config.output_dir);
    println!("  Verify candidates: {}", config.verify_candidates);
    println!();

    // Run search
    let result = parallel_search(&config)?;

    // Print results
    print_results(&result, &config);

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

    let results = verify_known_witnesses();
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
