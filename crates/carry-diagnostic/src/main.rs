mod analysis;
#[allow(dead_code)] // fused sieve functions used in tests and benchmarking
mod fast_sieve;
mod scanner;

use analysis::{analyze_near_miss, barrier_primes, NearMissReport};
use clap::Parser;
use erdos396::PrimeSieve;
use scanner::scan_near_misses;

#[derive(Parser)]
#[command(name = "carry-diagnostic")]
#[command(about = "Near-miss scanner for carry structure compensation analysis")]
struct Cli {
    /// Start of scan range
    #[arg(long)]
    start: u64,

    /// Number of integers to scan
    #[arg(long)]
    count: u64,

    /// Block width parameter (default 14)
    #[arg(long, default_value_t = 14)]
    k: u32,

    /// Number of threads (default: all cores)
    #[arg(long, short = 't')]
    threads: Option<usize>,

    /// Output as JSON instead of table
    #[arg(long)]
    json: bool,

    /// Print only the aggregate summary
    #[arg(long)]
    summary_only: bool,
}

fn main() {
    let cli = Cli::parse();

    if cli.k < 2 {
        eprintln!("Error: k must be >= 2");
        std::process::exit(1);
    }
    if cli.start < cli.k as u64 {
        eprintln!("Error: start must be >= k (to avoid underflow)");
        std::process::exit(1);
    }
    if cli.start.checked_add(cli.count).is_none() {
        eprintln!("Error: start + count overflows u64");
        std::process::exit(1);
    }

    // Configure thread pool
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }
    let num_threads = rayon::current_num_threads();

    let max_n = cli.start + cli.count + 100;
    eprintln!("Initializing sieve for range up to {}...", max_n);
    let sieve = PrimeSieve::for_range(max_n);

    eprintln!(
        "Scanning {} integers starting at {} (k={}, {} threads)...",
        cli.count, cli.start, cli.k, num_threads
    );

    let t_start = std::time::Instant::now();
    let near_misses = scan_near_misses(
        cli.start,
        cli.count,
        cli.k,
        sieve.primes(),
        &sieve,
        !cli.json,
    );
    let elapsed = t_start.elapsed().as_secs_f64();
    let rate = cli.count as f64 / elapsed / 1e6;
    eprintln!(
        "Scan complete in {:.1}s ({:.1}M/s, {} threads)",
        elapsed, rate, num_threads
    );

    let reports: Vec<NearMissReport> = near_misses
        .iter()
        .map(|nm| analyze_near_miss(nm, cli.k, &sieve))
        .collect();

    if cli.json {
        print_json(&reports);
    } else {
        if !cli.summary_only {
            for report in &reports {
                print_report(report);
            }
        }
        print_summary(&reports, cli.start, cli.count, cli.k);
    }
}

fn print_report(report: &NearMissReport) {
    let nm = &report.near_miss;
    let failing: Vec<u64> = report
        .barrier_primes
        .iter()
        .filter(|a| a.is_failing)
        .map(|a| a.p)
        .collect();
    let purity = if report.pure_barrier_only {
        " [PURE]"
    } else {
        ""
    };
    println!(
        "Near-miss at n={}, m={} (j={}), fails at {:?}{}",
        nm.n, nm.m, nm.j, failing, purity
    );

    println!("  Marginal (m vs n):");
    for a in &report.barrier_primes {
        let tag = if a.is_failing { " [FAILING]" } else { "" };
        println!(
            "    p={:<3} demand={:<3} supply={:<3} deficit={:<4} bonus={:<4} gap={:<4} cascade={}{}",
            a.p, a.demand, a.supply, a.deficit, a.bonus, a.gap, a.cascade_len, tag
        );
    }

    println!("  Witness condition (full block):");
    for w in &report.witness_condition {
        println!(
            "    p={:<3} total_demand={:<4} top_supply={:<4} witness_gap={}",
            w.p, w.total_demand, w.top_supply, w.witness_gap
        );
    }

    let status = if report.witness_gaps_all_nonpositive {
        "YES *** CANDIDATE ***"
    } else {
        "NO"
    };
    println!(
        "  Failing primes compensated: {}/{}  |  Witness condition met: {}",
        report.compensated_count, report.failing_prime_count, status
    );
    println!();
}

fn print_summary(reports: &[NearMissReport], start: u64, count: u64, k: u32) {
    let total = reports.len();
    let candidates = reports
        .iter()
        .filter(|r| r.witness_gaps_all_nonpositive)
        .count();
    let bp = barrier_primes(k);

    println!("---");
    println!("Scanned {} integers starting at {} (k={})", count, start, k);
    let pure_count = reports.iter().filter(|r| r.pure_barrier_only).count();
    println!(
        "Near-misses found: {} ({} pure barrier-only)",
        total, pure_count
    );

    if total > 0 {
        let worst = reports
            .iter()
            .flat_map(|r| r.witness_condition.iter())
            .map(|w| w.witness_gap)
            .max()
            .unwrap_or(0);

        let best_worst_gap: Option<(i64, u64)> = reports
            .iter()
            .filter(|r| r.failing_prime_count > 0)
            .map(|r| {
                let worst_failing = r
                    .barrier_primes
                    .iter()
                    .filter(|a| a.is_failing)
                    .map(|a| a.gap)
                    .max()
                    .unwrap_or(0);
                (worst_failing, r.near_miss.n)
            })
            .min_by_key(|(gap, _)| *gap);

        println!("  Worst witness_gap: {}", worst);
        if let Some((gap, n)) = best_worst_gap {
            println!(
                "  Best (closest to full compensation): gap={} (at n={})",
                gap, n
            );
        }

        // Per-prime compensation frequency
        let reports_with_failures: usize =
            reports.iter().filter(|r| r.failing_prime_count > 0).count();
        if reports_with_failures > 0 {
            let mut compensated_by_prime: Vec<(u64, usize)> = Vec::new();
            let mut never_compensated: Vec<u64> = Vec::new();
            for &p in &bp {
                let comp_count = reports
                    .iter()
                    .filter(|r| {
                        r.barrier_primes
                            .iter()
                            .any(|a| a.p == p && a.is_failing && a.gap <= 0)
                    })
                    .count();
                let fail_count = reports
                    .iter()
                    .filter(|r| r.barrier_primes.iter().any(|a| a.p == p && a.is_failing))
                    .count();
                if fail_count > 0 {
                    compensated_by_prime.push((p, comp_count));
                    if comp_count == 0 {
                        never_compensated.push(p);
                    }
                }
            }
            compensated_by_prime.sort_by(|a, b| b.1.cmp(&a.1));
            let top: Vec<String> = compensated_by_prime
                .iter()
                .take(5)
                .map(|(p, c)| format!("p={} ({}/{})", p, c, reports_with_failures))
                .collect();
            if !top.is_empty() {
                println!("  Primes most often compensated: {}", top.join(", "));
            }
            if !never_compensated.is_empty() {
                println!("  Primes never compensated: {:?}", never_compensated);
            }
        }
    }

    println!("  Near-misses with ALL witness gaps <= 0: {}", candidates);

    let pure_candidates = reports
        .iter()
        .filter(|r| r.witness_gaps_all_nonpositive && r.pure_barrier_only)
        .count();
    println!(
        "  Pure barrier candidates (no large-prime kill): {}",
        pure_candidates
    );

    if candidates > 0 {
        println!("  *** CANDIDATE NON-GOVERNOR WITNESSES FOUND! ***");
        for r in reports.iter().filter(|r| r.witness_gaps_all_nonpositive) {
            let tag = if r.pure_barrier_only {
                " [PURE - INVESTIGATE]"
            } else {
                " [has non-barrier failure]"
            };
            println!("    n={}, m={}{}", r.near_miss.n, r.near_miss.m, tag);
        }
    }
}

fn print_json(reports: &[NearMissReport]) {
    println!("{}", serde_json::to_string_pretty(reports).unwrap());
}
