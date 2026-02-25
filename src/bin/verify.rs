//! Standalone verification tool for Erdős 396 witnesses.
//!
//! Usage:
//!   verify -k `<k>` -n `<n>`     - Verify a single candidate
//!   verify --known              - Verify all known OEIS witnesses

use clap::Parser;
use erdos396::{verify::WitnessVerifier, BuildInfo, KNOWN_WITNESSES};
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(name = "verify")]
#[command(author = "Erdős 396 Research Team")]
#[command(version = "1.0.0")]
#[command(about = "Verify Erdős Problem #396 witnesses")]
struct Cli {
    /// Target k value
    #[arg(short = 'k', long)]
    target_k: Option<u32>,

    /// Candidate n value
    #[arg(short = 'n', long)]
    n: Option<u64>,

    /// Verify all known OEIS witnesses
    #[arg(long)]
    known: bool,

    /// Show detailed verification output
    #[arg(short = 'v', long)]
    verbose: bool,

    /// Output machine-readable JSON
    #[arg(long)]
    json: bool,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if cli.known {
        if cli.target_k.is_some() || cli.n.is_some() {
            anyhow::bail!("--known cannot be combined with -k/--target-k or -n/--n");
        }
        verify_all_known(cli.verbose, cli.json)?;
    } else if let (Some(k), Some(n)) = (cli.target_k, cli.n) {
        if k == 0 {
            anyhow::bail!("--target-k must be >= 1");
        }
        verify_single(k, n, cli.verbose, cli.json)?;
    } else {
        eprintln!("Error: Must provide either --known or both -k and -n");
        std::process::exit(1);
    }

    Ok(())
}

#[derive(Serialize)]
struct JsonVerification {
    build: BuildInfo,
    k: u32,
    n: u64,
    is_valid: bool,
    is_governor_run: bool,
    failing_prime: Option<u64>,
    demand: Vec<(u64, u64)>,
    supply: Vec<(u64, u64)>,
}

fn verify_single(k: u32, n: u64, verbose: bool, json: bool) -> anyhow::Result<()> {
    let verifier = WitnessVerifier::new(n);
    let result = verifier.verify(k, n)?;

    if json {
        let mut demand: Vec<(u64, u64)> = result.demand.iter().map(|(&p, &v)| (p, v)).collect();
        demand.sort_by_key(|(p, _)| *p);
        let mut supply: Vec<(u64, u64)> = result.supply.iter().map(|(&p, &v)| (p, v)).collect();
        supply.sort_by_key(|(p, _)| *p);
        let out = JsonVerification {
            build: BuildInfo::gather(),
            k,
            n,
            is_valid: result.is_valid,
            is_governor_run: result.is_governor_run,
            failing_prime: result.failing_prime,
            demand,
            supply,
        };
        println!("{}", serde_json::to_string_pretty(&out)?);
        std::process::exit(if result.is_valid { 0 } else { 1 });
    }

    println!("Verifying k={}, n={}", k, n);

    if verbose {
        println!("{}", result.summary(k, n));
    } else {
        if result.is_valid {
            println!("✓ VALID witness");
        } else {
            println!("✗ NOT a valid witness");
            if let Some(p) = result.failing_prime {
                println!(
                    "  Fails at p={} (demand={}, supply={})",
                    p,
                    result.demand.get(&p).unwrap_or(&0),
                    result.supply.get(&p).unwrap_or(&0)
                );
            }
        }
        println!("  Governor run: {}", result.is_governor_run);
    }

    if result.is_valid {
        std::process::exit(0);
    } else {
        std::process::exit(1);
    }
}

#[derive(Serialize)]
struct JsonKnown {
    build: BuildInfo,
    witnesses: Vec<(u32, u64, bool)>,
}

fn verify_all_known(verbose: bool, json: bool) -> anyhow::Result<()> {
    if !json {
        println!("Verifying all known OEIS A375077 witnesses...\n");
    }

    let max_n = KNOWN_WITNESSES
        .iter()
        .map(|(_, n)| *n)
        .max()
        .unwrap_or(1000);
    let verifier = WitnessVerifier::new(max_n);

    let mut all_pass = true;
    let mut out: Vec<(u32, u64, bool)> = Vec::with_capacity(KNOWN_WITNESSES.len());

    for &(k, n) in KNOWN_WITNESSES {
        let result = verifier.verify(k, n)?;
        out.push((k, n, result.is_valid));
        let status = if result.is_valid {
            "✓ PASS"
        } else {
            "✗ FAIL"
        };

        if !json {
            println!("k={:2}, n={:14}: {}", k, n, status);

            if verbose && !result.is_valid {
                println!("{}", result.summary(k, n));
            }
        }

        if !result.is_valid {
            all_pass = false;
        }
    }

    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&JsonKnown {
                build: BuildInfo::gather(),
                witnesses: out
            })?
        );
        if all_pass {
            return Ok(());
        }
        anyhow::bail!("Some witnesses failed verification!")
    }

    println!();
    if all_pass {
        println!(
            "All {} known witnesses verified successfully!",
            KNOWN_WITNESSES.len()
        );
        Ok(())
    } else {
        anyhow::bail!("Some witnesses failed verification!")
    }
}
