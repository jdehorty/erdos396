//! Standalone verification tool for Erdős 396 witnesses.
//!
//! Usage:
//!   verify `<k>` `<n>`           - Verify a single candidate
//!   verify --known           - Verify all known OEIS witnesses

use clap::Parser;
use erdos396::{verify::WitnessVerifier, KNOWN_WITNESSES};

#[derive(Parser, Debug)]
#[command(name = "verify")]
#[command(author = "Erdős 396 Research Team")]
#[command(version = "1.0.0")]
#[command(about = "Verify Erdős Problem #396 witnesses")]
struct Cli {
    /// Target k value
    #[arg(short = 'k', long)]
    k: Option<u32>,

    /// Candidate n value
    #[arg(short = 'n', long)]
    n: Option<u64>,

    /// Verify all known OEIS witnesses
    #[arg(long)]
    known: bool,

    /// Show detailed verification output
    #[arg(short = 'v', long)]
    verbose: bool,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if cli.known {
        verify_all_known(cli.verbose)?;
    } else if let (Some(k), Some(n)) = (cli.k, cli.n) {
        verify_single(k, n, cli.verbose)?;
    } else {
        eprintln!("Error: Must provide either --known or both -k and -n");
        std::process::exit(1);
    }

    Ok(())
}

fn verify_single(k: u32, n: u64, verbose: bool) -> anyhow::Result<()> {
    println!("Verifying k={}, n={}", k, n);

    let verifier = WitnessVerifier::new(n);
    let result = verifier.verify(k, n);

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

fn verify_all_known(verbose: bool) -> anyhow::Result<()> {
    println!("Verifying all known OEIS A375077 witnesses...\n");

    let max_n = KNOWN_WITNESSES.iter().map(|(_, n)| *n).max().unwrap_or(1000);
    let verifier = WitnessVerifier::new(max_n);

    let mut all_pass = true;

    for &(k, n) in KNOWN_WITNESSES {
        let result = verifier.verify(k, n);
        let status = if result.is_valid { "✓ PASS" } else { "✗ FAIL" };

        println!("k={:2}, n={:14}: {}", k, n, status);

        if verbose && !result.is_valid {
            println!("{}", result.summary(k, n));
        }

        if !result.is_valid {
            all_pass = false;
        }
    }

    println!();
    if all_pass {
        println!("All {} known witnesses verified successfully!", KNOWN_WITNESSES.len());
        Ok(())
    } else {
        anyhow::bail!("Some witnesses failed verification!")
    }
}
