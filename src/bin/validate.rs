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
//! For every integer `n` in `[start, end)`:
//!   1. Check the witness condition at each prime `p < 2k+1`:
//!
//!      ```text
//!      Σ v_p(n-i) for i=0..k  ≤  v_p(C(2n, n))
//!      ```
//!
//!   2. If ALL barrier primes pass → full p-adic verification (all primes)
//!   3. Report confirmed witnesses
//!
//! ## Why this finds all witnesses
//!
//! A number `n` is a `k`-witness iff for every prime `p`:
//!
//! `Σ v_p(n-i) (i=0..k) ≤ v_p(C(2n,n))`.
//!
//! This program enumerates every `n` in `[start, end)`. It first checks the
//! inequality for the barrier primes `p < 2k+1` (a necessary condition), then
//! fully verifies each screen-pass candidate at *all* primes. Therefore the
//! output contains exactly the witnesses in the range, assuming the arithmetic
//! kernels are correct.
//!
//! The Lean proof (Small Prime Barrier Theorem / Corollary 6) explains why the
//! barrier primes are the only place a witness can exhibit “non-governor”
//! behavior, making this screen effective for ruling out missed witnesses.
//!
//! ## Purpose
//!
//! Run this as a second pass over ranges already covered by the Governor Set
//! search. If the validator finds no new witnesses below the known minimum,
//! this constitutes a *computational* certificate of minimality — independent
//! of the Governor Set Completeness Conjecture.
//!
//! Implementation note: the underlying algorithm is justified by a Lean proof,
//! but this Rust binary is not itself formally verified. See `docs/trust.md`.

use clap::Parser;
use erdos396::governor::{vp_central_binom_dispatch, vp_factorial};
use erdos396::int_math::{isqrt_2n_u64, isqrt_u128};
use erdos396::sieve::PrimeSieve;
use erdos396::verify::{primes_up_to, WitnessVerifier};
use erdos396::BuildInfo;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::sync::Mutex;
use std::time::{Duration, Instant};

mod serde_u128_string {
    use serde::{Deserialize, Deserializer, Serializer};

    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Repr {
        Str(String),
        Num(u64),
    }

    pub fn serialize<S>(value: &u128, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&value.to_string())
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<u128, D::Error>
    where
        D: Deserializer<'de>,
    {
        match Repr::deserialize(deserializer)? {
            Repr::Str(s) => s.parse::<u128>().map_err(serde::de::Error::custom),
            Repr::Num(n) => Ok(n as u128),
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "validate")]
#[command(about = "Small-Prime Sieve Validator — provably complete witness search (Corollary 6)")]
#[command(long_about = r#"
Validates that no Erdős 396 k-witnesses were missed by the Governor Set search.

Uses the Small Prime Barrier Theorem: checks the witness condition at primes
p < 2k+1 for EVERY integer in the range (not just Governor Set members).
Any witness — governor run or not — must pass this screen.

This is a second-pass tool. Run it over ranges already covered by the main
search (erdos396) to produce a computational certificate of minimality
(see docs/trust.md).

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

    /// Sub-batch size for the fused sieve computation (default: 1M).
    /// Each worker processes its range in chunks of this size.
    #[arg(long, default_value = "1000000")]
    batch_size: u64,

    /// Verbose output
    #[arg(short = 'v', long)]
    verbose: bool,

    /// Startup cross-check samples (Kummer vs Legendre supply, and direct vs Legendre demand).
    ///
    /// This is a low-cost sanity check that the compiled binary behaves as expected on the
    /// target machine. Set to 0 to disable.
    #[arg(long, default_value = "0")]
    self_check_samples: u32,

    /// Periodically audit arithmetic kernels during the scan (0 disables).
    ///
    /// This runs the same checks as `--self-check-samples`, but on an in-range
    /// value encountered during the scan. It is intended as a low-rate guard
    /// against silent miscompilation or hardware faults in long runs.
    #[arg(long, default_value = "0")]
    audit_interval: u64,
}

#[derive(Serialize, Deserialize)]
struct ValidateCheckpoint {
    worker_id: usize,
    target_k: u32,
    start: u64,
    end: u64,
    current_pos: u64,
    checked: u64,
    /// Sum of all n values checked by this worker (coverage invariant).
    #[serde(default, with = "serde_u128_string")]
    sum_checked: u128,
    /// XOR of all n values checked by this worker (coverage invariant).
    #[serde(default)]
    xor_checked: u64,
    small_prime_candidates: u64,
    confirmed_witnesses: Vec<u64>,
    false_alarms: u64,
    timestamp: String,
}

/// Sum of integers in [a, b) computed in u128 (exact for all u64 inputs).
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

/// XOR of integers in [0, n] (inclusive).
#[inline]
fn xor_upto(n: u64) -> u64 {
    match n & 3 {
        0 => n,
        1 => 1,
        2 => n + 1,
        _ => 0,
    }
}

/// XOR of integers in [a, b).
#[inline]
fn xor_range(a: u64, b: u64) -> u64 {
    if b <= a {
        return 0;
    }
    let hi = xor_upto(b - 1);
    let lo = if a == 0 { 0 } else { xor_upto(a - 1) };
    hi ^ lo
}

#[cfg(test)]
mod tests {
    use super::*;
    use erdos396::sieve::PrimeSieve;

    fn naive_fails_at_large_prime(n: u64, barrier_threshold: u64, base_primes: &[u64]) -> bool {
        if n == 0 {
            return true;
        }
        if n <= 1 {
            return false;
        }

        let mut remaining = n;
        for &p in base_primes {
            if (p as u128) * (p as u128) > remaining as u128 {
                break;
            }
            if remaining % p == 0 {
                let mut exp = 0u32;
                while remaining % p == 0 {
                    exp += 1;
                    remaining /= p;
                }
                if p >= barrier_threshold {
                    let supply = vp_supply(n, p);
                    if (exp as u64) > supply {
                        return true;
                    }
                }
            }
        }

        if remaining > 1 && remaining >= barrier_threshold {
            if remaining > isqrt_2n_u64(n) {
                return true;
            }
            if 1 > vp_supply(n, remaining) {
                return true;
            }
        }

        false
    }

    #[test]
    fn test_sum_range_u128_small() {
        assert_eq!(sum_range_u128(0, 0), 0);
        assert_eq!(sum_range_u128(0, 1), 0);
        assert_eq!(sum_range_u128(0, 2), 1);
        assert_eq!(sum_range_u128(0, 10), 45);
        assert_eq!(sum_range_u128(5, 6), 5);
        assert_eq!(sum_range_u128(5, 10), 35);
    }

    #[test]
    fn test_xor_range_small() {
        assert_eq!(xor_range(0, 0), 0);
        assert_eq!(xor_range(0, 1), 0);
        assert_eq!(xor_range(0, 5), 4); // 0^1^2^3^4 = 4
        assert_eq!(xor_range(5, 10), 5 ^ 6 ^ 7 ^ 8 ^ 9);
        assert_eq!(xor_range(1, 2), 1);
    }

    #[test]
    fn test_small_prime_screen_matches_naive_small_ranges() {
        for k in 1u32..=6 {
            let barrier_primes = primes_up_to(k.saturating_mul(2));
            let start_n = k as u64 + 1;
            let mut demand_sums = init_demand_window(start_n, k, &barrier_primes);

            for n in start_n..=5000u64 {
                let screen = passes_small_prime_screen_window(n, &barrier_primes, &demand_sums);

                let two_n = n * 2;
                let mut naive_ok = true;
                for (idx, &p) in barrier_primes.iter().enumerate() {
                    let demand_direct: u64 = (0..=k).map(|i| vp_term(n - i as u64, p)).sum();
                    assert_eq!(
                        demand_sums[idx], demand_direct,
                        "window demand mismatch at k={}, n={}, p={}",
                        k, n, p
                    );
                    let supply_legendre = vp_factorial(two_n, p) - 2 * vp_factorial(n, p);
                    if demand_direct > supply_legendre {
                        naive_ok = false;
                    }
                }

                assert_eq!(screen, naive_ok, "screen mismatch at k={}, n={}", k, n);

                if n < 5000 {
                    advance_demand_window(n, k, &barrier_primes, &mut demand_sums).unwrap();
                }
            }
        }
    }

    #[test]
    fn test_validate_batch_matches_naive_large_prime_failures() {
        let cases = [
            (0u64, 512usize),
            (10_000u64, 1024usize),
            (1_000_000u64, 2048usize),
            (25_000_000_000_000u64 - 4096, 4096usize),
        ];

        for k in 1u32..=8 {
            let barrier_threshold = 2 * k as u64 + 1;
            for &(lo, len) in &cases {
                let max_n = lo + len as u64 - 1;
                let sieve_limit = isqrt_u128(2u128 * max_n as u128) as u64 + 64;
                let sieve = PrimeSieve::new(sieve_limit);
                let batch = ValidateBatch::compute(lo, len, sieve.primes(), barrier_threshold);

                for idx in 0..len {
                    let n = lo + idx as u64;
                    let expected = naive_fails_at_large_prime(n, barrier_threshold, sieve.primes());
                    assert_eq!(
                        batch.fails_at_large_prime[idx], expected,
                        "fails-at-large mismatch: k={}, n={}, barrier={}",
                        k, n, barrier_threshold
                    );
                }
            }
        }
    }
}

/// Compute v_p(C(2n, n)) using the best available method for each prime.
#[inline]
fn vp_supply(n: u64, p: u64) -> u64 {
    vp_central_binom_dispatch(n, p)
}

/// Initialize the sliding-window demand sums:
/// `demand[p] = Σ_{i=0..k} v_p(n-i)` for the current `n`.
#[inline]
#[allow(dead_code)]
fn init_demand_window(n: u64, k: u32, barrier_primes: &[u64]) -> Vec<u64> {
    debug_assert!(n > k as u64, "init_demand_window requires n > k");
    let start = n - k as u64;
    let mut sums = Vec::with_capacity(barrier_primes.len());
    for &p in barrier_primes {
        let mut s = 0u64;
        for t in start..=n {
            s += vp_term(t, p);
        }
        sums.push(s);
    }
    sums
}

/// Advance the sliding window from `n` to `n+1`.
///
/// Old window: `{n-k, ..., n}`, new window: `{n-k+1, ..., n+1}`.
#[inline]
#[allow(dead_code)]
fn advance_demand_window(
    n: u64,
    k: u32,
    barrier_primes: &[u64],
    sums: &mut [u64],
) -> anyhow::Result<()> {
    debug_assert_eq!(barrier_primes.len(), sums.len());
    debug_assert!(n > k as u64, "advance_demand_window requires n > k");

    let leaving = n - k as u64; // term leaving the window
    let entering = n + 1; // term entering the window (caller ensures n+1 is in range)

    for (idx, &p) in barrier_primes.iter().enumerate() {
        let leave_v = vp_term(leaving, p);
        let enter_v = vp_term(entering, p);
        let cur = sums[idx];
        if leave_v > cur {
            anyhow::bail!(
                "Internal error: demand window underflow at n={}, k={}, p={} (cur={}, leaving={})",
                n,
                k,
                p,
                cur,
                leave_v
            )
        }
        sums[idx] = cur - leave_v + enter_v;
    }
    Ok(())
}

/// Check the witness condition at all barrier primes (p < 2k+1) using a demand window.
#[inline]
#[allow(dead_code)]
fn passes_small_prime_screen_window(n: u64, barrier_primes: &[u64], demand_sums: &[u64]) -> bool {
    debug_assert_eq!(barrier_primes.len(), demand_sums.len());
    for (idx, &p) in barrier_primes.iter().enumerate() {
        let demand = demand_sums[idx];
        if demand == 0 {
            continue; // No block term divisible by p → always pass
        }
        let supply = vp_supply(n, p);
        if demand > supply {
            return false;
        }
    }
    true
}

#[inline]
fn lcg_next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

#[inline]
fn vp_term(mut n: u64, p: u64) -> u64 {
    let mut v = 0u64;
    while n % p == 0 {
        v += 1;
        n /= p;
    }
    v
}

/// Secondary screen: check if all block terms pass the governor condition at
/// primes p >= 2k+1.  (Legacy per-number approach — superseded by ValidateBatch.)
#[allow(dead_code)]
fn block_passes_large_prime_screen(n: u64, k: u32, screen_primes: &[u64], b_squared: u64) -> bool {
    let barrier_min = 2 * k as u64 + 1;
    let sqrt_2n = isqrt_2n_u64(n);

    for j in 0..=k {
        let term = n - j as u64;
        if term <= 1 {
            continue;
        }

        let mut remaining = term;
        for &p in screen_primes {
            if p * p > remaining {
                break;
            }
            if remaining % p == 0 {
                let mut exp = 0u64;
                while remaining % p == 0 {
                    exp += 1;
                    remaining /= p;
                }
                if p >= barrier_min {
                    let supply = vp_supply(n, p);
                    if exp > supply {
                        return false;
                    }
                }
            }
            if remaining == 1 {
                break;
            }
        }

        // Cofactor analysis.  remaining > 1 means there are prime factors > B.
        if remaining > 1 && remaining >= barrier_min && remaining <= b_squared {
            // All prime factors of remaining are > B, and remaining ≤ B²,
            // so remaining is prime.
            if remaining > sqrt_2n {
                // v_p(C(2n,n)) = 0 for primes p > √(2n), demand ≥ 1
                return false;
            }
            let supply = vp_supply(n, remaining);
            if 1 > supply {
                return false;
            }
        }
        // If remaining > B², it could be composite.
        // We cannot determine its prime structure cheaply, so we do NOT
        // reject — this preserves soundness. The candidate proceeds
        // to full verification.
    }

    true
}

/// Fused sieve for the validate binary: computes, for each position in a batch,
/// whether that number fails the governor condition at any prime >= `barrier_threshold`.
///
/// By the Small Prime Barrier Theorem (Corollary 6), a true k-witness cannot have
/// any block term that fails the governor test at a prime q >= 2k+1.  So if
/// `fails_at_large_prime[i]` is true for any block term of candidate n, then n is
/// not a witness and can be rejected immediately.
///
/// The sieve is a modified version of `FusedBatchResult::compute` from prefilter.rs,
/// with one key difference: instead of tracking exact governor membership, it tracks
/// whether the governor failure (if any) occurs at a prime >= barrier_threshold.
/// This enables a powerful skip optimisation: once a position is known to fail at a
/// large prime, all further division work for that position is skipped.
struct ValidateBatch {
    /// For each position i: true iff (lo + i) fails governor at ANY prime >= barrier_threshold.
    fails_at_large_prime: Vec<bool>,
}

impl ValidateBatch {
    /// Compute the `fails_at_large_prime` bitmap for every number in `[lo, lo + len)`.
    ///
    /// `base_primes` must contain all primes up to at least `sqrt(2 * (lo + len))`.
    /// `barrier_threshold` is `2k + 1` — primes at or above this value are "large".
    fn compute(lo: u64, len: usize, base_primes: &[u64], barrier_threshold: u64) -> Self {
        if len == 0 {
            return Self {
                fails_at_large_prime: vec![],
            };
        }

        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(1);
        let large_start = base_primes.partition_point(|&p| p < barrier_threshold);

        // remaining[i] = cofactor of (lo + i) after dividing out primes <= sieve_limit.
        let mut remaining = Vec::with_capacity(len);
        for i in 0..len {
            remaining.push(lo + i as u64);
        }
        let mut fails_large = vec![false; len];

        // Edge cases: n <= 1
        for i in 0..len {
            let n = lo + i as u64;
            if n <= 1 {
                remaining[i] = 1;
                if n == 0 {
                    // n=0 is not a governor; arguably fails at every prime.
                    // Mark it so any block containing n=0 is rejected.
                    fails_large[i] = true;
                }
            }
        }

        // === Phase 1: strip all "small-prime" factors (p < 2k+1) ===
        //
        // We only care about failures at large primes. For small primes we still
        // divide out factors to keep `remaining` accurate, but we skip any supply
        // checks entirely.
        for &p in &base_primes[..large_start] {
            if p > sieve_limit {
                break;
            }

            if p == 2 {
                // Strip powers of 2 in one pass using trailing-zero count.
                for value in &mut remaining {
                    let v = value.trailing_zeros();
                    if v > 0 {
                        *value >>= v;
                    }
                }
                continue;
            }

            // Find first multiple of p >= lo
            let first_multiple = if lo == 0 {
                p // skip n=0
            } else if lo % p == 0 {
                lo
            } else {
                lo + (p - lo % p)
            };

            let mut idx = (first_multiple - lo) as usize;
            while idx < len {
                while remaining[idx] % p == 0 {
                    remaining[idx] /= p;
                }
                idx += p as usize;
            }
        }

        // === Phase 2: strip large-prime factors and mark large-prime failures ===
        for &p in &base_primes[large_start..] {
            if p > sieve_limit {
                break;
            }

            // Find first multiple of p >= lo
            let first_multiple = if lo == 0 {
                p // skip n=0
            } else if lo % p == 0 {
                lo
            } else {
                lo + (p - lo % p)
            };

            let p_sq = (p as u128) * (p as u128);
            let mut idx = (first_multiple - lo) as usize;
            while idx < len {
                // Skip positions already known to fail at a large prime.
                if fails_large[idx] {
                    idx += p as usize;
                    continue;
                }

                if remaining[idx] % p == 0 {
                    let mut exp = 0u32;
                    while remaining[idx] % p == 0 {
                        exp += 1;
                        remaining[idx] /= p;
                    }

                    let n = lo + idx as u64;
                    let two_n = 2u128 * (n as u128);
                    if p_sq > two_n {
                        // If p > √(2n), then v_p(C(2n,n)) = 0.
                        // Since exp >= 1, this is an immediate large-prime failure.
                        fails_large[idx] = true;
                    } else {
                        let supply = vp_supply(n, p);
                        if (exp as u64) > supply {
                            fails_large[idx] = true;
                        }
                    }
                }

                idx += p as usize;
            }
        }

        // Post-sieve: remaining > 1 means a prime factor > sieve_limit >= sqrt(2n).
        // Such a prime q satisfies v_q(C(2n,n)) = 0 but v_q(n) >= 1 → governor failure.
        // Since sieve_limit >> barrier_threshold, this is always a "large prime" failure.
        for i in 0..len {
            if !fails_large[i] && remaining[i] > 1 {
                fails_large[i] = true;
            }
        }

        Self {
            fails_at_large_prime: fails_large,
        }
    }
}

/// Check the witness condition at all barrier primes from scratch (no sliding window).
///
/// This is only called for the ~1-in-10^7 candidates where every block term passes
/// the large-prime governor check.  The cost per call is O(k * |barrier_primes|)
/// Kummer evaluations, which is negligible at this call rate.
#[inline]
fn barrier_check_from_scratch(n: u64, k: u32, barrier_primes: &[u64]) -> bool {
    for &p in barrier_primes {
        let demand: u64 = (0..=k).map(|j| vp_term(n - j as u64, p)).sum();
        let supply = vp_supply(n, p);
        if demand > supply {
            return false;
        }
    }
    true
}

fn check_kernels_at(
    n: u64,
    k: u32,
    barrier_primes: &[u64],
    window_demand: Option<&[u64]>,
) -> anyhow::Result<()> {
    if n <= k as u64 {
        anyhow::bail!("kernel check requires n > k (k={}, n={})", k, n);
    }
    if let Some(w) = window_demand {
        if w.len() != barrier_primes.len() {
            anyhow::bail!(
                "kernel check window size mismatch: demand_sums has len {}, barrier_primes has len {}",
                w.len(),
                barrier_primes.len()
            );
        }
    }
    let two_n = n
        .checked_mul(2)
        .ok_or_else(|| anyhow::anyhow!("kernel check requires 2n to fit in u64 (n={})", n))?;

    for (idx, &p) in barrier_primes.iter().enumerate() {
        // Supply: Kummer vs Legendre.
        let kummer = vp_supply(n, p);
        let legendre = vp_factorial(two_n, p) - 2 * vp_factorial(n, p);
        if kummer != legendre {
            anyhow::bail!(
                "Supply mismatch at n={}, p={}: kummer={}, legendre={}",
                n,
                p,
                kummer,
                legendre
            );
        }

        // Demand: direct term valuation sum vs factorial-difference (Legendre).
        let block_bottom = n - k as u64;
        let demand_legendre = vp_factorial(n, p) - vp_factorial(block_bottom - 1, p);
        let demand_direct: u64 = (0..=k).map(|i| vp_term(n - i as u64, p)).sum();
        if demand_direct != demand_legendre {
            anyhow::bail!(
                "Demand mismatch at n={}, k={}, p={}: direct={}, legendre={}",
                n,
                k,
                p,
                demand_direct,
                demand_legendre
            );
        }
        if let Some(w) = window_demand {
            let window = w[idx];
            if demand_direct != window {
                anyhow::bail!(
                    "Window mismatch at n={}, k={}, p={}: window={}, direct={}",
                    n,
                    k,
                    p,
                    window,
                    demand_direct
                );
            }
        }
    }

    Ok(())
}

fn run_self_check(
    samples: u32,
    k: u32,
    scan_start: u64,
    end: u64,
    barrier_primes: &[u64],
) -> anyhow::Result<()> {
    let scan_size = end.saturating_sub(scan_start);
    if samples == 0 || scan_size == 0 {
        return Ok(());
    }

    let mut state = 0x9c9b_a3b2_4a3f_1d77u64 ^ (k as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    for _ in 0..samples {
        let n = scan_start + (lcg_next(&mut state) % scan_size);
        check_kernels_at(n, k, barrier_primes, None)?;
    }

    Ok(())
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
    sum_checked: u128,
    xor_checked: u64,
    candidates: u64,
    witnesses: Vec<u64>,
    false_alarms: u64,
}

#[derive(Serialize)]
struct ValidateReport {
    version: &'static str,
    build: BuildInfo,
    args: Vec<String>,
    target_k: u32,
    start: u64,
    end: u64,
    scan_start: u64,
    workers: usize,
    self_check_samples: u32,
    audit_interval: u64,
    barrier_primes: Vec<u64>,
    checked: u64,
    expected_checked: u64,
    #[serde(with = "serde_u128_string")]
    sum_checked: u128,
    #[serde(with = "serde_u128_string")]
    expected_sum_checked: u128,
    xor_checked: u64,
    expected_xor_checked: u64,
    screen_pass: u64,
    false_alarms: u64,
    witnesses: Vec<u64>,
    duration_secs: f64,
    rate_per_sec: f64,
    generated_at: String,
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
                    let worker_range = w_end.saturating_sub(w_start);
                    let pct = if worker_range > 0 {
                        (cp.current_pos.saturating_sub(w_start)) as f64 / worker_range as f64
                            * 100.0
                    } else {
                        100.0
                    };
                    log::info!(
                        "Worker {} resuming from {} ({:.2}% of worker range)",
                        worker_id,
                        cp.current_pos,
                        pct
                    );
                    let initial_scan_start = std::cmp::max(w_start, k as u64 + 1);
                    if cp.checked > 0
                        && cp.sum_checked == 0
                        && cp.xor_checked == 0
                        && cp.current_pos >= initial_scan_start
                        && cp.checked == cp.current_pos - initial_scan_start
                    {
                        let mut patched = cp;
                        patched.sum_checked =
                            sum_range_u128(initial_scan_start, patched.current_pos);
                        patched.xor_checked = xor_range(initial_scan_start, patched.current_pos);
                        return patched;
                    }
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
        sum_checked: 0,
        xor_checked: 0,
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
    if k == 0 {
        anyhow::bail!("--target-k must be >= 1");
    }
    let barrier_threshold = k.saturating_mul(2);
    let barrier_primes = primes_up_to(barrier_threshold);
    let num_workers = cli.workers.unwrap_or_else(num_cpus::get);
    let scan_start_global = std::cmp::max(cli.start, k as u64 + 1);
    let scan_size = cli.end.saturating_sub(scan_start_global);

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
    println!("  Batch size:        {}M", cli.batch_size / 1_000_000);
    println!("  Output dir:        {:?}", cli.output_dir);
    if cli.self_check_samples > 0 {
        println!(
            "  Self-check:        {} samples (supply + demand cross-check)",
            cli.self_check_samples
        );
    }
    if cli.audit_interval > 0 {
        println!(
            "  Audit interval:    {} scanned values (in-run kernel checks)",
            cli.audit_interval
        );
    }
    println!();
    println!("  Theorem (Corollary 6): every true k-witness passes this barrier-prime screen.");
    println!("  Note: this program is tested but not formally verified; see docs/trust.md.");
    println!("═══════════════════════════════════════════════════════════════════");
    println!();

    std::fs::create_dir_all(&cli.output_dir)?;

    if cli.self_check_samples > 0 {
        log::info!(
            "Running self-check: {} random samples (Kummer vs Legendre supply; direct vs Legendre demand)...",
            cli.self_check_samples
        );
        run_self_check(
            cli.self_check_samples,
            k,
            scan_start_global,
            cli.end,
            &barrier_primes,
        )?;
        log::info!("Self-check passed.");
    }

    // The fused sieve needs primes up to √(2·end) (not just √(end)) because the
    // governor check uses v_p(C(2n,n)), which requires sieving up to that bound.
    // This sieve is also sufficient for full verification (which only needs √(end)).
    log::info!("Building prime sieve (up to √(2·end))...");
    let sieve_start = Instant::now();
    let sieve_limit = isqrt_2n_u64(cli.end) + 1000;
    let sieve = PrimeSieve::new(sieve_limit);
    log::info!(
        "Sieve built: {} primes up to {} in {:?}",
        sieve.len(),
        sieve_limit,
        sieve_start.elapsed()
    );

    // The barrier threshold: primes >= this value are "large primes" in the
    // barrier theorem.  A block term failing governor at such a prime means
    // the candidate is not a witness.
    let sieve_barrier_threshold = 2 * k as u64 + 1;

    let verifier = WitnessVerifier::with_sieve(sieve.clone());

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

    // Sanity-check partitioning: contiguous, non-overlapping coverage.
    if let Some((_, first_start, _)) = ranges.first() {
        if *first_start != cli.start {
            anyhow::bail!(
                "Internal error: first worker start {} != cli.start {}",
                first_start,
                cli.start
            );
        }
    }
    if let Some((_, _, last_end)) = ranges.last() {
        if *last_end != cli.end {
            anyhow::bail!(
                "Internal error: last worker end {} != cli.end {}",
                last_end,
                cli.end
            );
        }
    }
    for w in ranges.windows(2) {
        let (_, _, a_end) = w[0];
        let (_, b_start, _) = w[1];
        if a_end != b_start {
            anyhow::bail!(
                "Internal error: worker ranges not contiguous ({} != {})",
                a_end,
                b_start
            );
        }
    }

    let total_checked = Arc::new(AtomicU64::new(0));
    let total_candidates = Arc::new(AtomicU64::new(0));
    let total_witnesses = Arc::new(AtomicU64::new(0));
    let stop_reporter = Arc::new(AtomicBool::new(false));
    let coverage_error = Arc::new(AtomicBool::new(false));
    let audit_error = Arc::new(AtomicBool::new(false));
    let audit_message: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));

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
                let pct = if scan_size > 0 {
                    c as f64 / scan_size as f64 * 100.0
                } else {
                    100.0
                };
                let eta_hours = if rate > 0.0 {
                    (scan_size.saturating_sub(c)) as f64 / rate / 3600.0
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
            let coverage_error = coverage_error.clone();
            let audit_error = audit_error.clone();
            let audit_message = audit_message.clone();
            let audit_interval = cli.audit_interval;
            let batch_size = cli.batch_size;
            let checkpoint_path = cli
                .output_dir
                .join(format!("validate_k{}_w{:02}.json", k, worker_id));

            let mut state = load_or_init(worker_id, k, w_start, w_end, &checkpoint_path);
            let initial_scan_start = std::cmp::max(w_start, k as u64 + 1);
            let scan_start = std::cmp::max(state.current_pos, k as u64 + 1);
            let mut checkpoint_counter = 0u64;

            if scan_start >= w_end {
                // This worker range is entirely below `k+1` (or already complete).
                // Mark it done so coverage checks are stable.
                state.current_pos = w_end;
            }

            if scan_start < w_end {
                // ── Batch-based processing ──────────────────────────────────
                //
                // For each sub-batch of `batch_size` numbers:
                //   1. Compute ValidateBatch (fused sieve) for the batch,
                //      extended by k positions to cover all block terms.
                //   2. For each n in the batch: check if ANY block term
                //      fails governor at a large prime (O(k) bitmap lookups).
                //   3. If all block terms pass → check barrier screen from
                //      scratch (extremely rare) → full verification.

                let mut n_pos = scan_start;

                while n_pos < w_end {
                    if audit_error.load(Ordering::Relaxed) {
                        break;
                    }

                    let batch_end = std::cmp::min(n_pos + batch_size, w_end);
                    // The sieve must cover block terms: for n in [n_pos, batch_end),
                    // block terms range from (n_pos - k) to (batch_end - 1).
                    let sieve_lo = n_pos.saturating_sub(k as u64);
                    let sieve_len = (batch_end - sieve_lo) as usize;

                    let batch = ValidateBatch::compute(
                        sieve_lo,
                        sieve_len,
                        sieve.primes(),
                        sieve_barrier_threshold,
                    );

                    for n in n_pos..batch_end {
                        if audit_error.load(Ordering::Relaxed) {
                            break;
                        }

                        state.checked += 1;
                        state.sum_checked += n as u128;
                        state.xor_checked ^= n;
                        state.current_pos = n + 1;
                        total_checked.fetch_add(1, Ordering::Relaxed);

                        // Periodic audit: cross-check Kummer vs Legendre supply
                        if audit_interval > 0 && state.checked % audit_interval == 0 {
                            if let Err(e) =
                                check_kernels_at(n, k, &barrier_primes, None)
                            {
                                audit_error.store(true, Ordering::Relaxed);
                                let mut msg = audit_message.lock().unwrap();
                                if msg.is_none() {
                                    *msg = Some(format!("Audit failed at n={}: {:#}", n, e));
                                }
                                break;
                            }
                        }

                        // Block term check: does ANY block term fail governor at
                        // a prime >= 2k+1?  By the barrier theorem, a true witness
                        // cannot have such a failure, so this is a sound rejection.
                        let mut any_fails_large = false;
                        for j in 0..=k {
                            let term = n - j as u64;
                            let idx = (term - sieve_lo) as usize;
                            if batch.fails_at_large_prime[idx] {
                                any_fails_large = true;
                                break;
                            }
                        }

                        if !any_fails_large {
                            // Extremely rare (~1 in 10^7): all block terms pass
                            // governor at large primes.  Check barrier screen.
                            state.small_prime_candidates += 1;
                            total_candidates.fetch_add(1, Ordering::Relaxed);

                            if barrier_check_from_scratch(n, k, &barrier_primes) {
                                // Cross-check with full p-adic verification
                                let result = match verifier.verify(k, n) {
                                    Ok(r) => r,
                                    Err(e) => {
                                        audit_error.store(true, Ordering::Relaxed);
                                        let mut msg = audit_message.lock().unwrap();
                                        if msg.is_none() {
                                            *msg = Some(format!(
                                                "Internal verification error at n={}: {}",
                                                n, e
                                            ));
                                        }
                                        break;
                                    }
                                };
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
                            } else {
                                state.false_alarms += 1;
                            }
                        }
                    }

                    // Update position and checkpoint at batch boundaries
                    state.current_pos = batch_end;
                    checkpoint_counter += batch_end - n_pos;
                    if checkpoint_counter >= cli.checkpoint_interval {
                        state.timestamp = chrono::Utc::now().to_rfc3339();
                        save_checkpoint(&checkpoint_path, &state);
                        checkpoint_counter = 0;
                    }

                    n_pos = batch_end;
                }
            }

            // Coverage sanity checks: this catches gaps/overlaps and resume errors.
            if !audit_error.load(Ordering::Relaxed) {
                let expected_checked = w_end.saturating_sub(initial_scan_start);
                let expected_sum = sum_range_u128(initial_scan_start, w_end);
                let expected_xor = xor_range(initial_scan_start, w_end);
                if state.checked != expected_checked || state.current_pos != w_end {
                    coverage_error.store(true, Ordering::Relaxed);
                    log::error!(
                        "Worker {} coverage mismatch: checked={} (expected {}), current_pos={} (expected {})",
                        worker_id,
                        state.checked,
                        expected_checked,
                        state.current_pos,
                        w_end
                    );
                }
                if state.sum_checked != expected_sum || state.xor_checked != expected_xor {
                    coverage_error.store(true, Ordering::Relaxed);
                    log::error!(
                        "Worker {} coverage invariant mismatch: sum={} (expected {}), xor={} (expected {})",
                        worker_id,
                        state.sum_checked,
                        expected_sum,
                        state.xor_checked,
                        expected_xor
                    );
                }
            }

            state.timestamp = chrono::Utc::now().to_rfc3339();
            save_checkpoint(&checkpoint_path, &state);

            log::info!(
                "Worker {} done: checked={}, block_pass={}, witnesses={}, false_alarms={}",
                worker_id,
                state.checked,
                state.small_prime_candidates,
                state.confirmed_witnesses.len(),
                state.false_alarms
            );

            WorkerResult {
                checked: state.checked,
                sum_checked: state.sum_checked,
                xor_checked: state.xor_checked,
                candidates: state.small_prime_candidates,
                witnesses: state.confirmed_witnesses,
                false_alarms: state.false_alarms,
            }
        })
        .collect();

    let duration = search_start.elapsed();

    stop_reporter.store(true, Ordering::Relaxed);
    let _ = reporter.join();

    if audit_error.load(Ordering::Relaxed) {
        if let Some(msg) = audit_message.lock().unwrap().as_ref() {
            anyhow::bail!("{}", msg);
        }
        anyhow::bail!("Audit failed (no details available)");
    }

    let agg_checked: u64 = results.iter().map(|r| r.checked).sum();
    let agg_sum_checked: u128 = results.iter().map(|r| r.sum_checked).sum();
    let agg_xor_checked: u64 = results.iter().map(|r| r.xor_checked).fold(0, |a, b| a ^ b);
    let agg_candidates: u64 = results.iter().map(|r| r.candidates).sum();
    let agg_false_alarms: u64 = results.iter().map(|r| r.false_alarms).sum();
    let mut agg_witnesses: Vec<u64> = results
        .iter()
        .flat_map(|r| r.witnesses.iter().copied())
        .collect();
    agg_witnesses.sort();
    agg_witnesses.dedup();

    let rate = agg_checked as f64 / duration.as_secs_f64();

    let expected_checked_total = scan_size;
    if agg_checked != expected_checked_total {
        anyhow::bail!(
            "Coverage mismatch: checked {} values but expected {} (range [{}, {}), k={})",
            agg_checked,
            expected_checked_total,
            cli.start,
            cli.end,
            k
        );
    }
    let expected_sum_total = sum_range_u128(scan_start_global, cli.end);
    let expected_xor_total = xor_range(scan_start_global, cli.end);
    if agg_sum_checked != expected_sum_total || agg_xor_checked != expected_xor_total {
        anyhow::bail!(
            "Coverage invariant mismatch: sum {} (expected {}), xor {} (expected {})",
            agg_sum_checked,
            expected_sum_total,
            agg_xor_checked,
            expected_xor_total
        );
    }
    if coverage_error.load(Ordering::Relaxed) {
        anyhow::bail!("One or more workers reported coverage inconsistencies; see logs");
    }

    let report = ValidateReport {
        version: env!("CARGO_PKG_VERSION"),
        build: BuildInfo::gather(),
        args: std::env::args().collect(),
        target_k: k,
        start: cli.start,
        end: cli.end,
        scan_start: scan_start_global,
        workers: num_workers,
        self_check_samples: cli.self_check_samples,
        audit_interval: cli.audit_interval,
        barrier_primes: barrier_primes.clone(),
        checked: agg_checked,
        expected_checked: expected_checked_total,
        sum_checked: agg_sum_checked,
        expected_sum_checked: expected_sum_total,
        xor_checked: agg_xor_checked,
        expected_xor_checked: expected_xor_total,
        screen_pass: agg_candidates,
        false_alarms: agg_false_alarms,
        witnesses: agg_witnesses.clone(),
        duration_secs: duration.as_secs_f64(),
        rate_per_sec: rate,
        generated_at: chrono::Utc::now().to_rfc3339(),
    };
    let report_path = cli.output_dir.join(format!(
        "validate_report_k{}_{}_{}.json",
        k, cli.start, cli.end
    ));
    let json = serde_json::to_string_pretty(&report)?;
    let tmp = report_path.with_extension("tmp");
    std::fs::write(&tmp, &json)?;
    std::fs::rename(&tmp, &report_path)?;

    println!();
    println!("═══════════════════════════════════════════════════════════════════");
    println!("  VALIDATION COMPLETE");
    println!("═══════════════════════════════════════════════════════════════════");
    println!("  Target k:          {}", k);
    println!("  Range:             [{}, {})", cli.start, cli.end);
    println!("  Report:            {:?}", report_path);
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
        println!(
            "  The range [{}, {}) is clear of k={} witnesses.",
            cli.start, cli.end, k
        );
        println!("  This conclusion relies on Corollary 6 (algorithmic completeness)");
        println!("  plus the standard trust assumptions for software/hardware.");
        println!("  See docs/trust.md for the exact scope and assumptions.");
        println!();
        println!("  Combined with the Governor Set search results, this");
        println!("  constitutes a computational certificate of minimality for");
        println!("  the known minimum witness.");
    } else {
        println!("  CONFIRMED WITNESSES:");
        for w in &agg_witnesses {
            println!("    k={}, n={}", k, w);
        }
    }

    println!("═══════════════════════════════════════════════════════════════════");

    Ok(())
}
