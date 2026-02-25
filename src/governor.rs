//! Governor Set membership checking.
//!
//! The Governor Set G is defined as G = {n : n | C(2n, n)}.
//! A number n is in G if and only if for all primes p dividing n:
//!     v_p(n) ≤ v_p(C(2n, n))
//!
//! We compute v_p(C(2n, n)) using the Legendre formula without
//! computing the actual central binomial coefficient.
//!
//! ## Optimizations
//!
//! - `is_governor_fast()`: Integrates p-adic check directly into trial division
//!   for early exit. Also includes the √(2n) barrier: if the remaining cofactor
//!   after trial division is an odd prime > √(2n), immediately return false.

use crate::factor::Factorization;
use crate::int_math::{isqrt_2n_u64, isqrt_u64};
use crate::sieve::PrimeSieve;

/// Checker for Governor Set membership with cached prime sieve.
#[derive(Debug, Clone)]
pub struct GovernorChecker {
    sieve: PrimeSieve,
}

impl GovernorChecker {
    /// Create a new checker that can handle numbers up to `max_n`.
    pub fn new(max_n: u64) -> Self {
        Self {
            sieve: PrimeSieve::for_range(max_n),
        }
    }

    /// Create a checker from an existing prime sieve.
    pub fn with_sieve(sieve: PrimeSieve) -> Self {
        Self { sieve }
    }

    /// Get the underlying prime sieve.
    #[inline]
    pub fn sieve(&self) -> &PrimeSieve {
        &self.sieve
    }

    /// Check if `n` is in the Governor Set (original method).
    ///
    /// n ∈ G iff n | C(2n, n)
    /// iff for all primes p: v_p(n) ≤ v_p(C(2n, n))
    #[inline]
    pub fn is_governor(&self, n: u64) -> bool {
        if n <= 1 {
            return n == 1;
        }

        let factors = Factorization::of(n, &self.sieve);

        for (p, exp) in factors.iter() {
            let supply = vp_central_binom(n, p);
            if (exp as u64) > supply {
                return false;
            }
        }

        true
    }

    /// Fast governor check with early-exit during factorization.
    ///
    /// Key optimizations over `is_governor()`:
    /// 1. Checks each prime factor's p-adic condition AS it's found during
    ///    trial division, bailing out immediately on failure.
    /// 2. After trial division, if the remaining cofactor is a prime > √(2n),
    ///    returns false immediately (Sanna's √(2n) barrier).
    /// 3. Avoids allocating a SmallVec for the full factorization.
    #[inline]
    pub fn is_governor_fast(&self, n: u64) -> bool {
        if n <= 1 {
            return n == 1;
        }

        let mut remaining = n;
        let mut sqrt_remaining = isqrt_u64(remaining);

        for &p in self.sieve.primes() {
            if p > sqrt_remaining {
                break;
            }

            if remaining % p == 0 {
                let mut exp = 0u32;
                while remaining % p == 0 {
                    exp += 1;
                    remaining /= p;
                }

                // Check supply immediately — early exit on failure
                // Use Kummer's theorem for all primes
                let supply = if p == 2 {
                    vp_central_binom_p2(n)
                } else if p == 3 {
                    vp_central_binom_p3(n)
                } else if p == 5 {
                    vp_central_binom_p5(n)
                } else {
                    vp_central_binom_kummer_fast(n, p)
                };
                if (exp as u64) > supply {
                    return false;
                }

                if remaining == 1 {
                    return true; // All factors checked and passed
                }
                sqrt_remaining = isqrt_u64(remaining);
            }
        }

        // Remaining factor > sqrt(n) is prime (single prime factor)
        if remaining > 1 {
            // √(2n) barrier: if this prime > √(2n), then v_p(C(2n,n)) = 0
            // but v_p(n) = 1, so n ∉ G. Skip the Legendre computation entirely.
            let sqrt_2n = isqrt_2n_u64(n);
            if remaining > sqrt_2n {
                return false;
            }

            // For primes ≤ √(2n), do the full check
            let supply = vp_central_binom_kummer_fast(n, remaining);
            if 1 > supply {
                return false;
            }
        }

        true
    }

    /// Check if all integers in [n-k, n] are in the Governor Set.
    ///
    /// This is equivalent to checking if the run starting at n-k has length ≥ k+1.
    #[inline]
    pub fn is_governor_run(&self, n: u64, k: u32) -> bool {
        (0..=k).all(|i| self.is_governor(n - i as u64))
    }

    /// Find the length of the Governor Set run ending at `n`.
    ///
    /// Returns 0 if n is not in G, otherwise returns the run length.
    pub fn run_length_ending_at(&self, n: u64) -> usize {
        if !self.is_governor(n) {
            return 0;
        }

        let mut len = 1;
        let mut current = n;

        while current > 1 {
            current -= 1;
            if self.is_governor(current) {
                len += 1;
            } else {
                break;
            }
        }

        len
    }
}

/// Compute v_p(n!) using Legendre's formula.
///
/// v_p(n!) = floor(n/p) + floor(n/p²) + floor(n/p³) + ...
///
/// This avoids computing n! directly, which would overflow for large n.
#[inline]
pub fn vp_factorial(n: u64, p: u64) -> u64 {
    debug_assert!(p >= 2, "p must be at least 2");

    let mut v = 0u64;
    let mut power = p;

    while power <= n {
        v += n / power;

        // Prevent overflow: if power * p > n, we're done
        if power > n / p {
            break;
        }
        power *= p;
    }

    v
}

/// Compute v_p(C(2n, n)) = v_p((2n)!) - 2*v_p(n!)
///
/// This is the p-adic "supply" from the central binomial coefficient.
#[inline]
pub fn vp_central_binom(n: u64, p: u64) -> u64 {
    let two_n = n
        .checked_mul(2)
        .expect("vp_central_binom requires 2n to fit in u64");
    vp_factorial(two_n, p) - 2 * vp_factorial(n, p)
}

/// Fast v_2(C(2n, n)) using Kummer's theorem: equals popcount(n).
///
/// Kummer's theorem: v_p(C(m+n, m)) = number of carries when adding m, n in base p.
/// For C(2n, n) with p=2: carry_out at each bit position equals the bit value of n,
/// so total carries = number of set bits = popcount(n).
///
/// This replaces ~34 loop iterations (two Legendre calls) with a single
/// hardware POPCNT instruction (~1 cycle on modern x86).
#[inline(always)]
pub fn vp_central_binom_p2(n: u64) -> u64 {
    n.count_ones() as u64
}

/// Fast v_3(C(2n, n)) using Kummer's theorem.
///
/// Single loop over base-3 digits of n, counting carries when adding n+n.
/// Since max digit is 2, max sum = 2*2+1 = 5, so carry is always 0 or 1.
/// Uses comparison `sum >= 3` instead of division for the carry.
///
/// At n=10B: ~21 iterations vs Legendre's ~42 (two vp_factorial calls × ~21 each).
#[inline(always)]
pub fn vp_central_binom_p3(n: u64) -> u64 {
    let mut carries = 0u64;
    let mut n_remaining = n;
    let mut carry = 0u32;
    while n_remaining > 0 {
        let digit = (n_remaining % 3) as u32;
        let sum = 2 * digit + carry;
        carry = (sum >= 3) as u32;
        carries += carry as u64;
        n_remaining /= 3;
    }
    carries
}

/// Fast v_5(C(2n, n)) using Kummer's theorem.
///
/// Single loop over base-5 digits of n, counting carries when adding n+n.
/// Since max digit is 4, max sum = 2*4+1 = 9, so carry is always 0 or 1.
/// Uses comparison `sum >= 5` instead of division for the carry.
///
/// At n=10B: ~14 iterations vs Legendre's ~28 (two vp_factorial calls × ~14 each).
#[inline(always)]
pub fn vp_central_binom_p5(n: u64) -> u64 {
    let mut carries = 0u64;
    let mut n_remaining = n;
    let mut carry = 0u32;
    while n_remaining > 0 {
        let digit = (n_remaining % 5) as u32;
        let sum = 2 * digit + carry;
        carry = (sum >= 5) as u32;
        carries += carry as u64;
        n_remaining /= 5;
    }
    carries
}

/// Fast v_p(C(2n, n)) using Kummer's theorem for any prime p.
///
/// Single loop counting carries when adding n+n in base p.
/// For any prime p ≥ 2, the carry is always 0 or 1 since
/// max sum = 2(p-1)+1 = 2p-1, and (2p-1)/p = 1.
///
/// ~2x fewer iterations than Legendre (one loop vs two vp_factorial calls).
/// Note: uses runtime division by p (not compile-time constant), so slightly
/// slower per-iteration than the specialized p=3,5 functions.
#[inline(always)]
pub fn vp_central_binom_kummer_fast(n: u64, p: u64) -> u64 {
    debug_assert!(p >= 2, "p must be at least 2");
    let mut carries = 0u64;
    let mut n_remaining = n;
    let mut carry = 0u64;
    // Avoid overflow in `2*digit + carry` by using threshold comparisons.
    // carry_out = 1  iff  2*digit + carry_in ≥ p.
    let half = p / 2;
    let half_up = half + (p & 1); // ceil(p/2)
    while n_remaining > 0 {
        let digit = n_remaining % p;
        let threshold = if carry == 0 { half_up } else { half };
        carry = (digit >= threshold) as u64;
        carries += carry;
        n_remaining /= p;
    }
    carries
}

/// Alternative computation using Kummer's theorem.
///
/// v_p(C(m+n, m)) = (number of carries when adding m and n in base p)
///
/// For C(2n, n), this is carries when adding n + n in base p.
#[allow(dead_code)]
pub fn vp_central_binom_kummer(n: u64, p: u64) -> u64 {
    let mut carries = 0u64;
    let mut n_remaining = n;
    let mut carry = 0u64;

    while n_remaining > 0 || carry > 0 {
        let digit = n_remaining % p;
        let sum = digit + digit + carry;
        carry = sum / p;
        carries += carry;
        n_remaining /= p;
    }

    carries
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg_next(state: &mut u64) -> u64 {
        // Deterministic pseudo-random generator (no external deps).
        // Constants from PCG-style LCGs; good enough for test sampling.
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state
    }

    #[test]
    fn test_vp_factorial() {
        // v_2(10!) = 8 (10/2 + 10/4 + 10/8 = 5 + 2 + 1)
        assert_eq!(vp_factorial(10, 2), 8);

        // v_5(100!) = 24 (100/5 + 100/25 = 20 + 4)
        assert_eq!(vp_factorial(100, 5), 24);

        // v_2(8!) = 7
        assert_eq!(vp_factorial(8, 2), 7);
    }

    #[test]
    fn test_vp_central_binom() {
        // v_2(C(4,2)) = v_2(6) = 1
        // C(4,2) = 6 = 2 * 3
        assert_eq!(vp_central_binom(2, 2), 1);

        // v_3(C(6,3)) = v_3(20) = 0
        // C(6,3) = 20 = 4 * 5
        assert_eq!(vp_central_binom(3, 3), 0);
    }

    #[test]
    fn test_vp_central_binom_p2_matches_legendre() {
        // Verify popcount formula matches Legendre for p=2 across many values
        for n in 1..=100_000u64 {
            let legendre = vp_central_binom(n, 2);
            let popcount = vp_central_binom_p2(n);
            assert_eq!(
                legendre, popcount,
                "p=2 mismatch at n={}: Legendre={}, popcount={}",
                n, legendre, popcount
            );
        }
    }

    #[test]
    fn test_vp_central_binom_p2_large_numbers() {
        // Verify at scales relevant to actual search (billions)
        for &n in &[
            10_000_000_000u64,
            10_000_000_001,
            339_949_252,      // k=8 witness
            17_609_764_993,   // k=9 witness
            u64::MAX / 2 - 1, // near overflow boundary
        ] {
            let legendre = vp_central_binom(n, 2);
            let popcount = vp_central_binom_p2(n);
            assert_eq!(
                legendre, popcount,
                "p=2 mismatch at n={}: Legendre={}, popcount={}",
                n, legendre, popcount
            );
        }
    }

    #[test]
    fn test_kummer_equivalence() {
        // Verify Legendre and Kummer give same results
        for n in [10, 100, 1000, 10000] {
            for p in [2, 3, 5, 7, 11] {
                assert_eq!(
                    vp_central_binom(n, p),
                    vp_central_binom_kummer(n, p),
                    "Mismatch for n={}, p={}",
                    n,
                    p
                );
            }
        }
    }

    #[test]
    fn test_vp_central_binom_p3_matches_legendre() {
        for n in 1..=100_000u64 {
            let legendre = vp_central_binom(n, 3);
            let kummer = vp_central_binom_p3(n);
            assert_eq!(
                legendre, kummer,
                "p=3 mismatch at n={}: Legendre={}, Kummer={}",
                n, legendre, kummer
            );
        }
    }

    #[test]
    fn test_vp_central_binom_p5_matches_legendre() {
        for n in 1..=100_000u64 {
            let legendre = vp_central_binom(n, 5);
            let kummer = vp_central_binom_p5(n);
            assert_eq!(
                legendre, kummer,
                "p=5 mismatch at n={}: Legendre={}, Kummer={}",
                n, legendre, kummer
            );
        }
    }

    #[test]
    fn test_vp_central_binom_p3_p5_large_numbers() {
        for &n in &[
            10_000_000_000u64,
            10_000_000_001,
            339_949_252,
            17_609_764_993,
            u64::MAX / 2 - 1,
        ] {
            let l3 = vp_central_binom(n, 3);
            let k3 = vp_central_binom_p3(n);
            assert_eq!(
                l3, k3,
                "p=3 mismatch at n={}: Legendre={}, Kummer={}",
                n, l3, k3
            );

            let l5 = vp_central_binom(n, 5);
            let k5 = vp_central_binom_p5(n);
            assert_eq!(
                l5, k5,
                "p=5 mismatch at n={}: Legendre={}, Kummer={}",
                n, l5, k5
            );
        }
    }

    #[test]
    fn test_kummer_fast_matches_legendre_all_small_primes() {
        // Verify the generic Kummer function matches Legendre for primes 7..23
        let primes = [7u64, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        for &p in &primes {
            for n in 1..=10_000u64 {
                let legendre = vp_central_binom(n, p);
                let kummer = vp_central_binom_kummer_fast(n, p);
                assert_eq!(
                    legendre, kummer,
                    "p={} mismatch at n={}: Legendre={}, Kummer={}",
                    p, n, legendre, kummer
                );
            }
        }
    }

    #[test]
    fn test_kummer_fast_large_numbers() {
        let primes = [7u64, 11, 13, 17, 19, 97, 997];
        for &p in &primes {
            for &n in &[10_000_000_000u64, 339_949_252, 17_609_764_993] {
                let legendre = vp_central_binom(n, p);
                let kummer = vp_central_binom_kummer_fast(n, p);
                assert_eq!(
                    legendre, kummer,
                    "p={} mismatch at n={}: Legendre={}, Kummer={}",
                    p, n, legendre, kummer
                );
            }
        }
    }

    #[test]
    fn test_kummer_fast_matches_legendre_random_large() {
        // Random sampling at the 10^13 scale (k=13 search range).
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 97, 997];
        let max_n = 25_000_000_000_000u64;

        let mut state = 0x3d2c_ba98_7c01_4e55u64;
        for _ in 0..2000 {
            let n = (lcg_next(&mut state) % max_n).saturating_add(1);
            let p = primes[(lcg_next(&mut state) as usize) % primes.len()];

            let legendre = vp_central_binom(n, p);
            let kummer = match p {
                2 => vp_central_binom_p2(n),
                3 => vp_central_binom_p3(n),
                5 => vp_central_binom_p5(n),
                _ => vp_central_binom_kummer_fast(n, p),
            };

            assert_eq!(
                legendre, kummer,
                "Random mismatch at n={}, p={}: Legendre={}, Kummer={}",
                n, p, legendre, kummer
            );
        }
    }

    #[test]
    fn test_governor_membership() {
        let checker = GovernorChecker::new(1000);

        // n ∈ G iff n | C(2n, n)
        let members = [1, 2, 6, 20];
        for &n in &members {
            assert!(checker.is_governor(n), "{} should be in Governor Set", n);
            assert!(
                checker.is_governor_fast(n),
                "{} should be in Governor Set (fast)",
                n
            );
        }

        let non_members = [3, 4, 5, 7, 8, 9, 10];
        for &n in &non_members {
            assert!(
                !checker.is_governor(n),
                "{} should NOT be in Governor Set",
                n
            );
            assert!(
                !checker.is_governor_fast(n),
                "{} should NOT be in Governor Set (fast)",
                n
            );
        }
    }

    #[test]
    fn test_fast_agrees_with_original() {
        // Verify is_governor_fast agrees with is_governor for a range of numbers
        let checker = GovernorChecker::new(100_000);

        for n in 1..=50_000u64 {
            let orig = checker.is_governor(n);
            let fast = checker.is_governor_fast(n);
            assert_eq!(
                orig, fast,
                "Mismatch at n={}: original={}, fast={}",
                n, orig, fast
            );
        }
    }

    #[test]
    fn test_governor_run() {
        let checker = GovernorChecker::new(1_000_000_000);

        // Known k=8 witness: n = 339,949,252
        assert!(checker.is_governor_run(339_949_252, 8));
    }

    #[test]
    fn test_verify_known_runs_of_9() {
        use crate::KNOWN_RUNS_OF_9;

        let max_n = *KNOWN_RUNS_OF_9.iter().max().unwrap_or(&1);
        let checker = GovernorChecker::new(max_n + 100);

        println!("\nVerifying run lengths for all known runs of 9:");
        for &end_pos in KNOWN_RUNS_OF_9 {
            let actual_len = checker.run_length_ending_at(end_pos);
            println!("  n={}: actual run length = {}", end_pos, actual_len);

            assert!(
                actual_len >= 9,
                "Run ending at {} should be at least 9, got {}",
                end_pos,
                actual_len
            );

            if actual_len > 9 {
                println!(
                    "  *** FOUND EXTENSION: run at {} is actually {} ***",
                    end_pos, actual_len
                );
            }

            let start_pos = end_pos - (actual_len as u64 - 1);
            if start_pos > 1 {
                let before_run = !checker.is_governor(start_pos - 1);
                assert!(
                    before_run,
                    "Position {} should NOT be in Governor Set (start of run)",
                    start_pos - 1
                );
            }
        }
    }
}
