//! Batch pre-filtering and fused sieve+governor computation.
//!
//! Implements optimizations based on Sanna/Pomerance and Ford-Konyagin:
//!
//! 1. **Odd prime exclusion**: Odd primes are never Governors.
//!    For odd prime p: v_p(C(2p,p)) = 0 but v_p(p) = 1, so p ∤ C(2p,p).
//!
//! 2. **Large prime factor barrier (√(2n) rule)**: If n has an odd prime
//!    factor p > √(2n), then v_p(C(2n,n)) = 0 but v_p(n) ≥ 1, so n ∉ G.
//!
//! 3. **Fused sieve+governor**: The segmented sieve divides out all small
//!    primes AND checks v_p(C(2n,n)) inline during the same pass. This
//!    eliminates redundant trial division — the sieve visits only multiples
//!    of each prime (amortized), while trial division tests every prime
//!    against every number.

use crate::governor::{
    vp_central_binom_kummer_fast, vp_central_binom_p2, vp_central_binom_p3, vp_central_binom_p5,
};
use crate::int_math::isqrt_u128;

/// Fused batch result: exact governor membership computed in a single sieve pass.
///
/// The segmented sieve divides out primes AND checks the p-adic condition
/// simultaneously. After the sieve, `is_governor[i]` is the exact answer —
/// no further factorization or governor checking needed.
pub struct FusedBatchResult {
    /// For each position i: true iff (lo + i) IS a Governor.
    pub is_governor: Vec<bool>,

    /// Number of governors found in this batch
    pub governor_count: usize,

    /// Rejected because n is an odd prime (cofactor == n after sieve)
    pub rejected_odd_prime: usize,

    /// Rejected because n has a large prime factor > √(2n)
    pub rejected_large_pf: usize,

    /// Rejected because v_p(n) > v_p(C(2n,n)) for some prime p
    pub rejected_vp_fail: usize,
}

impl FusedBatchResult {
    /// Compute exact governor membership for every number in [lo, lo + len).
    ///
    /// Algorithm:
    /// 1. Initialize remaining\[i\] = lo + i
    /// 2. For each prime p ≤ sieve_limit:
    ///    - Visit multiples of p in the batch (amortized O(batch/p) per prime)
    ///    - Divide out all factors of p, recording exponent
    ///    - Check v_p(C(2n,n)) ≥ exponent inline; reject on failure
    ///    - Skip already-rejected numbers entirely (saves division work)
    /// 3. After sieve: remaining\[i\] > 1 means large prime factor → reject
    ///
    /// `base_primes` must contain all primes up to at least √(2·(lo+len)).
    pub fn compute(lo: u64, len: usize, base_primes: &[u64]) -> Self {
        if len == 0 {
            return Self {
                is_governor: vec![],
                governor_count: 0,
                rejected_odd_prime: 0,
                rejected_large_pf: 0,
                rejected_vp_fail: 0,
            };
        }

        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(1);

        // remaining[i] tracks the cofactor of (lo+i) after dividing out small primes.
        // For non-rejected numbers, this is maintained accurately.
        // For rejected numbers, we skip divisions (remaining is invalid but unused).
        let mut remaining = Vec::with_capacity(len);
        for i in 0..len {
            remaining.push(lo + i as u64);
        }
        let mut is_governor = vec![true; len];
        let mut rejected_vp_fail = 0usize;

        // Handle edge cases: n = 0 or 1
        for i in 0..len {
            let n = lo + i as u64;
            if n <= 1 {
                is_governor[i] = n == 1;
                remaining[i] = 1; // prevent sieve from processing
            }
        }

        // === Fused sieve: divide out primes AND check v_p inline ===
        for &p in base_primes {
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

            let mut idx = (first_multiple - lo) as usize;
            while idx < len {
                // Skip already-rejected numbers entirely.
                // Their remaining[] is invalid (not fully reduced), but we
                // never read it again — they're already marked not-governor.
                if !is_governor[idx] {
                    idx += p as usize;
                    continue;
                }

                if remaining[idx] % p == 0 {
                    // Count and divide out all factors of p
                    let mut exp = 0u32;
                    while remaining[idx] % p == 0 {
                        exp += 1;
                        remaining[idx] /= p;
                    }

                    // Check v_p(C(2n,n)) ≥ v_p(n) inline
                    // Use Kummer's theorem: single loop counting carries
                    // when adding n+n in base p (~2x fewer iterations than
                    // Legendre). p=2 uses hardware POPCNT; p=3,5 use
                    // specialized functions with compile-time constant divisor.
                    let n = lo + idx as u64;
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
                        is_governor[idx] = false;
                        rejected_vp_fail += 1;
                    }
                }

                idx += p as usize;
            }
        }

        // === Post-sieve: check for large prime factors ===
        // For numbers still marked as governor, if remaining > 1, the cofactor
        // is a single prime > sieve_limit > √(2n) → not a governor.
        let mut rejected_odd_prime = 0usize;
        let mut rejected_large_pf = 0usize;

        for i in 0..len {
            if is_governor[i] && remaining[i] > 1 {
                is_governor[i] = false;
                if remaining[i] == lo + i as u64 {
                    // No small prime divided n → n is prime
                    rejected_odd_prime += 1;
                } else {
                    // n is composite with a large prime factor
                    rejected_large_pf += 1;
                }
            }
        }

        let governor_count = is_governor.iter().filter(|&&x| x).count();

        Self {
            is_governor,
            governor_count,
            rejected_odd_prime,
            rejected_large_pf,
            rejected_vp_fail,
        }
    }
}

/// Pre-filter results for a batch of numbers [lo, lo + len).
/// (Legacy two-step approach — kept for benchmarking comparison.)
pub struct BatchFilter {
    /// For each position i: true if (lo + i) *might* be a Governor.
    /// false guarantees (lo + i) is NOT a Governor.
    pub can_be_governor: Vec<bool>,

    /// Number of candidates that passed the filter
    pub passed: usize,

    /// Number of candidates rejected
    pub rejected: usize,

    /// Number rejected because they're odd primes
    pub rejected_odd_prime: usize,

    /// Number rejected because of large prime factor > √(2n)
    pub rejected_large_pf: usize,
}

impl BatchFilter {
    /// Run pre-filter on the range [lo, lo + len).
    pub fn compute(lo: u64, len: usize, base_primes: &[u64]) -> Self {
        if len == 0 {
            return Self {
                can_be_governor: vec![],
                passed: 0,
                rejected: 0,
                rejected_odd_prime: 0,
                rejected_large_pf: 0,
            };
        }

        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(1);

        let mut remaining = Vec::with_capacity(len);
        for i in 0..len {
            remaining.push(lo + i as u64);
        }

        for &p in base_primes {
            if p > sieve_limit {
                break;
            }

            let first_multiple = if lo == 0 {
                p
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

        let mut can_be_governor = vec![true; len];
        let mut rejected_odd_prime = 0usize;
        let mut rejected_large_pf = 0usize;

        for i in 0..len {
            let n = lo + i as u64;

            if n <= 1 {
                can_be_governor[i] = n == 1;
                continue;
            }

            if remaining[i] > 1 {
                can_be_governor[i] = false;
                if remaining[i] == n {
                    rejected_odd_prime += 1;
                } else {
                    rejected_large_pf += 1;
                }
            }
        }

        let rejected = rejected_odd_prime + rejected_large_pf;
        let passed = len - rejected;

        Self {
            can_be_governor,
            passed,
            rejected,
            rejected_odd_prime,
            rejected_large_pf,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::governor::GovernorChecker;
    use crate::sieve::PrimeSieve;

    #[test]
    fn test_odd_primes_rejected() {
        let lo = 100u64;
        let len = 100usize;
        let sieve = PrimeSieve::new(200);
        let filter = BatchFilter::compute(lo, len, sieve.primes());

        for &p in &[
            101u64, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
            191, 193, 197, 199,
        ] {
            if p >= lo && p < lo + len as u64 {
                let idx = (p - lo) as usize;
                assert!(
                    !filter.can_be_governor[idx],
                    "Odd prime {} should be rejected by prefilter",
                    p
                );
            }
        }
    }

    #[test]
    fn test_smooth_numbers_pass() {
        let lo = 100u64;
        let len = 100usize;
        let sieve = PrimeSieve::new(200);
        let filter = BatchFilter::compute(lo, len, sieve.primes());

        let idx = (120 - lo) as usize;
        assert!(filter.can_be_governor[idx], "120 should pass prefilter");

        let idx = (128 - lo) as usize;
        assert!(filter.can_be_governor[idx], "128 should pass prefilter");
    }

    #[test]
    fn test_large_pf_rejected() {
        let lo = 10000u64;
        let len = 100usize;
        let sieve = PrimeSieve::new(200);
        let filter = BatchFilter::compute(lo, len, sieve.primes());

        let idx = (10007 - lo) as usize;
        assert!(
            !filter.can_be_governor[idx],
            "Prime 10007 should be rejected"
        );
    }

    #[test]
    fn test_rejection_rate_large_range() {
        let lo = 1_000_000u64;
        let len = 100_000usize;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(1000);
        let sieve = PrimeSieve::new(sieve_limit);
        let filter = BatchFilter::compute(lo, len, sieve.primes());

        let rate = (filter.rejected as f64 / (filter.passed + filter.rejected) as f64) * 100.0;
        assert!(
            rate > 40.0,
            "Expected rejection rate > 40%, got {:.1}%",
            rate
        );
    }

    #[test]
    fn test_prefilter_agrees_with_governor() {
        use crate::governor::GovernorChecker;

        let lo = 10000u64;
        let len = 1000usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);
        let filter = BatchFilter::compute(lo, len, sieve.primes());

        let checker = GovernorChecker::new(hi);

        for i in 0..len {
            let n = lo + i as u64;
            if !filter.can_be_governor[i] {
                assert!(
                    !checker.is_governor(n),
                    "Prefilter rejected n={} but it IS a Governor!",
                    n
                );
            }
        }
    }

    // ==================== Fused batch tests ====================

    #[test]
    fn test_fused_exact_match_small_range() {
        // Verify fused results match original is_governor for a small range
        let lo = 2u64;
        let len = 200usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());
        let checker = GovernorChecker::new(hi);

        for i in 0..len {
            let n = lo + i as u64;
            let expected = checker.is_governor(n);
            assert_eq!(
                fused.is_governor[i], expected,
                "Fused disagrees with is_governor at n={}: fused={}, expected={}",
                n, fused.is_governor[i], expected
            );
        }
    }

    #[test]
    fn test_fused_exact_match_medium_range() {
        // Test at a larger range to catch edge cases
        let lo = 10_000u64;
        let len = 5_000usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());
        let checker = GovernorChecker::new(hi);

        let mut mismatches = 0;
        for i in 0..len {
            let n = lo + i as u64;
            let expected = checker.is_governor(n);
            if fused.is_governor[i] != expected {
                mismatches += 1;
                eprintln!(
                    "MISMATCH at n={}: fused={}, expected={}",
                    n, fused.is_governor[i], expected
                );
            }
        }
        assert_eq!(
            mismatches, 0,
            "Found {} mismatches in [{}..{})",
            mismatches, lo, hi
        );
    }

    #[test]
    fn test_fused_exact_match_large_numbers() {
        // Test at the scale of actual k=8 witness
        let lo = 339_949_240u64;
        let len = 20usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());
        let checker = GovernorChecker::new(hi);

        // The k=8 witness run: 339,949,244 through 339,949,252 should all be governors
        for i in 0..len {
            let n = lo + i as u64;
            let expected = checker.is_governor(n);
            assert_eq!(
                fused.is_governor[i], expected,
                "Fused disagrees at n={}: fused={}, expected={}",
                n, fused.is_governor[i], expected
            );
        }
    }

    #[test]
    fn test_fused_exact_match_25t_scale() {
        // Test at the top end of the k=13 exhaustive search range (~25T).
        let lo = 24_999_999_999_900u64;
        let len = 20usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());
        let checker = GovernorChecker::with_sieve(sieve);

        for i in 0..len {
            let n = lo + i as u64;
            let expected = checker.is_governor(n);
            assert_eq!(
                fused.is_governor[i], expected,
                "Fused disagrees at n={}: fused={}, expected={}",
                n, fused.is_governor[i], expected
            );
        }

        // Also sanity-check that the direct checker could factor the entire range.
        assert!(
            checker.sieve().can_factor(hi),
            "Test sieve should be sufficient to factor n up to {}",
            hi
        );
    }

    #[test]
    fn test_fused_governor_count_matches() {
        let lo = 100_000u64;
        let len = 10_000usize;
        let hi = lo + len as u64;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());
        let checker = GovernorChecker::new(hi);

        let expected_count = (lo..hi).filter(|&n| checker.is_governor(n)).count();
        assert_eq!(
            fused.governor_count, expected_count,
            "Governor count mismatch: fused={}, expected={}",
            fused.governor_count, expected_count
        );
    }

    #[test]
    fn test_fused_stats_add_up() {
        let lo = 50_000u64;
        let len = 10_000usize;
        let max_n: u128 = (lo as u128) + (len as u128) - 1;
        let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
        let sieve = PrimeSieve::new(sieve_limit);

        let fused = FusedBatchResult::compute(lo, len, sieve.primes());

        // governor_count + rejected should equal len (minus edge cases)
        let total_rejected =
            fused.rejected_odd_prime + fused.rejected_large_pf + fused.rejected_vp_fail;
        // Some numbers at the very start might be edge cases (n<=1), but lo=50000 so no edge cases
        assert_eq!(
            fused.governor_count + total_rejected,
            len,
            "Stats don't add up: governors={} + rejected({} + {} + {})={} != len={}",
            fused.governor_count,
            fused.rejected_odd_prime,
            fused.rejected_large_pf,
            fused.rejected_vp_fail,
            total_rejected,
            len
        );
    }
}
