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

use crate::governor::vp_central_binom_dispatch;
use crate::int_math::isqrt_u128;
use crate::sieve::{build_prime_data, PrimeData};

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

/// Sub-block size for L1 cache residency (32K elements × 8 bytes = 256KB).
const BLOCK_SIZE: usize = 32_768;

// ---------------------------------------------------------------------------
// Const-generic prime factor stripping (compile-time constant multiplication)
// ---------------------------------------------------------------------------

/// Strip all factors of compile-time-constant prime P from `remaining`,
/// returning the exponent. Returns 0 if P does not divide `*remaining`.
#[inline(always)]
fn strip_prime_const<const P: u64>(remaining: &mut u64) -> u32 {
    const fn inv_const(p: u64) -> u64 {
        let mut inv = p.wrapping_mul(3) ^ 2;
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv
    }
    const fn limit_const(p: u64) -> u64 {
        u64::MAX / p
    }

    let inv = inv_const(P);
    let limit = limit_const(P);

    let q = remaining.wrapping_mul(inv);
    if q > limit {
        return 0;
    }

    *remaining = q;
    let mut exp = 1u32;
    loop {
        let q2 = remaining.wrapping_mul(inv);
        if q2 > limit {
            break;
        }
        *remaining = q2;
        exp += 1;
    }
    exp
}

/// Strip all factors of runtime prime using precomputed PrimeData.
/// Returns exponent (0 if prime does not divide `*remaining`).
#[inline(always)]
fn strip_prime_dyn(remaining: &mut u64, inv_p: u64, max_quot: u64) -> u32 {
    let q = remaining.wrapping_mul(inv_p);
    if q > max_quot {
        return 0;
    }

    *remaining = q;
    let mut exp = 1u32;
    loop {
        let q2 = remaining.wrapping_mul(inv_p);
        if q2 > max_quot {
            break;
        }
        *remaining = q2;
        exp += 1;
    }
    exp
}

macro_rules! dispatch_strip {
    ($remaining:expr, $pd:expr, [$($prime:literal),* $(,)?]) => {
        match $pd.p as u64 {
            $($prime => strip_prime_const::<$prime>($remaining),)*
            _ => strip_prime_dyn($remaining, $pd.inv_p, $pd.max_quot),
        }
    };
}

// ---------------------------------------------------------------------------
// Bucketed large-prime processing
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct BucketItem {
    pd_idx: u32,
    offset: u32,
}

struct PrimeBucket {
    items: Vec<BucketItem>,
}

impl PrimeBucket {
    fn new() -> Self {
        PrimeBucket {
            items: Vec::with_capacity(256),
        }
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.items.clear();
    }

    #[inline(always)]
    fn push(&mut self, pd_idx: u32, offset: u32) {
        self.items.push(BucketItem { pd_idx, offset });
    }
}

impl FusedBatchResult {
    /// Compute exact governor membership for every number in [lo, lo + len).
    ///
    /// Uses modular-inverse arithmetic for divisibility testing and exact
    /// division (~3 cycles per wrapping_mul vs ~35 cycles for hardware div).
    /// Processes in 32K sub-blocks to keep the `remaining[]` working set in
    /// L1 cache.
    ///
    /// `prime_data` must contain all primes up to at least √(2·(lo+len)).
    pub fn compute(lo: u64, len: usize, prime_data: &[PrimeData]) -> Self {
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

        let total_primes = prime_data.partition_point(|pd| (pd.p as u64) <= sieve_limit);
        let first_odd = if total_primes > 0 && prime_data[0].p == 2 {
            1
        } else {
            0
        };

        let mut remaining = vec![0u64; len];
        let mut is_governor = vec![true; len];
        let mut rejected_vp_fail = 0usize;

        // Phase 1: Initialize remaining[] — strip factor of 2, check v_2 inline.
        for i in 0..len {
            let n = lo + i as u64;
            if n <= 1 {
                is_governor[i] = n == 1;
                remaining[i] = 1;
                continue;
            }
            let tz = n.trailing_zeros();
            remaining[i] = n >> tz;
            if (n.count_ones() as u32) < tz {
                is_governor[i] = false;
                rejected_vp_fail += 1;
            }
        }

        // Phase 2: Compute initial offsets for all odd primes (Barrett reduction).
        let mut prime_offsets = vec![0u32; total_primes];
        for idx in first_odd..total_primes {
            let pd = &prime_data[idx];
            let p = pd.p as u64;
            let first_multiple = if lo == 0 {
                p
            } else {
                let num = lo + p - 1;
                let start_c = (((num as u128).wrapping_mul(pd.magic as u128)) >> 64) as u64
                    >> pd.shift as u64;
                let fm = start_c * p;
                if fm < lo { fm + p } else { fm }
            };
            prime_offsets[idx] = (first_multiple - lo) as u32;
        }

        // Phase 3: Process in sub-blocks for L1 cache residency.
        let first_large = prime_data[..total_primes]
            .partition_point(|pd| (pd.p as usize) <= BLOCK_SIZE);

        // Pre-bucket large primes by target block index.
        let num_blocks = len.div_ceil(BLOCK_SIZE);
        let mut buckets: Vec<PrimeBucket> =
            (0..num_blocks).map(|_| PrimeBucket::new()).collect();
        for pidx in first_large..total_primes {
            let p = prime_data[pidx].p;
            let mut sj = prime_offsets[pidx];
            while (sj as usize) < len {
                let block_idx = (sj as usize) / BLOCK_SIZE;
                let offset = (sj as usize) % BLOCK_SIZE;
                buckets[block_idx].push(pidx as u32, offset as u32);
                sj += p;
            }
        }

        let mut block_start: usize = 0;
        while block_start < len {
            let block_end = (block_start + BLOCK_SIZE).min(len);
            let block_len = (block_end - block_start) as u32;

            // Small primes (p <= BLOCK_SIZE): stride through this block
            // Uses const-generic dispatch for primes 3-97 (compile-time inverse).
            for pidx in first_odd..first_large {
                let pd = &prime_data[pidx];
                let mut sj = prime_offsets[pidx];
                if sj < block_len {
                    while sj < block_len {
                        let abs_idx = block_start + sj as usize;
                        if is_governor[abs_idx] {
                            let exp = dispatch_strip!(
                                &mut remaining[abs_idx],
                                pd,
                                [
                                    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                                    61, 67, 71, 73, 79, 83, 89, 97
                                ]
                            );
                            if exp > 0 {
                                let n = lo + abs_idx as u64;
                                let supply = vp_central_binom_dispatch(n, pd.p as u64);
                                if (exp as u64) > supply {
                                    is_governor[abs_idx] = false;
                                    rejected_vp_fail += 1;
                                }
                            }
                        }
                        sj += pd.p;
                    }
                    prime_offsets[pidx] = sj - block_len;
                } else {
                    prime_offsets[pidx] = sj - block_len;
                }
            }

            // Large primes: process from pre-built buckets with prefetch
            {
                let block_idx = block_start / BLOCK_SIZE;
                let bucket = &buckets[block_idx];
                let pd_ptr = prime_data.as_ptr();
                for (bi, item) in bucket.items.iter().enumerate() {
                    // Prefetch upcoming prime data (8 items ahead)
                    if bi + 8 < bucket.items.len() {
                        let future_idx = bucket.items[bi + 8].pd_idx as usize;
                        let ptr = unsafe { pd_ptr.add(future_idx) as *const u8 };
                        unsafe {
                            #[cfg(target_arch = "aarch64")]
                            std::arch::asm!(
                                "prfm pldl2keep, [{ptr}]",
                                ptr = in(reg) ptr,
                                options(nostack, preserves_flags)
                            );
                            #[cfg(target_arch = "x86_64")]
                            std::arch::x86_64::_mm_prefetch(
                                ptr as *const i8,
                                std::arch::x86_64::_MM_HINT_T1,
                            );
                        }
                    }

                    let abs_idx = block_start + item.offset as usize;
                    if !is_governor[abs_idx] {
                        continue;
                    }
                    let pd = &prime_data[item.pd_idx as usize];
                    let exp = strip_prime_dyn(
                        &mut remaining[abs_idx],
                        pd.inv_p,
                        pd.max_quot,
                    );
                    if exp > 0 {
                        let n = lo + abs_idx as u64;
                        let supply = vp_central_binom_dispatch(n, pd.p as u64);
                        if (exp as u64) > supply {
                            is_governor[abs_idx] = false;
                            rejected_vp_fail += 1;
                        }
                    }
                }
            }

            block_start = block_end;
        }

        // === Post-sieve: check for large prime factors ===
        let mut rejected_odd_prime = 0usize;
        let mut rejected_large_pf = 0usize;

        for i in 0..len {
            if is_governor[i] && remaining[i] > 1 {
                is_governor[i] = false;
                let n = lo + i as u64;
                let n_odd = n >> n.trailing_zeros();
                if remaining[i] == n_odd {
                    rejected_odd_prime += 1;
                } else {
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

    /// Legacy entry point: builds PrimeData on the fly from raw primes.
    ///
    /// Prefer `compute()` with pre-built PrimeData for production use.
    pub fn compute_from_primes(lo: u64, len: usize, base_primes: &[u64]) -> Self {
        let prime_data = build_prime_data(base_primes);
        Self::compute(lo, len, &prime_data)
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
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

        let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());

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

    #[test]
    fn test_fused_modinv_matches_original_wide_range() {
        // Test across multiple ranges to catch edge cases
        for &(lo, len) in &[
            (2u64, 500),
            (10_000, 5_000),
            (1_000_000, 10_000),
            (339_949_240, 20),        // k=8 witness region
            (1_019_547_840, 20),      // k=9 witness region
            (24_999_999_999_900, 20), // 25T scale
        ] {
            let hi = lo + len as u64;
            let max_n: u128 = (lo as u128) + (len as u128) - 1;
            let sieve_limit = (isqrt_u128(2u128 * max_n) as u64).saturating_add(100);
            let sieve = PrimeSieve::new(sieve_limit);

            let fused = FusedBatchResult::compute_from_primes(lo, len, sieve.primes());
            let checker = GovernorChecker::with_sieve(sieve);

            for i in 0..len {
                let n = lo + i as u64;
                let expected = checker.is_governor(n);
                assert_eq!(
                    fused.is_governor[i], expected,
                    "Fused disagrees at n={} (range [{}, {})): fused={}, expected={}",
                    n, lo, hi, fused.is_governor[i], expected
                );
            }
        }
    }
}
