//! Prime sieve implementation using the Sieve of Eratosthenes.
//!
//! For searching up to n ~ 10^13–10^14, we need primes up to sqrt(n) ~ a few million.
//! This sieve is computed once at startup and shared across all workers.

use crate::int_math::isqrt_u64;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// PrimeData — precomputed per-prime data for fast modular arithmetic
// ---------------------------------------------------------------------------

/// Precomputed per-prime data for fast modular arithmetic.
///
/// Hot fields (inv_p, max_quot) are first so they share a cache line
/// during the inner strip loop.
#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct PrimeData {
    /// Modular inverse of p mod 2^64 (for odd p). Used for:
    /// - Divisibility test: `n.wrapping_mul(inv_p) <= max_quot` iff `p | n`
    /// - Exact division: `n / p == n.wrapping_mul(inv_p)` when `p | n`
    pub inv_p: u64,
    /// `u64::MAX / p` — divisibility threshold for the inverse test.
    pub max_quot: u64,
    /// Barrett magic constant: `ceil(2^(64+shift) / p)`. Used once per
    /// prime per batch to compute the first multiple of p in the batch
    /// without hardware division.
    pub magic: u64,
    /// The prime value.
    pub p: u32,
    /// `bit_width(p) - 1` — Barrett shift amount.
    pub shift: u8,
}

/// Compute modular inverse of odd `n` mod 2^64 via Newton's method.
///
/// After 5 Newton iterations starting from a 4-bit-accurate seed,
/// the result is exact mod 2^64.
#[inline(always)]
pub fn mod_inverse_u64(n: u64) -> u64 {
    debug_assert!(n & 1 == 1, "mod_inverse_u64 requires odd input");
    let mut inv = n.wrapping_mul(3) ^ 2; // correct mod 2^4
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv
}

/// Build `PrimeData` array from a slice of primes.
///
/// Precomputes modular inverse, Barrett magic, and divisibility
/// threshold for each prime. Called once at startup.
pub fn build_prime_data(primes: &[u64]) -> Vec<PrimeData> {
    primes
        .iter()
        .map(|&p| {
            let p32 = p as u32;
            let inv_p = if p % 2 != 0 { mod_inverse_u64(p) } else { 0 };
            let shift = (32 - p32.leading_zeros()).saturating_sub(1) as u8;
            let magic = if p > 1 {
                (((1u128 << (64 + shift as u32)) / p as u128) + 1) as u64
            } else {
                0
            };
            PrimeData {
                inv_p,
                max_quot: u64::MAX / p,
                magic,
                p: p32,
                shift,
            }
        })
        .collect()
}

/// A precomputed prime sieve for fast factorization.
///
/// The sieve stores all primes up to `limit` in a contiguous Vec for
/// cache-friendly sequential access during trial division.
#[derive(Debug, Clone)]
pub struct PrimeSieve {
    /// All primes up to the limit, sorted ascending
    primes: Arc<Vec<u64>>,
    /// The upper bound used when creating this sieve
    limit: u64,
}

impl PrimeSieve {
    /// Create a new prime sieve containing all primes up to `limit`.
    ///
    /// # Example
    /// ```
    /// use erdos396::PrimeSieve;
    /// let sieve = PrimeSieve::new(100);
    /// assert_eq!(sieve.primes(), &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]);
    /// ```
    pub fn new(limit: u64) -> Self {
        let primes = sieve_of_eratosthenes(limit as usize);
        Self {
            primes: Arc::new(primes),
            limit,
        }
    }

    /// Create a sieve suitable for factoring numbers up to `max_n`.
    ///
    /// This automatically computes the required limit as sqrt(max_n) + buffer.
    pub fn for_range(max_n: u64) -> Self {
        // sqrt(max_n) + 1000 for safety margin
        let limit = isqrt_u64(max_n).saturating_add(1000);
        Self::new(limit)
    }

    /// Get the slice of all primes in this sieve.
    #[inline]
    pub fn primes(&self) -> &[u64] {
        &self.primes
    }

    /// Get the number of primes in this sieve.
    #[inline]
    pub fn len(&self) -> usize {
        self.primes.len()
    }

    /// Check if the sieve is empty (should never be for limit >= 2).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.primes.is_empty()
    }

    /// Get the limit used when creating this sieve.
    #[inline]
    pub fn limit(&self) -> u64 {
        self.limit
    }

    /// Check if this sieve is sufficient for factoring `n`.
    #[inline]
    pub fn can_factor(&self, n: u64) -> bool {
        let sqrt_n = isqrt_u64(n);
        sqrt_n <= self.limit
    }
}

/// Classic Sieve of Eratosthenes.
///
/// Uses a bit-packed representation for memory efficiency.
fn sieve_of_eratosthenes(limit: usize) -> Vec<u64> {
    if limit < 2 {
        return vec![];
    }

    // is_prime[i] = true if i is prime
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    // Only need to sieve up to sqrt(limit)
    let sqrt_limit = isqrt_u64(limit as u64) as usize;

    for i in 2..=sqrt_limit {
        if is_prime[i] {
            // Mark all multiples of i as composite
            // Start from i*i since smaller multiples were already marked
            let mut j = i * i;
            while j <= limit {
                is_prime[j] = false;
                j += i;
            }
        }
    }

    // Collect primes
    is_prime
        .iter()
        .enumerate()
        .filter_map(|(i, &is_p)| if is_p { Some(i as u64) } else { None })
        .collect()
}

/// Segmented sieve for generating primes in a range [low, high).
///
/// Useful when only primes in a specific range are needed.
#[allow(dead_code)]
pub fn segmented_sieve(low: u64, high: u64, base_primes: &[u64]) -> Vec<u64> {
    if high <= low {
        return vec![];
    }

    let size = (high - low) as usize;
    let mut is_prime = vec![true; size];

    for &p in base_primes {
        if p * p >= high {
            break;
        }

        // Find first multiple of p >= low
        let mut start = low.div_ceil(p) * p;
        if start == p {
            start += p; // Don't mark p itself
        }

        let mut j = (start - low) as usize;
        while j < size {
            is_prime[j] = false;
            j += p as usize;
        }
    }

    is_prime
        .iter()
        .enumerate()
        .filter_map(|(i, &is_p)| if is_p { Some(low + i as u64) } else { None })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_sieve() {
        let sieve = PrimeSieve::new(30);
        assert_eq!(sieve.primes(), &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_prime_count() {
        // π(100) = 25
        let sieve = PrimeSieve::new(100);
        assert_eq!(sieve.len(), 25);

        // π(1000) = 168
        let sieve = PrimeSieve::new(1000);
        assert_eq!(sieve.len(), 168);

        // π(10000) = 1229
        let sieve = PrimeSieve::new(10000);
        assert_eq!(sieve.len(), 1229);

        // π(1_000_000) = 78498
        let sieve = PrimeSieve::new(1_000_000);
        assert_eq!(sieve.len(), 78_498);
    }

    #[test]
    fn test_for_range() {
        let sieve = PrimeSieve::for_range(1_000_000_000);
        assert!(sieve.can_factor(1_000_000_000));
        assert!(sieve.limit() >= 31623); // sqrt(10^9) ≈ 31623
    }

    #[test]
    fn test_mod_inverse_u64() {
        // For odd prime p, p * mod_inverse(p) ≡ 1 (mod 2^64)
        for &p in &[3u64, 5, 7, 11, 13, 97, 997, 7919, 104729] {
            let inv = super::mod_inverse_u64(p);
            assert_eq!(p.wrapping_mul(inv), 1u64, "inverse failed for p={}", p);
        }
    }

    #[test]
    fn test_prime_data_divisibility() {
        // n.wrapping_mul(inv_p) <= max_quot  iff  p | n
        let p = 7u64;
        let inv = super::mod_inverse_u64(p);
        let max_quot = u64::MAX / p;

        // Multiples of 7
        for m in 1..=1000u64 {
            let n = m * p;
            assert!(
                n.wrapping_mul(inv) <= max_quot,
                "false negative: {} should be divisible by {}",
                n,
                p
            );
        }
        // Non-multiples of 7
        for n in 1..=100u64 {
            if n % p != 0 {
                assert!(
                    n.wrapping_mul(inv) > max_quot,
                    "false positive: {} should NOT be divisible by {}",
                    n,
                    p
                );
            }
        }
    }

    #[test]
    fn test_prime_data_exact_division() {
        // When p | n: n / p == n.wrapping_mul(inv_p)
        let p = 13u64;
        let inv = super::mod_inverse_u64(p);
        for m in 1..=10_000u64 {
            let n = m * p;
            let quot = n.wrapping_mul(inv);
            assert_eq!(
                quot, m,
                "exact division failed: {} / {} should be {}, got {}",
                n, p, m, quot
            );
        }
    }

    #[test]
    fn test_build_prime_data_roundtrip() {
        let sieve = PrimeSieve::new(1000);
        let pd = super::build_prime_data(sieve.primes());
        assert_eq!(pd.len(), sieve.len());
        for (i, d) in pd.iter().enumerate() {
            assert_eq!(d.p as u64, sieve.primes()[i]);
            if d.p > 2 {
                assert_eq!(
                    (d.p as u64).wrapping_mul(d.inv_p),
                    1u64,
                    "inverse check failed for p={}",
                    d.p
                );
                assert_eq!(d.max_quot, u64::MAX / d.p as u64);
            }
        }
    }
}
