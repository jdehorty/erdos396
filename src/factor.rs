//! Integer factorization using trial division with precomputed primes.
//!
//! For numbers up to 10^10, trial division with a precomputed prime sieve
//! is highly efficient. Each factorization is O(sqrt(n) / ln(sqrt(n))).

use crate::sieve::PrimeSieve;

/// A prime factorization represented as (prime, exponent) pairs.
///
/// The factors are stored in ascending order by prime.
/// Uses a small-vector optimization since most numbers have few prime factors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Factorization {
    /// Factors stored as (prime, exponent) pairs
    factors: smallvec::SmallVec<[(u64, u32); 8]>,
}

// Use smallvec to avoid heap allocation for numbers with ≤8 distinct prime factors
use smallvec::SmallVec;

impl Factorization {
    /// Create an empty factorization (represents 1).
    #[inline]
    pub fn one() -> Self {
        Self {
            factors: SmallVec::new(),
        }
    }

    /// Factorize `n` using the provided prime sieve.
    ///
    /// # Panics
    /// Panics in debug mode if the sieve is insufficient for `n`.
    #[inline]
    pub fn of(mut n: u64, sieve: &PrimeSieve) -> Self {
        debug_assert!(
            sieve.can_factor(n),
            "Sieve limit {} is insufficient for n={}",
            sieve.limit(),
            n
        );

        if n <= 1 {
            return Self::one();
        }

        let mut factors: SmallVec<[(u64, u32); 8]> = SmallVec::new();
        let mut sqrt_n = (n as f64).sqrt() as u64;

        for &p in sieve.primes() {
            if p > sqrt_n {
                break;
            }

            if n % p == 0 {
                let mut exp = 0u32;
                while n % p == 0 {
                    exp += 1;
                    n /= p;
                }
                factors.push((p, exp));

                if n == 1 {
                    break;
                }
                sqrt_n = (n as f64).sqrt() as u64;
            }
        }

        // Any remaining n > 1 is prime
        if n > 1 {
            factors.push((n, 1));
        }

        Self { factors }
    }

    /// Get the prime factors as a slice of (prime, exponent) pairs.
    #[inline]
    pub fn factors(&self) -> &[(u64, u32)] {
        &self.factors
    }

    /// Check if this is the factorization of 1 (no prime factors).
    #[inline]
    pub fn is_one(&self) -> bool {
        self.factors.is_empty()
    }

    /// Get the number of distinct prime factors (ω function).
    #[inline]
    pub fn num_distinct_primes(&self) -> usize {
        self.factors.len()
    }

    /// Get the total number of prime factors with multiplicity (Ω function).
    #[inline]
    pub fn num_prime_factors(&self) -> u32 {
        self.factors.iter().map(|(_, e)| *e).sum()
    }

    /// Reconstruct the original number from its factorization.
    #[inline]
    pub fn to_number(&self) -> u64 {
        self.factors
            .iter()
            .fold(1u64, |acc, &(p, e)| acc * p.pow(e))
    }

    /// Get the exponent of prime `p` in this factorization.
    /// Returns 0 if `p` does not divide the number.
    #[inline]
    pub fn exponent(&self, p: u64) -> u32 {
        self.factors
            .iter()
            .find(|&&(prime, _)| prime == p)
            .map(|&(_, exp)| exp)
            .unwrap_or(0)
    }

    /// Iterate over (prime, exponent) pairs.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (u64, u32)> + '_ {
        self.factors.iter().copied()
    }
}

impl IntoIterator for Factorization {
    type Item = (u64, u32);
    type IntoIter = smallvec::IntoIter<[(u64, u32); 8]>;

    fn into_iter(self) -> Self::IntoIter {
        self.factors.into_iter()
    }
}

/// Compute the p-adic valuation v_p(n) directly without full factorization.
///
/// This is faster when only one prime's exponent is needed.
#[inline]
pub fn vp(n: u64, p: u64) -> u32 {
    if n == 0 || p <= 1 {
        return 0;
    }

    let mut n = n;
    let mut v = 0u32;
    while n % p == 0 {
        v += 1;
        n /= p;
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorization() {
        let sieve = PrimeSieve::new(100);

        // 1 = (empty)
        let f = Factorization::of(1, &sieve);
        assert!(f.is_one());
        assert_eq!(f.to_number(), 1);

        // 12 = 2^2 * 3
        let f = Factorization::of(12, &sieve);
        assert_eq!(f.factors(), &[(2, 2), (3, 1)]);
        assert_eq!(f.to_number(), 12);

        // 100 = 2^2 * 5^2
        let f = Factorization::of(100, &sieve);
        assert_eq!(f.factors(), &[(2, 2), (5, 2)]);

        // 97 is prime
        let f = Factorization::of(97, &sieve);
        assert_eq!(f.factors(), &[(97, 1)]);
    }

    #[test]
    fn test_vp() {
        assert_eq!(vp(12, 2), 2);
        assert_eq!(vp(12, 3), 1);
        assert_eq!(vp(12, 5), 0);
        assert_eq!(vp(1000, 2), 3);
        assert_eq!(vp(1000, 5), 3);
    }

    #[test]
    fn test_large_number() {
        let sieve = PrimeSieve::new(50000);
        let n = 339_949_252u64; // Known k=8 witness

        let f = Factorization::of(n, &sieve);
        assert_eq!(f.to_number(), n);
    }
}
