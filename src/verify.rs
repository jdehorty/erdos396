//! Witness verification for Erdős Problem #396.
//!
//! A k-witness is an integer n such that n(n-1)(n-2)···(n-k) | C(2n, n).
//!
//! Verification requires checking that for every prime p:
//!     sum of v_p(n-i) for i=0..k ≤ v_p(C(2n, n))

use crate::factor::Factorization;
use crate::governor::{
    vp_central_binom_kummer_fast, vp_central_binom_p2, vp_central_binom_p3, vp_central_binom_p5,
    vp_factorial, GovernorChecker,
};
use crate::sieve::PrimeSieve;
use crate::{Error, Result};
use std::collections::HashMap;

/// Verifier for Erdős 396 witnesses.
#[derive(Debug, Clone)]
pub struct WitnessVerifier {
    sieve: PrimeSieve,
}

impl WitnessVerifier {
    /// Create a new verifier for numbers up to `max_n`.
    pub fn new(max_n: u64) -> Self {
        Self {
            sieve: PrimeSieve::for_range(max_n),
        }
    }

    /// Create a verifier from an existing prime sieve.
    pub fn with_sieve(sieve: PrimeSieve) -> Self {
        Self { sieve }
    }

    /// Verify that n is a valid k-witness.
    ///
    /// This checks: n(n-1)(n-2)···(n-k) | C(2n, n)
    ///
    /// Returns detailed verification result.
    pub fn verify(&self, k: u32, n: u64) -> Result<VerificationResult> {
        // Check basic validity
        if n <= k as u64 {
            return Ok(VerificationResult {
                is_valid: false,
                is_governor_run: false,
                failing_prime: None,
                demand: HashMap::new(),
                supply: HashMap::new(),
            });
        }
        if n.checked_mul(2).is_none() {
            return Err(Error::InvalidParameter(format!(
                "n too large: 2n overflows u64 (n={})",
                n
            )));
        }

        // First check if it's a Governor run (necessary condition based on conjecture)
        let checker = GovernorChecker::with_sieve(self.sieve.clone());
        let is_governor_run = checker.is_governor_run(n, k);

        // Collect all prime factors from the block with their total exponents
        let mut demand: HashMap<u64, u64> = HashMap::new();

        for i in 0..=k {
            let term = n - i as u64;
            let factors = Factorization::of(term, &self.sieve);
            for (p, exp) in factors.iter() {
                *demand.entry(p).or_insert(0) += exp as u64;
            }
        }

        // Check supply for each prime
        let mut supply: HashMap<u64, u64> = HashMap::new();
        let mut failing_prime = None;

        for (&p, &required) in &demand {
            let available = vp_supply_checked(n, p)?;
            supply.insert(p, available);

            if required > available && failing_prime.is_none() {
                failing_prime = Some(p);
            }
        }

        Ok(VerificationResult {
            is_valid: failing_prime.is_none(),
            is_governor_run,
            failing_prime,
            demand,
            supply,
        })
    }

    /// Quick check if n is a valid k-witness (without detailed results).
    #[inline]
    pub fn is_witness(&self, k: u32, n: u64) -> Result<bool> {
        if n <= k as u64 {
            return Ok(false);
        }

        // Collect all prime factors from the block
        let mut demand: HashMap<u64, u64> = HashMap::new();

        for i in 0..=k {
            let term = n - i as u64;
            let factors = Factorization::of(term, &self.sieve);
            for (p, exp) in factors.iter() {
                *demand.entry(p).or_insert(0) += exp as u64;
            }
        }

        // Check supply >= demand for each prime
        for (p, required) in demand {
            let available = vp_supply_checked(n, p)?;
            if required > available {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Fast verification using Kummer carry-monotonicity theorem.
    ///
    /// Intended for *governor-run candidates* (i.e., all `k+1` block terms are
    /// in the Governor Set).
    ///
    /// For a governor run, primes `p ≥ 2k+1` are automatically satisfied:
    /// among `k+1` consecutive integers at most one is divisible by such a `p`,
    /// and carry-invariance gives `v_p(C(2n,n)) = v_p(C(2(n-j), n-j))` when
    /// `p ∣ (n-j)` and `j ≤ k ≤ (p-1)/2`. So it suffices to check only the
    /// **barrier primes** `p < 2k+1`.
    ///
    /// Demand is computed via Legendre's formula on the block product:
    ///   demand_p = v_p(n!) - v_p((n-k-1)!)
    /// Supply is computed via Kummer dispatch.
    pub fn verify_fast(&self, k: u32, n: u64) -> Result<VerificationResult> {
        if n <= k as u64 {
            return Ok(VerificationResult {
                is_valid: false,
                is_governor_run: false,
                failing_prime: None,
                demand: HashMap::new(),
                supply: HashMap::new(),
            });
        }
        if n.checked_mul(2).is_none() {
            return Err(Error::InvalidParameter(format!(
                "n too large: 2n overflows u64 (n={})",
                n
            )));
        }

        let checker = GovernorChecker::with_sieve(self.sieve.clone());
        let is_governor_run = checker.is_governor_run(n, k);
        if !is_governor_run {
            // `verify_fast` is only sound as a verifier for governor-run candidates.
            // Fall back to the full verifier to avoid false positives.
            return self.verify(k, n);
        }

        let barrier_threshold = k.saturating_mul(2);
        let barrier_primes = primes_up_to(barrier_threshold);
        let mut demand: HashMap<u64, u64> = HashMap::new();
        let mut supply: HashMap<u64, u64> = HashMap::new();
        let mut failing_prime = None;

        let block_bottom = n - k as u64; // n-k is lowest term, block is n(n-1)...(n-k)

        for p in &barrier_primes {
            let p = *p;
            // demand_p = v_p(n!) - v_p((n-k-1)!)
            // This equals sum of v_p(term) for term in [n-k, n]
            let d = vp_factorial(n, p) - vp_factorial(block_bottom.saturating_sub(1), p);
            let s = vp_supply_checked(n, p)?;

            demand.insert(p, d);
            supply.insert(p, s);

            if d > s && failing_prime.is_none() {
                failing_prime = Some(p);
            }
        }

        Ok(VerificationResult {
            is_valid: failing_prime.is_none(),
            is_governor_run,
            failing_prime,
            demand,
            supply,
        })
    }

    /// Quick boolean check for governor-run candidates (barrier primes only).
    #[inline]
    pub fn is_witness_fast(&self, k: u32, n: u64) -> Result<bool> {
        if n <= k as u64 {
            return Ok(false);
        }
        if n.checked_mul(2).is_none() {
            return Err(Error::InvalidParameter(format!(
                "n too large: 2n overflows u64 (n={})",
                n
            )));
        }

        let checker = GovernorChecker::with_sieve(self.sieve.clone());
        if !checker.is_governor_run(n, k) {
            return self.is_witness(k, n);
        }

        let barrier_threshold = k.saturating_mul(2);
        let barrier_primes = primes_up_to(barrier_threshold);
        let block_bottom = n - k as u64;

        for &p in &barrier_primes {
            let d = vp_factorial(n, p) - vp_factorial(block_bottom.saturating_sub(1), p);
            let s = vp_supply_checked(n, p)?;
            if d > s {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Check barrier primes `p < 2k+1` for a block starting at n with k+1 terms.
    ///
    /// Returns `None` if all barrier primes pass, `Some(p)` for the first
    /// failing prime. Used by safety-net mode to detect potential
    /// counter-examples.
    #[inline]
    pub fn check_small_primes(&self, k: u32, n: u64) -> Result<Option<u64>> {
        if n <= k as u64 {
            return Ok(Some(2)); // trivially fails
        }

        let barrier_threshold = k.saturating_mul(2);
        let barrier_primes = primes_up_to(barrier_threshold);
        let block_bottom = n - k as u64;

        for &p in &barrier_primes {
            let d = vp_factorial(n, p) - vp_factorial(block_bottom.saturating_sub(1), p);
            let s = vp_supply_checked(n, p)?;
            if d > s {
                return Ok(Some(p));
            }
        }

        Ok(None)
    }
}

/// Compute v_p(C(2n, n)) using the best available method for prime p.
#[inline]
fn vp_supply(n: u64, p: u64) -> u64 {
    match p {
        2 => vp_central_binom_p2(n),
        3 => vp_central_binom_p3(n),
        5 => vp_central_binom_p5(n),
        _ => vp_central_binom_kummer_fast(n, p),
    }
}

#[inline]
fn vp_central_binom_legendre_checked(n: u64, p: u64) -> Result<u64> {
    let two_n = n.checked_mul(2).ok_or_else(|| {
        Error::InvalidParameter(format!("n too large: 2n overflows u64 (n={})", n))
    })?;
    Ok(vp_factorial(two_n, p) - 2 * vp_factorial(n, p))
}

#[inline]
fn vp_supply_checked(n: u64, p: u64) -> Result<u64> {
    let legendre = vp_central_binom_legendre_checked(n, p)?;
    let kummer = vp_supply(n, p);
    if legendre != kummer {
        return Err(Error::Audit(format!(
            "internal inconsistency: supply mismatch at n={}, p={}: Legendre={}, Kummer={}",
            n, p, legendre, kummer
        )));
    }
    Ok(legendre)
}

/// Return all primes up to k (tiny sieve, only used for small k).
pub fn primes_up_to(k: u32) -> Vec<u64> {
    if k < 2 {
        return Vec::new();
    }
    let limit = k as usize + 1;
    let mut is_prime = vec![true; limit];
    is_prime[0] = false;
    if limit > 1 {
        is_prime[1] = false;
    }
    let mut i = 2;
    while i * i < limit {
        if is_prime[i] {
            let mut j = i * i;
            while j < limit {
                is_prime[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    is_prime
        .iter()
        .enumerate()
        .filter_map(|(i, &is_p)| if is_p { Some(i as u64) } else { None })
        .collect()
}

/// Detailed result of witness verification.
#[derive(Debug, Clone)]
pub struct VerificationResult {
    /// Whether n is a valid k-witness
    pub is_valid: bool,

    /// Whether all block terms are in the Governor Set
    pub is_governor_run: bool,

    /// The first prime that failed (if any)
    pub failing_prime: Option<u64>,

    /// p-adic demand: sum of v_p(n-i) for each prime p
    pub demand: HashMap<u64, u64>,

    /// p-adic supply: v_p(C(2n, n)) for each prime p
    pub supply: HashMap<u64, u64>,
}

impl VerificationResult {
    /// Get a human-readable summary of the verification.
    pub fn summary(&self, k: u32, n: u64) -> String {
        let mut s = String::new();

        s.push_str(&format!("Verification of k={}, n={}\n", k, n));
        s.push_str(&format!("  Valid witness: {}\n", self.is_valid));
        s.push_str(&format!("  Governor run: {}\n", self.is_governor_run));

        if let Some(p) = self.failing_prime {
            s.push_str(&format!(
                "  Failing prime: p={} (demand={}, supply={})\n",
                p,
                self.demand.get(&p).unwrap_or(&0),
                self.supply.get(&p).unwrap_or(&0)
            ));
        }

        // Show a few primes with their demand/supply
        s.push_str("  Prime breakdown (first few):\n");
        let mut primes: Vec<_> = self.demand.keys().copied().collect();
        primes.sort();
        for p in primes.iter().take(10) {
            let d = self.demand.get(p).unwrap_or(&0);
            let sup = self.supply.get(p).unwrap_or(&0);
            let status = if d <= sup { "✓" } else { "✗" };
            s.push_str(&format!(
                "    p={}: demand={}, supply={} {}\n",
                p, d, sup, status
            ));
        }

        if primes.len() > 10 {
            s.push_str(&format!("    ... and {} more primes\n", primes.len() - 10));
        }

        s
    }
}

/// Verify all known OEIS witnesses for correctness.
pub fn verify_known_witnesses() -> Result<Vec<(u32, u64, bool)>> {
    use crate::KNOWN_WITNESSES;

    let max_n = KNOWN_WITNESSES
        .iter()
        .map(|(_, n)| *n)
        .max()
        .unwrap_or(1000);
    let verifier = WitnessVerifier::new(max_n);

    let mut out = Vec::with_capacity(KNOWN_WITNESSES.len());
    for &(k, n) in KNOWN_WITNESSES {
        out.push((k, n, verifier.is_witness(k, n)?));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factor::vp as vp_term;
    use crate::KNOWN_WITNESSES;

    fn lcg_next(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state
    }

    #[test]
    fn test_known_witnesses() {
        let results = verify_known_witnesses().unwrap();

        for (k, n, is_valid) in results {
            assert!(is_valid, "Known witness k={}, n={} should be valid", k, n);
        }
    }

    #[test]
    fn test_verification_details() {
        let verifier = WitnessVerifier::new(1_000_000_000);

        // k=8 witness
        let result = verifier.verify(8, 339_949_252).unwrap();
        assert!(result.is_valid);
        assert!(result.is_governor_run);
        assert!(result.failing_prime.is_none());

        println!("{}", result.summary(8, 339_949_252));
    }

    #[test]
    fn test_false_positive_detection() {
        let verifier = WitnessVerifier::new(100_000);

        // n=30592 is mentioned as a Governor run of 3 but NOT a valid k=2 witness
        let result = verifier.verify(2, 30592).unwrap();

        // It should be a governor run but not a valid witness
        // (This is the expected false positive case)
        if result.is_governor_run && !result.is_valid {
            println!("Confirmed false positive at n=30592");
            println!("{}", result.summary(2, 30592));
        }
    }

    #[test]
    fn test_non_witness() {
        let verifier = WitnessVerifier::new(1000);

        // n=100 should not be a k=3 witness
        let result = verifier.verify(3, 100).unwrap();
        assert!(!result.is_valid);
    }

    #[test]
    fn test_primes_up_to() {
        assert_eq!(primes_up_to(0), Vec::<u64>::new());
        assert_eq!(primes_up_to(1), Vec::<u64>::new());
        assert_eq!(primes_up_to(2), vec![2]);
        assert_eq!(primes_up_to(3), vec![2, 3]);
        assert_eq!(primes_up_to(5), vec![2, 3, 5]);
        assert_eq!(primes_up_to(7), vec![2, 3, 5, 7]);
        assert_eq!(primes_up_to(9), vec![2, 3, 5, 7]);
        assert_eq!(primes_up_to(10), vec![2, 3, 5, 7]);
        assert_eq!(primes_up_to(11), vec![2, 3, 5, 7, 11]);
    }

    #[test]
    fn test_verify_fast_matches_full_for_known_witnesses() {
        let max_n = KNOWN_WITNESSES
            .iter()
            .map(|(_, n)| *n)
            .max()
            .unwrap_or(1000);
        let verifier = WitnessVerifier::new(max_n);

        for &(k, n) in KNOWN_WITNESSES {
            let full = verifier.verify(k, n).unwrap();
            let fast = verifier.verify_fast(k, n).unwrap();

            assert_eq!(
                full.is_valid, fast.is_valid,
                "verify vs verify_fast disagree for k={}, n={}: full={}, fast={}",
                k, n, full.is_valid, fast.is_valid
            );
            assert!(
                fast.is_valid,
                "Known witness k={}, n={} should pass fast verification",
                k, n
            );
        }
    }

    #[test]
    fn test_verify_fast_catches_false_positive() {
        // n=17,842,967,551 is a run of 10 governors but fails at p=3
        let n = 17_842_967_551u64;
        let k = 9u32;
        let verifier = WitnessVerifier::new(n);

        let fast = verifier.verify_fast(k, n).unwrap();
        assert!(
            !fast.is_valid,
            "Known false positive should fail fast verify"
        );
        assert_eq!(fast.failing_prime, Some(3), "Should fail at p=3");
    }

    #[test]
    fn test_is_witness_fast_matches_is_witness() {
        let verifier = WitnessVerifier::new(100_000);

        // Check a range for k=2
        let checker = GovernorChecker::with_sieve(PrimeSieve::for_range(100_000));
        for n in 3..=10_000u64 {
            if checker.is_governor_run(n, 2) {
                let full = verifier.is_witness(2, n).unwrap();
                let fast = verifier.is_witness_fast(2, n).unwrap();
                assert_eq!(
                    full, fast,
                    "is_witness vs is_witness_fast disagree at n={}: full={}, fast={}",
                    n, full, fast
                );
            }
        }
    }

    #[test]
    fn test_check_small_primes_known_witnesses() {
        let max_n = KNOWN_WITNESSES
            .iter()
            .map(|(_, n)| *n)
            .max()
            .unwrap_or(1000);
        let verifier = WitnessVerifier::new(max_n);

        for &(k, n) in KNOWN_WITNESSES {
            let result = verifier.check_small_primes(k, n).unwrap();
            assert!(
                result.is_none(),
                "Known witness k={}, n={} should pass check_small_primes, but failed at p={:?}",
                k,
                n,
                result
            );
        }
    }

    #[test]
    fn test_check_small_primes_false_positive() {
        let n = 17_842_967_551u64;
        let k = 9u32;
        let verifier = WitnessVerifier::new(n);

        let result = verifier.check_small_primes(k, n).unwrap();
        assert_eq!(result, Some(3), "Should detect failure at p=3");
    }

    #[test]
    fn test_verify_fast_matches_full_for_governor_runs_small_range() {
        let verifier = WitnessVerifier::new(100_000);
        let checker = GovernorChecker::with_sieve(PrimeSieve::for_range(100_000));

        // Check all governor runs of 3+ in [100, 10000] for k=2
        for n in 102..=10_000u64 {
            if checker.is_governor_run(n, 2) {
                let full = verifier.verify(2, n).unwrap();
                let fast = verifier.verify_fast(2, n).unwrap();
                assert_eq!(
                    full.is_valid, fast.is_valid,
                    "Disagree at n={}: full={}, fast={}",
                    n, full.is_valid, fast.is_valid
                );
            }
        }
    }

    #[test]
    fn test_block_demand_sum_matches_factorial_difference() {
        // Sanity check for the core demand identity used by `validate` and `verify_fast`:
        //   v_p(∏_{i=0..k} (n-i)) = v_p(n!) - v_p((n-k-1)!)
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 19, 23];
        let mut state = 0x8d9a_aa52_2f6d_c7b1u64;

        for _ in 0..2000 {
            let k = (lcg_next(&mut state) % 20) as u32;
            let n = (lcg_next(&mut state) % 1_000_000_000_000u64).saturating_add(k as u64 + 1);
            let p = primes[(lcg_next(&mut state) as usize) % primes.len()];

            let sum_demand: u64 = (0..=k).map(|i| vp_term(n - i as u64, p) as u64).sum();

            let block_bottom = n - k as u64;
            let legendre_demand = vp_factorial(n, p) - vp_factorial(block_bottom - 1, p);

            assert_eq!(
                sum_demand, legendre_demand,
                "Demand mismatch at n={}, k={}, p={}: sum={}, factorial_diff={}",
                n, k, p, sum_demand, legendre_demand
            );
        }
    }
}
