use erdos396::factor::vp;
use erdos396::governor::vp_central_binom_dispatch;
use erdos396::Factorization;
use erdos396::PrimeSieve;
use serde::Serialize;

/// A near-miss: a block of k+1 consecutive integers with exactly one non-governor.
#[derive(Debug, Clone, Serialize)]
pub struct NearMiss {
    pub n: u64,
    pub m: u64,
    pub j: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct PrimeAnalysis {
    pub p: u64,
    pub demand: u64,
    pub supply: u64,
    pub deficit: i64,
    pub is_failing: bool,
    pub top_supply: u64,
    pub bonus: i64,
    pub gap: i64,
    pub cascade_len: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct WitnessAnalysis {
    pub p: u64,
    pub total_demand: u64,
    pub top_supply: u64,
    pub witness_gap: i64,
}

#[derive(Debug, Clone, Serialize)]
pub struct NearMissReport {
    pub near_miss: NearMiss,
    pub barrier_primes: Vec<PrimeAnalysis>,
    pub witness_condition: Vec<WitnessAnalysis>,
    pub failing_prime_count: u32,
    pub compensated_count: u32,
    pub witness_gaps_all_nonpositive: bool,
    /// True if m fails ONLY at barrier primes (no large-prime-factor kill).
    /// These are the only near-misses that could theoretically be non-governor witnesses.
    pub pure_barrier_only: bool,
}

/// Return all primes p with 2 <= p < 2k+1.
pub fn barrier_primes(k: u32) -> Vec<u64> {
    let bound = 2 * k as u64 + 1;
    let mut primes = Vec::new();
    for candidate in 2..bound {
        let is_prime = (2..candidate).all(|d| candidate % d != 0);
        if is_prime {
            primes.push(candidate);
        }
    }
    primes
}

/// Count consecutive critical digits in the base-p representation of n,
/// starting at digit position `start_pos`.
///
/// A digit d is "critical" if d >= floor(p/2). Critical digits propagate
/// carries during self-addition (Kummer's theorem). The cascade length
/// determines the maximum supply bonus from carry propagation.
pub fn cascade_length(n: u64, p: u64, start_pos: u32) -> u32 {
    let threshold = p / 2; // floor(p/2)
    let mut remaining = n;

    // Skip to start_pos by dividing out p^start_pos
    for _ in 0..start_pos {
        remaining /= p;
    }

    let mut count = 0u32;
    while remaining > 0 {
        let digit = remaining % p;
        if digit >= threshold {
            count += 1;
        } else {
            break;
        }
        remaining /= p;
    }
    count
}

/// Compute per-prime marginal analysis for a near-miss.
pub fn analyze_prime(n: u64, m: u64, p: u64) -> PrimeAnalysis {
    let demand = vp(m, p) as u64;
    let supply = vp_central_binom_dispatch(m, p);
    let deficit = demand as i64 - supply as i64;
    let is_failing = deficit > 0;
    let top_supply = vp_central_binom_dispatch(n, p);
    let bonus = top_supply as i64 - supply as i64;
    let gap = deficit - bonus;

    let start_pos = if demand > 0 { demand as u32 } else { 0 };
    let cascade_len = cascade_length(n, p, start_pos);

    PrimeAnalysis {
        p,
        demand,
        supply,
        deficit,
        is_failing,
        top_supply,
        bonus,
        gap,
        cascade_len,
    }
}

/// Compute the full witness condition at barrier prime p for a block [n-k, n].
pub fn analyze_witness_condition(n: u64, k: u32, p: u64) -> WitnessAnalysis {
    let mut total_demand = 0u64;
    for i in 0..=k as u64 {
        total_demand += vp(n - i, p) as u64;
    }
    let top_supply = vp_central_binom_dispatch(n, p);
    let witness_gap = total_demand as i64 - top_supply as i64;

    WitnessAnalysis {
        p,
        total_demand,
        top_supply,
        witness_gap,
    }
}

/// Compute the full analysis of a near-miss block.
pub fn analyze_near_miss(nm: &NearMiss, k: u32, sieve: &PrimeSieve) -> NearMissReport {
    let bp = barrier_primes(k);

    let barrier_analysis: Vec<PrimeAnalysis> =
        bp.iter().map(|&p| analyze_prime(nm.n, nm.m, p)).collect();

    let witness_analysis: Vec<WitnessAnalysis> = bp
        .iter()
        .map(|&p| analyze_witness_condition(nm.n, k, p))
        .collect();

    let failing_prime_count = barrier_analysis.iter().filter(|a| a.is_failing).count() as u32;

    let compensated_count = barrier_analysis
        .iter()
        .filter(|a| a.is_failing && a.gap <= 0)
        .count() as u32;

    let witness_gaps_all_nonpositive = witness_analysis.iter().all(|w| w.witness_gap <= 0);

    let pure_barrier_only = is_pure_barrier_failure(nm.m, k, sieve);

    NearMissReport {
        near_miss: nm.clone(),
        barrier_primes: barrier_analysis,
        witness_condition: witness_analysis,
        failing_prime_count,
        compensated_count,
        witness_gaps_all_nonpositive,
        pure_barrier_only,
    }
}

/// Check whether m fails the governor test at any barrier prime for the given k.
pub fn has_barrier_prime_failure(m: u64, k: u32, sieve: &PrimeSieve) -> bool {
    let bp = barrier_primes(k);
    let factors = Factorization::of(m, sieve);
    for (p, exp) in factors.iter() {
        if bp.contains(&p) {
            let supply = vp_central_binom_dispatch(m, p);
            if (exp as u64) > supply {
                return true;
            }
        }
    }
    false
}

/// Check whether m fails the governor test ONLY at barrier primes.
///
/// Returns true if m has no governor failure at any non-barrier prime (p >= 2k+1).
/// This filters out near-misses killed by the sqrt(2n) barrier or other large-prime
/// failures, isolating cases where carry compensation is the sole question.
///
/// By Theorem 1, a non-governor witness can only exist if the non-governor term
/// passes at all primes >= 2k+1. So "pure barrier failures" are the only near-misses
/// that could theoretically be non-governor witnesses.
pub fn is_pure_barrier_failure(m: u64, k: u32, sieve: &PrimeSieve) -> bool {
    let bound = 2 * k as u64 + 1;
    let factors = Factorization::of(m, sieve);
    let mut has_barrier_fail = false;

    for (p, exp) in factors.iter() {
        let supply = vp_central_binom_dispatch(m, p);
        let fails = (exp as u64) > supply;

        if p < bound {
            // Barrier prime — allowed to fail
            if fails {
                has_barrier_fail = true;
            }
        } else {
            // Non-barrier prime — must NOT fail
            if fails {
                return false;
            }
        }
    }

    has_barrier_fail
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barrier_primes_k14() {
        let bp = barrier_primes(14);
        assert_eq!(bp, vec![2, 3, 5, 7, 11, 13, 17, 19, 23]);
    }

    #[test]
    fn test_barrier_primes_k2() {
        let bp = barrier_primes(2);
        assert_eq!(bp, vec![2, 3]);
    }

    #[test]
    fn test_barrier_primes_k8() {
        let bp = barrier_primes(8);
        assert_eq!(bp, vec![2, 3, 5, 7, 11, 13]);
    }

    #[test]
    fn test_cascade_length_p2_all_ones() {
        // n = 0b1111111 = 127. All digits are 1 in base 2.
        // floor(2/2) = 1, so critical digit threshold is 1.
        assert_eq!(cascade_length(127, 2, 0), 7);
    }

    #[test]
    fn test_cascade_length_p2_mixed() {
        // n = 0b1010110 = 86. Digits (LSB first): 0, 1, 1, 0, 1, 0, 1
        assert_eq!(cascade_length(86, 2, 0), 0);
        assert_eq!(cascade_length(86, 2, 1), 2);
    }

    #[test]
    fn test_cascade_length_p3() {
        // n = 26 in base 3: (222)_3. All digits = 2. floor(3/2) = 1.
        assert_eq!(cascade_length(26, 3, 0), 3);
    }

    #[test]
    fn test_cascade_length_p7() {
        // n = 73. Base 7: (1, 3, 3). Digits LSB first: 3, 3, 1.
        // floor(7/2) = 3. d=3 is critical, d=1 is not.
        assert_eq!(cascade_length(73, 7, 0), 2);
    }

    #[test]
    fn test_cascade_length_start_above_zero() {
        // n = 0b1011110 = 94. Digits LSB first: 0, 1, 1, 1, 1, 0, 1
        // Start at position 2: pos2=1 critical, pos3=1, pos4=1, pos5=0 stop. cascade = 3.
        assert_eq!(cascade_length(94, 2, 2), 3);
    }

    #[test]
    fn test_analyze_prime_proposition2() {
        // Proposition 2 tightness example: m=4, n=5, p=2
        let result = analyze_prime(5, 4, 2);
        assert_eq!(result.p, 2);
        assert_eq!(result.demand, 2); // v_2(4) = 2
        assert_eq!(result.supply, 1); // s_2(4) = popcount(4) = 1
        assert_eq!(result.deficit, 1);
        assert!(result.is_failing);
        assert_eq!(result.top_supply, 2); // s_2(5) = popcount(5) = 2
        assert_eq!(result.bonus, 1);
        assert_eq!(result.gap, 0); // exactly compensated!
        assert_eq!(result.cascade_len, 1);
    }

    #[test]
    fn test_analyze_prime_not_failing() {
        // m=6, p=2: v_2(6)=1, s_2(6)=popcount(6)=2. Not failing.
        let result = analyze_prime(7, 6, 2);
        assert_eq!(result.demand, 1);
        assert_eq!(result.supply, 2);
        assert_eq!(result.deficit, -1);
        assert!(!result.is_failing);
    }

    #[test]
    fn test_analyze_prime_p3() {
        // m=9, p=3: v_3(9)=2, s_3(9)=0 (carries in (100)_3+(100)_3 = 0).
        // n=10: s_3(10)=0. bonus=0. gap=2.
        let result = analyze_prime(10, 9, 3);
        assert_eq!(result.demand, 2);
        assert_eq!(result.supply, 0);
        assert_eq!(result.deficit, 2);
        assert!(result.is_failing);
        assert_eq!(result.top_supply, 0);
        assert_eq!(result.bonus, 0);
        assert_eq!(result.gap, 2);
    }

    #[test]
    fn test_witness_analysis_p2() {
        // Block [3, 4, 5] (k=2, n=5). p=2.
        // v_2(5)=0, v_2(4)=2, v_2(3)=0. total_demand=2.
        // top_supply = s_2(5) = popcount(5) = 2. witness_gap = 0.
        let result = analyze_witness_condition(5, 2, 2);
        assert_eq!(result.total_demand, 2);
        assert_eq!(result.top_supply, 2);
        assert_eq!(result.witness_gap, 0);
    }

    #[test]
    fn test_witness_analysis_p3_k2() {
        // Block [7, 8, 9] (k=2, n=9). p=3.
        // v_3(9)=2, v_3(8)=0, v_3(7)=0. total_demand=2.
        // top_supply = s_3(9) = 0. witness_gap = 2.
        let result = analyze_witness_condition(9, 2, 3);
        assert_eq!(result.total_demand, 2);
        assert_eq!(result.top_supply, 0);
        assert_eq!(result.witness_gap, 2);
    }

    #[test]
    fn test_analyze_near_miss_proposition2() {
        // m=4, n=5, k=1. Barrier primes = {2} (since 2*1+1=3).
        let sieve = erdos396::PrimeSieve::for_range(100);
        let nm = NearMiss { n: 5, m: 4, j: 1 };
        let report = analyze_near_miss(&nm, 1, &sieve);

        assert_eq!(report.barrier_primes.len(), 1);
        assert_eq!(report.barrier_primes[0].p, 2);
        assert!(report.barrier_primes[0].is_failing);
        assert_eq!(report.barrier_primes[0].gap, 0);

        assert_eq!(report.witness_condition.len(), 1);
        assert_eq!(report.witness_condition[0].witness_gap, 0);

        assert_eq!(report.failing_prime_count, 1);
        assert_eq!(report.compensated_count, 1);
        assert!(report.witness_gaps_all_nonpositive);
        assert!(report.pure_barrier_only); // m=4=2^2, only fails at p=2 (barrier)
    }

    #[test]
    fn test_has_barrier_prime_failure_true() {
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(has_barrier_prime_failure(4, 1, &sieve));
    }

    #[test]
    fn test_has_barrier_prime_failure_false_governor() {
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(!has_barrier_prime_failure(6, 2, &sieve));
    }

    #[test]
    fn test_has_barrier_prime_failure_false_large_prime_only() {
        // m=7 fails at p=7, but for k=2 barrier primes = {2, 3}. p=7 is NOT barrier.
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(!has_barrier_prime_failure(7, 2, &sieve));
    }

    #[test]
    fn test_is_pure_barrier_failure_true() {
        // m=4=2^2: fails only at p=2 (barrier for k>=1). No large prime factor.
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(is_pure_barrier_failure(4, 1, &sieve));
    }

    #[test]
    fn test_is_pure_barrier_failure_false_large_factor() {
        // m=14=2*7: v_2(14)=1, s_2(14)=popcount(14)=3 (passes at p=2).
        //           v_7(14)=1, s_7(14)=0 (fails at p=7).
        // For k=2, barrier primes={2,3}. p=7 is non-barrier. So not pure.
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(!is_pure_barrier_failure(14, 2, &sieve));
    }

    #[test]
    fn test_is_pure_barrier_failure_false_governor() {
        // m=6: is a governor (passes at all primes). No barrier failure at all.
        let sieve = erdos396::PrimeSieve::for_range(100);
        assert!(!is_pure_barrier_failure(6, 2, &sieve));
    }

    #[test]
    fn test_known_k2_witness_no_barrier_failures() {
        let checker = erdos396::GovernorChecker::new(2600);
        for i in 0..=2u64 {
            assert!(
                checker.is_governor_fast(2480 - i),
                "Block term {} should be a governor",
                2480 - i
            );
        }
    }

    #[test]
    fn test_known_k3_non_governor_witness_analysis() {
        // n=1441378 is a known k=3 non-governor witness.
        // The block [n-3, n] should contain exactly one non-governor whose
        // failures are confined to barrier primes (p < 2*3+1 = 7), i.e. {2, 3, 5}.
        let n = 1_441_378u64;
        let k = 3u32;
        let sieve = erdos396::PrimeSieve::for_range(n + 100);
        let checker = erdos396::GovernorChecker::with_sieve(sieve.clone());

        // Find the non-governor in the block
        let mut non_gov_pos = None;
        for j in 0..=k {
            let m = n - j as u64;
            if !checker.is_governor_fast(m) {
                assert!(
                    non_gov_pos.is_none(),
                    "Block should have exactly one non-governor"
                );
                non_gov_pos = Some(j);
            }
        }
        let j = non_gov_pos.expect("Block must contain a non-governor");
        let m = n - j as u64;

        // Verify it's a pure barrier failure
        assert!(
            is_pure_barrier_failure(m, k, &sieve),
            "k=3 non-governor witness m={m} should be a pure barrier failure"
        );

        // Full near-miss report
        let nm = NearMiss { n, m, j };
        let report = analyze_near_miss(&nm, k, &sieve);
        assert!(report.pure_barrier_only);
        assert!(
            report.failing_prime_count > 0,
            "Non-governor must fail at least one barrier prime"
        );
    }

    #[test]
    fn test_known_k4_non_governor_witness_analysis() {
        // n=2366563 is a known k=4 non-governor witness.
        // Barrier primes for k=4: {2, 3, 5, 7} (primes < 2*4+1 = 9).
        let n = 2_366_563u64;
        let k = 4u32;
        let sieve = erdos396::PrimeSieve::for_range(n + 100);
        let checker = erdos396::GovernorChecker::with_sieve(sieve.clone());

        // Find the non-governor in the block
        let mut non_gov_pos = None;
        for j in 0..=k {
            let m = n - j as u64;
            if !checker.is_governor_fast(m) {
                assert!(
                    non_gov_pos.is_none(),
                    "Block should have exactly one non-governor"
                );
                non_gov_pos = Some(j);
            }
        }
        let j = non_gov_pos.expect("Block must contain a non-governor");
        let m = n - j as u64;

        // Verify it's a pure barrier failure
        assert!(
            is_pure_barrier_failure(m, k, &sieve),
            "k=4 non-governor witness m={m} should be a pure barrier failure"
        );

        // Full near-miss report
        let nm = NearMiss { n, m, j };
        let report = analyze_near_miss(&nm, k, &sieve);
        assert!(report.pure_barrier_only);
        assert!(
            report.failing_prime_count > 0,
            "Non-governor must fail at least one barrier prime"
        );
    }

    #[test]
    fn test_proposition2_is_candidate() {
        let sieve = erdos396::PrimeSieve::for_range(100);
        let nm = NearMiss { n: 5, m: 4, j: 1 };
        let report = analyze_near_miss(&nm, 1, &sieve);
        assert!(
            report.witness_gaps_all_nonpositive,
            "Proposition 2 example should be a candidate"
        );
        // Verify with WitnessVerifier that it's NOT a true witness.
        // 5*4=20 does not divide C(10,5)=252.
        let verifier = erdos396::WitnessVerifier::new(5 + 100);
        let is_witness = verifier.is_witness(1, 5).unwrap();
        assert!(
            !is_witness,
            "Proposition 2 example passes barrier primes but is NOT a true witness"
        );
    }
}
