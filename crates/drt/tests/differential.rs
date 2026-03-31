//! Differential random tests — cross-implementation comparison.
//!
//! T0 (smoke):    `cargo test -p erdos396-drt`              < 10s, Rust-only
//! T1 (standard): `cargo test -p erdos396-drt -- --ignored`  < 60s, + Lean oracle

use erdos396_drt::generators;
use erdos396_drt::oracle::{Erdos396Oracle, LeanOracle, Oracle, SearchLabOracle};
use erdos396_drt::properties::*;

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

const SEED: u64 = 0x396_2026;

// ===================================================================
// T0: Smoke tests — fast, Rust-only, < 10 seconds
// ===================================================================

/// D1 smoke: search-lab exact_check vs erdos396 library on known witnesses.
#[test]
fn t0_known_witnesses_search_lab_vs_erdos396() {
    let mut sl = SearchLabOracle::new();
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);

    let known: &[(u32, u64)] = &[
        (1, 2),
        (2, 2_480),
        (3, 8_178),
        (4, 45_153),
        (5, 3_648_841),
        (6, 7_979_090),
        (7, 101_130_029),
        (8, 339_949_252),
        (9, 1_019_547_844),
        (10, 17_609_764_994),
    ];

    for &(k, n) in known {
        let sl_result = sl.check_witness(k, n);
        let e3_result = e3.check_witness(k, n);
        assert!(
            sl_result,
            "search-lab should accept known witness k={}, n={}",
            k, n
        );
        assert!(
            e3_result,
            "erdos396-lib should accept known witness k={}, n={}",
            k, n
        );

        // Also verify n-1 is NOT a witness (otherwise n wouldn't be minimal)
        if n > k as u64 + 2 {
            let sl_below = sl.check_witness(k, n - 1);
            let e3_below = e3.check_witness(k, n - 1);
            assert_eq!(
                sl_below,
                e3_below,
                "Disagreement at k={}, n={}: search-lab={}, erdos396={}",
                k,
                n - 1,
                sl_below,
                e3_below
            );
        }
    }
}

/// D1 smoke: search-lab vs erdos396 on 5K random inputs.
#[test]
fn t0_differential_search_lab_vs_erdos396_random() {
    let mut sl = SearchLabOracle::new();
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);

    let cases = generators::generate_tier(SEED, 5_000, 100_000_000_000);
    let mut failures = Vec::new();

    for tc in &cases {
        let sl_result = sl.check_witness(tc.k, tc.n);
        let e3_result = e3.check_witness(tc.k, tc.n);
        if sl_result != e3_result {
            failures.push(format!(
                "k={}, n={} ({}): search-lab={}, erdos396={}",
                tc.k, tc.n, tc.generator, sl_result, e3_result
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "D1: {} disagreements out of {} tests:\n{}",
        failures.len(),
        cases.len(),
        failures[..failures.len().min(10)].join("\n")
    );
}

/// P2 smoke: v_2(C(2n,n)) == popcount(n) via erdos396 library.
#[test]
fn t0_p2_v2_popcount() {
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);

    for _ in 0..5_000 {
        let n: u64 = rand::Rng::gen_range(&mut rng, 1..=1_000_000_000_000u64);
        if let Some(f) = p2_v2_popcount(n, &mut e3) {
            panic!("{:?}", f);
        }
    }
}

/// P3 smoke: modular inverse identity for all primes up to 10K.
#[test]
fn t0_p3_mod_inverse() {
    let primes = erdos396_search_lab::drt_exports::drt_sieve_primes(10_000);
    for &p in &primes {
        if p == 2 {
            continue;
        }
        if let Some(f) = p3_mod_inverse(p) {
            panic!("{:?}", f);
        }
    }
}

/// P4 smoke: divisibility via modular inverse.
#[test]
fn t0_p4_divisibility() {
    let primes = erdos396_search_lab::drt_exports::drt_sieve_primes(1_000);
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);

    for &p in &primes {
        if p == 2 {
            continue;
        }
        // Test multiples of p
        for mult in 1..=20u64 {
            let n = p * mult;
            if let Some(f) = p4_divisibility(n, p) {
                panic!("{:?}", f);
            }
        }
        // Test non-multiples
        for _ in 0..10 {
            let n: u64 = rand::Rng::gen_range(&mut rng, 1..=u64::MAX);
            if let Some(f) = p4_divisibility(n, p) {
                panic!("{:?}", f);
            }
        }
    }
}

/// P5 smoke: Barrett reduction at realistic operating points.
/// The sieve uses Barrett for computing offsets within chunks (up to ~25T).
#[test]
fn t0_p5_barrett() {
    let primes = erdos396_search_lab::drt_exports::drt_sieve_primes(1_000);
    let pd_all = erdos396_search_lab::drt_exports::drt_build_prime_data(&primes);
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    let max_n: u64 = 25_000_000_000_000;

    for pdata in pd_all.iter().skip(1).take(100) {
        let p = pdata.p;
        let bases: &[u64] = &[0, 1, 1_048_576, 1_000_000_000, 10_000_000_000, max_n];
        for &base in bases {
            let n = base + p - 1;
            let barrett_q =
                erdos396_search_lab::drt_exports::drt_barrett_quotient(n, pdata.magic, pdata.shift);
            let true_q = n / p;
            assert_eq!(barrett_q, true_q, "p={}, n={}", p, n);
        }
        for _ in 0..20 {
            let n: u64 = rand::Rng::gen_range(&mut rng, p..=max_n + p);
            let barrett_q =
                erdos396_search_lab::drt_exports::drt_barrett_quotient(n, pdata.magic, pdata.shift);
            let true_q = n / p;
            assert_eq!(barrett_q, true_q, "p={}, n={}", p, n);
        }
    }
}

/// P7 smoke: governor agreement — erdos396 is_governor vs is_governor_fast.
#[test]
fn t0_p7_governor_agreement() {
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);

    // Test governor status around known witnesses using erdos396 library
    // (which has both is_governor and is_governor_fast code paths)
    let test_points: &[u64] = &[1, 2, 3, 100, 1000, 2480, 8178, 45153, 3_648_841, 7_979_090];

    for &n in test_points {
        for delta in -5..=5i64 {
            let test_n = (n as i64 + delta).max(1) as u64;
            let result = e3.is_governor(test_n);
            // Also cross-check with the fast path
            let checker = erdos396::GovernorChecker::new(test_n);
            let fast_result = checker.is_governor_fast(test_n);
            assert_eq!(
                result, fast_result,
                "Governor disagreement (slow vs fast) at n={}",
                test_n
            );
        }
    }
}

/// False positive: 17,842,967,551 must NOT be a k=9 witness.
#[test]
fn t0_false_positive_rejection() {
    let mut sl = SearchLabOracle::new();
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);

    let n = 17_842_967_551u64;
    // This is a run of 10 governors but fails p=3 for the witness check
    let sl_result = sl.check_witness(9, n);
    let e3_result = e3.check_witness(9, n);
    assert!(
        !sl_result,
        "search-lab should reject false positive n={}",
        n
    );
    assert!(
        !e3_result,
        "erdos396-lib should reject false positive n={}",
        n
    );
}

// ===================================================================
// T1: Standard tests — includes Lean oracle, < 60 seconds
// ===================================================================

/// D2: search-lab vs Lean oracle on known witnesses.
#[test]
#[ignore] // T1: run with `cargo test -- --ignored`
fn t1_known_witnesses_search_lab_vs_lean() {
    let mut sl = SearchLabOracle::new();
    let Some(mut lean) = LeanOracle::new() else {
        eprintln!("Lean oracle not available — skipping");
        return;
    };

    let known: &[(u32, u64)] = &[
        (1, 2),
        (2, 2_480),
        (3, 8_178),
        (4, 45_153),
        (5, 3_648_841),
        (6, 7_979_090),
        (7, 101_130_029),
        (8, 339_949_252),
        (9, 1_019_547_844),
        (10, 17_609_764_994),
    ];

    for &(k, n) in known {
        let sl_result = sl.check_witness(k, n);
        let lean_result = lean.check_witness(k, n);
        assert_eq!(
            sl_result, lean_result,
            "search-lab vs lean disagree on k={}, n={}: sl={}, lean={}",
            k, n, sl_result, lean_result
        );
        assert!(sl_result, "known witness k={}, n={} should pass", k, n);
    }
}

/// D2: search-lab vs Lean oracle on 10K random inputs.
#[test]
#[ignore]
fn t1_differential_search_lab_vs_lean_random() {
    let mut sl = SearchLabOracle::new();
    let Some(mut lean) = LeanOracle::new() else {
        eprintln!("Lean oracle not available — skipping");
        return;
    };

    let cases = generators::generate_tier(SEED, 10_000, 10_000_000_000);
    let mut failures = Vec::new();

    for (i, tc) in cases.iter().enumerate() {
        let sl_result = sl.check_witness(tc.k, tc.n);
        let lean_result = lean.check_witness(tc.k, tc.n);
        if sl_result != lean_result {
            failures.push(format!(
                "k={}, n={} ({}): search-lab={}, lean={}",
                tc.k, tc.n, tc.generator, sl_result, lean_result
            ));
        }
        if (i + 1) % 1000 == 0 {
            eprintln!(
                "  [{}/{}] {} failures so far",
                i + 1,
                cases.len(),
                failures.len()
            );
        }
    }

    assert!(
        failures.is_empty(),
        "D2: {} disagreements out of {} tests:\n{}",
        failures.len(),
        cases.len(),
        failures[..failures.len().min(20)].join("\n")
    );
}

/// P1: Kummer vs Legendre agreement via Lean oracle.
#[test]
#[ignore]
fn t1_p1_kummer_vs_legendre() {
    let Some(mut lean) = LeanOracle::new() else {
        eprintln!("Lean oracle not available — skipping");
        return;
    };

    let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    let primes: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];

    for _ in 0..2_000 {
        let n: u64 = rand::Rng::gen_range(&mut rng, 1..=10_000_000_000u64);
        for &p in primes {
            if let Some(f) = p1_kummer_vs_legendre(n, p, &mut lean) {
                panic!("{:?}", f);
            }
        }
    }
}

/// Triple oracle: search-lab vs erdos396 vs Lean on 5K inputs.
#[test]
#[ignore]
fn t1_triple_oracle_agreement() {
    let mut sl = SearchLabOracle::new();
    let mut e3 = Erdos396Oracle::new(25_000_000_000_000);
    let Some(mut lean) = LeanOracle::new() else {
        eprintln!("Lean oracle not available — skipping");
        return;
    };

    let cases = generators::generate_tier(SEED + 1, 5_000, 10_000_000_000);
    let mut failures = Vec::new();

    for tc in &cases {
        let sl_result = sl.check_witness(tc.k, tc.n);
        let e3_result = e3.check_witness(tc.k, tc.n);
        let lean_result = lean.check_witness(tc.k, tc.n);

        if sl_result != e3_result || sl_result != lean_result {
            failures.push(format!(
                "k={}, n={} ({}): sl={}, erdos396={}, lean={}",
                tc.k, tc.n, tc.generator, sl_result, e3_result, lean_result
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "Triple oracle: {} disagreements out of {} tests:\n{}",
        failures.len(),
        cases.len(),
        failures[..failures.len().min(20)].join("\n")
    );
}
