//! Mathematical invariant property tests.
//!
//! These verify search-lab's primitives against independent implementations.
//! Each property is a mathematical identity that must hold for all valid inputs.

use crate::oracle::Oracle;

/// Failure details for a property test.
#[derive(Debug)]
pub struct PropertyFailure {
    pub property: &'static str,
    pub input: String,
    pub expected: String,
    pub actual: String,
}

// ---------------------------------------------------------------------------
// P1: Legendre self-consistency for v_p(C(2n, n))
// ---------------------------------------------------------------------------

/// `vp_central_binom(n, p) == vp_factorial(2n, p) - 2 * vp_factorial(n, p)`
///
/// Catches: vp_central_binom implementation bugs in the oracle.
pub fn p1_legendre_consistency(
    n: u64,
    p: u64,
    oracle: &mut impl Oracle,
) -> Option<PropertyFailure> {
    let combined = oracle.vp_central_binom(n, p);
    let vf_2n = oracle.vp_factorial(2 * n, p);
    let vf_n = oracle.vp_factorial(n, p);
    let expected = vf_2n - 2 * vf_n;
    if combined != expected {
        Some(PropertyFailure {
            property: "P1: vp_central_binom == vp_factorial(2n) - 2*vp_factorial(n)",
            input: format!("n={}, p={}", n, p),
            expected: format!("{}", expected),
            actual: format!("{}", combined),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P2: v_2(C(2n,n)) == popcount(n)
// ---------------------------------------------------------------------------

/// `v_2(C(2n, n)) == popcount(n)` (Kummer's theorem for p=2).
///
/// Catches: vp_central_binom bugs for p=2 in whichever oracle is passed.
pub fn p2_v2_popcount(n: u64, oracle: &mut impl Oracle) -> Option<PropertyFailure> {
    let supply = oracle.vp_central_binom(n, 2);
    let popcount = n.count_ones() as u64;
    if supply != popcount {
        Some(PropertyFailure {
            property: "P2: v_2(C(2n,n)) == popcount(n)",
            input: format!("n={}", n),
            expected: format!("{}", popcount),
            actual: format!("{}", supply),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P3: Modular inverse identity
// ---------------------------------------------------------------------------

/// `p * mod_inverse(p) == 1 (mod 2^64)` for all odd primes p.
///
/// Catches: modular inverse computation bugs.
pub fn p3_mod_inverse(p: u64) -> Option<PropertyFailure> {
    if p % 2 == 0 {
        return None;
    }
    let inv = erdos396_search_lab::drt_exports::drt_mod_inverse_u64(p);
    let product = p.wrapping_mul(inv);
    if product != 1 {
        Some(PropertyFailure {
            property: "P3: p * mod_inverse(p) == 1 (mod 2^64)",
            input: format!("p={}", p),
            expected: "1".to_string(),
            actual: format!("{}", product),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P4: Divisibility via modular inverse
// ---------------------------------------------------------------------------

/// `n * inv_p <= max_quot iff p divides n` for all (n, p).
///
/// Catches: divisibility test threshold bugs.
pub fn p4_divisibility(n: u64, p: u64) -> Option<PropertyFailure> {
    if p < 2 || p % 2 == 0 {
        return None;
    }
    let inv_p = erdos396_search_lab::drt_exports::drt_mod_inverse_u64(p);
    let max_quot = u64::MAX / p;
    let wrapped = n.wrapping_mul(inv_p);
    let inv_says_divisible = wrapped <= max_quot;
    let actually_divisible = n % p == 0;
    if inv_says_divisible != actually_divisible {
        Some(PropertyFailure {
            property: "P4: n*inv_p <= max_quot iff p | n",
            input: format!("n={}, p={}", n, p),
            expected: format!("divisible={}", actually_divisible),
            actual: format!("inv_says={}", inv_says_divisible),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P6: Factor stripping correctness
// ---------------------------------------------------------------------------

/// After `process_prime_dyn(p)` on a value divisible by p: result has no
/// factor of p, and `odd_n == result * p^(v_p(odd_n))`.
///
/// Catches: factor stripping bugs in the sieve hot loop.
pub fn p6_strip_correctness(n: u64, p: u64) -> Option<PropertyFailure> {
    if p < 3 || p % 2 == 0 || n == 0 {
        return None;
    }
    let odd_n = n >> n.trailing_zeros();
    // process_prime_dyn assumes it's called at p-stride offsets where the
    // value is divisible by p. Skip non-divisible inputs.
    if odd_n % p != 0 {
        return None;
    }

    let inv_p = erdos396_search_lab::drt_exports::drt_mod_inverse_u64(p);
    let max_quot = u64::MAX / p;

    let mut rem = vec![odd_n];
    let mut start_j = 0u64;
    erdos396_search_lab::drt_exports::drt_process_prime_dyn(
        p,
        inv_p,
        max_quot,
        &mut start_j,
        1,
        &mut rem,
    );
    let result = rem[0];

    // Verify: result has no factor of p
    if result % p == 0 && result != 0 {
        return Some(PropertyFailure {
            property: "P6: strip leaves no factor of p",
            input: format!("n={} (odd={}), p={}", n, odd_n, p),
            expected: format!("result % {} != 0", p),
            actual: format!("result={}", result),
        });
    }

    // Verify: odd_n == result * p^e for some e
    let mut temp = odd_n;
    while temp % p == 0 {
        temp /= p;
    }
    if temp != result && result != 0 {
        return Some(PropertyFailure {
            property: "P6: strip(n,p) == n / p^v_p(n)",
            input: format!("n={} (odd={}), p={}", n, odd_n, p),
            expected: format!("{}", temp),
            actual: format!("{}", result),
        });
    }

    None
}

// ---------------------------------------------------------------------------
// P7: Governor agreement
// ---------------------------------------------------------------------------

/// Governor agreement between two oracles (both must support `is_governor`).
///
/// Catches: governor membership computation disagreements between code paths.
pub fn p7_governor_agreement(
    n: u64,
    oracle_a: &mut impl Oracle,
    oracle_b: &mut impl Oracle,
) -> Option<PropertyFailure> {
    let a_result = oracle_a.is_governor(n);
    let b_result = oracle_b.is_governor(n);
    if a_result != b_result {
        Some(PropertyFailure {
            property: "P7: oracle_a governor == oracle_b governor",
            input: format!("n={}", n),
            expected: format!("oracle_a={}", a_result),
            actual: format!("oracle_b={}", b_result),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P8: Witness verification agreement
// ---------------------------------------------------------------------------

/// `oracle_a::check_witness(k, n) == oracle_b::check_witness(k, n)`.
///
/// Catches: any disagreement in the full witness verification path.
pub fn p8_witness_agreement(
    k: u32,
    n: u64,
    oracle_a: &mut impl Oracle,
    oracle_b: &mut impl Oracle,
) -> Option<PropertyFailure> {
    let a_result = oracle_a.check_witness(k, n);
    let b_result = oracle_b.check_witness(k, n);
    if a_result != b_result {
        Some(PropertyFailure {
            property: "P8: oracle_a witness == oracle_b witness",
            input: format!("k={}, n={}", k, n),
            expected: format!("oracle_a={}", a_result),
            actual: format!("oracle_b={}", b_result),
        })
    } else {
        None
    }
}
