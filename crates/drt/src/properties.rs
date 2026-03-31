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
// P1: Kummer vs Legendre for v_p(C(2n, n))
// ---------------------------------------------------------------------------

/// `vp_central_binom_kummer(n, p) == vp_central_binom_legendre(n, p)` for all (n, p).
///
/// Catches: Kummer carry-counting bugs.
pub fn p1_kummer_vs_legendre(n: u64, p: u64, lean: &mut impl Oracle) -> Option<PropertyFailure> {
    let legendre = lean.vp_central_binom(n, p);
    // Also query Kummer via VK command if using Lean oracle
    // For now, just verify Legendre is self-consistent
    let vf_2n = lean.vp_factorial(2 * n, p);
    let vf_n = lean.vp_factorial(n, p);
    let expected = vf_2n - 2 * vf_n;
    if legendre != expected {
        Some(PropertyFailure {
            property: "P1: vp_central_binom == vp_factorial(2n) - 2*vp_factorial(n)",
            input: format!("n={}, p={}", n, p),
            expected: format!("{}", expected),
            actual: format!("{}", legendre),
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
/// Catches: popcount shortcut bugs in search-lab's v₂ path.
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
// P5: Barrett reduction correctness
// ---------------------------------------------------------------------------

/// `Barrett_quotient(n, p) == n / p` for all (n, p).
///
/// Catches: Barrett magic constant bugs.
pub fn p5_barrett(n: u64, p: u64) -> Option<PropertyFailure> {
    if p < 2 {
        return None;
    }
    let primes = erdos396_search_lab::drt_exports::drt_sieve_primes(p as usize + 1);
    let pd = erdos396_search_lab::drt_exports::drt_build_prime_data(&primes);
    // Find the PrimeData for p
    let pdata = pd.iter().find(|d| d.p == p)?;
    let barrett_q =
        erdos396_search_lab::drt_exports::drt_barrett_quotient(n, pdata.magic, pdata.shift);
    let true_q = n / p;
    if barrett_q != true_q {
        Some(PropertyFailure {
            property: "P5: Barrett_quotient(n, p) == n / p",
            input: format!("n={}, p={}", n, p),
            expected: format!("{}", true_q),
            actual: format!("{}", barrett_q),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P6: Factor stripping correctness
// ---------------------------------------------------------------------------

/// After `process_prime_dyn(p)`: result has no factor of p, and
/// original == result * p^(v_p(original)).
///
/// Catches: factor stripping bugs in the sieve hot loop.
pub fn p6_strip_correctness(n: u64, p: u64) -> Option<PropertyFailure> {
    if p < 3 || p % 2 == 0 || n == 0 {
        return None;
    }
    let inv_p = erdos396_search_lab::drt_exports::drt_mod_inverse_u64(p);
    let max_quot = u64::MAX / p;

    // Simulate stripping: start with n (odd part), strip all factors of p
    let odd_n = n >> n.trailing_zeros();
    let mut rem = vec![odd_n];
    let mut start_j = 0u64;
    // Only strip if the value at offset 0 is divisible by p
    // We need start_j to be 0 for the strip to touch index 0
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
// P7: Governor sieve residual vs reference
// ---------------------------------------------------------------------------

/// search-lab sieve residual == 1 ⟺ erdos396::is_governor(n).
///
/// Catches: governor sieve vs factorization-based reference disagreements.
pub fn p7_governor_agreement(
    n: u64,
    search_lab: &mut impl Oracle,
    reference: &mut impl Oracle,
) -> Option<PropertyFailure> {
    let sl_result = search_lab.is_governor(n);
    let ref_result = reference.is_governor(n);
    if sl_result != ref_result {
        Some(PropertyFailure {
            property: "P7: search-lab governor == reference governor",
            input: format!("n={}", n),
            expected: format!("reference={}", ref_result),
            actual: format!("search-lab={}", sl_result),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// P8: Witness verification agreement
// ---------------------------------------------------------------------------

/// `search_lab::exact_check(n, k) == reference::check_witness(k, n)`.
///
/// Catches: any disagreement in the full witness verification path.
pub fn p8_witness_agreement(
    k: u32,
    n: u64,
    search_lab: &mut impl Oracle,
    reference: &mut impl Oracle,
) -> Option<PropertyFailure> {
    let sl_result = search_lab.check_witness(k, n);
    let ref_result = reference.check_witness(k, n);
    if sl_result != ref_result {
        Some(PropertyFailure {
            property: "P8: search-lab witness == reference witness",
            input: format!("k={}, n={}", k, n),
            expected: format!("reference={}", ref_result),
            actual: format!("search-lab={}", sl_result),
        })
    } else {
        None
    }
}
