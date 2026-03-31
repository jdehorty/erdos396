//! Targeted input generators for differential random testing.
//!
//! Not uniform random — biased toward zones where bugs hide.

use rand::Rng;

/// Maximum n for DRT test cases. Both the Lean oracle (5M prime sieve, sufficient
/// for factoring up to 25T since sqrt(25T) ≈ 5M) and search-lab (uses `2*n` which
/// overflows u64 for n > 2^63) have range limitations. We cap at the actual search
/// range ceiling.
pub const MAX_SAFE_N: u64 = 25_000_000_000_000;

/// A test case for DRT.
#[derive(Debug, Clone)]
pub struct TestCase {
    pub k: u32,
    pub n: u64,
    pub generator: &'static str,
}

/// Generate uniform random (k, n) pairs.
pub fn uniform_random(rng: &mut impl Rng, count: usize, max_n: u64) -> Vec<TestCase> {
    let safe_max = max_n.min(MAX_SAFE_N);
    (0..count)
        .map(|_| {
            let k = rng.gen_range(1..=14u32);
            let n = rng.gen_range((k as u64 + 1)..=safe_max);
            TestCase {
                k,
                n,
                generator: "uniform",
            }
        })
        .collect()
}

/// Generate n near block boundaries (multiples of BLOCK_SIZE ± small delta).
pub fn block_boundary(rng: &mut impl Rng, count: usize) -> Vec<TestCase> {
    let block_size = 32768u64;
    let deltas: &[i64] = &[0, 1, -1, 2, -2, 13, -13, 14, -14];
    (0..count)
        .map(|_| {
            let k = rng.gen_range(1..=14u32);
            let block_idx = rng.gen_range(1..=1_000_000u64);
            let delta = deltas[rng.gen_range(0..deltas.len())];
            let n = (block_idx * block_size).wrapping_add(delta as u64);
            let n = n.max(k as u64 + 1);
            TestCase {
                k,
                n,
                generator: "block_boundary",
            }
        })
        .collect()
}

/// Generate n near chunk boundaries (multiples of CHUNK_SIZE ± small delta).
pub fn chunk_boundary(rng: &mut impl Rng, count: usize) -> Vec<TestCase> {
    let chunk_size = 1_048_576u64;
    let deltas: &[i64] = &[0, 1, -1, 2, -2, 13, -13, 14, -14];
    (0..count)
        .map(|_| {
            let k = rng.gen_range(1..=14u32);
            let chunk_idx = rng.gen_range(1..=25_000_000u64);
            let delta = deltas[rng.gen_range(0..deltas.len())];
            let n = (chunk_idx * chunk_size).wrapping_add(delta as u64);
            let n = n.max(k as u64 + 1);
            TestCase {
                k,
                n,
                generator: "chunk_boundary",
            }
        })
        .collect()
}

/// Generate n near known witnesses (± small delta for regression testing).
pub fn near_witnesses(rng: &mut impl Rng) -> Vec<TestCase> {
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
    let mut cases = Vec::new();
    for &(k, witness_n) in known {
        for delta in -100..=100i64 {
            let n = (witness_n as i64 + delta).max(k as i64 + 1) as u64;
            // Test at the known k and also at k-1, k+1
            cases.push(TestCase {
                k,
                n,
                generator: "near_witness",
            });
            if k > 1 {
                cases.push(TestCase {
                    k: k - 1,
                    n,
                    generator: "near_witness",
                });
            }
            if k < 14 {
                cases.push(TestCase {
                    k: k + 1,
                    n,
                    generator: "near_witness",
                });
            }
        }
    }
    // Shuffle to avoid correlated runs
    use rand::seq::SliceRandom;
    cases.shuffle(rng);
    cases
}

/// Generate smooth numbers (all prime factors ≤ 47).
pub fn smooth_numbers(rng: &mut impl Rng, count: usize) -> Vec<TestCase> {
    let barrier_primes: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
    (0..count)
        .map(|_| {
            let k = rng.gen_range(1..=14u32);
            // Build a smooth number from random exponents
            let mut n = 1u64;
            for _ in 0..rng.gen_range(2..=20) {
                let p = barrier_primes[rng.gen_range(0..barrier_primes.len())];
                n = n.saturating_mul(p);
                if n > MAX_SAFE_N {
                    n = MAX_SAFE_N;
                    break;
                }
            }
            let n = n.max(k as u64 + 1);
            TestCase {
                k,
                n,
                generator: "smooth",
            }
        })
        .collect()
}

/// Generate prime powers (n = p^e).
pub fn prime_powers(rng: &mut impl Rng, count: usize) -> Vec<TestCase> {
    let primes: &[u64] = &[
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 127, 131, 251, 509, 1021, 2039, 4093, 8191, 65537,
    ];
    (0..count)
        .map(|_| {
            let k = rng.gen_range(1..=14u32);
            let p = primes[rng.gen_range(0..primes.len())];
            let max_exp = (64.0 / (p as f64).log2()).floor() as u32;
            let exp = rng.gen_range(1..=max_exp.max(1));
            let mut n = 1u64;
            for _ in 0..exp {
                n = n.saturating_mul(p);
                if n > MAX_SAFE_N {
                    n = MAX_SAFE_N;
                    break;
                }
            }
            let n = n.max(k as u64 + 1);
            TestCase {
                k,
                n,
                generator: "prime_power",
            }
        })
        .collect()
}

/// Generate false positive cases — known governor runs that fail witness verification.
pub fn false_positives() -> Vec<TestCase> {
    // FALSE_POSITIVE_RUNS: (end_n, run_length, failing_prime)
    let fps: &[(u64, u32, u32)] = &[(17_842_967_551, 10, 3)];
    let mut cases = Vec::new();
    for &(end_n, run_len, _) in fps {
        // Test as a witness with k = run_len - 1
        for k in 1..=run_len {
            cases.push(TestCase {
                k,
                n: end_n,
                generator: "false_positive",
            });
        }
        // Also test nearby values
        for delta in -5..=5i64 {
            let n = (end_n as i64 + delta) as u64;
            cases.push(TestCase {
                k: run_len - 1,
                n,
                generator: "false_positive",
            });
        }
    }
    cases
}

/// Generate all test cases for a given tier.
pub fn generate_tier(seed: u64, count: usize, max_n: u64) -> Vec<TestCase> {
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);

    let mut all = Vec::new();

    // Deterministic cases first
    all.extend(near_witnesses(&mut rng));
    all.extend(false_positives());

    // Targeted generators (each gets a proportional share)
    let remaining = count.saturating_sub(all.len());
    let per_gen = remaining / 6;

    all.extend(block_boundary(&mut rng, per_gen));
    all.extend(chunk_boundary(&mut rng, per_gen));
    all.extend(smooth_numbers(&mut rng, per_gen));
    all.extend(prime_powers(&mut rng, per_gen));
    // Fill remainder with uniform random
    let uniform_count = count.saturating_sub(all.len());
    all.extend(uniform_random(&mut rng, uniform_count, max_n));

    all
}
