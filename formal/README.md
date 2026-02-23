# Small Prime Barrier Theorem — Lean 4 Formalization

Machine-checked proof of the Small Prime Barrier Theorem for Erdős Problem #396.

## What is proven

This project formalizes three results from the paper:

1. **Carry Invariance** (`carry_invariance_supply`): For primes q ≥ 3 with q ∣ m,
   shifting m by j ≤ (q−1)/2 preserves the p-adic supply v_q(C(2m, m)).

2. **Small Prime Barrier** (`small_prime_barrier`): If n is a k-witness and a block
   term n−j fails the governor test at a prime q ≥ 2k+1, then n cannot be a
   k-witness — contradiction. Equivalently, non-governor failures in a witness
   block are confined to primes below 2k+1.

3. **Tightness** (`small_prime_barrier_tight`): The bound 2k+1 is optimal. For every
   k ≥ 1, there exist parameters where the carry invariance breaks at a prime
   q ≤ 2k (using Bertrand's postulate).

The corollary `governor_test_passes_at_large_primes` gives the contrapositive form
used by the `validate` binary: every block term of a k-witness passes the governor
test at all primes q ≥ 2k+1.

## Building

Requires [elan](https://github.com/leanprover/elan) (the Lean version manager).

```bash
cd formal
lake build
```

The first build downloads Lean 4.15.0 and Mathlib, which takes several minutes.
Subsequent builds are cached and complete in seconds.

## Dependencies

- **Lean 4**: v4.15.0
- **Mathlib**: v4.15.0 (via `leanprover-community/mathlib`)

Key Mathlib imports:
- `Mathlib.Data.Nat.Choose.Central` — central binomial coefficients
- `Mathlib.NumberTheory.Padics.PadicVal` — p-adic valuations
- `Mathlib.NumberTheory.Bertrand` — Bertrand's postulate

## Structure

```
formal/
├── SmallPrimeBarrier/
│   ├── Defs.lean    — supply, IsWitness definitions
│   └── Main.lean    — all theorems and proofs
├── lakefile.toml
├── lean-toolchain
└── README.md
```

## Relationship to the paper

This formalization machine-checks the theoretical foundation (Theorem 1 and
Proposition 2) that guarantees the `validate` binary's completeness: checking
only primes p < 2k+1 is sufficient to find every k-witness in any range,
including non-governor witnesses invisible to the standard Governor Set sieve.
