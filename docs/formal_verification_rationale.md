# Formal Verification Rationale

This document explains why this project uses machine-checked Lean 4 proofs
alongside its Rust implementation, what the formal proofs cover, and why this
level of rigor is warranted for the claims we make.

---

## The problem

For each k, we search for the smallest n such that

    n(n−1)(n−2)⋯(n−k) ∣ C(2n, n).

Our computational results include two kinds of claims:

1. **Existence:** "n = 339,949,252 is a k=8 witness." This is easy to check
   independently — evaluate the divisibility condition at (k, n) with
   arbitrary-precision arithmetic. A single computation settles it.

2. **Minimality:** "No k=8 witness exists below 339,949,252." This requires
   checking every integer in a range. There is no single computation that
   settles it — it is an implicit universal quantifier over billions or
   trillions of integers.

Claim (2) is where formal verification matters.

---

## Why testing alone is insufficient

The search algorithm is a sieve: it processes every integer in a range, applies
a fast rejection filter, and fully verifies the candidates that survive. The
filter is designed to have no false negatives — every true witness must survive
the filter.

If the filter has a bug that incorrectly rejects some true witness, the search
will report "no witness found" in a range where one exists. This is a silent
failure. No test case will catch it unless the test happens to cover the exact
input where the bug manifests.

Our test suite (see `crates/search-lab/src/main.rs`) verifies correctness on
known witnesses for k=1 through k=11 and tests each algorithmic component in
isolation. This provides strong empirical confidence. But the search ranges for
k=12 and k=13 extend to 25 trillion integers. Testing a representative sample
is not the same as proving a property holds universally.

The core question is: **why is the filter safe?** That is a mathematical
question, and it has a mathematical answer.

---

## What we prove in Lean

The formal Lean 4 project in `formal/` proves the theorems that justify the
algorithmic structure of the search — specifically, why the filter used by our
`validate` binary is complete.

### 1. The witness criterion is equivalent to a prime-by-prime condition

A positive integer n is a k-witness if and only if, for every prime p,

    v_p(n(n−1)⋯(n−k)) ≤ v_p(C(2n, n))

where v_p denotes the p-adic valuation. This is a consequence of the
fundamental theorem of arithmetic: a product divides another if and only if it
does so at every prime. The Lean formalization defines `IsWitness` via the
divisibility condition and `supply` as v_p(C(2n, n)), then proves the
connection.

**Why this matters:** The entire search strategy rests on checking individual
primes rather than computing a product of potentially hundreds of digits. If
this equivalence were wrong, the search would be checking the wrong condition.

### 2. The Small Prime Barrier Theorem

> **Theorem.** If n is a k-witness and a block term (n−j) fails the governor
> test at a prime q, then q < 2k+1.

> **Corollary.** Restricting the validation screen to primes p < 2k+1 yields
> zero false negatives.

This is the theorem that makes the `validate` binary's completeness claim
possible. For k=13, instead of checking all primes up to √(2n) ≈ 6 million,
the screen need only check 9 primes: {2, 3, 5, 7, 11, 13, 17, 19, 23}. The
Lean proof covers:

- **Carry invariance** (`carry_invariance_supply`): For primes q ≥ 3 with
  q ∣ m, shifting m by j ≤ (q−1)/2 preserves v_q(C(2m, m)). Proved via the
  digit-sum characterization of Kummer's theorem.

- **Unique divisibility** (`at_most_one_divisible`): Among k+1 consecutive
  integers, at most one is divisible by a prime q > k.

- **The barrier theorem** (`small_prime_barrier`): Combining these, any
  governor-test failure at a large prime leads to a contradiction with the
  witness hypothesis.

- **Tightness** (`small_prime_barrier_tight`): The bound 2k+1 cannot be
  improved, using Bertrand's postulate.

**Why this matters:** Without this theorem, the `validate` binary would need
to check every prime up to √(2n) for every integer in the range — a factor of
~700,000× more work for k=13. The theorem is what makes exhaustive validation
computationally feasible, and a bug in the theorem would silently invalidate
every minimality claim.

### The theorem addresses a real phenomenon

The Small Prime Barrier Theorem is not merely theoretical insurance. Non-governor
witnesses — witnesses where some block term is *not* in the Governor Set — exist
and have been computationally verified:

| k | Smallest governor-run witness | Smallest non-governor witness | Ratio |
|---|---:|---:|---:|
| 2 | 2,480 | 40,664 | 16× |
| 3 | 8,178 | 1,441,378 | 176× |
| 4 | 45,153 | 2,366,563 | 52× |

At k=2, 137 non-governor witnesses exist below 10 million. At k=4, 7,198 exist
below 100 billion. In every verified case, the non-governor block term fails the
governor test at a barrier prime (p < 2k+1), exactly as the theorem predicts.
The governor sieve would miss all of these witnesses. Without the Small Prime
Barrier Theorem justifying the `validate` binary's barrier-prime screen, there
would be no proven method to find them.

The k=14 search (targeting runs of 15 consecutive governors across [25T, 200T))
provides additional confirmation at scale. Among thousands of governor-run
near-misses of length 14, all failures occur at primes above the 2k+1 = 29
threshold — precisely where the theorem guarantees compensation is impossible.
Two detailed case studies (at n ≈ 93T and n ≈ 188T) show the characteristic
failure mode: a single prime factor above the barrier threshold creates
irrecoverable demand, while all barrier-prime conditions are satisfied with
comfortable surplus.

The `carry-diagnostic` crate (`crates/carry-diagnostic/`) is the tool that
produced this dataset. It scans ranges for near-miss blocks (windows with
exactly one non-governor), analyzes carry compensation at each barrier prime,
and identifies candidates where all witness gaps are non-positive. For
reproduction instructions and the full dataset, see `docs/non_governor_witnesses.md`.

### 3. Algorithmic identities used by the implementation

The file `formal/SmallPrimeBarrier/Algorithm.lean` proves identities that
correspond directly to the arithmetic performed by the Rust code:

- **Kummer's theorem** for C(2n, n): v_p(C(2n, n)) equals the number of
  carries when adding n + n in base p. This justifies the carry-counting
  implementation used in the governor test.

- **Demand as a factorial difference:** v_p(n(n−1)⋯(n−k)) = v_p(n!) −
  v_p((n−k−1)!). This justifies computing the block demand via Legendre's
  formula.

- **Sliding-window update rule:** The demand sum at position n+1 equals the
  demand sum at position n, minus the term leaving the window, plus the term
  entering. This justifies the incremental update in the `validate` binary's
  inner loop.

**Why this matters:** These are not deep theorems, but they are the exact
equations implemented in hot loops that process trillions of integers. An error
in any of them would propagate silently through the entire search.

---

## What we do *not* prove in Lean

The Lean proofs verify the mathematical foundations. They do not verify:

- The Rust implementation itself (that the code faithfully implements the
  proven identities)
- The compiler, standard library, or hardware
- The absence of race conditions in the parallel search

This gap between theorem and implementation is the standard trust model for
computational mathematics. We document it explicitly in `docs/trust.md` and
mitigate it through:

- An independent Python witness verifier (`certificates/scripts/verify_witness.py`)
- Cross-checks between the governor sieve and the small-prime sieve
- Unit tests comparing Kummer carry-counting against Legendre's formula
- Coverage invariants (count, sum, XOR) in search and validation reports
- Optional runtime self-checks (`--self-check-samples`, `--audit-interval`)

---

## The verification stack, summarized

```
┌─────────────────────────────────────────────────────────┐
│  Mathematical claim                                     │
│  "No k-witness exists in [a, b]"                        │
└───────────────┬─────────────────────────────────────────┘
                │ justified by
┌───────────────▼─────────────────────────────────────────┐
│  Small Prime Barrier Theorem (Lean 4, machine-checked)  │
│  "Screening at primes p < 2k+1 has zero false           │
│   negatives for any k-witness"                          │
└───────────────┬─────────────────────────────────────────┘
                │ implemented by
┌───────────────▼─────────────────────────────────────────┐
│  validate binary (Rust)                                 │
│  Checks every integer in [a, b] against barrier primes, │
│  fully verifies screen-pass candidates                  │
└───────────────┬─────────────────────────────────────────┘
                │ verified by
┌───────────────▼─────────────────────────────────────────┐
│  Empirical safeguards                                   │
│  Coverage invariants, cross-checks, independent Python  │
│  verifier, unit tests, optional runtime audits          │
└─────────────────────────────────────────────────────────┘
```

The Lean layer is the load-bearing element. If the Small Prime Barrier Theorem
were wrong, the validate binary's screen would have false negatives, and our
minimality claims would be unsound — regardless of how many tests pass. The
machine-checked proof removes this risk entirely: the theorem is as trustworthy
as the Lean kernel itself.

---

## Future formalization targets

The current Lean project proves the Small Prime Barrier Theorem, which is the
load-bearing result for the `validate` binary. Three conjectures emerging from
the non-governor witness research are natural candidates for extending the
formal verification:

**Unit Deficit Conjecture.** In every verified non-governor witness (hundreds
of examples across k=2,3,4), the non-governor block term fails the governor
test at exactly one barrier prime p, and the deficit v_p(m) − v_p(C(2m, m))
is exactly 1. The exponent v_p(m) varies widely (1 through 10), but the gap
is always 1. A Lean proof would formalize the structural limitation of the
carry compensation mechanism: a supply bonus of exactly 1 (from a single
extra carry at digit position 0) is the common case, and a deficit of 2
would require a carry cascade that is too constrained to coexist with the
remaining witness conditions.

**Minimum Witnesses Are Governor Runs.** For every k through k=13, the minimum
k-witness has all block terms in the Governor Set. Non-governor witnesses exist
but are always larger — sometimes by orders of magnitude. A Lean proof would
show that at the minimum witness n*, the supply-demand balance at some barrier
prime is at exact equality, leaving no slack for a governor deficit to be
compensated by carry bonuses. This would formally establish that the governor
sieve is sufficient for finding minimum witnesses, even though it is not
complete in general.

**Governor Set Completeness Threshold.** Non-governor witnesses exist at
k=2,3,4 but appear to vanish at higher k. The 100-billion scan at k=5 found
zero candidates. A Lean proof would use the Chinese Remainder Theorem to show
that for k above some threshold K*, the digit constraints required for
simultaneous carry compensation at all barrier primes are contradictory — the
probability of joint compensation decays multiplicatively with |P_k| ≈ k/ln(k).

These conjectures are not required for the current minimality certificates
(which rest entirely on the proven Small Prime Barrier Theorem). They represent
opportunities to strengthen the formal foundations and potentially simplify
future searches by proving that certain failure modes cannot occur.

---

## For reviewers

To verify the formal proofs:

```bash
cd formal && lake build
```

This typechecks all Lean files against Mathlib with no `sorry` or `admit`.
The build requires [elan](https://github.com/leanprover/elan) and downloads
Lean 4.15.0 + Mathlib on first run (several minutes); subsequent builds are
cached.

For the informal mapping between Lean statements and Rust code, see
`docs/lean_rust_bridge.md`. For the full trust model, see `docs/trust.md`.
