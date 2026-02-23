# Validation Architecture

**Context:** k=13 search complete, 25 trillion integers of continuous coverage

---

## 1. False Positives vs. False Negatives

### False Positives (Governor Run ≠ True Witness)

A **false positive** is a run of k+1 consecutive Governor Set members where the
product n(n-1)...(n-k) does NOT divide C(2n, n).

The governor condition checks each term individually against its own central
binomial coefficient. The witness condition checks the combined product against
a single target at n. A false positive occurs when individual demands are each
within their own supply, but the combined demand overflows the supply at n for
some prime p.

False positives are rare (~0.3% for runs of length 10) and are caught by the
full p-adic verification step that follows every candidate detection. **The risk
to witness claims is zero.**

### False Negatives (Missed Witnesses)

A **false negative** would be a true witness where some block term is NOT a
Governor Set member — invisible to the governor-run sieve.

The Governor Set Completeness Conjecture asserts no such witnesses exist. This
conjecture is supported by 25 trillion integers of empirical evidence but remains
unproven.

The **Small Prime Barrier Theorem** (see `small_prime_barrier.md`) proves that
any non-governor witness must have governor failures confined to primes p < 2k+1.
For k=13, this means only 9 primes: {2, 3, 5, 7, 11, 13, 17, 19, 23}.

---

## 2. Three-Phase Validation

### Phase 1: Verify All Recorded Runs (False Positive Catalog)

For every run of length L ≥ 6 in the search logs (~12M entries across 25T),
compute the exact p-adic verification for k = L-1. Produces a complete false
positive catalog with failing primes, demand, and supply.

**Output:** Whether each run is a true k-witness or a false positive (and why).

### Phase 2: Non-Governor Witness Search (Completeness Check)

The `validate` binary implements this phase. For each integer n in a range, it
checks the witness condition at the **barrier primes** p < 2k+1 using a sliding
window over the block {n-k, ..., n}.

By the Small Prime Barrier Theorem (Corollary 6), this is a **provably complete**
method: any k-witness — governor run or not — satisfies the witness condition at
all primes, including all barrier primes. Therefore every witness passes the
small-prime screen and is detected.

The sliding window maintains per-prime demand sums incrementally. As n advances
by 1, one term leaves and one enters the window. For k=13 with 9 barrier primes,
this costs ~60 arithmetic operations per integer — comparable to the governor sieve.

Candidates passing the screen undergo full all-prime verification.

**How the `validate` binary works:**

```bash
# Prove no k=13 witnesses exist below the known minimum
cargo run --release --bin validate -- -k 13 --start 0 --end 18253129921842
```

| | `erdos396` (Governor Search) | `validate` (Small-Prime Sieve) |
|---|---|---|
| **Tests** | Is each n a Governor Set member? | Does n satisfy the witness condition? |
| **Skips** | Non-governors (immediately) | Nothing — checks every integer |
| **Filter** | Governor run detection | Barrier primes p < 2k+1 only |
| **Purpose** | Find witnesses efficiently | Prove none were missed |
| **Completeness** | Assumes Completeness Conjecture | Provably complete (Theorem 1) |

### Phase 3: Cross-Validation

Run both the governor sieve and the small-prime sieve on a representative
subrange and confirm they find identical witness sets. If the sets match:
strong evidence for the Completeness Conjecture. If the direct set has extras:
non-governor witnesses exist.

---

## 3. Certificate of Minimality

Running all three phases produces a **Certificate of Minimality** for each k:

1. **Governor Search** found the witness at position N.
2. **Small-Prime Sieve** proved no witness exists in [0, N).
3. **Cross-Validation** confirmed both methods agree.

Together, these constitute a proof that N is the smallest k-witness —
independent of the Governor Set Completeness Conjecture.

---

## 4. Theoretical Foundation

The completeness of Phase 2 rests on the Small Prime Barrier Theorem:

> **Theorem 1.** If n is a k-witness and a block term n−j fails the governor
> test at prime q, then q < 2k+1.

> **Corollary 6.** Checking the witness condition at primes p < 2k+1 for every
> integer in a range is a provably complete method for finding all k-witnesses.

The formal Lean 4 proof is in `formal/SmallPrimeBarrier/`. See
`docs/small_prime_barrier.md` for the full mathematical exposition.
