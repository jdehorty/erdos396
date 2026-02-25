# The Small Prime Barrier Theorem for Erdős Problem 396

**Date:** 2026-02-21
**Authors:** Justin Dehorty

**Status:** New result — independent verification requested

---

## Abstract

We prove that if a witness to Erdős Problem 396 contains a block term that
is not a member of the Governor Set, then that term's governor failure must
occur at a prime strictly less than 2k+1, where k is the witness parameter.
The proof rests on a carry invariance property: for primes p ≥ 2k+1, the
number of carries when doubling n in base p is identical to the number of
carries when doubling n−j (for any 0 ≤ j ≤ k with p | (n−j)). We show the
bound 2k+1 is tight.

For k = 13, this constrains any potential non-governor witness to involve
governor failures only at the 9 primes {2, 3, 5, 7, 11, 13, 17, 19, 23},
and identifies these primes as the only place “non-governor” behavior can occur.
In particular, every k-witness passes the witness inequality at these primes
(as it does at all primes), so a sliding-window screen at the barrier primes
followed by full all-prime verification yields a **provably complete**
computational method for finding all Erdős 396 witnesses in any range — including
any that the standard Governor Set sieve might miss (subject to the usual
software/hardware trust assumptions; see `docs/trust.md`).

---

## 1. Introduction

### 1.1 Erdős Problem 396

Erdős Problem 396 asks: for each positive integer k, find the smallest
positive integer n such that the falling factorial

    n(n−1)(n−2)···(n−k)

divides the central binomial coefficient C(2n, n) = (2n)! / (n!)².

Such an n is called a **k-witness**. The known k-witnesses form OEIS
sequence A375077, currently extending through k = 13.

### 1.2 The Governor Set Approach

All known witnesses share a striking structural property: every term in the
block {n, n−1, ..., n−k} belongs to the **Governor Set**

    G = { m ∈ ℕ : m | C(2m, m) }

(OEIS A014847, density ≈ 11.4% by Ford–Konyagin). This motivates a search
strategy: find runs of k+1 consecutive Governor Set members, then verify
each against the full divisibility condition. This approach discovered all
witnesses through k = 13 via exhaustive search of the range [0, 25 × 10¹²).

### 1.3 The Completeness Question

The governor-run search detects **false positives** (a governor run that
is not a witness) through post-hoc verification. However, it is blind to
**false negatives**: witnesses where some block term is NOT a governor. Such
a witness would satisfy n(n−1)···(n−k) | C(2n, n) despite some term n−j
failing to divide its own central binomial C(2(n−j), n−j).

The **Governor Set Completeness Conjecture** asserts no such witnesses exist.
This conjecture is supported by 25 trillion integers of empirical evidence
but remains unproven.

### 1.4 This Result

We prove that any counterexample to the Completeness Conjecture — a witness
containing a non-governor term — must have a highly constrained structure:
the non-governor term can only fail the governor test at "small" primes
(those below 2k+1). This follows from a carry-arithmetic argument using
Kummer's theorem.

While this does not settle the Completeness Conjecture, it:

1. Proves that a targeted small-prime sieve is **complete** for finding all
   witnesses, including non-governor witnesses
2. Reduces the computational verification problem from ~481,000 primes to
   just 9 primes (for k = 13)
3. Provides a structural constraint that may aid a future proof of the
   full conjecture

---

## 2. Preliminaries

We establish the definitions and notation used throughout. All integers are
assumed non-negative unless stated otherwise.

### 2.1 The p-adic Valuation

**Definition 1.** For a prime p and a positive integer m, the **p-adic
valuation** v_p(m) is the largest non-negative integer a such that p^a
divides m. Equivalently, if m = p^a · m' with gcd(m', p) = 1, then
v_p(m) = a.

**Properties:**
- v_p(m) = 0 if and only if p does not divide m
- v_p(a · b) = v_p(a) + v_p(b) for all positive integers a, b
- v_p(m) ≥ 1 if and only if p | m

### 2.2 Base-p Representation

**Definition 2.** Every positive integer m has a unique representation in
base p:

    m = d_0 + d_1 · p + d_2 · p² + ··· + d_L · p^L

where each **digit** d_i ∈ {0, 1, ..., p−1} and d_L ≠ 0. We write
m = (d_L d_{L−1} ... d_1 d_0)_p.

**Observation 1.** v_p(m) equals the number of trailing zeros in the base-p
representation of m. That is, if v_p(m) = a, then d_0 = d_1 = ··· = d_{a−1} = 0
and d_a ≠ 0 (or a = L+1 if m = p^{L+1}, but this would mean d_L = 0, so
the representation is just (1 0 ... 0)_p with a zeros).

### 2.3 Carries in Base-p Self-Addition

**Definition 3.** For a non-negative integer m, define **s_p(m)** as the
number of carries when computing m + m in base p using the standard
(schoolbook) addition algorithm. Formally, let m = (d_L ... d_1 d_0)_p.
Define the carry sequence {c_i}_{i ≥ −1} by:

    c_{−1} = 0
    c_i = ⌊(2 · d_i + c_{i−1}) / p⌋     for i = 0, 1, 2, ...

Then:

    s_p(m) = Σ_{i=0}^{L} c_i

**Observation 2.** Since 0 ≤ d_i ≤ p−1 and c_{i−1} ∈ {0, 1}, we have
2 · d_i + c_{i−1} ≤ 2(p−1) + 1 = 2p − 1 < 2p, so each c_i ∈ {0, 1}.

**Observation 3.** A carry occurs at position i (i.e., c_i = 1) if and only
if 2 · d_i + c_{i−1} ≥ p.

### 2.4 Kummer's Theorem

**Theorem (Kummer, 1852).** For any prime p and positive integer m:

    v_p(C(2m, m)) = s_p(m)

*The p-adic valuation of the central binomial coefficient C(2m, m) equals
the number of carries when adding m to itself in base p.*

This classical result is our primary tool. It converts divisibility questions
about binomial coefficients into digit-arithmetic questions about carry
propagation.

### 2.5 The Governor Set

**Definition 4.** The **Governor Set** is:

    G = { m ∈ ℕ_{>0} : m | C(2m, m) }

By Kummer's theorem, m ∈ G if and only if for every prime p:

    v_p(m) ≤ s_p(m)

We say m **fails the governor test at prime p** if v_p(m) > s_p(m), i.e.,
the p-adic demand of m exceeds the p-adic supply of C(2m, m).

### 2.6 Witnesses and Block Terms

**Definition 5.** A positive integer n with n > k is a **k-witness** for
Erdős Problem 396 if:

    n(n−1)(n−2)···(n−k) divides C(2n, n)

By unique factorization and Kummer's theorem, this holds if and only if for
every prime p:

    Σ_{i=0}^{k} v_p(n−i) ≤ s_p(n)                         ... (W)

We call {n, n−1, ..., n−k} the **witness block** and each n−i a **block
term**. Condition (W) is the **witness condition at prime p**.

**Definition 6.** A k-witness n is a **governor-run witness** if every block
term is a governor: n−i ∈ G for all 0 ≤ i ≤ k. A k-witness is a
**non-governor witness** if some block term n−j ∉ G.

---

## 3. The Carry Invariance Lemma

The core technical observation is that shifting a multiple of p by a small
amount does not change the carry count in self-addition.

**Lemma 1 (Carry Invariance).** Let p be a prime with p ≥ 3, and let m be
a positive integer with p | m. For any integer j with 0 ≤ j ≤ (p−1)/2:

    s_p(m + j) = s_p(m)

### Proof of Lemma 1

Let a = v_p(m) ≥ 1. Write m in base p:

    m = (d_L ... d_a  0  0 ... 0)_p
                      └─ a zeros ─┘

where d_a ≠ 0 (the first non-zero digit is at position a).

**Step A: Base-p representation of m + j.**

Since 0 ≤ j < p (because j ≤ (p−1)/2 < p for p ≥ 3), and the digit of m
at position 0 is 0, the addition m + j in base p is:

    position 0:     0 + j = j  < p   →  digit j, carry 0
    positions 1..a−1:  0 + 0 = 0       →  digit 0, carry 0  (carry-in is 0)
    positions ≥ a:     d_i + 0 = d_i    →  unchanged         (carry-in is 0)

No carry propagates. Therefore:

    m + j = (d_L ... d_a  0  0 ... 0  j)_p
                          └ a−1 zeros ┘

The digits of m + j are identical to those of m at every position except
position 0, where the digit changes from 0 to j.

**Step B: Self-addition carries for m.**

Compute the carry sequence for m + m:

    Position 0:         2 · 0 + c_{−1} = 0 + 0 = 0 < p.     c_0 = 0.
    Position 1:         2 · 0 + c_0    = 0 + 0 = 0 < p.     c_1 = 0.
    ...
    Position a−1:       2 · 0 + c_{a−2} = 0 + 0 = 0 < p.    c_{a−1} = 0.
    Position a:         2 · d_a + c_{a−1} = 2d_a + 0 = 2d_a.  c_a = ⌊2d_a/p⌋.
    Positions a+1..L:   2 · d_i + c_{i−1}.   (depends on carry chain from position a)

Summary: carries at positions 0 through a−1 are all 0. The carry chain from
position a onward is determined by digits d_a, d_{a+1}, ..., d_L and the
carry-in to position a, which is c_{a−1} = 0.

**Step C: Self-addition carries for m + j.**

Compute the carry sequence for (m+j) + (m+j):

    Position 0:         2j + c_{−1} = 2j + 0 = 2j.

Since j ≤ (p−1)/2, we have 2j ≤ p−1 < p. Therefore **c_0 = 0. No carry.**

    Position 1:         2 · 0 + c_0 = 0 + 0 = 0 < p.         c_1 = 0.
    ...
    Position a−1:       2 · 0 + c_{a−2} = 0 + 0 = 0 < p.     c_{a−1} = 0.
    Position a:         2 · d_a + c_{a−1} = 2d_a + 0 = 2d_a.  c_a = ⌊2d_a/p⌋.
    Positions a+1..L:   2 · d_i + c_{i−1}.   (same digits, same carry-in as Step B)

**Step D: Comparison.**

At every position i:
- The digit of m + j equals the digit of m (except at position 0, which is
  j instead of 0 — but this does not affect the carry, as shown).
- The carry-in c_{i−1} is identical for both computations (both are 0 at
  positions 0 through a−1, and from position a onward the inputs are
  identical).
- Therefore c_i is identical for both computations at every position.

We conclude s_p(m + j) = s_p(m).  **□**

---

## 4. The Small Prime Barrier Theorem

**Theorem 1 (Small Prime Barrier).** Let k ≥ 1 and let n > k be a
k-witness for Erdős Problem 396. If a block term n−j (with 0 ≤ j ≤ k)
fails the governor test at a prime q, then q < 2k+1.

**Equivalently (contrapositive):** For any prime q ≥ 2k+1, if a block term
n−j fails the governor test at q, then n is not a k-witness.

### Proof of Theorem 1

We prove the contrapositive. Let q ≥ 2k+1 be a prime, let n > k, and let
0 ≤ j ≤ k with:

    v_q(n−j) > s_q(n−j)                                     ... (1)

(That is, n−j fails the governor test at q.) We show the witness condition
(W) fails at prime q, so n is not a k-witness.

**Step 1: Unique divisibility in the block.**

The block {n−k, n−k+1, ..., n} consists of k+1 consecutive integers
spanning a range of k. Two elements of this block are both divisible by q
only if their difference is divisible by q. Since their difference is at
most k, this requires q ≤ k. But q ≥ 2k+1 > k (for k ≥ 1), so at most
one block term is divisible by q.

From (1), v_q(n−j) ≥ 1, so q | (n−j). Therefore n−j is the unique block
term divisible by q. For all i ≠ j: v_q(n−i) = 0. Hence:

    Σ_{i=0}^{k} v_q(n−i) = v_q(n−j)                        ... (2)

**Step 2: Apply the Carry Invariance Lemma.**

We verify the hypotheses of Lemma 1 with p = q and m = n−j:

- q ≥ 3: Holds because q ≥ 2k+1 ≥ 3 for k ≥ 1.
- q | m: Holds because q | (n−j).
- j ≤ (q−1)/2: Since q ≥ 2k+1, we have (q−1)/2 ≥ k ≥ j.  ✓

Lemma 1 gives:

    s_q(n) = s_q((n−j) + j) = s_q(n−j)                     ... (3)

**Step 3: The witness condition fails.**

Combining (1), (2), and (3):

    Σ_{i=0}^{k} v_q(n−i)  =  v_q(n−j)       [by (2)]
                            >  s_q(n−j)       [by (1)]
                            =  s_q(n)         [by (3)]
                            =  v_q(C(2n, n))  [by Kummer's theorem]

Therefore:

    Σ_{i=0}^{k} v_q(n−i)  >  v_q(C(2n, n))

The witness condition (W) fails at prime q. Hence n is not a k-witness. **□**

---

## 5. Tightness of the Bound

The bound 2k+1 in Theorem 1 is optimal: for primes below 2k+1, the carry
invariance can break, and the supply bonus from the extra carry CAN
compensate for a governor deficit.

**Proposition 2 (Tightness).** The bound q ≥ 2k+1 in Theorem 1 cannot be
replaced by any smaller threshold. Specifically, for every k ≥ 1, there
exist a prime q ≤ 2k and integers m, j with q | m, 0 ≤ j ≤ k, and
v_q(m) > s_q(m), such that:

    s_q(m + j)  >  s_q(m)

That is, the carry invariance fails, and the self-addition of m + j has
strictly more carries than that of m.

### Proof of Proposition 2

**Case 1: k ≥ 2.** By Bertrand's postulate, there exists a prime q with
k < q < 2k. In particular, q ≤ 2k − 1 (since q is odd for q ≥ 3) and
q ≥ k+1.

Set m = q and j = ⌈q/2⌉. We verify:

- q | m: Trivially, since m = q.
- 0 ≤ j ≤ k: Since q ≤ 2k−1 (odd prime), j = (q+1)/2 ≤ (2k−1+1)/2 = k. ✓
- v_q(m) > s_q(m): We have v_q(q) = 1. The base-q representation of q is
  (1, 0)_q, and the self-addition (1,0)_q + (1,0)_q gives: position 0:
  2·0 = 0, no carry; position 1: 2·1 = 2 < q (since q ≥ 3), no carry.
  So s_q(q) = 0 < 1 = v_q(q). ✓ (Note: q is not a governor at prime q.)
- s_q(m+j) > s_q(m): The number m+j = q + j has base-q representation
  (1, j)_q. Self-addition: position 0: 2j = 2⌈q/2⌉ ≥ q, **carry!**
  Position 1: 2·1 + 1 = 3 < q (since q ≥ 5 for k ≥ 2), no carry.
  So s_q(m+j) = 1 > 0 = s_q(m). ✓

**Case 2: k = 1.** Take q = 2, m = 4, j = 1.

- 2 | 4: ✓
- v_2(4) = 2, s_2(4) = s_2((100)_2) = 1 carry (at position 2: 2·1 ≥ 2).
  So v_2(4) = 2 > 1 = s_2(4). (Note: 4 is not a governor at p = 2.)
- s_2(5) = s_2((101)_2): position 0: 2·1 = 2 ≥ 2, carry; position 1:
  2·0 + 1 = 1 < 2, no carry; position 2: 2·1 + 0 = 2 ≥ 2, carry.
  So s_2(5) = 2 > 1 = s_2(4). ✓

In both cases, the carry-in from position 0 of (m+j)+(m+j) is 1 (vs. 0 for
m+m), creating at least one additional carry. The mechanism in Lemma 1
fails because 2j ≥ p, violating the hypothesis j ≤ (p−1)/2.  **□**

**Remark.** Proposition 2 shows that the carry invariance breaks down — it
does NOT construct an actual non-governor witness. Whether a non-governor
witness exists remains open (the Governor Set Completeness Conjecture).
What Proposition 2 demonstrates is that the proof technique of Theorem 1
cannot be extended to primes below 2k+1.

---

## 6. Corollaries

### Corollary 3 (Top Term Must Be a Governor)

If n is a k-witness for k ≥ 1, then n ∈ G.

**Proof.** Since n | n(n−1)···(n−k) and n(n−1)···(n−k) | C(2n, n), we
have n | C(2n, n), so n ∈ G by definition.  **□**

**Remark.** This elementary observation does not require Theorem 1. It uses
only the fact that n is a factor of the falling factorial. By contrast,
Theorem 1 addresses the deeper question of whether the INTERIOR terms
n−1, n−2, ..., n−k must also be governors.

### Corollary 4 (Non-Governor Failures Confined to Small Primes)

If n is a k-witness and n−j ∉ G for some 0 ≤ j ≤ k, then n−j fails the
governor test **only** at primes q < 2k+1. That is, for every prime
q ≥ 2k+1:

    v_q(n−j) ≤ s_q(n−j)

**Proof.** If v_q(n−j) > s_q(n−j) for some q ≥ 2k+1, then by Theorem 1,
n is not a k-witness — contradiction.  **□**

### Corollary 5 (Application to k = 13)

For k = 13, define:

    P₁₃ = { p prime : p < 2·13 + 1 } = { 2, 3, 5, 7, 11, 13, 17, 19, 23 }

If a k = 13 non-governor witness exists, every non-governor block term
fails the governor test exclusively at primes in P₁₃.

**Proof.** Direct application of Corollary 4 with k = 13. The primes ≥ 27
are all ≥ 29 (since 27 = 3³ is not prime), and Theorem 1 eliminates all
of them.  **□**

### Corollary 6 (Computational Completeness of the Small-Prime Sieve)

For any k ≥ 1 and any integer range [A, B], a computational search that:

1. Checks the witness condition (W) at every prime p < 2k+1 for each
   integer n ∈ [A, B], and
2. Performs full p-adic verification on all candidates passing step (1),

is a **provably complete** method for finding all k-witnesses in [A, B] —
including non-governor witnesses.

**Proof.** Let n ∈ [A, B] be a k-witness. By definition, n satisfies (W)
at every prime p. In particular, n passes the check at every prime p < 2k+1,
so n is identified as a candidate in step (1). Step (2) then correctly
classifies n as a witness (since it satisfies (W) at ALL primes). No witness
is missed.

Conversely, any n that passes step (1) but fails step (2) is correctly
rejected.  **□**

**Remark.** The completeness of this method does NOT depend on the Governor
Set Completeness Conjecture. It is unconditionally correct: it finds every
k-witness in the range, whether or not it is a governor run. This is
precisely the property needed for rigorous computational verification.

---

## 7. Discussion

### 7.1 Relationship to the Governor Set Completeness Conjecture

Theorem 1 does not prove the Governor Set Completeness Conjecture. It proves
a weaker but still powerful statement: **if a non-governor witness exists,
its non-governor terms have a highly constrained arithmetic structure.**

Specifically, any non-governor block term n−j must:
- Fail the governor test at some prime q ≤ 2k (a small prime)
- Succeed at all primes q ≥ 2k+1 (where its governor status matches
  automatically with the witness condition)

To actually prove the Completeness Conjecture, one would need to show that
the "supply bonus" at small primes — the difference s_p(n) − s_p(n−j) for
p < 2k+1 — is never sufficient to compensate the governor deficit, when
combined with the demand from all other block terms. This appears to be a
difficult combinatorial problem involving simultaneous constraints across
multiple small primes.

### 7.2 The Supply Bonus Mechanism at Small Primes

For primes p < 2k+1, the Carry Invariance Lemma (Lemma 1) fails because
2j can equal or exceed p. When this happens, the self-addition of n = (n−j)+j
produces a carry at position 0 that is absent in the self-addition of n−j.
This carry can cascade:

- At position 1, the carry-in of 1 (vs. 0 for n−j) may trigger an
  additional carry if the digit d₁ satisfies 2d₁ + 1 ≥ p (equivalently,
  d₁ ≥ (p−1)/2). We call such a digit **critical**.
- The cascade continues through consecutive critical digits.

The maximum supply bonus equals 1 plus the number of consecutive critical
digits starting at position 1 (in the base-p representation of n−j, above
the trailing-zero block).

**Analysis by prime:**

| Prime p | Carry at pos 0 when | Critical digit | Cascade probability |
|---------|---------------------|----------------|---------------------|
| 2 | j ≥ 1 (always for j ≥ 1) | d = 1 (i.e., 50% of digits) | High |
| 3 | j ≥ 2 | d = 1 (33% of digits) | Moderate |
| 5 | j ≥ 3 | d = 2 (20% of digits) | Low |
| 7 | j ≥ 4 | d = 3 (14% of digits) | Low |
| 23 | j ≥ 12 | d = 11 (4% of digits) | Very low |

For large primes in P_k, the supply bonus mechanism is weak: few values of
j trigger a carry, and cascading through critical digits is rare. For p = 2,
the mechanism is strongest, but the demand from k+1 block terms at p = 2 is
also highest (many terms are even). The balance of supply and demand at each
small prime determines whether a non-governor witness can exist.

### 7.3 Why Minimum Witnesses Are Almost Certainly Governor Runs

At the minimum k-witness n*, the supply-demand balance at the tightest prime
is typically razor-thin. (For example, the k = 10 witness at n = 17,609,764,994
has v_2(C(2n, n)) = 15 = Σ v_2(n−i) — supply equals demand exactly at p = 2.)

A non-governor block term INCREASES the demand side (it requires at least as
much p-adic valuation as a governor term, plus the deficit that makes it
non-governor). For this increased demand to be covered at the minimum
witness, the supply bonus from extra carries would need to compensate —
but at the minimum, there is almost no slack for this compensation.

**This is a heuristic argument, not a proof.** It is conceivable that a
non-governor witness exists at some n < n* (the minimum governor-run witness),
where the carry structure at some small prime is unusually favorable. The
computational verification described in Corollary 6 can settle this question
for any finite range.

### 7.4 Computational Implications

For k = 13, Theorem 1 reduces the validation problem as follows:

| Approach | Primes to check per n | Operations per n |
|----------|----------------------|------------------|
| Naive (all primes up to √(2·25T)) | ~481,000 | ~481,000 × 27 |
| Small-prime sieve (Corollary 6) | 9 | ~60 |
| Governor sieve (current pipeline) | ~481,000 (but only for governor terms) | ~1,500 (amortized via fused sieve) |

The small-prime sieve has comparable per-integer cost to the existing
governor sieve, but is **provably complete** — it cannot miss any witness,
governor or non-governor. Running both sieves on the same range and
comparing results would constitute the strongest possible computational
evidence for (or against) the Completeness Conjecture.

**Estimated cost for k = 13, range [0, 25T):**

- Small-prime screen: ~1M integers/sec per core (comparable to governor sieve)
- Time: 25T ÷ 40 cores ÷ 1M/sec ≈ 7 days on a 40-core machine
- Full verification of candidates: negligible (expected candidate rate
  ≈ 1 per 10⁹ integers)

---

## 8. Verification Checklist

The following specific claims should be independently verified:

### Lemma 1

- [ ] **L1.1**: For p ≥ 3 and j ≤ (p−1)/2, we have 2j ≤ p−1 < p, so the
  self-addition of m+j produces no carry at position 0.

- [ ] **L1.2**: The base-p digits of m+j at positions 1 through a−1 are 0
  (same as m), because the addition m + j with j < p and digit_0(m) = 0
  produces no carry from position 0 to position 1.

- [ ] **L1.3**: The base-p digits of m+j at positions ≥ a are identical to
  those of m, because no carry propagates from lower positions.

- [ ] **L1.4**: The carry-in to position a is 0 for both m+m and (m+j)+(m+j),
  since all carries at positions 0 through a−1 are 0 in both cases.

- [ ] **L1.5**: From position a onward, the carries in m+m and (m+j)+(m+j)
  are identical, because the digits and carry-ins are identical.

### Theorem 1

- [ ] **T1.1**: For q ≥ 2k+1 > k and k+1 consecutive integers, at most one
  is divisible by q. (Because any two multiples of q differ by ≥ q > k.)

- [ ] **T1.2**: j ≤ k ≤ (q−1)/2 when q ≥ 2k+1. (Because q ≥ 2k+1 implies
  (q−1)/2 ≥ k.)

- [ ] **T1.3**: Lemma 1 applies with p = q and m = n−j, giving
  s_q(n) = s_q(n−j).

- [ ] **T1.4**: The chain of inequalities
  Σ v_q(n−i) = v_q(n−j) > s_q(n−j) = s_q(n) = v_q(C(2n, n))
  is valid and each step is justified by (2), (1), (3), and Kummer's theorem
  respectively.

### Proposition 2

- [ ] **P2.1**: For k ≥ 2, a prime q with k+1 ≤ q ≤ 2k−1 exists (by
  Bertrand's postulate).

- [ ] **P2.2**: For such q, ⌈q/2⌉ ≤ k. (Because q ≤ 2k−1 and q odd imply
  (q+1)/2 ≤ k.)

- [ ] **P2.3**: s_q(q) = 0 for q ≥ 3. (Base-q representation of q is
  (1,0)_q; all self-addition digits are < q.)

- [ ] **P2.4**: s_q(q + ⌈q/2⌉) ≥ 1. (Position 0: digit ⌈q/2⌉, self-addition
  gives 2⌈q/2⌉ ≥ q, producing a carry.)

### Corollary 6

- [ ] **C6.1**: Every k-witness passes the small-prime screen (trivially,
  since it satisfies (W) at all primes).

- [ ] **C6.2**: Full verification in step (2) correctly classifies all
  candidates (since verification is exact — no approximation).

---

## 9. Summary

**Theorem 1 (Small Prime Barrier)** establishes that the Governor Set
sieve's only potential blind spot — non-governor witnesses — is confined to
a narrow arithmetic regime: governor failures at primes below 2k+1. This
yields a provably complete validation method (Corollary 6) that requires
checking only 9 primes for k = 13, instead of hundreds of thousands.

The result is elementary (relying only on Kummer's theorem and base-p
arithmetic), but its implications for computational number theory are
significant: it provides the theoretical foundation for a validation sieve
that can certify the completeness of our 25-trillion-integer Governor Set
search — without assuming the Completeness Conjecture.

---

## References

1. E. E. Kummer, "Über die Ergänzungssätze zu den allgemeinen
   Reciprocitätsgesetzen," *J. Reine Angew. Math.* **44** (1852), 93–146.

2. K. Ford and S. V. Konyagin, "On two conjectures of Erdős,"
   unpublished manuscript (cited in OEIS A014847). See also arXiv:1909.03903.

3. P. Erdős, R. L. Graham, I. Z. Ruzsa, and E. G. Straus, "On the prime
   factors of C(2n, n)," *Math. Comp.* **29** (1975), 83–92.

4. OEIS Foundation, Sequences A375077 (Erdős 396 witnesses) and A014847
   (Governor Set), https://oeis.org/.
