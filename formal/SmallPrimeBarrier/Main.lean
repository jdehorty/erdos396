import SmallPrimeBarrier.Defs
import Mathlib.NumberTheory.Bertrand

/-!
# The Small Prime Barrier Theorem

*Reference:* `compute/docs/small_prime_barrier_theorem.md`

This file formalizes the **Small Prime Barrier Theorem** for Erdős Problem #396:
if a k-witness contains a block term that fails the governor test at some prime q,
then q < 2k + 1.

## Main Results

- `carry_invariance_supply`: For primes q ≥ 3 with q ∣ m, shifting m by a small
  amount j ≤ (q−1)/2 preserves the p-adic supply (Kummer carry count).
- `small_prime_barrier`: If n is a k-witness and n−j fails the governor test
  at prime q ≥ 2k+1, contradiction — n cannot be a k-witness.
- `governor_test_passes_at_large_primes`: Any block term of a k-witness passes
  the governor test at all primes q ≥ 2k+1.
- `small_prime_barrier_tight`: The bound 2k+1 is optimal.

## Proof Strategy

The carry invariance proof uses the digit-sum characterization of Kummer's theorem:
`(q−1) · v_q(C(2m, m)) = 2 · digitSum_q(m) − digitSum_q(2m)`. When q ∣ m, shifting
m by j ≤ (q−1)/2 only changes digit 0 (from 0 to j), adding j to each digit sum
on the right and 2j to the one on the left, which cancel.
-/

namespace Erdos396

open Nat

-- ============================================================
-- § 1. Digit Sum Shift Lemma
-- ============================================================

/-- When `q ∣ m`, the base-q digit sum of `m + j` equals the digit sum of `m`
    plus `j`, provided `j < q`. This is because `q ∣ m` means digit 0 of `m`
    is 0, and adding `j < q` just sets it to `j`. -/
private lemma digitSum_shift (q m j : ℕ) (hq1 : 1 < q) (hdvd : q ∣ m)
    (hm : 0 < m) (hj : j < q) :
    (Nat.digits q (m + j)).sum = (Nat.digits q m).sum + j := by
  obtain ⟨c, hc⟩ := hdvd
  subst hc
  have hc_pos : 0 < c := by
    rcases c with _ | c
    · simp at hm
    · omega
  -- digits q (q * c) = 0 :: digits q c  [since (q*c) % q = 0 and (q*c) / q = c]
  have hm_digits : Nat.digits q (q * c) = 0 :: Nat.digits q c := by
    rw [Nat.digits_def' hq1 (by omega : 0 < q * c)]
    congr 1
    · exact Nat.mul_mod_right q c
    · exact congrArg (Nat.digits q) (Nat.mul_div_cancel_left c (by omega : 0 < q))
  -- q * c + j = j + q * c
  -- digits q (j + q * c) = j :: digits q c  [by digits_add]
  have hmj_digits : Nat.digits q (q * c + j) = j :: Nat.digits q c := by
    rw [show q * c + j = j + q * c from by omega]
    exact Nat.digits_add q hq1 j c hj (Or.inr (by omega : c ≠ 0))
  rw [hm_digits, hmj_digits, List.sum_cons, List.sum_cons]
  omega

-- ============================================================
-- § 2. Carry Invariance (Lemma 1)
-- ============================================================

/--
**Carry Invariance Lemma (Supply Form).**

For a prime `q ≥ 3` with `q ∣ m` and `0 ≤ j ≤ (q − 1) / 2`:

    `supply q (m + j) = supply q m`

The proof uses the digit-sum characterization of Kummer's theorem:
`(q−1) · v_q(C(2x, x)) = 2 · Σdigits_q(x) − Σdigits_q(2x)`.

Since `q ∣ m`, digit 0 of `m` (base q) is 0. Adding `j` changes it to `j`,
increasing the digit sum by `j`. Similarly `q ∣ 2m`, and adding `2j` to `2m`
increases its digit sum by `2j`. The `2j` terms cancel:
`2(A+j) − (B+2j) = 2A − B`.
-/
theorem carry_invariance_supply (q m j : ℕ)
    (hq : Nat.Prime q) (hq3 : 3 ≤ q)
    (hdvd : q ∣ m) (hm : 0 < m)
    (hj : j ≤ (q - 1) / 2) :
    supply q (m + j) = supply q m := by
  haveI : Fact q.Prime := ⟨hq⟩
  have hq1 : 1 < q := hq.one_lt
  have hj_lt : j < q := by omega
  have h2j_lt : 2 * j < q := by omega
  -- Digit sum shifts: adding j to m increases digit sum by j
  have h_ds : (Nat.digits q (m + j)).sum = (Nat.digits q m).sum + j :=
    digitSum_shift q m j hq1 hdvd hm hj_lt
  -- Adding 2j to 2m increases digit sum by 2j (since q ∣ 2m)
  have h_ds2 : (Nat.digits q (2 * m + 2 * j)).sum =
      (Nat.digits q (2 * m)).sum + 2 * j :=
    digitSum_shift q (2 * m) (2 * j) hq1 (dvd_mul_of_dvd_right hdvd 2) (by omega) h2j_lt
  -- Kummer digit-sum formula for supply q m:
  -- (q-1) * v_q(C(2m, m)) = digitSum(m) + digitSum(m) - digitSum(2m)
  have h_m := sub_one_mul_padicValNat_choose_eq_sub_sum_digits (p := q)
    (k := m) (n := 2 * m) (by omega : m ≤ 2 * m)
  simp only [show 2 * m - m = m from by omega] at h_m
  -- Same formula for supply q (m+j):
  have h_mj := sub_one_mul_padicValNat_choose_eq_sub_sum_digits (p := q)
    (k := m + j) (n := 2 * (m + j)) (by omega : m + j ≤ 2 * (m + j))
  simp only [show 2 * (m + j) - (m + j) = m + j from by omega] at h_mj
  simp only [show 2 * (m + j) = 2 * m + 2 * j from by ring] at h_mj
  -- Cancel (q - 1): suffices (q-1) * supply q (m+j) = (q-1) * supply q m
  suffices hsuff : (q - 1) * supply q (m + j) = (q - 1) * supply q m by
    have : 0 < q - 1 := by omega
    exact mul_left_cancel₀ (by omega : (q : ℕ) - 1 ≠ 0) hsuff
  -- Unfold supply to padicValNat q (centralBinom _)
  simp only [supply, centralBinom]
  -- Now goal is about padicValNat q (choose (2*(m+j)) (m+j)) vs choose (2*m) m
  simp only [show 2 * (m + j) = 2 * m + 2 * j from by ring]
  -- Both sides equal digit-sum expressions that are equal after substitution
  rw [h_mj, h_m, h_ds, h_ds2]
  -- Goal: (A+j)+(A+j)-(B+2j) = A+A-B  where A = digitSum m, B = digitSum (2m)
  omega

-- ============================================================
-- § 3. Unique Prime Divisibility in Consecutive Blocks
-- ============================================================

/--
Among `k + 1` consecutive integers `{n − k, …, n}`, at most one is
divisible by a prime `q` with `q > k`.
-/
theorem at_most_one_divisible (q n k i₁ i₂ : ℕ)
    (_hq : Nat.Prime q) (hqk : k < q)
    (hn : k ≤ n)
    (hi₁ : i₁ ≤ k) (hi₂ : i₂ ≤ k)
    (h₁ : q ∣ (n - i₁)) (h₂ : q ∣ (n - i₂)) :
    i₁ = i₂ := by
  -- Two multiples of q in the block differ by at most k < q, so they coincide.
  rcases Nat.le_total i₁ i₂ with h_le | h_le
  · -- i₁ ≤ i₂: the difference i₂ - i₁ is divisible by q but < q, hence 0
    have hdiff : q ∣ (i₂ - i₁) := by
      have h := Nat.dvd_sub' h₁ h₂
      have heq : (n - i₁) - (n - i₂) = i₂ - i₁ := by omega
      rwa [heq] at h
    rcases Nat.eq_zero_or_pos (i₂ - i₁) with h0 | hpos
    · omega
    · exact absurd (Nat.le_of_dvd hpos hdiff) (by omega)
  · -- i₂ ≤ i₁: symmetric
    have hdiff : q ∣ (i₁ - i₂) := by
      have h := Nat.dvd_sub' h₂ h₁
      have heq : (n - i₂) - (n - i₁) = i₁ - i₂ := by omega
      rwa [heq] at h
    rcases Nat.eq_zero_or_pos (i₁ - i₂) with h0 | hpos
    · omega
    · exact absurd (Nat.le_of_dvd hpos hdiff) (by omega)

-- ============================================================
-- § 4. Helpers for the Main Theorem
-- ============================================================

/-- The descending factorial equals a product over `Finset.range`. -/
private theorem descFactorial_eq_prod_range (n k : ℕ) :
    descFactorial n k = ∏ i ∈ Finset.range k, (n - i) := by
  induction k with
  | zero => simp [descFactorial]
  | succ k ih => rw [descFactorial_succ, Finset.prod_range_succ, ih, mul_comm]

/-- The block term `n - j` divides `descFactorial n (k + 1)` when `j ≤ k`. -/
private theorem factor_dvd_descFactorial (n k j : ℕ) (hj : j ≤ k) :
    (n - j) ∣ descFactorial n (k + 1) := by
  rw [descFactorial_eq_prod_range]
  exact Finset.dvd_prod_of_mem _ (Finset.mem_range.mpr (by omega))

-- ============================================================
-- § 5. The Small Prime Barrier Theorem
-- ============================================================

/--
**Small Prime Barrier Theorem.**

If `n` is a k-witness for Erdős Problem 396, and a block term `n − j`
(with `j ≤ k`) fails the governor test at a prime `q ≥ 2k + 1`, then
we reach a contradiction — `n` cannot actually be a k-witness.

**Proof:**
1. By carry invariance: `supply q n = supply q (n − j)`.
2. By transitivity: `(n − j) ∣ descFactorial n (k+1) ∣ centralBinom n`.
3. By valuation monotonicity: `padicValNat q (n−j) ≤ supply q n`.
4. Chain: `supply q n = supply q (n−j) < padicValNat q (n−j) ≤ supply q n`.
   Contradiction.
-/
theorem small_prime_barrier (k n j q : ℕ)
    (hq : Nat.Prime q) (hq_large : 2 * k + 1 ≤ q)
    (hn : k < n) (hj : j ≤ k) (hjn : j ≤ n)
    (hdvd : q ∣ (n - j))
    (hfail : supply q (n - j) < padicValNat q (n - j))
    (hw : IsWitness k n) :
    False := by
  haveI : Fact q.Prime := ⟨hq⟩
  -- Step 1: supply q n = supply q (n - j)  [carry invariance]
  have h_supply_eq : supply q n = supply q (n - j) := by
    rcases Nat.eq_zero_or_pos j with hj0 | hj_pos
    · simp [hj0]
    · -- j ≥ 1 implies k ≥ 1 implies q ≥ 3
      conv_lhs => rw [show n = (n - j) + j from by omega]
      exact carry_invariance_supply q (n - j) j hq (by omega) hdvd (by omega) (by omega)
  -- Step 2: (n - j) ∣ centralBinom n  [via descFactorial and IsWitness]
  have h_nj_dvd : (n - j) ∣ centralBinom n :=
    dvd_trans (factor_dvd_descFactorial n k j hj) hw
  -- Step 3: padicValNat q (n - j) ≤ supply q n  [valuation monotonicity]
  have hcb_ne : centralBinom n ≠ 0 := Nat.centralBinom_ne_zero n
  have h_val_le : padicValNat q (n - j) ≤ supply q n := by
    -- p^v ∣ (n-j)  and  (n-j) ∣ centralBinom n  imply  p^v ∣ centralBinom n
    -- Hence v ≤ padicValNat q (centralBinom n) = supply q n
    show padicValNat q (n - j) ≤ padicValNat q (centralBinom n)
    have hpow_dvd : q ^ padicValNat q (n - j) ∣ centralBinom n :=
      dvd_trans pow_padicValNat_dvd h_nj_dvd
    exact (padicValNat_dvd_iff_le hcb_ne).1 hpow_dvd
  -- Step 4: contradiction
  -- supply q n = supply q (n-j) < padicValNat q (n-j) ≤ supply q n
  linarith

-- ============================================================
-- § 6. Corollaries
-- ============================================================

/--
**Corollary: Non-governor failures are confined to small primes.**

If `n` is a k-witness, then for every prime `q ≥ 2k + 1` and every
block index `j ≤ k` with `q ∣ (n − j)`, the block term `n − j` passes
the governor test at `q`:

    `padicValNat q (n − j) ≤ supply q (n − j)`

This is the contrapositive of `small_prime_barrier`.
-/
theorem governor_test_passes_at_large_primes (k n j q : ℕ)
    (hq : Nat.Prime q) (hq_large : 2 * k + 1 ≤ q)
    (hn : k < n) (hj : j ≤ k) (hjn : j ≤ n)
    (hdvd : q ∣ (n - j))
    (hw : IsWitness k n) :
    padicValNat q (n - j) ≤ supply q (n - j) := by
  by_contra h
  push_neg at h
  exact small_prime_barrier k n j q hq hq_large hn hj hjn hdvd h hw

-- ============================================================
-- § 7. Helpers for the Tightness Theorem
-- ============================================================

/-- Base-q digits of q itself: `[0, 1]`. -/
private lemma digits_self_base (q : ℕ) (hq1 : 1 < q) : Nat.digits q q = [0, 1] := by
  rw [Nat.digits_def' hq1 (by omega : 0 < q)]
  congr 1
  · exact Nat.mod_self q
  · rw [Nat.div_self (by omega : 0 < q)]
    exact Nat.digits_of_lt q 1 (by omega) hq1

/-- Base-q digits of `2 * q`: `[0, 2]`. Requires `q ≥ 3`. -/
private lemma digits_two_mul_self (q : ℕ) (hq1 : 1 < q) (hq3 : 3 ≤ q) :
    Nat.digits q (2 * q) = [0, 2] := by
  rw [Nat.digits_def' hq1 (by omega : 0 < 2 * q)]
  congr 1
  · rw [show 2 * q = q * 2 from by ring]; exact Nat.mul_mod_right q 2
  · rw [show 2 * q = q * 2 from by ring, Nat.mul_div_cancel_left 2 (by omega : 0 < q)]
    exact Nat.digits_of_lt q 2 (by omega) (by omega)

/-- `v_q(C(2q, q)) = 0` for primes `q ≥ 3`: the central binomial has no q-factor.
    Proof via Kummer digit-sum: digits of q in base q are `[0,1]` (sum=1),
    digits of 2q are `[0,2]` (sum=2), so `(q-1) · supply = 2·1 - 2 = 0`. -/
private lemma supply_self_eq_zero (q : ℕ) (hq : Nat.Prime q) (hq3 : 3 ≤ q) :
    supply q q = 0 := by
  haveI : Fact q.Prime := ⟨hq⟩
  have hq1 : 1 < q := hq.one_lt
  have h_kummer := sub_one_mul_padicValNat_choose_eq_sub_sum_digits (p := q)
    (k := q) (n := 2 * q) (by omega : q ≤ 2 * q)
  simp only [show 2 * q - q = q from by omega] at h_kummer
  rw [digits_self_base q hq1, digits_two_mul_self q hq1 hq3] at h_kummer
  simp only [List.sum_cons, List.sum_nil] at h_kummer
  -- h_kummer : (q - 1) * padicValNat q (Nat.choose (2 * q) q) = 0
  exact (Nat.mul_eq_zero.mp h_kummer).resolve_left (by omega)

-- ============================================================
-- § 8. Tightness of the Bound
-- ============================================================

/--
**Corollary: Tightness of the bound.**

The bound `2k + 1` in the Small Prime Barrier Theorem cannot be improved.
For every `k ≥ 1`, there exist a prime `q ≤ 2k` and integers `m, j` with
`q ∣ m`, `j ≤ k`, and `v_q(m) > supply q m`, such that
`supply q (m + j) > supply q m` — i.e., the carry invariance fails and
the supply bonus mechanism is active.

See `small_prime_barrier_theorem.md`, Proposition 2, for the construction
(using Bertrand's postulate to find q with k < q < 2k).
-/
theorem small_prime_barrier_tight :
    ∀ k : ℕ, 1 ≤ k →
      ∃ q m j : ℕ, Nat.Prime q ∧ q ≤ 2 * k ∧ q ∣ m ∧ j ≤ k ∧
        supply q m < padicValNat q m ∧
        supply q m < supply q (m + j) := by
  intro k hk
  rcases hk.eq_or_gt with hk1 | hk2
  · -- Case k = 1: witnesses q = 2, m = 4, j = 1
    subst hk1
    refine ⟨2, 4, 1, by decide, by omega, ⟨2, rfl⟩, le_refl 1, ?_, ?_⟩
    · -- supply 2 4 < padicValNat 2 4
      show padicValNat 2 (centralBinom 4) < padicValNat 2 4
      native_decide
    · -- supply 2 4 < supply 2 5
      show padicValNat 2 (centralBinom 4) < padicValNat 2 (centralBinom 5)
      native_decide
  · -- Case k ≥ 2: use Bertrand's postulate
    obtain ⟨q, hq_prime, hq_gt, hq_le⟩ :=
      Nat.exists_prime_lt_and_le_two_mul k (by omega)
    haveI : Fact q.Prime := ⟨hq_prime⟩
    have hq3 : 3 ≤ q := by omega
    have hq1 : 1 < q := hq_prime.one_lt
    set j := (q + 1) / 2 with hj_def
    refine ⟨q, q, j, hq_prime, hq_le, dvd_refl q, ?_, ?_, ?_⟩
    · -- j ≤ k: q is odd prime ≥ 3, so (q+1)/2 ≤ k
      have hq_odd_mod : q % 2 = 1 := by
        have : ¬ (q % 2 = 0) := by
          intro h0
          have h2dvd : 2 ∣ q := Nat.dvd_of_mod_eq_zero h0
          have := hq_prime.eq_one_or_self_of_dvd 2 h2dvd
          omega
        omega
      omega
    · -- supply q q = 0 < 1 = padicValNat q q
      rw [supply_self_eq_zero q hq_prime hq3, padicValNat_self]
      omega
    · -- supply q q = 0 < supply q (q + j)
      rw [supply_self_eq_zero q hq_prime hq3]
      suffices h : 1 ≤ supply q (q + j) by omega
      set n := q + j with hn_def
      -- Step 1: digits of n in base q = [j, 1]
      have hj_lt : j < q := by omega
      have hn_eq : n = j + q * 1 := by omega
      have h_digits_n : Nat.digits q n = [j, 1] := by
        rw [hn_eq, Nat.digits_add q hq1 j 1 hj_lt (Or.inr (by omega))]
        congr 1; exact Nat.digits_of_lt q 1 (by omega) hq1
      -- Step 2: Kummer formula for C(2n, n)
      have h_kummer := sub_one_mul_padicValNat_choose_eq_sub_sum_digits (p := q)
        (k := n) (n := 2 * n) (by omega : n ≤ 2 * n)
      simp only [show 2 * n - n = n from by omega] at h_kummer
      rw [h_digits_n] at h_kummer
      simp only [List.sum_cons, List.sum_nil] at h_kummer
      -- Step 3: Compute 2n = 3q + 1 (using q odd)
      have hq_odd_mod : q % 2 = 1 := by
        have : ¬ (q % 2 = 0) := by
          intro h0
          have := hq_prime.eq_one_or_self_of_dvd 2 (Nat.dvd_of_mod_eq_zero h0)
          omega
        omega
      have h2j_eq : 2 * j = q + 1 := by
        have := Nat.div_add_mod (q + 1) 2
        have : (q + 1) % 2 = 0 := by omega
        omega
      have h2n_eq : 2 * n = 3 * q + 1 := by omega
      -- Step 4: Case split q = 3 vs q ≥ 5
      rcases hq3.eq_or_gt with hq_eq_3 | hq_ge_4
      · -- q = 3: supply 3 5 ≥ 1 by direct computation
        subst hq_eq_3; simp only [hn_def, hj_def]
        show 1 ≤ supply 3 5
        native_decide
      · -- q ≥ 5: digits of 3q+1 in base q are [1, 3]
        have hq5 : 5 ≤ q := by omega
        have h_digits_2n : Nat.digits q (2 * n) = [1, 3] := by
          rw [h2n_eq, show 3 * q + 1 = 1 + q * 3 from by ring]
          rw [Nat.digits_add q hq1 1 3 (by omega) (Or.inr (by omega))]
          congr 1; exact Nat.digits_of_lt q 3 (by omega) (by omega)
        rw [h_digits_2n] at h_kummer
        simp only [List.sum_cons, List.sum_nil] at h_kummer
        -- (q-1) * supply = q - 1, so supply = 1
        have h_eq_1 : padicValNat q (Nat.choose (2 * n) n) = 1 :=
          mul_left_cancel₀ (show (q : ℕ) - 1 ≠ 0 by omega) (by omega)
        show 1 ≤ supply q n
        simp only [supply, centralBinom]
        omega

end Erdos396
