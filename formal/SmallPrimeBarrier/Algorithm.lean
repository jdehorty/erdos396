import SmallPrimeBarrier.Defs

import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.NumberTheory.Padics.PadicVal.Basic

/-!
# Algorithmic specifications used by `validate`

This file records a few key *mathematical* identities that justify the arithmetic
performed by the Rust `validate` binary:

* `supply p n = v_p(C(2n,n))` (definition).
* Kummer's theorem specialized to `C(2n,n)`: `supply p n` equals the number of
  carries when adding `n + n` in base `p` (`supply_eq_carry_count`).
* The block demand identity used by the sliding window:
  `v_p(∏_{i=0..k} (n-i)) = v_p(n!) - v_p((n-k-1)!)` (`padicValNat_descFactorial_eq_sub_factorial`).

These statements are independent of any particular implementation language.
-/

namespace Erdos396

open Nat
open scoped BigOperators

section Supply

/-- `supply p n` is the `p`-adic valuation of `C(2n,n)` (central binomial). -/
@[simp] theorem supply_eq_padicValNat_centralBinom (p n : ℕ) :
    supply p n = padicValNat p (centralBinom n) := by
  rfl

/--
Kummer's theorem specialized to the central binomial coefficient:
`v_p(C(2n,n))` equals the number of carries when adding `n + n` in base `p`,
expressed using the standard `mod p^i` predicate from Mathlib.
-/
theorem supply_eq_carry_count (p n b : ℕ) [Fact p.Prime] (hnb : log p (n + n) < b) :
    supply p n =
      ((Finset.Ico 1 b).filter fun i => p ^ i ≤ n % p ^ i + n % p ^ i).card := by
  -- `centralBinom n = choose (2n) n` by definition.
  -- Apply Kummer's theorem for `(n+n).choose n`.
  simpa [supply, centralBinom, two_mul, add_comm, add_left_comm, add_assoc] using
    (padicValNat_choose' (p := p) (n := n) (k := n) (b := b) hnb)

end Supply

section Demand

private theorem padicValNat_le_of_dvd {p a b : ℕ} [Fact p.Prime] (hb : b ≠ 0) (hab : a ∣ b) :
    padicValNat p a ≤ padicValNat p b := by
  -- `p^(v_p(a)) ∣ a` and `a ∣ b` implies `p^(v_p(a)) ∣ b`.
  have hpow : p ^ padicValNat p a ∣ b := dvd_trans (pow_padicValNat_dvd (p := p) (n := a)) hab
  exact (padicValNat_dvd_iff_le (p := p) hb).1 hpow

/--
The block product `n(n-1)…(n-k)` is `descFactorial n (k+1)`.

If `n` is a `k`-witness, then for every prime `p` the `p`-adic demand of the
block is bounded by the `p`-adic supply of `C(2n,n)`.
-/
theorem padicValNat_descFactorial_le_supply {k n p : ℕ} (hp : Nat.Prime p) (hw : IsWitness k n) :
    padicValNat p (descFactorial n (k + 1)) ≤ supply p n := by
  haveI : Fact p.Prime := ⟨hp⟩
  -- Witness condition is a divisibility.
  have hcb_ne : centralBinom n ≠ 0 := Nat.centralBinom_ne_zero n
  -- Monotonicity of `padicValNat` under divisibility.
  have := padicValNat_le_of_dvd (p := p) (a := descFactorial n (k + 1)) (b := centralBinom n)
    hcb_ne hw
  simpa [supply] using this

/-- A witness must satisfy `k+1 ≤ n` (otherwise the block product is `0`). -/
theorem isWitness_le {k n : ℕ} (hw : IsWitness k n) : k + 1 ≤ n := by
  by_contra hkn
  have hlt : n < k + 1 := Nat.lt_of_not_ge hkn
  have hdf0 : descFactorial n (k + 1) = 0 :=
    (Nat.descFactorial_eq_zero_iff_lt (n := n) (k := k + 1)).2 hlt
  have h0 : centralBinom n = 0 := by
    have hw' : descFactorial n (k + 1) ∣ centralBinom n := hw
    rw [hdf0] at hw'
    exact (zero_dvd_iff.mp hw')
  exact (Nat.centralBinom_ne_zero n) h0

/--
Demand identity used by `validate`:

`v_p(descFactorial n (k+1)) = v_p(n!) - v_p((n-k-1)!)` (for `k+1 ≤ n`).
-/
theorem padicValNat_descFactorial_eq_sub_factorial {p n k : ℕ} [Fact p.Prime]
    (hk : k + 1 ≤ n) :
    padicValNat p (descFactorial n (k + 1)) =
      padicValNat p (n !) - padicValNat p ((n - (k + 1))!) := by
  -- From `factorial_mul_descFactorial`: (n-(k+1))! * descFactorial n (k+1) = n!
  have hmul :
      (n - (k + 1))! * descFactorial n (k + 1) = n ! := by
    -- Avoid `simp` here; we want the exact `descFactorial` term.
    exact Nat.factorial_mul_descFactorial (n := n) (k := k + 1) hk
  -- Take `padicValNat` of both sides.
  have hleft0 : (n - (k + 1))! ≠ 0 := Nat.factorial_ne_zero (n - (k + 1))
  have hdf0 : descFactorial n (k + 1) ≠ 0 := by
    intro h0
    have hmul' := hmul
    rw [h0] at hmul'
    have : n ! = 0 := by simpa using hmul'.symm
    exact Nat.factorial_ne_zero n this
  have hval := congrArg (padicValNat p) hmul
  -- `v_p(a*b) = v_p(a)+v_p(b)` for nonzero.
  have hsum :
      padicValNat p ((n - (k + 1))!) + padicValNat p (descFactorial n (k + 1)) =
        padicValNat p (n !) := by
    -- Rewrite `v_p((n-(k+1))! * descFactorial n (k+1))` and use `hval`.
    have hmulVal :
        padicValNat p ((n - (k + 1))! * descFactorial n (k + 1)) =
          padicValNat p ((n - (k + 1))!) + padicValNat p (descFactorial n (k + 1)) :=
      padicValNat.mul (p := p) (a := (n - (k + 1))!) (b := descFactorial n (k + 1)) hleft0 hdf0
    exact hmulVal.symm.trans hval
  -- Rearrange.
  set a := padicValNat p ((n - (k + 1))!) with ha
  set d := padicValNat p (descFactorial n (k + 1)) with hd
  have had : a + d = padicValNat p (n !) := by simpa [ha, hd] using hsum
  calc
    padicValNat p (descFactorial n (k + 1)) = d := by simp [hd]
    _ = (a + d) - a := by simp
    _ = padicValNat p (n !) - a := by simp [had]
    _ = padicValNat p (n !) - padicValNat p ((n - (k + 1))!) := by simp [ha]

/--
For a witness `n`, the block demand written in factorial valuations is bounded by the supply:

`v_p(n!) - v_p((n-k-1)!) ≤ v_p(C(2n,n))`.

This is exactly the inequality checked by the `validate` small-prime screen.
-/
theorem padicValNat_factorial_sub_le_supply {k n p : ℕ} (hp : Nat.Prime p) (hw : IsWitness k n) :
    padicValNat p (n !) - padicValNat p ((n - (k + 1))!) ≤ supply p n := by
  have hk : k + 1 ≤ n := isWitness_le (k := k) (n := n) hw
  haveI : Fact p.Prime := ⟨hp⟩
  have hle :
      padicValNat p (descFactorial n (k + 1)) ≤ supply p n :=
    padicValNat_descFactorial_le_supply (k := k) (n := n) (p := p) hp hw
  have heq :
      padicValNat p (descFactorial n (k + 1)) =
        padicValNat p (n !) - padicValNat p ((n - (k + 1))!) :=
    padicValNat_descFactorial_eq_sub_factorial (p := p) (n := n) (k := k) hk
  have hle' := hle
  rw [heq] at hle'
  exact hle'

end Demand

section Window

/--
`demandSum p k n` is the p-adic demand of the witness block `{n-k, ..., n}` written as
an explicit sum:

`Σ_{i=0..k} v_p(n-i)`.

This is the quantity maintained incrementally by the Rust `validate` binary's sliding window.
-/
def demandSum (p k n : ℕ) : ℕ :=
  ∑ i in Finset.range (k + 1), padicValNat p (n - i)

/--
Sliding-window update rule for `demandSum` (algebraic form):

`demandSum p k (n+1) + v_p(n-k) = demandSum p k n + v_p(n+1)`.

This corresponds exactly to updating a running sum by removing the term leaving the window
(`n-k`) and adding the new term entering the window (`n+1`).
-/
theorem demandSum_succ_add (p k n : ℕ) [Fact p.Prime] :
    demandSum p k (n + 1) + padicValNat p (n - k) =
      demandSum p k n + padicValNat p (n + 1) := by
  have h_succ :
      demandSum p k (n + 1) =
        padicValNat p (n + 1) + ∑ i in Finset.range k, padicValNat p (n - i) := by
    simp [demandSum, Finset.sum_range_succ', Nat.succ_sub_succ_eq_sub, add_assoc, add_comm,
      add_left_comm]

  have h_here :
      demandSum p k n = (∑ i in Finset.range k, padicValNat p (n - i)) + padicValNat p (n - k) := by
    simp [demandSum, Finset.sum_range_succ, add_assoc, add_comm, add_left_comm]

  calc
    demandSum p k (n + 1) + padicValNat p (n - k)
        = (padicValNat p (n + 1) + ∑ i in Finset.range k, padicValNat p (n - i)) +
            padicValNat p (n - k) := by simp [h_succ]
    _ = padicValNat p (n + 1) +
          ((∑ i in Finset.range k, padicValNat p (n - i)) + padicValNat p (n - k)) := by
          simp [add_assoc, add_comm, add_left_comm]
    _ = padicValNat p (n + 1) + demandSum p k n := by simp [h_here, add_assoc]
    _ = demandSum p k n + padicValNat p (n + 1) := by ac_rfl

end Window

end Erdos396
