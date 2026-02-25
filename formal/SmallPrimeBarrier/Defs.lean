import Mathlib.Data.Nat.Choose.Central
import Mathlib.NumberTheory.Padics.PadicVal.Basic

namespace Erdos396

open Nat

/-- `supply p n = v_p(C(2n, n))` -/
def supply (p n : ℕ) : ℕ := padicValNat p (centralBinom n)

/-- A k-witness: n(n-1)...(n-k) | C(2n, n) -/
def IsWitness (k n : ℕ) : Prop := descFactorial n (k + 1) ∣ centralBinom n

end Erdos396
