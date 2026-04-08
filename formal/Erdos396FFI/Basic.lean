/-!
# Erdős 396 — Executable Oracle (no Mathlib dependency)

Efficient, standalone implementations of the core mathematical functions used by
the Rust search-lab. These serve as the reference oracle for differential random
testing (DRT): generate millions of inputs, run both Lean and Rust, assert agreement.

The formal proofs in `SmallPrimeBarrier/` (which depend on Mathlib) establish that
these algorithms correctly compute the mathematical definitions. This module
implements the same algorithms without Mathlib, so it compiles fast and produces
a small binary.

## Functions

* `sievePrimes`          — Sieve of Eratosthenes
* `vpFactorial`          — v_p(n!) via Legendre's formula
* `vpCentralBinom`       — v_p(C(2n,n)) via Legendre
* `vpCentralBinomKummer` — v_p(C(2n,n)) via carry counting (Kummer's theorem)
* `isqrt`                — Integer square root
* `vp`                   — v_p(m): p-adic valuation of m
* `isGovernor`           — Governor Set membership: n | C(2n, n)
* `checkWitness`         — Witness check: n(n-1)...(n-k) | C(2n, n)
-/

namespace Erdos396FFI

/-- Sieve of Eratosthenes: all primes p with 2 ≤ p ≤ `limit`. -/
def sievePrimes (limit : Nat) : Array Nat := Id.run do
  if limit < 2 then return #[]
  let mut isPrime := Array.mkArray (limit + 1) true
  isPrime := isPrime.set! 0 false
  isPrime := isPrime.set! 1 false
  let mut p := 2
  while p * p <= limit do
    if isPrime[p]! then
      let mut m := p * p
      while m <= limit do
        isPrime := isPrime.set! m false
        m := m + p
    p := p + 1
  let mut result : Array Nat := #[]
  for i in [2:limit + 1] do
    if isPrime[i]! then
      result := result.push i
  return result

/-- `v_p(n!)` via Legendre's formula: `Σ_{i≥1} ⌊n / p^i⌋`. -/
def vpFactorial (n p : Nat) : Nat := Id.run do
  if p < 2 then return 0
  let mut v := 0
  let mut power := p
  while power <= n do
    v := v + n / power
    -- Overflow guard: if power > n / p, next power would exceed n anyway
    if power > n / p then break
    power := power * p
  return v

/-- `v_p(C(2n, n))` via Legendre: `v_p((2n)!) - 2 · v_p(n!)`. -/
def vpCentralBinom (n p : Nat) : Nat :=
  vpFactorial (2 * n) p - 2 * vpFactorial n p

/-- `v_p(C(2n, n))` via Kummer's theorem: count carries when adding `n + n` in base `p`. -/
def vpCentralBinomKummer (n p : Nat) : Nat := Id.run do
  if p < 2 then return 0
  let mut carries := 0
  let mut remaining := n
  let mut carry : Nat := 0
  while remaining > 0 do
    let digit := remaining % p
    let sum := digit + digit + carry
    if sum >= p then
      carry := 1
      carries := carries + 1
    else
      carry := 0
    remaining := remaining / p
  return carries

/-- Integer square root via Newton's method. Returns `⌊√n⌋`. -/
def isqrt (n : Nat) : Nat := Id.run do
  if n <= 1 then return n
  let mut x := n
  let mut y := (x + 1) / 2
  while y < x do
    x := y
    y := (x + n / x) / 2
  return x

/-- `v_p(m)`: how many times `p` divides `m`. -/
def vp (m p : Nat) : Nat := Id.run do
  if p < 2 || m == 0 then return 0
  let mut v := 0
  let mut x := m
  while x % p == 0 do
    v := v + 1
    x := x / p
  return v

/--
Governor Set membership: `n ∈ G ⟺ n | C(2n, n)`.

Equivalent to: for all primes `p` dividing `n`, `v_p(n) ≤ v_p(C(2n, n))`.

Uses trial division with the supplied prime table.
-/
def isGovernor (n : Nat) (primes : Array Nat) : Bool := Id.run do
  if n <= 1 then return (n == 1)
  let mut remaining := n
  for p in primes do
    if p * p > remaining then break
    if remaining % p == 0 then
      let mut exp := 0
      while remaining % p == 0 do
        exp := exp + 1
        remaining := remaining / p
      let supply := vpCentralBinom n p
      if exp > supply then return false
  -- Remaining factor > sqrt(n) is prime
  if remaining > 1 then
    let supply := vpCentralBinom n remaining
    if 1 > supply then return false
  return true

/--
Witness check: `n(n-1)···(n-k) | C(2n, n)`.

For each prime `p` appearing in any block term `{n-k, ..., n}`:
  demand_p = Σ v_p(n-i), for i = 0..k
  supply_p = v_p(C(2n, n))
  Require: demand_p ≤ supply_p

This uses direct trial-division factoring of each block term — algorithmically
independent from search-lab's modular-inverse sieve approach.
-/
def checkWitness (k n : Nat) (primes : Array Nat) : Bool := Id.run do
  if k == 0 || n <= k then return false

  -- Step 1: Factor each block term, accumulate per-prime demand
  -- Use a flat array of (prime, demand) pairs
  let mut demandPairs : Array (Nat × Nat) := #[]

  for i in [0:k + 1] do
    let ni := n - i
    let mut remaining := ni
    -- Trial division by small primes
    for p in primes do
      if p * p > remaining then break
      if remaining % p == 0 then
        let mut e := 0
        while remaining % p == 0 do
          e := e + 1
          remaining := remaining / p
        demandPairs := addDemand demandPairs p e
    -- Large prime factor (residual after trial division)
    if remaining > 1 then
      demandPairs := addDemand demandPairs remaining 1

  -- Step 2: Check supply ≥ demand for every prime
  for (p, demand) in demandPairs do
    let supply := vpCentralBinom n p
    if demand > supply then return false

  return true
where
  /-- Add `e` to the demand for prime `p` in the association list. -/
  addDemand (pairs : Array (Nat × Nat)) (p e : Nat) : Array (Nat × Nat) := Id.run do
    let mut pairs := pairs
    for idx in [0:pairs.size] do
      if pairs[idx]!.1 == p then
        pairs := pairs.set! idx (p, pairs[idx]!.2 + e)
        return pairs
    return pairs.push (p, e)

end Erdos396FFI
