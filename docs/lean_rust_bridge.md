# Lean ↔ Rust bridge (informal)

The Lean project in `formal/` proves the *mathematical statements* used to justify the
structure of the Rust algorithms, but it does **not** formally verify the Rust code.

This note records the intended correspondence for reviewers.

## `validate` (barrier-prime screen + full verification)

Rust code:
- `src/bin/validate.rs`
  - Maintains a sliding window `demand_sums[p] = Σ v_p(n-i)` (i=0..k) via
    `init_demand_window` / `advance_demand_window`, then checks
    `demand_sums[p] ≤ v_p(C(2n,n))` at each barrier prime (`passes_small_prime_screen_window`).
    This is equivalent to checking
    `v_p(n!) - v_p((n-k-1)!) ≤ v_p(C(2n,n))`.

Lean statements:
- `formal/SmallPrimeBarrier/Algorithm.lean`
  - `Erdos396.isWitness_le` (witness ⇒ `k+1 ≤ n`)
  - `Erdos396.padicValNat_descFactorial_eq_sub_factorial`
  - `Erdos396.padicValNat_factorial_sub_le_supply`
  - `Erdos396.demandSum` / `Erdos396.demandSum_succ_add` (sliding-window demand sum + update rule)

These show that the *same inequality* holds for every prime `p` whenever `IsWitness k n` holds.

## Kummer / carry counting for `v_p(C(2n,n))`

Rust code:
- `src/governor.rs`
  - `vp_central_binom_kummer_fast` and prime-specialized variants `vp_central_binom_p2/p3/p5`.

Lean statement:
- `formal/SmallPrimeBarrier/Algorithm.lean`
  - `Erdos396.supply_eq_carry_count` (Kummer, specialized to central binomial).

The Rust repo includes runtime cross-checks (Kummer vs Legendre) in unit tests and (optionally)
in `validate` via `--self-check-samples` and `--audit-interval`. See `docs/trust.md`.

## Small Prime Barrier theorem (why “non-governor behavior” is confined to small primes)

Rust code:
- `src/governor.rs` defines the governor test.
- `src/bin/validate.rs` uses only primes `p < 2k+1` for the *screen*, then fully verifies
  screen-pass candidates.

Lean statement:
- `formal/SmallPrimeBarrier/Main.lean`
  - `Erdos396.small_prime_barrier` and related corollaries (see `formal/README.md`).

## Scope reminder

The items above are an *informal mapping* between a verified theorem (Lean) and an implementation
(Rust). The remaining assumptions and mitigations are stated in `docs/trust.md`.
