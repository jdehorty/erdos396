# Trust model and verification scope

This repository contains both a **machine-checked theorem** (Lean 4 / Mathlib)
and **high-performance search / validation code** (Rust). They support each
other, but they are not formally connected by verified compilation or a proof
of correctness of the Rust implementation.

This document clarifies what is proven, what is computed, and what assumptions
remain.

---

## What is *mathematically* proven (Lean)

The Lean project in `formal/` typechecks (no `sorry`/`admit`) and proves the
Small Prime Barrier results used to justify the *algorithmic* structure of
validation.

See:
- `formal/SmallPrimeBarrier/Main.lean`
- `formal/SmallPrimeBarrier/Algorithm.lean`
- `formal/README.md`
- `docs/lean_rust_bridge.md` (informal mapping to Rust)

---

## What is *computationally* checked (Rust)

There are two distinct computational claims:

1. **Witness validity (Claim A)**  
   A specific `(k, n)` is a witness iff the divisibility condition holds:
   `n(n-1)…(n-k) ∣ C(2n, n)`.
   This is independently checkable; the `verify` binary runs a full p-adic check.

2. **Minimality / exhaustiveness in a range (Claim B)**  
   The `validate` binary checks *every* `n` in a range against the barrier-prime
   screen (primes `p < 2k+1`), then fully verifies each screen-pass candidate.
   The Lean theorem justifies *why the screen is safe to use*, but the execution
   is still conventional software running on conventional hardware.

---

## Remaining assumptions (“air gap”)

Even with a machine-checked theorem, a computational minimality result depends on:

- Correctness of the Rust implementation
- Correctness of the compiler, standard library, and CPU execution
- No silent data corruption during long runs (memory / storage / hardware)

This is the standard trust base for computational mathematics.

---

## Mitigations implemented in this repo

To reduce practical risk in the Rust layer, we do the following:

- **Independent witness verifier (Claim A)**  
  `scripts/verify_witness.py` is a stdlib-only Python implementation of the
  witness divisibility test. It provides an independent cross-check of witness
  validity separate from the Rust implementation.

- **No floating point in correctness-critical paths**  
  All square-root bounds used for trial division and the `√(2n)` barrier are
  computed with exact integer square roots (see `src/int_math.rs`).

- **Optional fused-sieve audits in `erdos396`**  
  The main `erdos396` search uses a high-performance fused sieve. For additional
  confidence on a given machine, enable `erdos396 --fused-self-check-samples N`
  (startup samples per worker) and/or `erdos396 --fused-audit-interval M`
  (periodic in-run cross-checks) to compare fused-sieve membership against the
  direct governor test on sampled values.

- **Range coverage invariants in `erdos396`**  
  Each worker checkpoint tracks simple coverage invariants (sum and XOR of the
  scanned `n` values), and `erdos396` asserts they match the worker’s assigned
  range on successful completion. This helps detect accidental gaps/overlaps in
  the scan loop or resume logic.

- **Auditable run logs**  
  `erdos396` writes append-only `runs_k*_w*.jsonl` files recording significant governor runs.
  The `audit_runs` binary can re-verify these run log entries by full p-adic verification and
  produce a false-positive catalog. If a long run was stopped and resumed, the run logs may
  contain a small number of duplicate entries; `audit_runs` skips out-of-order duplicates by
  enforcing monotone-increasing run starts per file.

- **Cross-check tests for arithmetic kernels**  
  Unit tests compare Kummer carry-counting against Legendre’s formula on both
  systematic small ranges and random large samples.

- **Range coverage invariants in `validate`**  
  `validate` asserts that the total number of checked integers equals the
  expected scan size (excluding the trivial `n ≤ k` prefix), and also checks
  simple coverage invariants (sum and XOR of all scanned `n` values) to detect
  accidental gaps/overlaps. It writes a `validate_report_*.json` with the
  measured and expected invariants. The stdlib-only script
  `scripts/check_validate_report.py` independently recomputes and checks these
  invariants from the report, and can also verify that multiple chunk reports
  form a contiguous partition (`--check-partition`).
  Reports also embed build metadata (`git_hash`, toolchain) so results can be
  tied to a specific source revision.

If you need even higher confidence, rerun `validate` with a different worker
count and compare results and checkpoints, and consider enabling the startup
cross-check (`validate --self-check-samples N`) to sanity-check the arithmetic
kernels on the target machine and/or the periodic audit (`validate --audit-interval M`)
for long runs.
