# Reproducibility checklist

This repo contains:

- A Lean 4 / Mathlib proof (`formal/`) establishing the Small Prime Barrier results used to
  justify the validation *algorithm*.
- Rust binaries implementing the search (`erdos396`), witness verification (`verify`), and the
  provably complete range validator (`validate`).

This checklist is aimed at reviewers who want to reproduce the main claims from scratch.

## 1) Build + sanity checks (minutes)

This repo pins the Rust toolchain via `rust-toolchain.toml` and the Lean toolchain via
`formal/lean-toolchain`.

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo build --release
```

## 1b) Optional runtime audits (recommended for long runs)

- `erdos396` supports optional cross-checks of the fused sieve against the direct governor test:

  ```bash
  cargo run --release --bin erdos396 -- --fused-self-check-samples 100 --fused-audit-interval 100000000 ...
  ```

- `validate` supports optional kernel audits during range scans:

  ```bash
  cargo run --release --bin validate -- --self-check-samples 100 --audit-interval 100000000 ...
  ```

## 2) Machine-check the Lean formalization (minutes)

```bash
cd formal
lake build
```

The Lean project is expected to build with no `sorry`/`admit`.

## 3) Check witness validity (seconds)

```bash
# Validity of the listed witnesses (divisibility condition)
cargo run --release --bin verify -- --known

# Independent cross-check (stdlib-only, no Rust involved)
python3 scripts/verify_witness.py --known

# Optional: machine-readable output
cargo run --release --bin verify -- --known --json > verify_known.json
```

## 4) Check minimality in a range (days)

Run the provably complete validator (Corollary 6 screen + full verification):

```bash
cargo run --release --bin validate -- -k 13 --start 0 --end 18253129921842 --workers 40
python3 scripts/check_validate_report.py validate_checkpoints/validate_report_k13_0_18253129921842.json
```

The `validate_report_*.json` includes embedded build metadata (`git_hash`, toolchain) so reports can
be tied to a specific source revision. For published / archival runs, we recommend building from a
clean checkout so `build.git_dirty=false` in the report.

Recommended robustness checks:

- Rerun with a different `--workers` count and confirm the witness list and report invariants match.
- If you run multiple contiguous chunks, verify the combined coverage invariants with
  `python3 scripts/check_validate_report.py --check-partition validate_checkpoints/validate_report_k13_*.json`.
- For long runs, consider enabling `--self-check-samples` (startup kernel checks) and/or
  `--audit-interval` (periodic in-run kernel checks).

## 5) Trust model

The Small Prime Barrier theorem is machine-checked, but the Rust code is not formally verified.
For the exact scope of what is proven vs. computed (and the remaining assumptions), see
`docs/trust.md`.

## 6) Optional: reproduce the Governor-run search (days to weeks)

The `validate` pass is the *completeness/minimality* certificate. The `erdos396` binary is the
high-performance Governor Set run search used to discover candidates.

If you want to reproduce the Governor-run scan over a range, run:

```bash
cargo run --release --bin erdos396 -- -k 13 --start 0 --end 25000000000000 --workers 40
```

On successful completion, `erdos396` writes a `search_report_*.json` into the output directory
(default `checkpoints/`). Independently check the report invariants with:

```bash
python3 scripts/check_search_report.py checkpoints/search_report_k13_0_25000000000000.json
```

If you also want to audit the run logs (`runs_k13_w*.jsonl`) to produce a false-positive catalog,
run:

```bash
cargo run --release --bin audit_runs -- -k 13 -d checkpoints -o audit_runs_k13_report.json
```
