# Erdős Problem #396 — Witness Search

High-performance parallel search for witnesses to [Erdős Problem #396](https://www.erdosproblems.com/396) (OEIS [A375077](https://oeis.org/A375077)):

> For each *k*, find the smallest *n* such that *n*(*n*−1)(*n*−2)⋯(*n*−*k*) divides C(2*n*, *n*).

## Key Insight: Governor Set Approach

All known witnesses have every block term in the **Governor Set** G = {n : n | C(2n, n)}.
This reduces the search to finding runs of consecutive Governor Set members, then verifying
each candidate with a full p-adic valuation test.

## Building

```bash
cargo build --release
```

## Usage

### Search

```bash
# Search for k=9 witnesses in a range
cargo run --release -- --start 1000000 --end 2000000 -k 9

# Use a specific number of workers
cargo run --release -- --start 1000000 --end 2000000 -k 9 --workers 8

# List known witnesses
cargo run --release -- --list-known

# Verify all known witnesses
cargo run --release -- --verify-known
```

### Verification

```bash
# Verify a specific candidate
cargo run --release --bin verify -- -k 8 -n 339949252

# Verify all known witnesses
cargo run --release --bin verify -- --known
```

## Known Witnesses (OEIS A375077)

| k  | Smallest witness n | Notes |
|----|-------------------|-------|
| 1  | 2                 |       |
| 2  | 2,480             |       |
| 3  | 8,178             |       |
| 4  | 45,153            |       |
| 5  | 3,648,841         |       |
| 6  | 7,979,090         |       |
| 7  | 101,130,029       |       |
| 8  | 339,949,252       |       |
| 9  | 17,609,764,993    | Run of 11 consecutive governors |
| 10 | 17,609,764,994    | Same run of 11 |

## Architecture

| File | Purpose |
|------|---------|
| `sieve.rs` | Prime sieve (Sieve of Eratosthenes) |
| `factor.rs` | Integer factorization via trial division |
| `governor.rs` | Governor Set membership checking |
| `verify.rs` | Full witness verification with p-adic analysis |
| `search.rs` | Parallel search for Governor runs |
| `checkpoint.rs` | Checkpoint management for long-running searches |

## Checkpoints

Search progress is automatically saved to the `checkpoints/` directory.
Each worker maintains its own checkpoint file (`checkpoint_k{k}_w{id}.json`),
allowing searches to be paused and resumed without losing progress.

## License

MIT — see [LICENSE](LICENSE).
