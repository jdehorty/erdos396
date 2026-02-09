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

## Known Witnesses

The official sequence is OEIS [A375077](https://oeis.org/A375077) ([b-file](https://oeis.org/A375077/b375077.txt)), which currently lists k=1 through k=7:

| k | Smallest witness n | Source |
|---|-------------------|--------|
| 1 | 2                 | OEIS   |
| 2 | 2,480             | OEIS   |
| 3 | 8,178             | OEIS   |
| 4 | 45,153            | OEIS   |
| 5 | 3,648,841         | OEIS   |
| 6 | 7,979,090         | OEIS   |
| 7 | 101,130,029       | OEIS   |

The following witnesses were discovered using this tool and have not yet been added to the OEIS:

| k  | Smallest witness n | Source | Notes |
|----|-------------------|--------|-------|
| 8  | 339,949,252       | This project | Found 2025-01-17 |
| 9  | 17,609,764,993    | This project | Found 2025-01-20; run of 11 consecutive governors |
| 10 | 17,609,764,994    | This project | Found 2025-01-20; same run of 11 |

## How It Works

The search pipeline uses a fused sieve+governor computation to reject ~87% of candidates
in a single pass, then detects consecutive governor runs and verifies them:

```mermaid
graph TD
    A["Input: batch of numbers<br/>[n, n+1, ..., n+batch_size)"] --> B

    subgraph FUSED["Fused Sieve Pass"]
        direction TB
        B{"n is odd prime?"}
        B -- "Yes (4.3%)" --> R1["Reject"]
        B -- No --> C["Divide out small primes p,<br/>check v_p via Kummer's theorem"]
        C --> D{"v_p(C(2n,n)) &ge; v_p(n)<br/>for all small p?"}
        D -- "No (25.8%)" --> R2["Reject"]
        D -- Yes --> E{"Cofactor after sieve > 1?<br/>(large prime factor > &radic;2n)"}
        E -- "Yes (57.2%)" --> R3["Reject"]
        E -- No --> G["Governor &check;"]
    end

    G --> H{"Run of k+1<br/>consecutive<br/>governors?"}
    H -- No --> I["Continue scanning"]
    H -- Yes --> J["Full p-adic<br/>verification"]
    J --> K{"All primes<br/>satisfied?"}
    K -- No --> L["False positive"]
    K -- Yes --> M["Witness found"]

    style FUSED fill:#1a1a2e,stroke:#e94560,stroke-width:2px,color:#eee
    style G fill:#0f3460,stroke:#16c79a,color:#eee
    style M fill:#0f3460,stroke:#16c79a,stroke-width:3px,color:#eee
    style R1 fill:#1a1a2e,stroke:#e94560,color:#999
    style R2 fill:#1a1a2e,stroke:#e94560,color:#999
    style R3 fill:#1a1a2e,stroke:#e94560,color:#999
    style L fill:#1a1a2e,stroke:#e94560,color:#999
```

Only ~12.6% of candidates survive the fused sieve (matching the Ford–Konyagin governor density).
Kummer's theorem replaces Legendre's formula for computing v_p(C(2n,n)), halving the number of
iterations per prime — and for p=2 it reduces to a single `POPCNT` instruction.

## Architecture

| File | Purpose |
|------|---------|
| `sieve.rs` | Prime sieve (Sieve of Eratosthenes) |
| `factor.rs` | Integer factorization via trial division |
| `governor.rs` | Governor Set membership checking via Kummer's theorem |
| `prefilter.rs` | Fused sieve+governor batch computation (performance-critical hot loop) |
| `verify.rs` | Full witness verification with p-adic analysis |
| `search.rs` | Parallel search with rayon, run detection, checkpointing |
| `checkpoint.rs` | Checkpoint save/resume for long-running searches |

## Checkpoints

Search progress is automatically saved to the `checkpoints/` directory.
Each worker maintains its own checkpoint file (`checkpoint_k{k}_w{id}.json`),
allowing searches to be paused and resumed without losing progress.

## License

MIT — see [LICENSE](LICENSE).
