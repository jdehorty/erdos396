# search-lab Benchmark (Rust)

Measures per-k throughput of the witness search (sieve time excluded).
k=8, k=9, k=10 run over full witness ranges. k=13 runs a 10B-candidate
slice at production scale (sieve_limit ~ 6M primes).

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **rustc** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **Commit** | `fb274e8` |
| **Date** | 2026-03-29 00:51:34 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 8 | `[1, 339949253)` | 339,949,252 | 0.0728 | 4669.91 |
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.1928 | 3525.16 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 5.6256 | 2949.06 |
| 13 | `[18243129921842, 18253129921843)` | 10,000,000,001 | 6.7576 | 1479.81 |

