# search-lab Benchmark (Rust)

Per-k throughput (sieve time excluded). k=9,10 use full witness
ranges. k=11 is capped at 50B candidates. k=13 matches the k=9-11
total (~67.3B candidates) at production sieve depth.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **Compiler** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **Commit** | `7548c86` |
| **Date** | 2026-03-29 02:12:08 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.2254 | 3015.1 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 5.9754 | 2776.4 |
| 11 | `[17609764994, 67609764994)` | 50,000,000,000 | 19.2883 | 2592.2 |
| 13 | `[18185829921842, 18253129921842)` | 67,300,000,000 | 41.9431 | 1604.6 |

