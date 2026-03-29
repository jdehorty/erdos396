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
| **rustc** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **Commit** | `42b43c1` |
| **Date** | 2026-03-29 01:08:12 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.168 | 4042.8 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 5.086 | 3261.7 |
| 11 | `[17609764994, 67609764994)` | 50,000,000,000 | 19.629 | 2547.3 |
| 13 | `[18185829921842, 18253129921842)` | 67,300,000,000 | 38.233 | 1760.2 |

