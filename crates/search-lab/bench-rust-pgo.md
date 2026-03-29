# search-lab Benchmark (Rust + PGO)

Steady-state throughput with Profile-Guided Optimization. Trained on k=9,10
full ranges + k=13 10B slice before optimized rebuild. Fixed n-ranges
tuned for ~6s runtime; candidate counts and witness counts serve as
the deterministic accuracy signature. 8s cool-down between.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **Compiler** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **PGO** | yes (profile-generate + profile-use) |
| **Commit** | `ea9c6b3` |
| **Date** | 2026-03-29 02:58:32 UTC |

## Results

| k | Range | Candidates | Witnesses | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|----------:|---------:|-----------------:|
| 8 | `[339,949,252, 18,339,949,252)` | 18,000,000,000 | 39 | 5.5470 | 3245.0 |
| 9 | `[1,019,547,844, 19,019,547,844)` | 18,000,000,000 | 2 | 4.9441 | 3640.7 |
| 10 | `[17,609,764,994, 34,609,764,994)` | 17,000,000,000 | 0 | 5.0622 | 3358.3 |
| 11 | `[1,070,858,041,585, 1,086,858,041,585)` | 16,000,000,000 | 0 | 5.5643 | 2875.5 |
| 12 | `[5,048,891,644,646, 5,061,891,644,646)` | 13,000,000,000 | 0 | 5.4689 | 2377.1 |
| 13 | `[18,253,129,921,842, 18,264,129,921,842)` | 11,000,000,000 | 0 | 5.1125 | 2151.6 |

