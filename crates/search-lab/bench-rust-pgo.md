# search-lab Benchmark (Rust + PGO)

Per-k throughput with Profile-Guided Optimization. Trained on k=9,10
full ranges + k=13 10B slice before optimized rebuild.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **Compiler** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **PGO** | yes (profile-generate + profile-use) |
| **Commit** | `7548c86` |
| **Date** | 2026-03-29 02:14:21 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.1541 | 4410.1 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 4.9274 | 3366.9 |
| 11 | `[17609764994, 67609764994)` | 50,000,000,000 | 17.3843 | 2876.2 |
| 13 | `[18185829921842, 18253129921842)` | 67,300,000,000 | 43.5797 | 1544.3 |

