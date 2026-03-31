# search-lab Benchmark (Rust)

Steady-state throughput at each k's witness position. Fixed n-ranges
tuned for ~6s runtime; candidate counts and witness counts serve as
the deterministic accuracy signature.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **Compiler** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **PGO** | yes (profile-generate + profile-use) |
| **Commit** | `84d2a62` |
| **Date** | 2026-03-31 08:02:46 UTC |

## Results

| k | Range | Candidates | Witnesses | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|----------:|---------:|-----------------:|
| 8 | `[169,974,626, 18,509,923,878)` | 18,339,949,252 | 40 | 5.9733 | 3070.3 |
| 9 | `[509,773,922, 19,029,321,766)` | 18,519,547,844 | 3 | 5.5627 | 3329.2 |
| 10 | `[8,804,882,497, 25,749,999,991)` | 16,945,117,494 | 1 | 4.0951 | 4137.9 |
| 11 | `[1,062,858,041,585, 1,078,999,999,999)` | 16,141,958,414 | 1 | 4.6942 | 3438.7 |
| 12 | `[5,042,391,644,646, 5,055,499,999,999)` | 13,108,355,353 | 1 | 4.4434 | 2950.1 |
| 13 | `[18,247,629,921,842, 18,258,749,999,999)` | 11,120,078,157 | 1 | 4.6988 | 2366.6 |

