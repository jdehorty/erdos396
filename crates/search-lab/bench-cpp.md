# search-lab Benchmark (C++ reference, Sharvil Kesarwani)

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
| **Compiler** | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| **Commit** | `84d2a62` |
| **Date** | 2026-03-31 08:04:56 UTC |

## Results

| k | Range | Candidates | Witnesses | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|----------:|---------:|-----------------:|
| 8 | `[169,974,626, 18,509,923,878)` | 18,339,949,252 | 40 | 4.9581 | 3699.0 |
| 9 | `[509,773,922, 19,029,321,766)` | 18,519,547,844 | 3 | 5.1458 | 3599.0 |
| 10 | `[8,804,882,497, 25,749,999,991)` | 16,945,117,494 | 1 | 4.6561 | 3639.3 |
| 11 | `[1,062,858,041,585, 1,078,999,999,999)` | 16,141,958,414 | 1 | 5.0269 | 3211.1 |
| 12 | `[5,042,391,644,646, 5,055,499,999,999)` | 13,108,355,353 | 1 | 4.8653 | 2694.2 |
| 13 | `[18,247,629,921,842, 18,258,749,999,999)` | 11,120,078,157 | 1 | 4.6706 | 2380.9 |

