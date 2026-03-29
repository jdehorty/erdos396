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
| **Commit** | `16620b7` |
| **Date** | 2026-03-29 03:18:11 UTC |

## Results

| k | Range | Candidates | Witnesses | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|----------:|---------:|-----------------:|
| 8 | `[169,974,626, 18,509,923,878)` | 18,339,949,252 | 40 | 5.8231 | 3149.5 |
| 9 | `[509,773,922, 19,029,321,766)` | 18,519,547,844 | 3 | 7.3276 | 2527.4 |
| 10 | `[8,804,882,497, 25,749,999,991)` | 16,945,117,494 | 1 | 5.4188 | 3127.1 |
| 11 | `[1,062,858,041,585, 1,078,999,999,999)` | 16,141,958,414 | 1 | 6.7559 | 2389.3 |
| 12 | `[5,042,391,644,646, 5,055,499,999,999)` | 13,108,355,353 | 1 | 5.5174 | 2375.8 |
| 13 | `[18,247,629,921,842, 18,258,749,999,999)` | 11,120,078,157 | 1 | 5.4677 | 2033.8 |

