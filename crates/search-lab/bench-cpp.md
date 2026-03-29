# search-lab Benchmark (C++ reference, Sharvil Kesarwani)

Measures wall-clock throughput of the witness search over the known
witness ranges for k=8, k=9, and k=10. These ranges reproduce the
same search space traversed by the C++ reference implementation,
providing a standardized baseline for comparing optimizations.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 (hardcoded to hardware_concurrency) |
| **C++** | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| **Commit** | `a00061d` |
| **Date** | 2026-03-29 00:01:05 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 8 | `[1, 339949253)` | 339,949,252 | 0.1120 | 3034.35 |
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.2385 | 2849.84 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 5.9057 | 2809.21 |

