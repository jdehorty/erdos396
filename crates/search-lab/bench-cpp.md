# search-lab Benchmark (C++ reference, Sharvil Kesarwani)

Measures wall-clock throughput of the witness search. k=8, k=9, k=10
run over the full known witness ranges. k=13 runs a 10B-candidate
slice at production scale (sieve_limit ~ 6M primes).

Note: thread count is hardcoded to hardware_concurrency() in the binary.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **C++** | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| **Commit** | `fb274e8` |
| **Date** | 2026-03-29 00:52:10 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 8 | `[1, 339949253)` | 339,949,252 | 0.0784 | 4337.85 |
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.1603 | 4239.20 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 5.2623 | 3152.64 |
| 13 | `[18243129921842, 18253129921843)` | 10,000,000,001 | 7.3914 | 1352.92 |

