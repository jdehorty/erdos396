# search-lab Benchmark (C++ reference, Sharvil Kesarwani)

Per-k throughput (sieve time excluded). k=9,10 use full witness
ranges. k=11 is capped at 50B candidates. k=13 matches the k=9-11
total (~67.3B candidates) at production sieve depth.

Note: thread count is hardcoded to hardware_concurrency() in the binary.

## Environment

| | |
|---|---|
| **CPU** | Apple M3 Max |
| **Memory** | 128 GB |
| **OS** | macOS 15.6.1 |
| **Threads** | 16 |
| **C++** | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| **Commit** | `4dd7054` |
| **Date** | 2026-03-29 01:54:37 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 9 | `[339949252, 1019547845)` | 679,598,593 | 0.2007 | 3386.1 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 6.7796 | 2447.1 |
| 11 | `[17609764994, 67609764994)` | 50,000,000,000 | 20.9288 | 2389.1 |
| 13 | `[18185829921842, 18253129921842)` | 67,300,000,000 | 44.8077 | 1502.0 |

