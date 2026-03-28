# search-lab Benchmark (C++ port)

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
| **Threads** | 16 |
| **C++** | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| **Commit** | `1ea2a1f` |
| **Date** | 2026-03-28 23:54:55 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 8 | `[1, 339949253)` | 339,949,252 | 1.545 | 220.0 |
| 9 | `[339949252, 1019547845)` | 679,598,593 | 1.587 | 428.3 |
