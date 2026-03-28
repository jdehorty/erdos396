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
| **Commit** | `1ea2a1f` |
| **Date** | 2026-03-28 23:51:16 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
