# search-lab Benchmark (Rust)

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
| **rustc** | `rustc 1.87.0 (17067e9ac 2025-05-09) (Homebrew)` |
| **Commit** | `e2cb440` |
| **Date** | 2026-03-28 23:39:06 UTC |

## Results

| k | Range | Candidates | Time (s) | Throughput (M/s) |
|--:|:------|-----------:|---------:|-----------------:|
| 8 | `[1, 339949253)` | 339,949,252 | 1.433 | 237.3 |
| 9 | `[339949252, 1019547845)` | 679,598,593 | 1.533 | 443.3 |
| 10 | `[1019547844, 17609764995)` | 16,590,217,151 | 15.743 | 1053.8 |
