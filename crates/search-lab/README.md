# search-lab

A minimal, bare-bones implementation of the Erdos Problem #396 witness search,
focused on maximum throughput with zero unnecessary abstractions.

This crate serves as a lightweight testing ground for search algorithm
improvements. New optimization ideas — modular-inverse exact division,
const-generic Kummer kernels, cache-resident blocking strategies, etc. — are
prototyped and benchmarked here before being integrated into the main `erdos396`
crate.

## Structure

```
search-lab/
├── justfile            Build, benchmark, and clean recipes
├── src/main.rs         Rust implementation
└── cpp-reference/
    ├── main.cpp                  C++ port
    ├── erdos_problem_396.cpp     Reference implementation by Sharvil Kesarwani
    ├── justfile
    └── README.md                 C++ source attribution and notes
```

## Building

```sh
just build        # both Rust and C++
just build-rust   # Rust only
just build-cpp    # C++ only
```

## Benchmarking

Benchmarks run k=8, k=9, k=10 over the known witness ranges (matching what the
C++ reference implementation traverses):

```sh
just bench        # both implementations side by side
just bench-rust   # Rust only
just bench-cpp    # C++ only

# Override thread count
just threads=8 bench
```

## Usage

```sh
# Run witness search starting at k=9
cargo run -p erdos396-search-lab --release -- -k 9

# With custom thread count and range
cargo run -p erdos396-search-lab --release -- -k 9 --threads 8 --start 1000000 --end 1000000000
```
