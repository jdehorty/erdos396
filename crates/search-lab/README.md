# search-lab

An experimental workspace providing a common ground for comparing search
algorithm optimizations across languages and implementations. Each approach
runs against the same standardized benchmarks (k=8, k=9, k=10 over the known
witness ranges), making it straightforward to determine what actually improves
throughput and what does not.

Optimizations are prototyped and measured here before being integrated into the
main `erdos396` crate.

## Structure

```
search-lab/
├── justfile            Build, benchmark, and clean recipes
├── bench-rust.md       Baseline Rust benchmark results
├── src/main.rs         Rust implementation
└── cpp-reference/
    ├── erdos_problem_396.cpp     C++ reference implementation by Sharvil Kesarwani
    ├── justfile
    └── README.md                 Attribution and notes
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
