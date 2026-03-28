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
├── src/main.rs     Rust implementation (cargo build --release)
└── cpp-reference/
    ├── main.cpp                  C++ port (just build)
    ├── erdos_problem_396.cpp     Reference implementation by Sharvil Kesarwani
    ├── justfile
    └── README.md                 C++ source attribution and notes
```

## Building

**Rust:**
```sh
cargo build -p erdos396-search-lab --release
```

**C++:**
```sh
cd cpp-reference && just build
```

## Usage

```sh
# Run witness search starting at k=9
cargo run -p erdos396-search-lab --release -- -k 9

# With custom thread count and range
cargo run -p erdos396-search-lab --release -- -k 9 --threads 8 --start 1000000 --end 1000000000
```
