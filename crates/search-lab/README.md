# search-lab

An experimental workspace providing a common ground for comparing search
algorithm optimizations across languages and implementations. All
implementations target the same Universal Benchmark Interface (below),
making it straightforward to determine what actually improves throughput.

Optimizations are prototyped and measured here before being integrated into the
main `erdos396` crate.

## Structure

```
search-lab/
├── justfile                Build, benchmark, and clean recipes
├── bench-rust.md           Rust benchmark results
├── bench-cpp.md            C++ benchmark results
├── parse_bench.py          R-line → markdown table converter
├── .cargo/config.toml      Rust build flags (target-cpu=native)
├── src/main.rs             Rust implementation
└── cpp-reference/
    ├── erdos_problem_396.cpp     C++ reference by Sharvil Kesarwani
    ├── justfile
    └── README.md
```

## Building

```sh
just build        # both Rust and C++
just build-rust   # Rust only
just build-cpp    # C++ only
```

## Benchmarking

```sh
just bench-rust   # Rust
just bench-cpp    # C++

# Override thread count
just threads=8 bench-rust
```

Results are written to `bench-rust.md` / `bench-cpp.md`. The benchmark
sweeps k=8 through k=13, each over a fixed n-range that includes the
known minimum witness. Candidate counts and witness counts serve as a
deterministic accuracy signature — if the sieve is correct, these are
identical across implementations.

## Usage

```sh
cargo run -p erdos396-search-lab --release -- -k 9
cargo run -p erdos396-search-lab --release -- -k 13 --threads 8 --start 1000000 --end 25000000000000
```

---

## Universal Benchmark Interface

### How the benchmark works

The benchmark measures steady-state sieve throughput by sweeping a fixed
range of the number line for each k value (k=8 through k=13). Each range
is chosen so that:

1. **It includes the known minimum witness** for that k (OEIS A375077),
   placing the sweep at the sieve depth characteristic of real search work.
2. **It runs for roughly 6 seconds**, long enough to reach thermal
   steady state but short enough to avoid excessive thermal throttling.
3. **Candidates and witnesses are deterministic.** For a given range
   `[start, end)` and k, the number of candidates checked (`end − start`)
   and the number of witnesses found are fixed by the mathematics — they
   do not depend on the implementation, language, or hardware. These two
   numbers form the **accuracy signature** of the run.

The accuracy signature is the key to fidelity. If an implementation reports
the correct candidate count and witness count, its sieve and verification
logic are producing the same results as every other correct implementation.
Only then is the throughput number (candidates per second) a meaningful
comparison. If either signature value differs, the implementation has a bug
and its throughput is not comparable.

### Interface specification

Any implementation targeting search-lab benchmarks **must** conform to this
interface. The `just bench` recipe invokes the binary **once per k value**
(each with `--kmax` equal to `-k`) and pipes stdout through
`parse_bench.py` to produce the results table.

### CLI contract

The benchmark harness invokes the binary as:

```
<binary> -k K --kmax K --start N --end N --threads T --bench 999
```

| Flag | Required | Description |
|------|----------|-------------|
| `-k K` | yes | The k value for this run |
| `--kmax K` | yes | Set equal to `-k` (one k per invocation) |
| `--start N` | yes | First candidate position, inclusive |
| `--end N` | yes | Last candidate position, exclusive. The range `[start, end)` defines the sweep. |
| `--threads T` | yes | Number of worker threads |
| `--bench SECS` | yes (bench) | **Benchmark mode.** The binary must: **(1)** sieve the **entire** range `[start, end)` without terminating early when a witness is found, **(2)** count every witness discovered, and **(3)** disable all file I/O (checkpointing, result logging) so disk writes do not perturb timing. The `SECS` argument is a time-limit safety valve (the harness passes `999`; position-based termination via `--end` takes effect first). |

Without `--bench`, the binary may operate in normal search mode (stop at
the first witness, checkpoint progress, chain k values, etc.). The
benchmark harness always passes `--bench`.

### Output contract

**Stdout** must contain exactly one R-line per k value. Additional lines
starting with `#` are permitted (informational; ignored by the parser).
No other output may appear on stdout.

**R-line format** — seven tab-separated fields, prefixed with the literal
character `R`:

```
R<TAB>k<TAB>witness<TAB>elapsed_s<TAB>candidates<TAB>speed_mcps<TAB>witnesses_found
```

| Field | Type | Format | Description |
|-------|------|--------|-------------|
| `R` | literal | `R` | Line prefix (allows the parser to skip `#` comment lines) |
| `k` | int | `%d` | The k value |
| `witness` | string | | `"bench"` in benchmark mode. In search mode: the witness n, or `"none"`. |
| `elapsed_s` | float | `%.4f` | Wall-clock sieve time in seconds. Must **exclude** prime-table generation. |
| `candidates` | int | `%d` | Exact number of candidates in the sweep: `end − start`. This is a deterministic value for a given range — not an approximation. |
| `speed_mcps` | float | `%.2f` | `candidates / elapsed_s / 1e6` (millions of candidates per second) |
| `witnesses_found` | int | `%d` | Total witnesses verified in the range. Together with `candidates`, this forms the **deterministic accuracy signature** for the run. |

**Example R-line** (k=8 sweep):

```
R	8	bench	5.9359	18339949252	3089.74	40
```

**Stderr** is free for progress output, diagnostics, and warnings.

### Benchmark sweep ranges and expected signature

The `just bench` recipe sweeps the following predetermined ranges. Each
range covers a fixed segment of the number line at that k's sieve depth;
the known minimum OEIS A375077 witness for each k falls within its range.

| k | Range | Candidates | Witnesses |
|--:|:------|-----------:|----------:|
| 8 | `[169,974,626, 18,509,923,878)` | 18,339,949,252 | 40 |
| 9 | `[509,773,922, 19,029,321,766)` | 18,519,547,844 | 3 |
| 10 | `[8,804,882,497, 25,749,999,991)` | 16,945,117,494 | 1 |
| 11 | `[1,062,858,041,585, 1,078,999,999,999)` | 16,141,958,414 | 1 |
| 12 | `[5,042,391,644,646, 5,055,499,999,999)` | 13,108,355,353 | 1 |
| 13 | `[18,247,629,921,842, 18,258,749,999,999)` | 11,120,078,157 | 1 |

**Both columns — Candidates and Witnesses — are the deterministic accuracy
signature.** Any correct implementation must produce these exact values for
the above ranges. A mismatch in either column indicates a bug in the sieve
or witness verification logic. Only after the signature is verified should
throughput numbers be compared.

### Adding a new implementation

1. Write a binary that conforms to the CLI and output contracts above.
2. Add a `build-<impl>` and `bench-<impl>` recipe to the justfile,
   following the existing Rust/C++ patterns.
3. Run `just bench-<impl>` and verify that candidates and witness counts
   match the expected signature table **exactly** before comparing
   throughput.

---

## Algorithm and Optimization Reference

This section documents every optimization technique used in the Rust and C++
implementations, organized top-down from architecture to micro-optimizations.

### 1. Overall Architecture: Bucketed Segmented Sieve

Both implementations use a three-level hierarchy:

**Chunk** (1M candidates) > **Block** (32K candidates) > **Element** (1 candidate)

```
for each chunk [L, L + 1M):
    precompute offsets for all primes (once per chunk)
    bucket large primes by target block

    for each block [B, B + 32K):
        init rem[j] = odd part of (B + j)
        strip small primes (p <= 32K, every block)
        strip bucketed large primes (only those with multiples here)
        scan for k+1 consecutive rem <= 1
        carry overlap buffer for cross-block runs

    if candidate found: exact_check (full governor verification)
```

**Why this structure matters:**
- **Chunking** amortizes prime offset computation across ~32 blocks
- **Blocking** keeps the `rem[]` working set in L1 cache (~256 KB for 32K u64s)
- **Bucketing** ensures large primes (p > 32K) only touch blocks where they
  have multiples, eliminating O(num_primes) wasted iteration per block

### 2. Prime Data Precomputation

Each prime has precomputed constants stored once:

| Field | Purpose | Used |
|-------|---------|------|
| `inv_p` | Modular inverse of p mod 2^64 | Every strip iteration |
| `max_quot` | u64::MAX / p (divisibility threshold) | Every strip iteration |
| `magic` | Barrett magic constant: ceil(2^(64+shift) / p) | Once per prime per chunk |
| `p` | Prime value | Stepping, bucketing |
| `shift` | bit_width(p) - 1 | Barrett reduction |

**Rust:** `#[repr(C)]` with hot fields (`inv_p`, `max_quot`) first so they share
a cache line. Two PrimeData entries (29 bytes each, padded to 32) fit in one
64-byte cache line.

**C++:** `PrimeData` struct with the same fields, compiler-ordered.

### 3. Modular-Inverse Exact Division

The core operation in the strip loop is: "does p divide x? If so, divide."

**Hardware division** (`x / p`, `x % p`): ~30-40 cycles on x86-64, ~20 on ARM64.

**Modular inverse** (`x * inv_p`): ~3-4 cycles (single multiply + compare).

The identity: for odd prime p with modular inverse `inv_p = p^(-1) mod 2^64`:
- `p | x` iff `x * inv_p <= u64::MAX / p`
- When `p | x`: `x / p == x * inv_p` (exact, no remainder needed)

Both implementations precompute `inv_p` via Newton's method (4-5 iterations
from a ~4-bit seed, converging to 64-bit precision):

```
inv = p * 3 ^ 2                        // ~4 correct bits
inv = inv * (2 - p * inv)              // 8 bits
inv = inv * (2 - p * inv)              // 16 bits
inv = inv * (2 - p * inv)              // 32 bits
inv = inv * (2 - p * inv)              // 64 bits
```

### 4. Barrett Reduction for Offset Computation

At the start of each chunk, we compute the first multiple of each prime p
that falls at or after the chunk start L:

```
start_c = ceil(L / p)
offset  = start_c * p - L
```

The naive `L / p` is a hardware division (~20-35 cycles). With ~10,000 primes
at k=13 scale, this adds ~200-350K cycles per chunk.

**Barrett reduction** replaces this with a 128-bit multiply + shift (~3-4 cycles):

```
start_c = (L * magic) >> (64 + shift)
```

where `magic = ceil(2^(64+shift) / p)` is precomputed.

**Rust:** Uses `u128` widening multiply, which compiles to `umulh` on ARM64.

**C++:** Uses `unsigned __int128` with the same codegen.

**Impact:** 10-15x faster offset computation per prime. At k=8 (where offset
computation dominates), this was the single biggest optimization: 260 to
3,800+ M/s.

### 5. Const-Generic / Template Monomorphization

For the 24 smallest primes (3 through 97), the strip function is specialized
at compile time. The compiler replaces modular-inverse constants with immediate
values, eliminating memory loads and enabling instruction scheduling.

**Rust:** `process_prime<const P: u32>()` with const-evaluated `inv_p` and
`max_quot`. The `dispatch_strip!` macro generates a match statement:

```rust
match p {
    3 => process_prime::<3>(&mut sj, w, rem),
    5 => process_prime::<5>(&mut sj, w, rem),
    ...
    97 => process_prime::<97>(&mut sj, w, rem),
    _ => process_prime_dyn(p, inv_p, max_quot, &mut sj, w, rem),
}
```

**C++:** `process_prime_p<p>()` templates with `constexpr` inverse/limit.
A preprocessor macro `PROCESS_PRIME(p)` generates the switch cases.

Primes > 97 fall through to the dynamic path, which uses precomputed values
from `PrimeData` (no recomputation).

### 6. Bucketed Large-Prime Sieve

For primes p > BLOCK_SIZE (32,768), a given block of 32K elements has at most
one multiple of p (since p > block size). Iterating all large primes for every
block wastes time on primes that have no multiples there.

**Solution:** At the start of each chunk, pre-assign each large prime to the
bucket corresponding to the block where its next multiple falls:

```
for each large prime p:
    sj = offset[p]
    while sj < chunk_width:
        buckets[sj >> BLOCK_SHIFT].push(p_index, sj & BLOCK_MASK)
        sj += p
```

Each block then processes only its own bucket. This reduces per-block work from
O(total_primes) to O(bucket_size), which is typically 100-500 items instead of
10,000+.

**Prefetch:** Both implementations prefetch the PrimeData for 8 items ahead in
the bucket loop, since consecutive bucket items may reference non-contiguous
PrimeData entries (cache-unfriendly access pattern).

- **Rust (ARM64):** `prfm pldl2keep, [ptr]` via inline assembly
- **Rust (x86-64):** `_mm_prefetch(ptr, _MM_HINT_T1)`
- **C++:** `__builtin_prefetch(ptr, 0, 1)`

### 7. Factor-of-2 Folded into Initialization

Rather than stripping factors of 2 in the sieve loop, both implementations fold
this into the `rem[]` initialization:

```
rem[j] = x >> trailing_zeros(x)    // = odd part of x
```

This uses a single `ctz` (count trailing zeros) + shift per element, replacing
what would otherwise be a strip loop for p=2. On ARM64, `ctz` compiles to the
`rbit` + `clz` instruction pair (~2 cycles).

### 8. Strip-Only Sieve with Deferred Governor Check

**Key architectural decision:** The sieve loop only strips prime factors from
each element. It does NOT check the governor condition (v_p(n) <= v_p(C(2n,n)))
inline. Instead, the scan phase identifies candidate runs (k+1 consecutive
elements with `rem <= 1`), and the cold-path `exact_check` function performs
the full governor verification.

**Why this is faster:** The governor check (Kummer/Legendre) involves per-element
branching and additional arithmetic that pollutes the instruction pipeline. Since
>99.99% of elements are non-governors (rejected by `rem > 1`), deferring the
check to the rare candidate path avoids this overhead on the hot path.

### 9. Sliding Window Scan with Skip-Ahead

The consecutive-run detection uses a backward scan with skip-ahead:

```
while scan_j + k < w_search:
    i = k
    while i >= 0 and rem[scan_j + i] <= 1:
        i -= 1
    if i < 0:
        exact_check(candidate)       // all k+1 elements passed
        scan_j += 1
    else:
        scan_j += i + 1              // skip past the failing position
```

The skip-ahead (`scan_j += i + 1`) avoids re-checking positions that are known
to fail, reducing average scan work from O(k * n) to O(n).

### 10. Overlap Buffer for Cross-Block Runs

A k+1 consecutive governor run could span two blocks. To handle this, the last
k elements of each block are carried over to the start of the next block's
`rem[]` buffer:

```
rem[0..k] = rem[w_search-k..w_search]    // overlap from previous block
rem[k..k+BLOCK_SIZE] = new block data     // fresh data
```

The scan position (`scan_j`) is adjusted to account for the overlap. This
ensures no cross-boundary runs are missed.

**Rust:** Uses `ptr::copy` (equivalent to `memmove`) for the overlap transfer.
**C++:** Uses a simple element-by-element copy loop.

### 11. Exact Check: Three-Part Governor Verification

Called only on rare candidates (a few dozen times per k). Three phases:

**Part 1 (v_2 check):** Uses popcount (hardware `popcnt` / `cnt` instruction)
to verify the 2-adic condition via Kummer's theorem:
`v_2(C(2n,n)) = popcount(n)`.

**Part 2 (Legendre for p <= 2k):** For each small prime p, computes
v_p(n * (n-1) * ... * (n-k)) and v_p(C(2n,n)) via Legendre's formula and
checks the divisibility condition.

**Part 3 (per-element large prime check):** For each element n-i in the run:
1. Strip factors of 2 (trailing_zeros)
2. Strip all primes p <= 2k using modular-inverse exact division
3. Factor the residual: for each prime p > 2k dividing the residual, compute
   v_p(n-i) and check against v_p(C(2n,n)) via Legendre
4. If the residual is a single large prime > 2k, check it directly

This part was added to support the strip-only sieve architecture (the sieve
no longer checks the governor condition inline).

### 12. Threading Model

**Rust:** `std::thread::scope` with an atomic chunk counter. Workers grab
chunks via `fetch_add` and terminate when `l_chunk > global_min_n` or
`l_chunk >= end_l`. No work-stealing; simple atomic counter distribution.

**C++:** `std::thread` pool with the same atomic counter pattern. Workers
persist across the entire solve call (no per-k thread creation overhead).

Both use `Relaxed` memory ordering for the chunk counter and `global_min_n`
since exact ordering is not required (worst case: a worker processes one
extra chunk before seeing the updated minimum).

### 13. Compiler and Build Optimizations

**Rust (`Cargo.toml` + `.cargo/config.toml`):**
- `opt-level = 3` (maximum optimization)
- `lto = "fat"` (full cross-crate link-time optimization)
- `codegen-units = 1` (single codegen unit for better whole-program optimization)
- `panic = "abort"` (no unwinding overhead)
- `target-cpu = native` (enables M3-specific instruction scheduling, NEON, etc.)

**C++ (Makefile/justfile flags):**
- `-O3 -march=native -mtune=native` (equivalent to Rust's opt-level + target-cpu)
- `-flto` (link-time optimization)
- `-funroll-loops` (aggressive loop unrolling)
- `-fomit-frame-pointer` (frees a register)
- `-fno-exceptions -fno-rtti` (removes C++ runtime overhead)

### 14. Cache Layout Optimization

**PrimeData struct** (Rust): `#[repr(C)]` with fields ordered by access frequency.
`inv_p` and `max_quot` (accessed every strip iteration) come first, sharing the
first 16 bytes of the struct. `magic`, `p`, and `shift` (accessed once per chunk)
follow. At 32 bytes per entry, two PrimeData structs fit in a single 64-byte
cache line.

**BucketItem struct:** 8 bytes (two u32s), 8 items per cache line. Sequential
access pattern in the bucket loop is cache-friendly.

**rem[] buffer:** 32K u64s = 256 KB, sized to fit in L1/L2 cache. The
strip loop accesses `rem[j]` with stride p (one access per p elements),
which has good spatial locality for small primes.

### 15. Platform-Specific Optimizations

**ARM64 / Apple Silicon (Rust):**
- `prfm pldl2keep` inline assembly for bucket prefetch (x86 `_mm_prefetch` is
  dead code on ARM)
- `trailing_zeros()` compiles to `rbit` + `clz` (2 cycles)
- `count_ones()` compiles to `cnt` + horizontal add (~3 cycles)
- `u128` multiply compiles to `umulh` (upper 64 bits, 1 cycle)

**x86-64 (Rust):**
- `_mm_prefetch` with `_MM_HINT_T1` for L2 prefetch
- `trailing_zeros()` compiles to `tzcnt` (1 cycle with BMI)
- `count_ones()` compiles to `popcnt` (1 cycle)

**C++ (both architectures):**
- `__builtin_prefetch` maps to the platform-appropriate prefetch instruction
- `std::countr_zero`, `std::popcount` map to hardware instructions
- `unsigned __int128` for Barrett reduction

---

## Performance History

Relative improvement of each optimization, measured against the preceding
baseline. The starting point was a naive Rust port with runtime hardware
division and no bucketing.

| # | Optimization | Δ k=8 | Δ k=13 | Key insight |
|--:|:-------------|------:|-------:|:------------|
| 0 | Baseline (runtime div, no bucketing) | — | — | |
| 1 | Bucketed sieve architecture | +5% | +0% | Only process primes with multiples in block |
| 2 | Barrett reduction | +14.6× | +1.13× | Replace hardware div with 128-bit mul+shift |
| 3 | ARM64 prefetch + target-cpu=native | +5.6% | +8.5% | `prfm` instruction was entirely missing on ARM |
| 4 | repr(C) cache-aligned PrimeData | +16.8% | −1.6% | Hot fields first in cache line |
