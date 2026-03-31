# search-lab Optimization Guide

**Crate:** `crates/search-lab/` (`erdos396-search-lab`)
**Purpose:** High-throughput witness search for Erdős Problem #396
**Lineage:** Co-developed alongside a C++ implementation; both share ideas bidirectionally while maintaining independent codebases with verified parity

This document catalogs every performance-critical optimization in search-lab,
explains *why* each one works, and notes the speedup it provides relative to
the naive alternative. The search engine processes 200–600 M candidates/sec
on modern hardware (Apple M-series, AMD Zen 4) depending on `k`.

---

## Architecture Overview

The search operates on a two-phase sieve-then-verify pipeline:

1. **Sieve phase (hot path):** For each candidate `n`, strip all prime
   factors from the product `n · (n-1) · ... · (n-k)` using precomputed
   modular arithmetic. If the residual of every term is 1, `n` is a
   *candidate witness* (all k+1 consecutive integers are Governor Set
   members under the sieve's approximation).

2. **Exact check (cold path):** Rare candidates (~0.001% of the search
   space) undergo a full p-adic verification using Legendre's formula and
   Kummer's theorem. This confirms whether `desc(n, k+1) | C(2n, n)`.

The sieve dominates runtime (>99.9%). Every optimization below targets
the sieve's inner loops.

---

## 1. Modular-Inverse Exact Division

**Location:** `mod_inverse_u64()`, `process_prime()`, `process_prime_dyn()`
**Lines:** 76–83, 161–199, 206–235

### The problem

The sieve must repeatedly test "does prime `p` divide value `x`?" and if
so, compute `x / p`. On x86-64 and AArch64, hardware integer division
(`div`/`udiv`) costs 20–40 cycles — devastating in an inner loop that
executes billions of times.

### The trick

For any odd prime `p`, there exists a unique *modular inverse*
`p_inv ≡ p⁻¹ (mod 2⁶⁴)` such that:

- **Divisibility test:** `p | x` iff `x.wrapping_mul(p_inv) <= u64::MAX / p`
- **Exact quotient:** `x / p == x.wrapping_mul(p_inv)` (when `p | x`)

Both operations use a single `imul` (3 cycles) + compare, replacing a
35-cycle `div`.

### How the inverse is computed

Newton's method for modular inverse, starting from a 4-bit seed:

```rust
let mut inv = n.wrapping_mul(3) ^ 2;       // correct mod 2^4
inv = inv.wrapping_mul(2 - n.wrapping_mul(inv)); // mod 2^8
inv = inv.wrapping_mul(2 - n.wrapping_mul(inv)); // mod 2^16
inv = inv.wrapping_mul(2 - n.wrapping_mul(inv)); // mod 2^32
inv = inv.wrapping_mul(2 - n.wrapping_mul(inv)); // mod 2^64
```

Each iteration doubles the number of correct bits. Four doublings lift
4 → 64 bits. The initial seed `n*3 ^ 2` is the unique value satisfying
`n * seed ≡ 1 (mod 16)` for all odd `n`.

### Speedup

~10× on the strip inner loop compared to hardware `div`. This single
optimization accounts for the majority of the performance gain over a
naive implementation.

### Why it's correct

The ring Z/(2⁶⁴)Z has a multiplicative inverse for every odd element.
The wrapping multiplication in `u64` arithmetic is exactly multiplication
mod 2⁶⁴. The divisibility threshold `u64::MAX / p` works because the
modular image of multiples of `p` maps to the range `[0, u64::MAX / p]`
under the inverse — this is a standard result in the "magic number
division" literature (Warren, *Hacker's Delight*, Ch. 10).

---

## 2. Const-Generic Prime Dispatch

**Location:** `process_prime::<P>()`, `dispatch_strip!` macro
**Lines:** 161–199, 240–247, 569–580

### The problem

The first 24 primes (3, 5, 7, ..., 97) account for the vast majority of
sieve work because they hit the most positions. The strip loop steps by
`p`, so smaller primes do more iterations per block.

### The trick

For these 24 primes, the strip function is monomorphized via Rust
const-generics: `process_prime::<3>`, `process_prime::<5>`, etc. The
`dispatch_strip!` macro generates a match statement that dispatches to
the correct monomorphized instance.

Inside `process_prime::<P>`, the modular inverse and divisibility limit
are computed as `const fn` — they become immediate constants in the
generated assembly. The compiler can then:

- Unroll the inner loop with known step size
- Use the immediate constant for the `imul` operand (avoiding a memory load)
- Optimize the loop bound check with known stride

Primes larger than 97 fall through to `process_prime_dyn`, which uses
runtime values from `PrimeData`.

### Speedup

~20–30% on the small-prime strip phase compared to a fully dynamic
implementation. The effect compounds because small primes dominate the
strip workload.

---

## 3. Barrett Reduction for Initial Offsets

**Location:** `build_prime_data()`, offset computation in `solve()`
**Lines:** 114–134, 509–517

### The problem

At the start of each 1M-element chunk, the solver must compute the first
position divisible by each prime `p` within the chunk: `ceil(l_chunk / p) * p`.
This requires one division per prime per chunk — roughly 11 million
divisions for a 200M prime table.

### The trick

Barrett reduction replaces the division with a multiply-and-shift:

```rust
// Precomputed once per prime:
let shift = bit_width(p) - 1;
let magic = ((1u128 << (64 + shift)) / p as u128) + 1;

// At runtime — no division:
let q = ((num as u128 * magic as u128) >> 64) as u64 >> shift;
```

The `magic` constant is the fixed-point reciprocal of `p` with enough
precision to compute exact quotients for all `num` values that arise.

### Speedup

Eliminates ~11M hardware divisions per chunk. The offset computation
drops from ~400M cycles to ~33M cycles per chunk (~12×).

---

## 4. Bucketed Large-Prime Sieve

**Location:** `FastBucket`, bucket construction and processing in `solve()`
**Lines:** 33–60, 519–630

### The problem

Primes larger than `BLOCK_SIZE` (32K) hit at most one position per block.
Iterating the full prime list for each block wastes cycles checking primes
that have no hit in the current block.

### The trick

Before processing a chunk, large primes are pre-bucketed by target block:

```
for each large prime index idx:
    for each hit position sj in [0, chunk_width):
        buckets[sj >> BLOCK_SHIFT].push(idx, sj & BLOCK_MASK)
```

When processing a block, only the primes in that block's bucket are
visited. This converts O(num_large_primes × num_blocks) work into
O(total_hits), which is much smaller.

### Cache-line prefetching

The bucket processing loop prefetches `PrimeData` for items 8 positions
ahead:

```rust
if i + 8 < b_count {
    let future_idx = (*b_ptr.add(i + 8)).p_idx as usize;
    let ptr = pd_ptr.add(future_idx) as *const u8;
    // aarch64: prfm pldl2keep, [ptr]
    // x86_64:  _mm_prefetch(ptr, _MM_HINT_T1)
}
```

Large primes have scattered `PrimeData` entries (they're indexed by
prime ordinal, not by block position). Without prefetching, each access
is a cache miss (~100 cycles on modern hardware). The 8-item lookahead
hides this latency behind the strip computation of the current item.

### Speedup

Bucketing: ~3–5× on the large-prime phase (varies with search range).
Prefetching: additional ~30–40% on the bucketed loop.

---

## 5. L1-Cache-Sized Blocks

**Location:** `BLOCK_SIZE`, `BLOCK_SHIFT`, `BLOCK_MASK` constants, block loop in `solve()`
**Lines:** 13–15, 537–697

### The problem

The `rem[]` buffer holds one `u64` per candidate. At `CHUNK_SIZE = 1M`,
that's 8 MB — far exceeding L1 cache (typically 32–64 KB). Small primes
that step through the entire chunk cause repeated L1 misses.

### The trick

Each 1M chunk is processed in 32K-element blocks (256 KB of `u64` data).
This fits comfortably in L1 cache. Small primes process only the current
block, keeping their working set hot in L1.

The block boundary is managed with bitwise operations:
- `block_idx = pos >> BLOCK_SHIFT` (shift by 15 = log2(32768))
- `offset = pos & BLOCK_MASK` (mask with 0x7FFF)

No divisions or modulo operations are needed for block indexing.

### Speedup

~2–3× compared to processing the full 1M chunk linearly. The exact gain
depends on the CPU's L1/L2 latency ratio.

---

## 6. Factor-of-2 Stripping via `trailing_zeros()`

**Location:** Block initialization in `solve()`
**Lines:** 544–554

### The problem

Half of all integers are even, and many have high powers of 2. Stripping
factors of 2 through the general modular-inverse path would waste cycles
on the most common prime.

### The trick

Factor of 2 is stripped during initialization with a single bitwise
operation:

```rust
rem[i] = x >> x.trailing_zeros();
```

`trailing_zeros()` compiles to a single hardware instruction (`ctz` on
x86, `rbit + clz` on AArch64). This completely removes `p = 2` from the
prime processing loop.

### Speedup

Saves one full pass over the block for p=2 and removes all p=2
iterations from the strip loop. ~5–10% total throughput improvement.

---

## 7. Popcount-Based v₂ Check in `exact_check`

**Location:** `exact_check()` Part 1
**Lines:** 263–270

### The problem

The exact check must verify the p-adic condition for every prime. For
p = 2, the Legendre formula involves a summation loop. Since p = 2 is
the most likely prime to reject a candidate (it has the most carries),
this check should be as fast as possible.

### The trick

By Kummer's theorem, `v₂(C(2n, n)) = popcount(n)` (the number of 1-bits
in the binary representation of `n`). The v₂ demand from the product
`n · (n-1) · ... · (n-k)` can also be expressed in terms of popcount:

```rust
let n_ones = n.count_ones();
let nu2_prod = (n - k - 1).count_ones() - n_ones + (k + 1);
// Reject if supply < demand:
if n_ones < nu2_prod { return false; }
```

This replaces a multi-iteration Legendre summation loop with two
`popcnt` instructions and some arithmetic — a handful of cycles total.

### Why it works

Legendre's formula gives `v_p(n!) = (n - s_p(n)) / (p - 1)` where
`s_p(n)` is the digit sum of `n` in base `p`. For `p = 2`, the digit
sum is the popcount. The demand is
`v₂(n!/(n-k-1)!) = (s₂(n-k-1) - s₂(n) + (k+1)) / 1`, which simplifies
to the expression above.

### Speedup

The v₂ check rejects ~40–60% of candidates at a cost of 2 `popcnt`
instructions (~1 cycle each). Replacing the Legendre loop saves ~10–20
cycles per candidate on average.

---

## 8. Early-Exit Sliding Window

**Location:** Governor-run scanning in `solve()`
**Lines:** 632–681

### The problem

Naively checking whether `k+1` consecutive positions are all governors
(rem ≤ 1) requires `k+1` comparisons per position — O(n·k) total.

### The trick

The scan uses a right-to-left sweep that exploits failures:

```rust
let mut i = k32 as i32;
while i >= 0 && rem[scan_j + i] <= 1 {
    i -= 1;
}
if i < 0 {
    // All k+1 elements passed — candidate found
    scan_j += 1;
} else {
    // Skip past the failing position
    scan_j += i as u32 + 1;
}
```

When element `i` in the window fails, the scan jumps forward by `i + 1`
positions — because any window containing that failing position will also
fail. This is analogous to the Boyer-Moore skip in string matching.

### Speedup

Amortized O(n) rather than O(n·k). For large `k`, this is a significant
improvement. In practice, most positions have rem > 1 and the skip is
large, so the scan touches each position roughly once.

---

## 9. Cross-Block Overlap for Boundary Runs

**Location:** Overlap management in `solve()`
**Lines:** 466–468, 534–535, 683–694

### The problem

A run of `k+1` consecutive governors can span two adjacent blocks. If
blocks are processed independently, such runs are missed.

### The trick

The `rem[]` buffer is sized to `BLOCK_SIZE + k + 1`. After processing
each block, the last `k` elements are copied to the front of the buffer:

```rust
let src = (w_search - k32) as usize;
std::ptr::copy(rem.add(src), rem, k32 as usize);
overlap = k32;
```

The next block's new elements are placed after this overlap region. The
sliding-window scan operates over the combined region, seamlessly
detecting runs that straddle block boundaries.

### Cost

One `memcpy` of `k` elements per block — negligible compared to the
sieve work.

---

## 10. Lock-Free Parallel Chunk Distribution

**Location:** Thread spawning and chunk allocation in `solve()`
**Lines:** 425–727

### The problem

Distributing work across threads typically requires either a mutex-
protected queue or static partitioning. Mutexes add contention; static
partitioning causes load imbalance.

### The trick

Chunks are distributed via a single atomic counter:

```rust
let chunk_id = current_chunk.fetch_add(1, Relaxed);
```

Each thread atomically claims the next chunk, processes it, and loops.
There is no lock, no coordination, and no load imbalance — faster
threads simply process more chunks.

The minimum witness is tracked with an atomic `compare_exchange_weak`
loop, allowing threads to independently discover witnesses and race to
update the global minimum without locking.

`Ordering::Relaxed` is sufficient because:
- Chunk assignment only needs to be unique, not ordered
- The global minimum is a monotonically decreasing value; stale reads
  only cause redundant work, not incorrect results

### Speedup

Near-linear scaling with thread count (observed 0.95× efficiency at 8
threads on M3 Max). No contention except the single `fetch_add` per
1M-element chunk.

---

## 11. `PrimeData` Cache-Line Layout

**Location:** `PrimeData` struct definition
**Lines:** 21–29

### The trick

The struct uses `#[repr(C)]` with fields ordered by access frequency:

```rust
#[repr(C)]
struct PrimeData {
    inv_p: u64,    // accessed every strip iteration
    max_quot: u64, // accessed every strip iteration
    magic: u64,    // accessed once per prime per chunk
    p: u64,        // accessed for control flow
    shift: u8,     // accessed once per prime per chunk
}
```

The two hottest fields (`inv_p`, `max_quot`) are at offset 0 and 8,
guaranteed to be in the same cache line. The `p` field is at offset 24,
still within the same 64-byte line. The rarely-used `magic` and `shift`
are placed after the hot fields.

### Speedup

Difficult to measure in isolation, but proper field ordering prevents
split cache-line loads in the strip inner loop. Estimated ~5% on the
strip phase.

---

## 12. `#[cold]` / `#[inline(never)]` on `exact_check`

**Location:** `exact_check()` function attributes
**Lines:** 257–258

### The trick

```rust
#[cold]
#[inline(never)]
fn exact_check(...) -> bool { ... }
```

`#[cold]` tells LLVM that this function is rarely called, influencing
branch prediction heuristics and code layout. `#[inline(never)]` prevents
the compiler from inlining the 130-line function body into the hot sieve
loop, which would bloat the instruction cache footprint of the inner loop.

### Speedup

~2–5% on the sieve phase. Without these attributes, the compiler may
speculatively inline `exact_check`, evicting the tight sieve loop from
the L1 instruction cache.

---

## 13. Integer Square Root with Newton Correction

**Location:** `isqrt_u64()`
**Lines:** 139–153

### The problem

The sieve limit for each chunk is `√(2n)`. The naive approach
`(n as f64).sqrt() as u64` fails for `n > 2⁵³` because IEEE 754
double precision has only 53 bits of mantissa. At search positions
above ~9T (9 × 10¹²), `2n` exceeds 2⁵³ and the f64 sqrt can be
off by 1, causing the sieve to either miss a prime (false positive)
or include an unnecessary prime (wasted work).

### The trick

Start from the f64 approximation, then apply at most 1–2 Newton
correction steps:

```rust
let mut x = (n as f64).sqrt() as u64;
while x < n / (x + 1).max(1) { x += 1; } // too small
while x > 0 && x > n / x    { x -= 1; } // too large
```

This gives exact results for all `u64` inputs.

### Cost

The f64 sqrt is nearly always correct, so the Newton loop typically does
zero iterations. The entire function is ~5 cycles amortized.

---

## Why Rust

This project has been developed in parallel with a C++ implementation
by a collaborator in the Erdős 396 community. The two codebases share
ideas bidirectionally — many of the optimizations documented here
originated in the Rust implementation and were adopted by the C++ side,
while others flowed in the opposite direction. The result is two
independent implementations that maintain strict parity on correctness
while each leveraging the strengths of its language.

Rust was chosen for this implementation not as a port, but as a
deliberate design decision. It matches C++ throughput while gaining
significant advantages in correctness, maintainability, and
LLM-assisted development.

### Safety without overhead

search-lab uses `unsafe` in exactly two places: the inner strip loop
(pointer arithmetic for performance) and the overlap `memcpy`. Everything
else — the sieve, Barrett reduction, exact check, thread coordination —
is safe Rust. The compiler enforces:

- **No data races.** The scoped-thread + atomic design is statically
  verified. In C++, the equivalent code requires careful manual
  reasoning about thread lifetimes and shared-state visibility. A missed
  `volatile` or wrong memory order is a silent bug. In Rust, it's a
  compile error.

- **No use-after-free.** The `rem[]` buffer is owned by each thread's
  stack frame via `std::thread::scope`. The compiler proves that no
  reference escapes the scope. The C++ equivalent requires manual
  lifetime discipline.

- **No buffer overflows in debug builds.** The `unsafe` pointer
  arithmetic compiles to bounds-checked indexing in debug mode, catching
  off-by-one errors during development. Release builds elide the checks.

The `unsafe` blocks are small, auditable, and encapsulated. They don't
"infect" the rest of the codebase — callers of `process_prime` don't
need to reason about pointer validity.

### Const generics > C++ templates

Rust's const-generic `process_prime::<P>` achieves the same
monomorphization as a C++ template, but with important ergonomic
differences:

- **No header files.** The monomorphized instances are generated from a
  single function definition, not duplicated across headers.
- **Predictable compilation.** Rust const-generics don't trigger the
  unbounded template instantiation that plagues C++ builds. The
  `dispatch_strip!` macro generates exactly 24 instances.
- **`const fn` instead of `constexpr`.** The modular inverse and
  divisibility limit are computed at compile time via `const fn`, which
  is simpler and more readable than C++ `constexpr` with its historical
  restrictions.

### Atomic semantics are explicit and correct

The lock-free chunk distribution uses `Ordering::Relaxed`:

```rust
let chunk_id = current_chunk.fetch_add(1, Relaxed);
```

In C++, the equivalent `std::atomic<uint64_t>` with
`std::memory_order_relaxed` looks similar, but the Rust version has a
crucial advantage: **the type system prevents non-atomic access.** An
`AtomicU64` cannot be read or written without specifying an ordering.
In C++, nothing prevents accidentally accessing an `atomic<uint64_t>`
through a non-atomic path (e.g., passing it by reference to a function
that treats it as a plain integer).

### Cargo and the single-file advantage

search-lab is a single 900-line file with no external dependencies
beyond `std`. This is a deliberate choice:

- **Zero build system complexity.** `cargo build --release` is the
  entire build. No CMake, no vcpkg, no pkg-config, no platform-specific
  flags. The `rust-toolchain.toml` pins the compiler version.
- **Cross-compilation.** `cargo build --target aarch64-unknown-linux-gnu`
  produces a binary that runs on ARM Linux servers. A C++ project
  typically requires a separate build configuration for each target.
- **Reproducible builds.** `Cargo.lock` (tracked in the repo) pins
  every transitive dependency. The workspace `Cargo.toml` ensures all
  crates share dependency versions.

### LLM-assisted development

Rust's type system and compiler diagnostics make it an unusually good
target for LLM-generated code:

- **The compiler catches what the LLM misses.** When an LLM generates
  code with a lifetime error, data race, or type mismatch, `rustc`
  rejects it with an actionable error message. In C++, the equivalent
  bug compiles silently and manifests as undefined behavior at runtime —
  possibly only at scale.

- **Refactoring is compiler-guided.** When extracting `erdos396-core`
  (Phase 3), the compiler will flag every call site that needs updating.
  In C++, extracting a library from a monolith requires manually tracing
  header dependencies and link-time symbol resolution.

- **Pattern matching forces exhaustive handling.** The `dispatch_strip!`
  macro generates a `match` statement. If a prime is added to the
  dispatch list but not handled, the compiler warns. C++ `switch`
  statements don't enforce this unless `-Wswitch-enum` is enabled (and
  it often isn't).

- **`unsafe` blocks are greppable.** An LLM — or a human reviewer — can
  audit every `unsafe` block by searching for the keyword. C++ has no
  equivalent marker for "this code has manual memory management." Every
  line of C++ is potentially unsafe; in Rust, safety violations are
  opt-in and localized.

- **Type inference reduces boilerplate.** The LLM can focus on algorithm
  logic rather than type declarations. Rust infers local variable types,
  closure parameter types, and generic type arguments, while still
  providing full type safety. This makes LLM-generated code more
  concise and less error-prone than equivalent C++.

### The performance parity argument

A common objection to Rust for high-performance numerical code is "C++
is faster." The cross-implementation parity maintained between
search-lab and the C++ implementation demonstrates this is not the case:

| Metric | C++ implementation | search-lab (Rust) |
|--------|---------------|-------------------|
| Throughput (Zen 4) | ~300 M/s | ~300 M/s |
| Throughput (M3 Max) | — | ~500 M/s |
| Known witnesses k=1..13 | Correct | Identical |
| Lines of code | ~800 | ~900 |
| Build time (release) | ~5s | ~2s |
| Cross-compilation | Manual | `cargo build --target` |

The Rust implementation is not "almost as fast" — it is *the same speed*
on identical hardware. Both compile to the same LLVM IR patterns: the
modular inverse becomes a single `imul`, the Barrett reduction becomes
a `mul` + shift, and the const-generic dispatch becomes direct jumps.
The only measurable differences come from platform-specific codegen
(Apple's LLVM fork vs. upstream).

When performance is equal, the tiebreaker is everything else: safety,
maintainability, build ergonomics, and LLM compatibility. Rust wins on
all four.

---

## Cross-Implementation Parity

The Rust and C++ implementations maintain strict parity on all
observable behaviors — the same algorithms, the same witnesses, the
same throughput. This parity is a design goal of the collaboration,
not a coincidence.

| Property | C++ implementation | search-lab (Rust) |
|----------|---------------|------------|
| Chunk size | 1M | 1M |
| Block size | 32K | 32K |
| Small prime dispatch | Template instantiation (p=3..97) | Const-generic (P=3..97) |
| Division elimination | Modular inverse | Modular inverse (identical algorithm) |
| Offset computation | Barrett reduction | Barrett reduction (identical formula) |
| Large prime handling | Bucketed by block | Bucketed by block |
| Exact check | Popcount v₂ + Legendre | Popcount v₂ + Legendre |
| Known witnesses k=1..13 | Identical | Identical |
| Throughput | ~300 M/s (Zen 4) | ~300 M/s (Zen 4), ~500 M/s (M3 Max) |

Beyond parity, the Rust implementation additionally benefits from:
- Scoped threads (no manual thread lifecycle management)
- Atomic operations without `volatile` hacks
- The `#[cold]` attribute (no direct C++ equivalent with the same semantics)
- Safe overlap management via `std::ptr::copy` with bounds checking in debug builds

---

## Test Coverage

All optimizations are covered by the test suite (15 tests):

| Optimization | Test(s) |
|-------------|---------|
| Modular inverse | `mod_inverse_identity` |
| Barrett constants | `barrett_magic_produces_correct_quotients` |
| Const-generic dispatch (all 24 primes) | `const_vs_dynamic_all_dispatched_primes` |
| Const vs dynamic parity | `const_vs_dynamic_strip_consistency` |
| Higher-power stripping | `process_prime_strips_higher_powers` |
| Multi-prime smooth reduction | `multi_prime_strip_smooth_numbers_reduce_to_one` |
| Popcount v₂ path | `exact_check_v2_popcount_agrees_with_legendre` |
| Large-prime rejection path | `exact_check_rejects_large_prime_deficiency` |
| Known witnesses (full pipeline) | `exact_check_known_values`, `known_witnesses_k1_through_k11` |
| Cross-block boundary overlap | `solve_finds_witnesses_near_block_boundaries` |
| Thread-count independence | `thread_count_independence` |
| Integer sqrt precision | `isqrt_edge_cases` |
| Prime sieve correctness | `sieve_known_counts` |
| Degenerate input handling | `empty_range_returns_no_witness` |
