# Search Algorithm Optimizations for A375077

In the past couple of months, the community has made considerable progress on [Erdős Problem #396](https://www.erdosproblems.com/396) and the associated OEIS sequence [A375077](https://oeis.org/A375077). Terms $a(8)$ through $a(14)$ have been identified:

| $k$ | $a(k)$ |
|-----|--------|
| 1 | 2 |
| 2 | 2,480 |
| 3 | 8,178 |
| 4 | 45,153 |
| 5 | 3,648,841 |
| 6 | 7,979,090 |
| 7 | 101,130,029 |
| 8 | 339,949,252 |
| 9 | 1,019,547,844 |
| 10 | 17,609,764,994 |
| 11 | 1,070,858,041,585 |
| 12 | 5,048,891,644,646 |
| 13 | 18,253,129,921,842 |
| 14 | 359,503,904,702,190 |

It is entirely possible that $a(15)$ is another order of magnitude beyond $a(14)$. A formal verification of minimality for such a candidate would require immense computational effort with the current search algorithms. To help organize future progress, this document catalogs every significant optimization that has been discussed in the [forum thread](https://www.erdosproblems.com/forum/thread/396) or implemented in either reference codebase ([jdehorty/erdos396](https://github.com/jdehorty/erdos396) in Rust, [sharky564/ErdosProblems](https://github.com/sharky564/ErdosProblems) in C++), with attribution and measured impact where available.

---

## Catalogue of Significant Optimizations

Impact is measured single-threaded at $n \approx 10^{10}$ unless noted otherwise.

| # | Optimization | Origin | Impact |
|---|---|---|---|
| 1 | Kummer carry check | Justin Dehorty (Feb); Max Alekseyev indep. (Mar) | Foundational |
| 2 | Governor Set sieve | Justin Dehorty (Feb) | Foundational |
| 3 | $\sqrt{2n}$ small-prime barrier | Justin Dehorty (Feb); Max Alekseyev indep. (Mar) | Foundational |
| 4 | Fused segmented sieve | Justin Dehorty (Feb) | Major |
| 5 | $v_2\binom{2n}{n} = \text{popcount}(n)$ | Justin Dehorty (Feb) | Major |
| 6 | Compile-time constant divisors | Justin Dehorty (Mar); Sharvil Kesarwani indep. (Mar) | Major |
| 7 | Branchless Kummer (all primes) | Justin Dehorty (Feb) | Major |
| 8 | Modular inverse exact division | Sharvil Kesarwani (Mar); Justin Dehorty indep. (Mar) | Major |
| 9 | Boyer-Moore smoothness search | Justin Dehorty (Mar); Sharvil Kesarwani indep. (Mar) | Major |
| 10 | Skip primes with no multiples | Sharvil Kesarwani (Mar) | Moderate |
| 11 | Barrett reduction | Sharvil Kesarwani (Mar) | Moderate |

---

## Descriptions

1. **[Kummer carry check](https://github.com/jdehorty/erdos396/commit/a997d3b63d3eed76a93d82f35329d8f9aaa695fd).** Replace direct computation of $\binom{2n}{n}$ with per-prime Kummer carry counting: $v_p(\binom{2n}{n}) = s_p(n)$, the number of carries when adding $n + n$ in base $p$. As a corollary, no odd prime is a governor ($v_p(\binom{2p}{p}) = 0$ for odd $p$), which provides an additional ~4.3% free rejection rate.

2. **[Governor Set sieve](https://github.com/jdehorty/erdos396/commit/a997d3b63d3eed76a93d82f35329d8f9aaa695fd).** Restrict search to runs of $k+1$ consecutive members of $G = \{n : n \mid \binom{2n}{n}\}$. Governor density $\approx 11.4\%$ (Ford-Konyagin), so this rejects ~88% of integers before any heavy work.

3. **[Small-prime barrier](https://github.com/jdehorty/erdos396/commit/c219a1b29c2d4a92b73c00e705ad2cbc6d84c855)** ($\sqrt{2n}$ rejection). If the largest prime factor of $n$ exceeds $\sqrt{2n}$, then $n \notin G$. Falls out naturally from the sieve: after stripping all primes up to the sieve limit, any remaining factor $> 1$ triggers rejection.

4. **[Fused segmented sieve](https://github.com/jdehorty/erdos396/commit/a997d3b63d3eed76a93d82f35329d8f9aaa695fd).** Single pass strips small prime factors AND checks the Kummer condition $v_p(n) \leq s_p(n)$ inline, rejecting ~87% of integers. Eliminates all redundant trial division between prefilter and governor check.

5. **[$v_2(\binom{2n}{n}) = \text{popcount}(n)$](https://github.com/jdehorty/erdos396/commit/a997d3b63d3eed76a93d82f35329d8f9aaa695fd).** Kummer's theorem at $p = 2$ reduces to a single POPCNT hardware instruction, replacing ~34 iterations of the carry loop.

6. **[Compile-time constant divisors](https://github.com/jdehorty/erdos396/commit/f0dfb4169c5fb726f29753612ca9bd0fea98abad).** Replace runtime div instructions (~30-40 cycles) with multiply+shift (~3-4 cycles) for known small primes. The Rust implementation uses const generics for primes 2-47; the C++ implementation ([`85e3a88`](https://github.com/sharky564/ErdosProblems/commit/85e3a8859150cba763849f9c595aab328869d1d9)) uses a bespoke libdivide-style approach, reporting a 33% speedup.

7. **[Branchless Kummer carry for all primes](https://github.com/jdehorty/erdos396/commit/a997d3b63d3eed76a93d82f35329d8f9aaa695fd).** Generic single-loop Kummer with branchless carry: for any prime, the carry is always 0 or 1. Runs in $\sim\log_p(n)$ iterations vs Legendre's $\sim 2\log_p(n)$.

8. **[Modular inverse exact division](https://github.com/sharky564/ErdosProblems/commit/23a22eea62232d732441591d1353bc686668daf2).** For the factor-stripping loop: use the multiplicative inverse $p^{-1} \pmod{2^{64}}$ for exact division (single imul, no shift). Newton-Raphson iteration computes the inverse. The Rust implementation ([`4956323`](https://github.com/jdehorty/erdos396/commit/49563233031c872411bc7d2e0f6d159d967ee1a2)) was developed independently.

9. **[Boyer-Moore smoothness search](https://github.com/sharky564/ErdosProblems/commit/3a10d4d0ed7ff0c0574317677cbc35aee809666d).** Apply Boyer-Moore-style skipping to the smoothness array, reducing the number of elements inspected when scanning for fully-factored candidates.

10. **[Skip primes with no multiples in interval](https://github.com/sharky564/ErdosProblems/commit/95c469ee379d59061b245a830c2577195ced3332).** For primes near the sieve limit, check whether the current batch contains any multiples before entering the sieve loop. Avoids touching memory for primes that contribute nothing.

11. **[Barrett reduction](https://github.com/sharky564/ErdosProblems/commit/85e3a8859150cba763849f9c595aab328869d1d9).** Pre-computed Barrett constants for modular reduction in the sieve inner loop, avoiding hardware div for medium-sized primes.

---

## Combined Effect

Starting from a naive per-candidate governor check (~1M candidates/sec), these optimizations collectively achieve ~3 billion candidates/sec on modern hardware (Rust, 16 threads, ARM64). The C++ implementation independently achieved ~2 billion candidates/sec on an AMD Ryzen 7 with a partially overlapping set of techniques. wibble132 ran the C++ implementation on a 9800x3d processor, achieving ~1,366 M candidates/sec and finding $a(14)$ after approximately 2.5 days of compute. The convergence of two independent implementations on the same core ideas gives confidence in the approach.

Both implementations are organized side by side in the [`crates/search-lab/`](https://github.com/jdehorty/erdos396/tree/main/crates/search-lab) crate under a shared benchmark harness. The full codebase, validation certificates, and Lean 4 formal proofs are available at [github.com/jdehorty/erdos396](https://github.com/jdehorty/erdos396).

---

If anyone notices an error in the above or if I missed a major optimization technique, please let me know and I will adjust the list accordingly. Compiling this involved going through a large number of commits across multiple repositories as well as the lengthy [forum comment section](https://www.erdosproblems.com/forum/thread/396) for this problem, with the help of AI and manual spot-checking. I did my best to give attribution accurately, but mistakes are possible. Thanks to Sharvil Kesarwani, Max Alekseyev, and wibble132 for their contributions to this effort.
