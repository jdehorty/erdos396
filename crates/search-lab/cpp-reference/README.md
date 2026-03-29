# C++ Reference Implementation

**Author:** Sharvil Kesarwani ([@sharky564](https://github.com/sharky564))
**Source:** [sharky564/ErdosProblems](https://github.com/sharky564/ErdosProblems/blob/a82d3380e3010b1016883b7f780ef25e2c32e39a/erdos_problem_396.cpp)

The original high-performance C++ implementation of the Erdos Problem #396 witness
search. Included here as a convenient reference, as I have been collaborating
extensively with Sharvil on improvements to this algorithm.

Key techniques in this implementation:
- Precomputed `PrimeData` structs with magic constants and modular inverses
- Bucket-based sieve for cache-efficient prime processing
- Three-phase Kummer kernel (reject / accept / full check)
- Modular-inverse exact division (replacing hardware `div`)

## Building

```sh
just build
```
