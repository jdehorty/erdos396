//! # Cedar-Style Differential Random Testing for Erdős Problem #396
//!
//! Validates that search-lab's optimized code faithfully implements the proven
//! theorems by comparing against independent oracles:
//!
//! 1. **Lean oracle** (primary) — compiled Lean 4 executable reimplementing the
//!    same algorithms (Legendre, Kummer) that the formal proofs in
//!    `SmallPrimeBarrier/` verify. Not formally linked to the proofs, but shared
//!    algorithmic design + DRT self-tests provide high confidence.
//! 2. **erdos396 library** (secondary) — separate Rust code path using factorization
//!    rather than modular-inverse sieve tricks
//!
//! ## Trust Chain
//!
//! - Lean proofs (`SmallPrimeBarrier/`): algorithms = mathematical definitions
//!   (formally verified, zero sorry's)
//! - Lean oracle (`Erdos396FFI/`): reimplements same algorithms (code-reviewed,
//!   not formally linked to proofs)
//! - DRT: Rust search-lab = Lean oracle (empirically tested on thousands of inputs)
//! - Therefore: Rust code = mathematical definitions (high confidence, transitive)

pub mod generators;
pub mod oracle;
pub mod properties;

/// Known OEIS A375077 witnesses, shared across all tests.
pub const KNOWN_WITNESSES: &[(u32, u64)] = &[
    (1, 2),
    (2, 2_480),
    (3, 8_178),
    (4, 45_153),
    (5, 3_648_841),
    (6, 7_979_090),
    (7, 101_130_029),
    (8, 339_949_252),
    (9, 1_019_547_844),
    (10, 17_609_764_994),
];

/// Fixed seed for reproducibility.
pub const DRT_SEED: u64 = 0x0E4D_0539_6202_6000;
