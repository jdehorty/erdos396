//! # Cedar-Style Differential Random Testing for Erdős Problem #396
//!
//! Validates that search-lab's optimized code faithfully implements the proven
//! theorems by comparing against independent oracles:
//!
//! 1. **Lean oracle** (primary) — compiled Lean 4 executable implementing the same
//!    algorithms the formal proofs verify, but in a different language
//! 2. **erdos396 library** (secondary) — separate Rust code path using factorization
//!    rather than modular-inverse sieve tricks
//!
//! ## Trust Chain
//!
//! - Lean proofs: efficient algorithms = mathematical definitions (formally verified)
//! - DRT: Rust code = Lean oracle (empirically tested on millions of inputs)
//! - Therefore: Rust code = mathematical definitions (transitive)

pub mod generators;
pub mod oracle;
pub mod properties;

/// Fixed seed for reproducibility. Printed at test start.
pub const DRT_SEED: u64 = 0x0E4D_0539_6202_6000;

/// Test tier configuration.
pub struct TierConfig {
    pub name: &'static str,
    pub total_inputs: usize,
    pub use_lean: bool,
}

impl TierConfig {
    /// T0: Smoke — < 10s, 10K inputs, Rust-only oracles.
    pub fn smoke() -> Self {
        Self {
            name: "T0-smoke",
            total_inputs: 10_000,
            use_lean: false,
        }
    }

    /// T1: Standard — < 60s, 100K inputs, includes Lean oracle.
    pub fn standard() -> Self {
        Self {
            name: "T1-standard",
            total_inputs: 100_000,
            use_lean: true,
        }
    }

    /// T2: Full — < 10min, 1M inputs, Lean oracle + high volume.
    pub fn full() -> Self {
        Self {
            name: "T2-full",
            total_inputs: 1_000_000,
            use_lean: true,
        }
    }
}
