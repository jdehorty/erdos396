//! # Erdős Problem #396 — High-Performance Solver
//!
//! This crate provides optimized algorithms for finding witnesses to Erdős Problem #396:
//! For each k, find the smallest n such that n(n-1)(n-2)···(n-k) | C(2n, n).
//!
//! ## Key Discovery: Governor Set Approach
//!
//! All known witnesses have every block term in the **Governor Set** G = {n : n | C(2n, n)}.
//! This reduces the search to finding runs of consecutive Governor Set members.
//!
//! ## Optimizations
//!
//! - **Batch prefilter**: Segmented sieve rejects ~69% of candidates cheaply
//!   (odd primes + large prime factor barrier from Sanna/Pomerance)
//! - **Early-exit governor check**: Interleaves p-adic check with factorization
//! - **√(2n) barrier**: Immediate rejection when largest prime factor > √(2n)

pub mod build_info;
pub mod checkpoint;
pub mod factor;
pub mod false_positive;
pub mod governor;
pub mod int_math;
pub mod prefilter;
pub mod run_collection;
pub mod search;
pub mod sieve;
pub mod sieve_solver;
pub mod verify;

// Re-export main types for convenience
pub use build_info::BuildInfo;
pub use checkpoint::{Checkpoint, CheckpointManager, CounterExampleInfo, RunLogger};
pub use factor::Factorization;
pub use false_positive::{
    classify_run_windows, ClassificationConfig, ClassificationStats, RunWindowClassification,
    VerificationAudit,
};
pub use governor::GovernorChecker;
pub use run_collection::{
    build_run_corpus, BuildRunCorpusConfig, BuildRunCorpusStats, RunRecord, RunWindowRecord,
};
pub use search::{SearchConfig, SearchResult, SearchWorker};
pub use sieve::{build_prime_data, PrimeData, PrimeSieve};
pub use verify::WitnessVerifier;

/// Known OEIS A375077 witnesses for verification
/// These are the smallest n such that n(n-1)...(n-k) | C(2n, n)
pub const KNOWN_WITNESSES: &[(u32, u64)] = &[
    (1, 2),
    (2, 2_480),
    (3, 8_178),
    (4, 45_153),
    (5, 3_648_841),
    (6, 7_979_090),
    (7, 101_130_029),
    (8, 339_949_252),
    (9, 1_019_547_844),       // Corrected 2026-03-16 by full validate pass
    (10, 17_609_764_994),     // Found 2025-01-20! (same run of 11!)
    (11, 1_070_858_041_585),  // Found 2026-02-09 (run of 12)
    (12, 5_048_891_644_646),  // Found 2026-02-11 (run of 13)
    (13, 18_253_129_921_842), // Found 2026-02-16, confirmed 2026-02-21 (run of 14, exhaustive 15T-25T)
];

/// Known runs of 9 consecutive Governor Set members found during k=9 search
/// These are the ENDING positions of each run (i.e., [n-8, n] are all governors)
pub const KNOWN_RUNS_OF_9: &[u64] = &[
    2_583_896_932,  // run [2_583_896_924 .. 2_583_896_932]
    3_908_897_710,  // run [3_908_897_702 .. 3_908_897_710]
    4_041_886_053,  // run [4_041_886_045 .. 4_041_886_053]
    6_156_311_784,  // run [6_156_311_776 .. 6_156_311_784]
    6_289_399_016,  // run [6_289_399_008 .. 6_289_399_016]
    7_281_939_718,  // run [7_281_939_710 .. 7_281_939_718]
    7_692_223_001,  // run [7_692_222_993 .. 7_692_223_001]
    7_977_858_245,  // run [7_977_858_237 .. 7_977_858_245]
    8_479_362_120,  // run [8_479_362_112 .. 8_479_362_120]
    11_575_727_913, // run [11_575_727_905 .. 11_575_727_913]
];

/// Known runs of 10+ consecutive Governor Set members (k=9 witness candidates)
pub const KNOWN_RUNS_OF_10_PLUS: &[(u64, u32)] = &[
    (17_609_764_993, 11), // run of 11! [17_609_764_984 .. 17_609_764_994] - valid k=9 witness, but not minimal
    (17_842_967_551, 10), // run of 10 - FALSE POSITIVE (p=3 fails)
    (19_295_451_546, 10), // run of 10 - valid k=9 witness
    (29_661_303_231, 10), // run of 10 - valid k=9 witness
];

/// False positive runs of 10+ - Governor runs that fail p-adic verification
pub const FALSE_POSITIVE_RUNS: &[(u64, u32, u32)] = &[
    (17_842_967_551, 10, 3), // run of 10, fails at p=3: demand=10, supply=8
];

/// Approximate density of the Governor Set (Ford-Konyagin)
pub const GOVERNOR_DENSITY: f64 = 0.114;

/// Error types for the crate
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("Checkpoint error: {0}")]
    Checkpoint(String),

    #[error("Audit failure: {0}")]
    Audit(String),
}

pub type Result<T> = std::result::Result<T, Error>;
