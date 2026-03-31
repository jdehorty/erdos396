//! Oracle implementations for differential testing.
//!
//! Each oracle independently computes the same mathematical functions, using
//! different code paths. Any disagreement between oracles is a bug.

use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};

// ---------------------------------------------------------------------------
// Oracle trait
// ---------------------------------------------------------------------------

/// A reference oracle that can answer mathematical queries about Erdős 396.
pub trait Oracle {
    fn name(&self) -> &str;

    /// v_p(n!) via Legendre's formula.
    fn vp_factorial(&mut self, n: u64, p: u64) -> u64;

    /// v_p(C(2n, n)).
    fn vp_central_binom(&mut self, n: u64, p: u64) -> u64;

    /// Is n in the Governor Set? (n | C(2n, n))
    fn is_governor(&mut self, n: u64) -> bool;

    /// Is n a k-witness? (n(n-1)...(n-k) | C(2n, n))
    fn check_witness(&mut self, k: u32, n: u64) -> bool;
}

// ---------------------------------------------------------------------------
// Erdos396 library oracle (separate Rust code path)
// ---------------------------------------------------------------------------

/// Uses the erdos396 library crate — factorization-based, completely separate
/// from search-lab's modular-inverse sieve.
pub struct Erdos396Oracle {
    checker: erdos396::GovernorChecker,
    verifier: erdos396::WitnessVerifier,
}

impl Erdos396Oracle {
    pub fn new(max_n: u64) -> Self {
        let sieve = erdos396::PrimeSieve::for_range(max_n);
        Self {
            checker: erdos396::GovernorChecker::with_sieve(sieve.clone()),
            verifier: erdos396::WitnessVerifier::with_sieve(sieve),
        }
    }
}

impl Oracle for Erdos396Oracle {
    fn name(&self) -> &str {
        "erdos396-lib"
    }

    fn vp_factorial(&mut self, n: u64, p: u64) -> u64 {
        erdos396::governor::vp_factorial(n, p)
    }

    fn vp_central_binom(&mut self, n: u64, p: u64) -> u64 {
        erdos396::governor::vp_central_binom_dispatch(n, p)
    }

    fn is_governor(&mut self, n: u64) -> bool {
        self.checker.is_governor(n)
    }

    fn check_witness(&mut self, k: u32, n: u64) -> bool {
        self.verifier.is_witness(k, n).unwrap_or(false)
    }
}

// ---------------------------------------------------------------------------
// Search-lab oracle (the code under test)
// ---------------------------------------------------------------------------

/// Uses search-lab's optimized functions via drt_exports.
pub struct SearchLabOracle {
    ctx: erdos396_search_lab::drt_exports::DrtContext,
}

impl Default for SearchLabOracle {
    fn default() -> Self {
        Self::new()
    }
}

impl SearchLabOracle {
    pub fn new() -> Self {
        let ctx = erdos396_search_lab::drt_exports::DrtContext::new(200_000_000);
        Self { ctx }
    }
}

impl Oracle for SearchLabOracle {
    fn name(&self) -> &str {
        "search-lab"
    }

    fn vp_factorial(&mut self, n: u64, p: u64) -> u64 {
        // search-lab doesn't export a standalone vp_factorial — it's inlined
        // in exact_check. Use Legendre formula directly.
        let (mut v, mut pw) = (0u64, p);
        while pw <= n {
            v += n / pw;
            if pw > n / p {
                break;
            }
            pw *= p;
        }
        v
    }

    fn vp_central_binom(&mut self, n: u64, p: u64) -> u64 {
        let vf_2n = self.vp_factorial(2 * n, p);
        let vf_n = self.vp_factorial(n, p);
        vf_2n - 2 * vf_n
    }

    fn is_governor(&mut self, _n: u64) -> bool {
        // search-lab doesn't export a standalone governor check.
        // For DRT governor comparisons, we compare erdos396-lib vs Lean instead.
        // This method panics to catch accidental usage.
        unimplemented!("Use Erdos396Oracle or LeanOracle for governor checks")
    }

    fn check_witness(&mut self, k: u32, n: u64) -> bool {
        self.ctx.exact_check(n, k as u64)
    }
}

// ---------------------------------------------------------------------------
// Lean oracle (subprocess — compiled Lean 4 executable)
// ---------------------------------------------------------------------------

/// Communicates with the compiled Lean 4 oracle via stdin/stdout.
/// The Lean oracle implements the same algorithms the formal proofs verify.
pub struct LeanOracle {
    child: Child,
    stdin: std::process::ChildStdin,
    reader: BufReader<std::process::ChildStdout>,
}

impl LeanOracle {
    /// Find the Lean oracle binary. Checks:
    /// 1. `LEAN_ORACLE_PATH` environment variable
    /// 2. `formal/.lake/build/bin/erdos396ffi` relative to workspace root
    pub fn find_binary() -> Option<PathBuf> {
        if let Ok(path) = std::env::var("LEAN_ORACLE_PATH") {
            let p = PathBuf::from(path);
            if p.exists() {
                return Some(p);
            }
        }
        // Walk up from the crate to the workspace root
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        for _ in 0..5 {
            let candidate = dir.join("formal/.lake/build/bin/erdos396ffi");
            if candidate.exists() {
                return Some(candidate);
            }
            if !dir.pop() {
                break;
            }
        }
        None
    }

    /// Start the Lean oracle subprocess.
    pub fn new() -> Option<Self> {
        let binary = Self::find_binary()?;
        let mut child = Command::new(&binary)
            .args(["--batch", "5000000"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .ok()?;

        let stdin = child.stdin.take()?;
        let stdout = child.stdout.take()?;
        let stderr = child.stderr.take()?;

        // Read the READY line from stderr
        let mut stderr_reader = BufReader::new(stderr);
        let mut ready_line = String::new();
        stderr_reader.read_line(&mut ready_line).ok()?;
        if !ready_line.contains("READY") {
            eprintln!("Lean oracle failed to start: {}", ready_line.trim());
            return None;
        }
        eprintln!("Lean oracle: {}", ready_line.trim());

        Some(Self {
            child,
            stdin,
            reader: BufReader::new(stdout),
        })
    }

    fn query(&mut self, cmd: &str) -> Option<String> {
        writeln!(self.stdin, "{}", cmd).ok()?;
        self.stdin.flush().ok()?;
        let mut line = String::new();
        self.reader.read_line(&mut line).ok()?;
        let trimmed = line.trim().to_string();
        if trimmed == "ERR" {
            None
        } else {
            Some(trimmed)
        }
    }

    fn query_u64(&mut self, cmd: &str) -> u64 {
        self.query(cmd)
            .and_then(|s| s.parse().ok())
            .expect("Lean oracle returned invalid response")
    }

    fn query_bool(&mut self, cmd: &str) -> bool {
        self.query_u64(cmd) == 1
    }
}

impl Drop for LeanOracle {
    fn drop(&mut self) {
        let _ = writeln!(self.stdin, "Q");
        let _ = self.stdin.flush();
        let _ = self.child.wait();
    }
}

impl Oracle for LeanOracle {
    fn name(&self) -> &str {
        "lean4"
    }

    fn vp_factorial(&mut self, n: u64, p: u64) -> u64 {
        self.query_u64(&format!("VF {} {}", n, p))
    }

    fn vp_central_binom(&mut self, n: u64, p: u64) -> u64 {
        self.query_u64(&format!("VC {} {}", n, p))
    }

    fn is_governor(&mut self, n: u64) -> bool {
        self.query_bool(&format!("G {}", n))
    }

    fn check_witness(&mut self, k: u32, n: u64) -> bool {
        self.query_bool(&format!("W {} {}", k, n))
    }
}
