use serde::{Deserialize, Serialize};

/// Build metadata embedded at compile time (best effort).
///
/// This is intended for reproducibility and review: reports produced by long-running
/// computations can record exactly which source version and toolchain built them.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuildInfo {
    pub crate_version: String,
    pub git_hash: Option<String>,
    pub git_describe: Option<String>,
    pub git_dirty: Option<bool>,
    pub rustc_version: Option<String>,
    pub profile: String,
}

impl BuildInfo {
    /// Gather build metadata from compile-time environment variables.
    ///
    /// Values may be `None` if the binary was built outside of a Git checkout or the
    /// build environment lacked the relevant tools.
    pub fn gather() -> Self {
        Self {
            crate_version: env!("CARGO_PKG_VERSION").to_string(),
            git_hash: option_env!("ERDOS396_GIT_HASH").map(|s| s.to_string()),
            git_describe: option_env!("ERDOS396_GIT_DESCRIBE").map(|s| s.to_string()),
            git_dirty: option_env!("ERDOS396_GIT_DIRTY").map(|s| s == "1"),
            rustc_version: option_env!("ERDOS396_RUSTC_VERSION").map(|s| s.to_string()),
            profile: if cfg!(debug_assertions) {
                "debug".to_string()
            } else {
                "release".to_string()
            },
        }
    }
}
