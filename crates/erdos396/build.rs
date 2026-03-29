use std::path::Path;
use std::process::Command;

fn git_output(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?;
    let s = s.trim().to_string();
    if s.is_empty() {
        None
    } else {
        Some(s)
    }
}

fn git_dirty() -> Option<bool> {
    let out = Command::new("git")
        .args(["status", "--porcelain"])
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    Some(!out.stdout.is_empty())
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    if Path::new(".git/HEAD").exists() {
        println!("cargo:rerun-if-changed=.git/HEAD");
    }
    if Path::new(".git/index").exists() {
        println!("cargo:rerun-if-changed=.git/index");
    }

    if let Some(hash) = git_output(&["rev-parse", "HEAD"]) {
        println!("cargo:rustc-env=ERDOS396_GIT_HASH={hash}");
    }
    if let Some(describe) = git_output(&["describe", "--always", "--dirty"]) {
        println!("cargo:rustc-env=ERDOS396_GIT_DESCRIBE={describe}");
    }
    if let Some(is_dirty) = git_dirty() {
        let dirty = if is_dirty { "1" } else { "0" };
        println!("cargo:rustc-env=ERDOS396_GIT_DIRTY={dirty}");
    }

    let rustc = Command::new("rustc").arg("--version").output();
    if let Ok(out) = rustc {
        if out.status.success() {
            if let Ok(s) = String::from_utf8(out.stdout) {
                let s = s.trim().to_string();
                if !s.is_empty() {
                    println!("cargo:rustc-env=ERDOS396_RUSTC_VERSION={s}");
                }
            }
        }
    }
}
