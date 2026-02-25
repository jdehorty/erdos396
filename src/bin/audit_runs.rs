//! Audit tool for Governor run logs (`runs_k{K}_w{NN}.jsonl`).
//!
//! This is intended for reviewer-facing reproducibility:
//! given the run logs emitted by `erdos396`, verify each recorded run
//! as a candidate witness for `k = length - 1` at `n = start + length - 1`.
//!
//! This does *not* prove completeness of witness discovery in a range; that role
//! is served by the `validate` binary (Small Prime Barrier screen + full verification).
//! See `docs/validation.md` and `docs/trust.md`.

use anyhow::Context;
use clap::Parser;
use erdos396::verify::WitnessVerifier;
use erdos396::BuildInfo;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const FNV1A_64_OFFSET: u64 = 14695981039346656037;
const FNV1A_64_PRIME: u64 = 1099511628211;

#[inline]
fn fnv1a64_update(mut hash: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(FNV1A_64_PRIME);
    }
    hash
}

#[derive(Parser, Debug)]
#[command(name = "audit_runs")]
#[command(about = "Audit Governor run logs by full witness verification")]
#[command(long_about = r#"
Reads `runs_k{K}_w{NN}.jsonl` files produced by `erdos396` and verifies each run.

For each JSONL record {"start":S,"length":L}, we verify the witness condition at:
  k = L - 1,  n = S + L - 1

This is useful for producing a false-positive catalog or cross-checking run logs.
"#)]
struct Cli {
    /// The k used in the `erdos396` run (selects files `runs_k{K}_w*.jsonl`).
    #[arg(short = 'k', long)]
    target_k: u32,

    /// Directory containing the run logs (and usually checkpoints).
    #[arg(short = 'd', long, default_value = "checkpoints")]
    input_dir: PathBuf,

    /// Minimum run length to audit (default matches `significant_run_threshold`).
    #[arg(long, default_value = "6")]
    min_length: usize,

    /// Maximum number of runs to process (0 means no limit).
    #[arg(long, default_value = "0")]
    limit: u64,

    /// Maximum n to size the prime sieve (optional).
    ///
    /// If omitted, the tool will try to infer it from checkpoints in `--input-dir`,
    /// falling back to a quick scan of the run logs.
    #[arg(long)]
    max_n: Option<u64>,

    /// Write false positives as JSONL (one record per line).
    #[arg(long)]
    false_positive_jsonl: Option<PathBuf>,

    /// Write verified witnesses as JSONL (one record per line).
    #[arg(long)]
    witness_jsonl: Option<PathBuf>,

    /// Output summary report JSON path.
    #[arg(short = 'o', long, default_value = "audit_runs_report.json")]
    output: PathBuf,

    /// Verbose logging.
    #[arg(short = 'v', long)]
    verbose: bool,
}

#[derive(Debug, Deserialize)]
struct RunLine {
    start: u64,
    length: usize,
}

#[derive(Debug, Clone)]
struct RunCandidate {
    start: u64,
    length: usize,
}

#[derive(Debug, Serialize)]
struct RunOutcome {
    start: u64,
    length: usize,
    k: u32,
    n: u64,
    is_witness: bool,
    failing_prime: Option<u64>,
    demand: Option<u64>,
    supply: Option<u64>,
}

#[derive(Debug, Default, Serialize)]
struct LengthStats {
    runs: u64,
    witnesses: u64,
    false_positives: u64,
}

#[derive(Debug, Serialize)]
struct RunLogFileStats {
    file: String,
    fnv1a64: String,
    total_lines: u64,
    parsed_runs: u64,
    kept_runs: u64,
    skipped_out_of_order: u64,
    first_start: Option<u64>,
    last_start: Option<u64>,
}

#[derive(Debug, Serialize)]
struct AuditReport {
    version: &'static str,
    build: BuildInfo,
    args: Vec<String>,
    target_k: u32,
    input_dir: String,
    min_length: usize,
    limit: u64,
    max_n: u64,
    files: Vec<String>,
    file_stats: Vec<RunLogFileStats>,
    total_runs: u64,
    witnesses: u64,
    false_positives: u64,
    by_length: HashMap<usize, LengthStats>,
    failing_prime_counts: HashMap<u64, u64>,
    generated_at: String,
}

fn list_run_logs(dir: &Path, target_k: u32) -> anyhow::Result<Vec<PathBuf>> {
    let prefix = format!("runs_k{target_k}_w");
    let mut out = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let name = match path.file_name().and_then(|s| s.to_str()) {
            Some(s) => s,
            None => continue,
        };
        if name.starts_with(&prefix) && name.ends_with(".jsonl") {
            out.push(path);
        }
    }
    out.sort();
    Ok(out)
}

fn infer_max_n_from_checkpoints(dir: &Path, target_k: u32) -> anyhow::Result<Option<u64>> {
    let prefix = format!("checkpoint_k{target_k}_w");
    let mut max_n = 0u64;
    let mut found = false;
    for entry in std::fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let name = match path.file_name().and_then(|s| s.to_str()) {
            Some(s) => s,
            None => continue,
        };
        if !name.starts_with(&prefix) || !name.ends_with(".json") {
            continue;
        }
        let cp = erdos396::checkpoint::Checkpoint::load(&path)
            .with_context(|| format!("load checkpoint {}", path.display()))?;
        if cp.target_k != target_k {
            continue;
        }
        found = true;
        max_n = max_n.max(cp.end);
    }
    Ok(if found { Some(max_n) } else { None })
}

fn infer_max_n_from_run_logs(
    files: &[PathBuf],
    min_length: usize,
    limit: u64,
) -> anyhow::Result<u64> {
    let mut max_n = 0u64;
    let mut seen = 0u64;

    for path in files {
        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let reader = BufReader::new(file);
        for (lineno, line) in reader.lines().enumerate() {
            let line =
                line.with_context(|| format!("read line {}:{}", path.display(), lineno + 1))?;
            let line = line.trim_end_matches('\r');
            if line.trim().is_empty() {
                continue;
            }
            let rec: RunLine = serde_json::from_str(line).with_context(|| {
                format!("parse JSON {}:{}: {}", path.display(), lineno + 1, line)
            })?;
            if rec.length < min_length || rec.length == 0 {
                continue;
            }
            let n = rec
                .start
                .checked_add(rec.length as u64)
                .and_then(|x| x.checked_sub(1))
                .context("overflow computing run end n")?;
            max_n = max_n.max(n);
            seen += 1;
            if limit > 0 && seen >= limit {
                return Ok(max_n);
            }
        }
    }

    Ok(max_n)
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    if cli.target_k == 0 {
        anyhow::bail!("--target-k must be >= 1");
    }

    let files = list_run_logs(&cli.input_dir, cli.target_k)?;
    if files.is_empty() {
        anyhow::bail!(
            "No run logs found in {} for target_k={}",
            cli.input_dir.display(),
            cli.target_k
        );
    }
    log::info!(
        "Found {} run log files for target_k={}",
        files.len(),
        cli.target_k
    );

    let max_n = if let Some(n) = cli.max_n {
        n
    } else if let Some(n) = infer_max_n_from_checkpoints(&cli.input_dir, cli.target_k)? {
        log::info!("Inferred max_n={} from checkpoints.", n);
        n
    } else {
        log::info!(
            "Inferring max_n from run logs (min_length={})...",
            cli.min_length
        );
        infer_max_n_from_run_logs(&files, cli.min_length, cli.limit)?
    };
    log::info!(
        "Using max_n={} for prime sieve sizing (override with --max-n)",
        max_n
    );

    let verifier = WitnessVerifier::new(max_n.max(1000));

    let mut fp_writer = match &cli.false_positive_jsonl {
        Some(path) => Some(BufWriter::new(
            File::create(path).with_context(|| format!("create {}", path.display()))?,
        )),
        None => None,
    };
    let mut w_writer = match &cli.witness_jsonl {
        Some(path) => Some(BufWriter::new(
            File::create(path).with_context(|| format!("create {}", path.display()))?,
        )),
        None => None,
    };

    let mut by_length: HashMap<usize, LengthStats> = HashMap::new();
    let mut failing_prime_counts: HashMap<u64, u64> = HashMap::new();
    let mut file_stats: Vec<RunLogFileStats> = Vec::new();

    let mut total_runs = 0u64;
    let mut witnesses = 0u64;
    let mut false_positives = 0u64;

    // Chunked parallel verification.
    const CHUNK: usize = 50_000;
    let mut chunk: Vec<RunCandidate> = Vec::with_capacity(CHUNK);
    let mut seen = 0u64;
    for path in &files {
        let mut last_start_kept: Option<u64> = None;
        let mut h = FNV1A_64_OFFSET;
        let mut total_lines = 0u64;
        let mut parsed_runs = 0u64;
        let mut kept_runs = 0u64;
        let mut skipped_out_of_order = 0u64;
        let mut first_start: Option<u64> = None;
        let mut last_start: Option<u64> = None;

        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let reader = BufReader::new(file);
        for (lineno, line) in reader.lines().enumerate() {
            let line =
                line.with_context(|| format!("read line {}:{}", path.display(), lineno + 1))?;
            let line = line.trim_end_matches('\r');
            if line.trim().is_empty() {
                continue;
            }
            total_lines += 1;
            h = fnv1a64_update(h, line.as_bytes());
            h = fnv1a64_update(h, b"\n");
            let rec: RunLine = serde_json::from_str(line).with_context(|| {
                format!("parse JSON {}:{}: {}", path.display(), lineno + 1, line)
            })?;
            parsed_runs += 1;
            if rec.length < cli.min_length || rec.length == 0 {
                continue;
            }
            if let Some(prev) = last_start_kept {
                if rec.start <= prev {
                    skipped_out_of_order += 1;
                    continue;
                }
            }
            last_start_kept = Some(rec.start);
            kept_runs += 1;
            if first_start.is_none() {
                first_start = Some(rec.start);
            }
            last_start = Some(rec.start);
            chunk.push(RunCandidate {
                start: rec.start,
                length: rec.length,
            });
            seen += 1;
            if cli.limit > 0 && seen >= cli.limit {
                break;
            }
            if chunk.len() >= CHUNK {
                let outcomes: Vec<RunOutcome> = chunk
                    .par_iter()
                    .map(|c| -> anyhow::Result<RunOutcome> {
                        let length = c.length;
                        let k = (length - 1) as u32;
                        let n = c.start + length as u64 - 1;
                        let r = verifier.verify(k, n)?;
                        let mut out = RunOutcome {
                            start: c.start,
                            length,
                            k,
                            n,
                            is_witness: r.is_valid,
                            failing_prime: r.failing_prime,
                            demand: None,
                            supply: None,
                        };
                        if let Some(p) = r.failing_prime {
                            out.demand = Some(*r.demand.get(&p).unwrap_or(&0));
                            out.supply = Some(*r.supply.get(&p).unwrap_or(&0));
                        }
                        Ok(out)
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                chunk.clear();

                for out in outcomes {
                    total_runs += 1;
                    let stats = by_length.entry(out.length).or_default();
                    stats.runs += 1;

                    if out.is_witness {
                        witnesses += 1;
                        stats.witnesses += 1;
                        if let Some(w) = w_writer.as_mut() {
                            serde_json::to_writer(&mut *w, &out)?;
                            w.write_all(b"\n")?;
                        }
                    } else {
                        false_positives += 1;
                        stats.false_positives += 1;
                        if let Some(p) = out.failing_prime {
                            *failing_prime_counts.entry(p).or_insert(0) += 1;
                        }
                        if let Some(w) = fp_writer.as_mut() {
                            serde_json::to_writer(&mut *w, &out)?;
                            w.write_all(b"\n")?;
                        }
                    }
                }
            }
        }
        file_stats.push(RunLogFileStats {
            file: path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string(),
            fnv1a64: format!("{:016x}", h),
            total_lines,
            parsed_runs,
            kept_runs,
            skipped_out_of_order,
            first_start,
            last_start,
        });
        if cli.limit > 0 && seen >= cli.limit {
            break;
        }
    }

    if !chunk.is_empty() {
        let outcomes: Vec<RunOutcome> = chunk
            .par_iter()
            .map(|c| -> anyhow::Result<RunOutcome> {
                let length = c.length;
                let k = (length - 1) as u32;
                let n = c.start + length as u64 - 1;
                let r = verifier.verify(k, n)?;
                let mut out = RunOutcome {
                    start: c.start,
                    length,
                    k,
                    n,
                    is_witness: r.is_valid,
                    failing_prime: r.failing_prime,
                    demand: None,
                    supply: None,
                };
                if let Some(p) = r.failing_prime {
                    out.demand = Some(*r.demand.get(&p).unwrap_or(&0));
                    out.supply = Some(*r.supply.get(&p).unwrap_or(&0));
                }
                Ok(out)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        for out in outcomes {
            total_runs += 1;
            let stats = by_length.entry(out.length).or_default();
            stats.runs += 1;

            if out.is_witness {
                witnesses += 1;
                stats.witnesses += 1;
                if let Some(w) = w_writer.as_mut() {
                    serde_json::to_writer(&mut *w, &out)?;
                    w.write_all(b"\n")?;
                }
            } else {
                false_positives += 1;
                stats.false_positives += 1;
                if let Some(p) = out.failing_prime {
                    *failing_prime_counts.entry(p).or_insert(0) += 1;
                }
                if let Some(w) = fp_writer.as_mut() {
                    serde_json::to_writer(&mut *w, &out)?;
                    w.write_all(b"\n")?;
                }
            }
        }
    }

    if let Some(w) = fp_writer.as_mut() {
        w.flush()?;
    }
    if let Some(w) = w_writer.as_mut() {
        w.flush()?;
    }

    let report = AuditReport {
        version: "audit_runs_report_v1",
        build: BuildInfo::gather(),
        args: std::env::args().collect(),
        target_k: cli.target_k,
        input_dir: cli.input_dir.display().to_string(),
        min_length: cli.min_length,
        limit: cli.limit,
        max_n,
        files: files
            .iter()
            .map(|p| {
                p.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string()
            })
            .collect(),
        file_stats,
        total_runs,
        witnesses,
        false_positives,
        by_length,
        failing_prime_counts,
        generated_at: chrono::Utc::now().to_rfc3339(),
    };

    let json = serde_json::to_string_pretty(&report)?;
    let tmp = cli.output.with_extension("tmp");
    std::fs::write(&tmp, json)?;
    std::fs::rename(&tmp, &cli.output)?;

    println!(
        "Audit complete: {} runs, {} witnesses, {} false positives (report: {})",
        total_runs,
        witnesses,
        false_positives,
        cli.output.display()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_list_run_logs_filters_by_k() {
        let dir = tempdir().unwrap();
        std::fs::write(dir.path().join("runs_k1_w00.jsonl"), b"").unwrap();
        std::fs::write(dir.path().join("runs_k2_w00.jsonl"), b"").unwrap();
        std::fs::write(dir.path().join("not_a_run_log.txt"), b"").unwrap();

        let k1 = list_run_logs(dir.path(), 1).unwrap();
        assert_eq!(k1.len(), 1);
        assert!(k1[0]
            .file_name()
            .unwrap()
            .to_string_lossy()
            .contains("runs_k1_w"));
    }

    #[test]
    fn test_infer_max_n_from_run_logs_respects_limit() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("runs_k1_w00.jsonl");
        let content = br#"{"start":100,"length":6}
{"start":1000,"length":7}
"#;
        std::fs::write(&path, content).unwrap();

        let files = vec![path];
        let max1 = infer_max_n_from_run_logs(&files, 6, 1).unwrap();
        assert_eq!(max1, 105); // 100 + 6 - 1

        let max2 = infer_max_n_from_run_logs(&files, 6, 0).unwrap();
        assert_eq!(max2, 1006); // 1000 + 7 - 1
    }

    #[test]
    fn test_infer_max_n_from_checkpoints() {
        let dir = tempdir().unwrap();
        let mut cp0 = erdos396::checkpoint::Checkpoint::new_worker(13, 0, 100, 0);
        cp0.current_pos = 100;
        cp0.checked = 100;
        cp0.sum_checked = 4950;
        cp0.xor_checked = 0;
        cp0.save(dir.path().join("checkpoint_k13_w00.json"))
            .unwrap();

        let mut cp1 = erdos396::checkpoint::Checkpoint::new_worker(13, 100, 250, 1);
        cp1.current_pos = 250;
        cp1.checked = 150;
        cp1.sum_checked = 26175;
        cp1.xor_checked = 0;
        cp1.save(dir.path().join("checkpoint_k13_w01.json"))
            .unwrap();

        let inferred = infer_max_n_from_checkpoints(dir.path(), 13).unwrap();
        assert_eq!(inferred, Some(250));
    }
}
