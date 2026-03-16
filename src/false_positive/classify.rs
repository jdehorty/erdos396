use crate::build_info::BuildInfo;
use crate::false_positive::{
    ClassificationConfig, ClassificationStatusSource, RunWindowClassification, VerificationAudit,
};
use crate::run_collection::parquet::{
    classification_writer, latest_classifications, read_run_windows, ClassificationFilter,
    RunWindowFilter,
};
use crate::run_collection::SCHEMA_VERSION;
use crate::verify::WitnessVerifier;
use anyhow::Result;
use chrono::Utc;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use serde::Serialize;

#[derive(Debug, Clone, Default, Serialize)]
pub struct ClassificationStats {
    pub selected_windows: usize,
    pub witnesses: usize,
    pub false_positives: usize,
    pub errors: usize,
}

pub fn verify_window(
    verifier: &WitnessVerifier,
    window_length: usize,
    window_n: u64,
) -> Result<VerificationAudit> {
    let result = verifier.verify(window_length as u32 - 1, window_n)?;
    let (demand, supply) = match result.failing_prime {
        Some(p) => (
            Some(*result.demand.get(&p).unwrap_or(&0)),
            Some(*result.supply.get(&p).unwrap_or(&0)),
        ),
        None => (None, None),
    };

    Ok(VerificationAudit {
        is_witness: result.is_valid,
        failing_prime: result.failing_prime,
        demand,
        supply,
    })
}

fn verify_governor_window(
    verifier: &WitnessVerifier,
    window_length: usize,
    window_n: u64,
) -> Result<VerificationAudit> {
    let result = verifier.verify_fast(window_length as u32 - 1, window_n)?;
    let (demand, supply) = match result.failing_prime {
        Some(p) => (
            Some(*result.demand.get(&p).unwrap_or(&0)),
            Some(*result.supply.get(&p).unwrap_or(&0)),
        ),
        None => (None, None),
    };

    Ok(VerificationAudit {
        is_witness: result.is_valid,
        failing_prime: result.failing_prime,
        demand,
        supply,
    })
}

pub fn classify_run_windows(config: &ClassificationConfig) -> Result<ClassificationStats> {
    let latest = latest_classifications(
        &config.output_root,
        &ClassificationFilter {
            window_length: config.window_length,
        },
    )?;
    let mut windows = read_run_windows(
        &config.output_root,
        &RunWindowFilter {
            window_length: config.window_length,
            start_n: config.start_n,
            end_n: config.end_n,
            limit: None,
        },
    )?;

    windows.retain(|window| match config.status_source {
        ClassificationStatusSource::Recorded => window.recorded_status.is_some(),
        ClassificationStatusSource::Unclassified => !latest.contains_key(&window.window_id),
        ClassificationStatusSource::All => true,
    });

    if let Some(limit) = config.limit {
        windows.truncate(limit);
    }

    let mut stats = ClassificationStats {
        selected_windows: windows.len(),
        ..ClassificationStats::default()
    };
    if windows.is_empty() {
        return Ok(stats);
    }

    let max_n = windows
        .iter()
        .map(|window| window.window_n)
        .max()
        .unwrap_or(1000);
    let verifier = WitnessVerifier::new(max_n.max(1000));
    let verified_at_utc = Utc::now().to_rfc3339();
    let invocation_id = verified_at_utc
        .replace(':', "")
        .replace('-', "")
        .replace('.', "")
        .replace('+', "plus");
    let build_git_hash = BuildInfo::gather()
        .git_hash
        .unwrap_or_else(|| "unknown".to_string());

    let pool = ThreadPoolBuilder::new()
        .num_threads(config.workers.max(1))
        .build()?;
    let rows = pool.install(|| {
        windows
            .par_iter()
            .map(|window| {
                let audit = verify_governor_window(
                    &verifier,
                    window.window_length as usize,
                    window.window_n,
                );
                let (audit_status, failing_prime, demand, supply) = match audit {
                    Ok(audit) if audit.is_witness => ("witness".to_string(), None, None, None),
                    Ok(audit) => (
                        "false_positive".to_string(),
                        audit.failing_prime,
                        audit.demand,
                        audit.supply,
                    ),
                    Err(_) => ("error".to_string(), None, None, None),
                };

                RunWindowClassification {
                    schema_version: SCHEMA_VERSION.to_string(),
                    classification_id: format!("{}|{}", window.window_id, verified_at_utc),
                    window_id: window.window_id.clone(),
                    window_length: window.window_length,
                    window_k: window.window_k,
                    window_start: window.window_start,
                    window_end: window.window_end,
                    window_n: window.window_n,
                    audit_status,
                    failing_prime,
                    demand,
                    supply,
                    verifier_mode: "fast_governor_barrier".to_string(),
                    verified_at_utc: verified_at_utc.clone(),
                    build_git_hash: build_git_hash.clone(),
                }
            })
            .collect::<Vec<_>>()
    });

    let mut writer = classification_writer(&config.output_root, &invocation_id);
    for row in rows {
        match row.audit_status.as_str() {
            "witness" => stats.witnesses += 1,
            "false_positive" => stats.false_positives += 1,
            _ => stats.errors += 1,
        }
        writer.push(row)?;
    }
    writer.finish()?;

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::checkpoint::{Checkpoint, RunInfo};
    use crate::run_collection::parquet::{read_run_windows, RunWindowFilter};
    use crate::run_collection::{build_run_corpus, BuildRunCorpusConfig};
    use tempfile::tempdir;

    #[test]
    fn persists_false_positive_classification() {
        let fixture = tempdir().unwrap();
        let archive_root = fixture.path().join("archives/server_a_k12_10T-15T");
        std::fs::create_dir_all(&archive_root).unwrap();

        let false_positive_n = 5_431_258_726_986u64;
        let run_length = 13usize;
        let run_start = false_positive_n + 1 - run_length as u64;

        let mut checkpoint = Checkpoint::new_worker(12, run_start, false_positive_n + 1, 0);
        checkpoint.version = 2;
        checkpoint.significant_runs.push(RunInfo {
            start: run_start,
            length: run_length,
            is_witness: None,
            failing_prime: None,
        });
        checkpoint.candidates.push(false_positive_n);
        checkpoint.false_positives.push(false_positive_n);
        checkpoint.record_false_positive_detail(false_positive_n, run_length, 3, 8, 7);
        checkpoint
            .save(archive_root.join("checkpoint_k12_w00.json"))
            .unwrap();

        let output = tempdir().unwrap();
        build_run_corpus(&BuildRunCorpusConfig {
            source_root: fixture.path().to_path_buf(),
            output_root: output.path().to_path_buf(),
            min_length: 6,
            max_length: 14,
            include_overlaps: false,
        })
        .unwrap();

        let stats = classify_run_windows(&ClassificationConfig {
            output_root: output.path().to_path_buf(),
            window_length: Some(13),
            start_n: Some(false_positive_n),
            end_n: Some(false_positive_n),
            status_source: ClassificationStatusSource::All,
            workers: 1,
            limit: None,
        })
        .unwrap();

        assert_eq!(stats.false_positives, 1);

        let windows = read_run_windows(
            output.path(),
            &RunWindowFilter {
                window_length: Some(13),
                start_n: Some(false_positive_n),
                end_n: Some(false_positive_n),
                limit: None,
            },
        )
        .unwrap();
        assert_eq!(windows.len(), 1);
        assert_eq!(
            windows[0].recorded_status.as_deref(),
            Some("false_positive")
        );

        let latest = latest_classifications(
            output.path(),
            &ClassificationFilter {
                window_length: Some(13),
            },
        )
        .unwrap();
        let classification = latest.get(&windows[0].window_id).unwrap();
        assert_eq!(classification.audit_status, "false_positive");
        assert_eq!(classification.failing_prime, Some(3));
    }
}
