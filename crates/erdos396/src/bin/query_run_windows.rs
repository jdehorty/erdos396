use anyhow::{anyhow, Result};
use clap::Parser;
use erdos396::false_positive::RunWindowClassification;
use erdos396::run_collection::parquet::{
    latest_classifications, read_run_windows, ClassificationFilter, QuerySummary, RunWindowFilter,
};
use erdos396::run_collection::{RunWindowRecord, DEFAULT_OUTPUT_ROOT};
use serde::Serialize;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "query_run_windows")]
#[command(about = "Query Parquet run-window corpus as JSON")]
struct Cli {
    #[arg(long, default_value = DEFAULT_OUTPUT_ROOT)]
    output_root: PathBuf,

    #[arg(long)]
    window_length: u32,

    #[arg(long, default_value = "any")]
    status: String,

    #[arg(long)]
    start_n: Option<u64>,

    #[arg(long)]
    end_n: Option<u64>,

    #[arg(long)]
    limit: Option<usize>,

    #[arg(long, default_value_t = false)]
    summary: bool,
}

#[derive(Debug, Clone, Copy)]
enum QueryStatus {
    Any,
    RecordedFalsePositive,
    RecordedWitness,
    AuditFalsePositive,
    AuditWitness,
    CoalescedFalsePositive,
    CoalescedWitness,
}

#[derive(Debug, Serialize)]
struct QueryRow {
    #[serde(flatten)]
    window: RunWindowRecord,
    audit_status: Option<String>,
    audit_failing_prime: Option<u64>,
    audit_demand: Option<u64>,
    audit_supply: Option<u64>,
    verified_at_utc: Option<String>,
    coalesced_status: Option<String>,
}

impl QueryStatus {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "any" => Ok(Self::Any),
            "recorded_false_positive" => Ok(Self::RecordedFalsePositive),
            "recorded_witness" => Ok(Self::RecordedWitness),
            "audit_false_positive" => Ok(Self::AuditFalsePositive),
            "audit_witness" => Ok(Self::AuditWitness),
            "coalesced_false_positive" => Ok(Self::CoalescedFalsePositive),
            "coalesced_witness" => Ok(Self::CoalescedWitness),
            _ => Err(anyhow!(
                "invalid status {value}; expected any|recorded_false_positive|recorded_witness|audit_false_positive|audit_witness|coalesced_false_positive|coalesced_witness"
            )),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let status = QueryStatus::parse(&cli.status)?;
    let latest = latest_classifications(
        &cli.output_root,
        &ClassificationFilter {
            window_length: Some(cli.window_length),
        },
    )?;
    let mut rows = Vec::new();
    for window in read_run_windows(
        &cli.output_root,
        &RunWindowFilter {
            window_length: Some(cli.window_length),
            start_n: cli.start_n,
            end_n: cli.end_n,
            limit: None,
        },
    )? {
        let latest_audit = latest.get(&window.window_id);
        if !matches_status(&window, latest_audit, status) {
            continue;
        }
        rows.push(QueryRow {
            coalesced_status: coalesced_status(&window, latest_audit),
            audit_status: latest_audit.map(|audit| audit.audit_status.clone()),
            audit_failing_prime: latest_audit.and_then(|audit| audit.failing_prime),
            audit_demand: latest_audit.and_then(|audit| audit.demand),
            audit_supply: latest_audit.and_then(|audit| audit.supply),
            verified_at_utc: latest_audit.map(|audit| audit.verified_at_utc.clone()),
            window,
        });
        if let Some(limit) = cli.limit {
            if rows.len() >= limit {
                break;
            }
        }
    }

    if cli.summary {
        println!(
            "{}",
            serde_json::to_string_pretty(&QuerySummary { count: rows.len() })?
        );
    } else {
        println!("{}", serde_json::to_string_pretty(&rows)?);
    }
    Ok(())
}

fn matches_status(
    window: &RunWindowRecord,
    latest_audit: Option<&RunWindowClassification>,
    status: QueryStatus,
) -> bool {
    match status {
        QueryStatus::Any => true,
        QueryStatus::RecordedFalsePositive => {
            window.recorded_status.as_deref() == Some("false_positive")
        }
        QueryStatus::RecordedWitness => window.recorded_status.as_deref() == Some("witness"),
        QueryStatus::AuditFalsePositive => {
            latest_audit.map(|audit| audit.audit_status.as_str()) == Some("false_positive")
        }
        QueryStatus::AuditWitness => {
            latest_audit.map(|audit| audit.audit_status.as_str()) == Some("witness")
        }
        QueryStatus::CoalescedFalsePositive => {
            coalesced_status(window, latest_audit).as_deref() == Some("false_positive")
        }
        QueryStatus::CoalescedWitness => {
            coalesced_status(window, latest_audit).as_deref() == Some("witness")
        }
    }
}

fn coalesced_status(
    window: &RunWindowRecord,
    latest_audit: Option<&RunWindowClassification>,
) -> Option<String> {
    latest_audit
        .map(|audit| audit.audit_status.clone())
        .or_else(|| window.recorded_status.clone())
}
