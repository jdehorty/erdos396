use anyhow::Result;
use clap::Parser;
use erdos396::false_positive::{
    classify_run_windows, ClassificationConfig, ClassificationStatusSource,
};
use erdos396::run_collection::DEFAULT_OUTPUT_ROOT;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "classify_run_windows")]
#[command(about = "Fully verify run windows and persist witness / false-positive classifications")]
struct Cli {
    #[arg(long, default_value = DEFAULT_OUTPUT_ROOT)]
    output_root: PathBuf,

    #[arg(long)]
    window_length: Option<u32>,

    #[arg(long)]
    start_n: Option<u64>,

    #[arg(long)]
    end_n: Option<u64>,

    #[arg(long, default_value = "all")]
    status_source: String,

    #[arg(long, default_value_t = 1)]
    workers: usize,

    #[arg(long)]
    limit: Option<usize>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let stats = classify_run_windows(&ClassificationConfig {
        output_root: cli.output_root,
        window_length: cli.window_length,
        start_n: cli.start_n,
        end_n: cli.end_n,
        status_source: ClassificationStatusSource::parse(&cli.status_source)?,
        workers: cli.workers,
        limit: cli.limit,
    })?;

    println!("{}", serde_json::to_string_pretty(&stats)?);
    Ok(())
}
