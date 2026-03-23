use anyhow::Result;
use clap::Parser;
use erdos396::run_collection::{build_run_corpus, BuildRunCorpusConfig, DEFAULT_OUTPUT_ROOT};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "build_run_corpus")]
#[command(about = "Build Parquet datasets for maximal runs and embedded run windows")]
struct Cli {
    #[arg(long)]
    source_root: PathBuf,

    #[arg(long, default_value = DEFAULT_OUTPUT_ROOT)]
    output_root: PathBuf,

    #[arg(long, default_value_t = 6)]
    min_length: usize,

    #[arg(long, default_value_t = 14)]
    max_length: usize,

    #[arg(long, default_value_t = false)]
    include_overlaps: bool,

    #[arg(long = "extra-v3-dir")]
    extra_v3_dirs: Vec<PathBuf>,

    /// Append to existing corpus instead of replacing it
    #[arg(long, default_value_t = false)]
    append: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let stats = build_run_corpus(&BuildRunCorpusConfig {
        source_root: cli.source_root,
        output_root: cli.output_root,
        min_length: cli.min_length,
        max_length: cli.max_length,
        include_overlaps: cli.include_overlaps,
        extra_v3_dirs: cli.extra_v3_dirs,
        append: cli.append,
    })?;

    println!("{}", serde_json::to_string_pretty(&stats)?);
    Ok(())
}
