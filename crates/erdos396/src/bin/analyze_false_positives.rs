use anyhow::{Context, Result};
use arrow::array::{Array, StringArray, UInt32Array, UInt64Array};
use clap::Parser;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

#[derive(Parser, Debug)]
#[command(name = "analyze_false_positives")]
#[command(about = "Stream false-positive analytics from run_window_classifications parquet")]
struct Cli {
    #[arg(long, default_value = "data/run_corpus/v1")]
    output_root: PathBuf,
}

#[derive(Debug, Default)]
struct MutableLengthStats {
    audited: u64,
    false_positives: u64,
    prime_counts: HashMap<u64, u64>,
    deficit_counts: HashMap<u64, u64>,
    trillion_counts: HashMap<u64, u64>,
    audited_mod9: [u64; 9],
    p3_false_positive_mod9: [u64; 9],
    audited_mod25: [u64; 25],
    p5_false_positive_mod25: [u64; 25],
    examples: Vec<Example>,
}

#[derive(Debug, Clone, Serialize)]
struct Example {
    n: u64,
    failing_prime: Option<u64>,
    deficit: Option<u64>,
}

#[derive(Debug, Clone, Serialize)]
struct CountPair {
    key: u64,
    count: u64,
}

#[derive(Debug, Serialize)]
struct LengthReport {
    window_length: u32,
    audited: u64,
    false_positives: u64,
    false_positive_rate: f64,
    failing_primes: Vec<CountPair>,
    deficits: Vec<CountPair>,
    trillion_buckets: Vec<CountPair>,
    p3_mod9_hotspots: Vec<ResidueRate>,
    p5_mod25_hotspots: Vec<ResidueRate>,
    examples: Vec<Example>,
}

#[derive(Debug, Serialize)]
struct Report {
    total_audited: u64,
    total_false_positives: u64,
    false_positive_rate: f64,
    unique_false_positive_n: u64,
    repeated_false_positive_n: u64,
    repeated_length_patterns: Vec<LengthPatternCount>,
    repeated_examples: Vec<RepeatedExample>,
    by_length: Vec<LengthReport>,
    overall_failing_primes: Vec<CountPair>,
    overall_deficits: Vec<CountPair>,
    overall_trillion_buckets: Vec<CountPair>,
    overall_p3_mod9_hotspots: Vec<ResidueRate>,
    overall_p5_mod25_hotspots: Vec<ResidueRate>,
}

#[derive(Debug, Serialize)]
struct LengthPatternCount {
    lengths: Vec<u32>,
    count: u64,
}

#[derive(Debug, Serialize)]
struct RepeatedExample {
    n: u64,
    lengths: Vec<u32>,
}

#[derive(Debug, Serialize)]
struct ResidueRate {
    modulus: u64,
    residue: u64,
    audited: u64,
    false_positives: u64,
    rate: f64,
    lift: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut total_audited = 0u64;
    let mut total_false_positives = 0u64;
    let mut overall_prime_counts = HashMap::new();
    let mut overall_deficit_counts = HashMap::new();
    let mut overall_trillion_counts = HashMap::new();
    let mut overall_audited_mod9 = [0u64; 9];
    let mut overall_p3_false_positive_mod9 = [0u64; 9];
    let mut overall_audited_mod25 = [0u64; 25];
    let mut overall_p5_false_positive_mod25 = [0u64; 25];
    let mut by_length: BTreeMap<u32, MutableLengthStats> = BTreeMap::new();
    let mut false_positive_patterns: HashMap<u64, u16> = HashMap::new();

    for path in parquet_files_under(&cli.output_root.join("run_window_classifications"))? {
        let file = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .with_context(|| format!("read parquet {}", path.display()))?;
        let reader = builder.build()?;
        for batch in reader {
            let batch = batch?;
            let window_length = u32_array(&batch, "window_length")?;
            let window_n = u64_array(&batch, "window_n")?;
            let audit_status = string_array(&batch, "audit_status")?;
            let failing_prime = u64_array(&batch, "failing_prime")?;
            let demand = u64_array(&batch, "demand")?;
            let supply = u64_array(&batch, "supply")?;

            for row in 0..batch.num_rows() {
                let len = window_length.value(row);
                let stats = by_length.entry(len).or_default();
                stats.audited += 1;
                total_audited += 1;
                let mod9 = (window_n.value(row) % 9) as usize;
                let mod25 = (window_n.value(row) % 25) as usize;
                stats.audited_mod9[mod9] += 1;
                stats.audited_mod25[mod25] += 1;
                overall_audited_mod9[mod9] += 1;
                overall_audited_mod25[mod25] += 1;

                if audit_status.value(row) != "false_positive" {
                    continue;
                }

                let n = window_n.value(row);
                let p = nullable_u64(failing_prime, row);
                let d = nullable_u64(demand, row);
                let s = nullable_u64(supply, row);
                let deficit = d.zip(s).map(|(d, s)| d.saturating_sub(s));
                let trillion_bucket = n / 1_000_000_000_000;

                stats.false_positives += 1;
                total_false_positives += 1;
                let bit = length_bit(len);
                false_positive_patterns
                    .entry(n)
                    .and_modify(|mask| *mask |= bit)
                    .or_insert(bit);
                *stats.trillion_counts.entry(trillion_bucket).or_insert(0) += 1;
                *overall_trillion_counts.entry(trillion_bucket).or_insert(0) += 1;
                if let Some(p) = p {
                    *stats.prime_counts.entry(p).or_insert(0) += 1;
                    *overall_prime_counts.entry(p).or_insert(0) += 1;
                    if p == 3 {
                        stats.p3_false_positive_mod9[mod9] += 1;
                        overall_p3_false_positive_mod9[mod9] += 1;
                    }
                    if p == 5 {
                        stats.p5_false_positive_mod25[mod25] += 1;
                        overall_p5_false_positive_mod25[mod25] += 1;
                    }
                }
                if let Some(deficit) = deficit {
                    *stats.deficit_counts.entry(deficit).or_insert(0) += 1;
                    *overall_deficit_counts.entry(deficit).or_insert(0) += 1;
                }
                if stats.examples.len() < 5 {
                    stats.examples.push(Example {
                        n,
                        failing_prime: p,
                        deficit,
                    });
                }
            }
        }
    }

    let by_length = by_length
        .into_iter()
        .map(|(window_length, stats)| LengthReport {
            window_length,
            audited: stats.audited,
            false_positives: stats.false_positives,
            false_positive_rate: ratio(stats.false_positives, stats.audited),
            failing_primes: top_counts(stats.prime_counts, 8),
            deficits: top_counts(stats.deficit_counts, 8),
            trillion_buckets: top_counts(stats.trillion_counts, 8),
            p3_mod9_hotspots: residue_hotspots(
                9,
                &stats.audited_mod9,
                &stats.p3_false_positive_mod9,
                4,
            ),
            p5_mod25_hotspots: residue_hotspots(
                25,
                &stats.audited_mod25,
                &stats.p5_false_positive_mod25,
                6,
            ),
            examples: stats.examples,
        })
        .collect();

    let repeated_false_positive_n = false_positive_patterns
        .values()
        .filter(|&&mask| mask.count_ones() > 1)
        .count() as u64;
    let mut repeated_pattern_counts: HashMap<u16, u64> = HashMap::new();
    let mut repeated_examples: Vec<RepeatedExample> = false_positive_patterns
        .iter()
        .filter_map(|(&n, &mask)| {
            if mask.count_ones() > 1 {
                *repeated_pattern_counts.entry(mask).or_insert(0) += 1;
                Some(RepeatedExample {
                    n,
                    lengths: mask_to_lengths(mask),
                })
            } else {
                None
            }
        })
        .collect();
    repeated_examples.sort_by_key(|example| example.n);
    repeated_examples.truncate(12);

    let report = Report {
        total_audited,
        total_false_positives,
        false_positive_rate: ratio(total_false_positives, total_audited),
        unique_false_positive_n: false_positive_patterns.len() as u64,
        repeated_false_positive_n,
        repeated_length_patterns: top_patterns(repeated_pattern_counts, 12),
        repeated_examples,
        by_length,
        overall_failing_primes: top_counts(overall_prime_counts, 12),
        overall_deficits: top_counts(overall_deficit_counts, 12),
        overall_trillion_buckets: top_counts(overall_trillion_counts, 12),
        overall_p3_mod9_hotspots: residue_hotspots(
            9,
            &overall_audited_mod9,
            &overall_p3_false_positive_mod9,
            9,
        ),
        overall_p5_mod25_hotspots: residue_hotspots(
            25,
            &overall_audited_mod25,
            &overall_p5_false_positive_mod25,
            10,
        ),
    };

    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}

fn parquet_files_under(root: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in WalkDir::new(root) {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() && path.extension().and_then(|ext| ext.to_str()) == Some("parquet") {
            files.push(path.to_path_buf());
        }
    }
    files.sort();
    Ok(files)
}

fn top_counts(counts: HashMap<u64, u64>, n: usize) -> Vec<CountPair> {
    let mut items: Vec<_> = counts
        .into_iter()
        .map(|(key, count)| CountPair { key, count })
        .collect();
    items.sort_by(|a, b| b.count.cmp(&a.count).then_with(|| a.key.cmp(&b.key)));
    items.truncate(n);
    items
}

fn ratio(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        num as f64 / den as f64
    }
}

fn length_bit(length: u32) -> u16 {
    1u16 << (length - 6)
}

fn mask_to_lengths(mask: u16) -> Vec<u32> {
    (6u32..=14u32)
        .filter(|length| (mask & length_bit(*length)) != 0)
        .collect()
}

fn top_patterns(counts: HashMap<u16, u64>, n: usize) -> Vec<LengthPatternCount> {
    let mut items: Vec<_> = counts
        .into_iter()
        .map(|(mask, count)| LengthPatternCount {
            lengths: mask_to_lengths(mask),
            count,
        })
        .collect();
    items.sort_by(|a, b| {
        b.count
            .cmp(&a.count)
            .then_with(|| a.lengths.cmp(&b.lengths))
    });
    items.truncate(n);
    items
}

fn residue_hotspots<const N: usize>(
    modulus: u64,
    audited: &[u64; N],
    false_positives: &[u64; N],
    limit: usize,
) -> Vec<ResidueRate> {
    let total_audited: u64 = audited.iter().sum();
    let total_false_positives: u64 = false_positives.iter().sum();
    let baseline = ratio(total_false_positives, total_audited);

    let mut rows: Vec<_> = (0..N)
        .map(|residue| {
            let audited_count = audited[residue];
            let false_positive_count = false_positives[residue];
            let rate = ratio(false_positive_count, audited_count);
            let lift = if baseline == 0.0 {
                0.0
            } else {
                rate / baseline
            };
            ResidueRate {
                modulus,
                residue: residue as u64,
                audited: audited_count,
                false_positives: false_positive_count,
                rate,
                lift,
            }
        })
        .collect();

    rows.sort_by(|a, b| {
        b.rate
            .partial_cmp(&a.rate)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.false_positives.cmp(&a.false_positives))
            .then_with(|| a.residue.cmp(&b.residue))
    });
    rows.truncate(limit);
    rows
}

fn string_array<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a StringArray> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<StringArray>()
        .with_context(|| format!("column {name} is not StringArray"))
}

fn u32_array<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a UInt32Array> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .with_context(|| format!("column {name} is not UInt32Array"))
}

fn u64_array<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a UInt64Array> {
    let index = batch.schema().index_of(name)?;
    batch
        .column(index)
        .as_any()
        .downcast_ref::<UInt64Array>()
        .with_context(|| format!("column {name} is not UInt64Array"))
}

fn nullable_u64(array: &UInt64Array, row: usize) -> Option<u64> {
    if array.is_null(row) {
        None
    } else {
        Some(array.value(row))
    }
}
