use crate::run_collection::{RunRecord, RunWindowRecord, RANGE_BUCKET_SIZE, SCHEMA_VERSION};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordedEvent {
    pub status: String,
    pub failing_prime: Option<u64>,
    pub demand: Option<u64>,
    pub supply: Option<u64>,
}

pub fn expand_run_windows<F>(
    run: &RunRecord,
    min_length: usize,
    max_length: usize,
    lookup_event: F,
) -> Vec<RunWindowRecord>
where
    F: Fn(u64, u32, u64) -> Option<RecordedEvent>,
{
    let mut out = Vec::new();
    let run_length = run.run_length as usize;
    let upper = max_length.min(run_length);

    if min_length > upper {
        return out;
    }

    for window_length in min_length..=upper {
        let span_count = run_length - window_length;
        for offset in 0..=span_count {
            let window_start = run.run_start + offset as u64;
            let window_end = window_start + window_length as u64 - 1;
            let event = lookup_event(window_start, window_length as u32, window_end);

            out.push(RunWindowRecord {
                schema_version: SCHEMA_VERSION.to_string(),
                window_id: format!("{}|{}|{}", run.run_id, window_length, offset),
                parent_run_id: run.run_id.clone(),
                checkpoint_id: run.checkpoint_id.clone(),
                source_label: run.source_label.clone(),
                archive_role: run.archive_role.clone(),
                server: run.server.clone(),
                source_path: run.source_path.clone(),
                campaign_target_k: run.campaign_target_k,
                worker_id: run.worker_id,
                max_run_start: run.run_start,
                max_run_end: run.run_end,
                max_run_length: run.run_length,
                window_length: window_length as u32,
                window_k: window_length as u32 - 1,
                window_offset: offset as u32,
                window_start,
                window_end,
                window_n: window_end,
                range_bucket_t: window_end / RANGE_BUCKET_SIZE,
                recorded_status: event.as_ref().map(|e| e.status.clone()),
                recorded_failing_prime: event.as_ref().and_then(|e| e.failing_prime),
                recorded_demand: event.as_ref().and_then(|e| e.demand),
                recorded_supply: event.as_ref().and_then(|e| e.supply),
                is_unique_coverage: run.is_unique_coverage,
            });
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn sample_run(length: u32) -> RunRecord {
        RunRecord {
            schema_version: SCHEMA_VERSION.to_string(),
            run_id: "run-1".to_string(),
            checkpoint_id: "checkpoint".to_string(),
            source_label: "fixture".to_string(),
            archive_role: "canonical".to_string(),
            source_kind: "checkpoint_v2".to_string(),
            server: "fixture".to_string(),
            source_path: "fixture/checkpoint.json".to_string(),
            source_record_ordinal: 0,
            campaign_target_k: 13,
            worker_id: 0,
            run_start: 1000,
            run_end: 1000 + length as u64 - 1,
            run_length: length,
            range_bucket_t: 0,
            is_unique_coverage: true,
        }
    }

    #[test]
    fn expands_length_fourteen_counts() {
        let run = sample_run(14);
        let windows = expand_run_windows(&run, 6, 14, |_, _, _| None);
        let mut counts = HashMap::new();
        for window in windows {
            *counts.entry(window.window_length).or_insert(0usize) += 1;
        }

        assert_eq!(counts.get(&6), Some(&9));
        assert_eq!(counts.get(&7), Some(&8));
        assert_eq!(counts.get(&8), Some(&7));
        assert_eq!(counts.get(&9), Some(&6));
        assert_eq!(counts.get(&10), Some(&5));
        assert_eq!(counts.get(&11), Some(&4));
        assert_eq!(counts.get(&12), Some(&3));
        assert_eq!(counts.get(&13), Some(&2));
        assert_eq!(counts.get(&14), Some(&1));
    }
}
