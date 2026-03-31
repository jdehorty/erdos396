use crate::analysis::{has_barrier_prime_failure, NearMiss};
use crate::fast_sieve::{build_prime_data, fast_governor_sieve, is_governor_cold, strip_only_sieve, PrimeData};
use erdos396::PrimeSieve;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const CHUNK_SIZE: u64 = 1_048_576; // 1M integers per chunk

/// Scan range [start, start+count) for near-miss blocks of width k+1.
///
/// Uses fast_governor_sieve (modular-inverse division, L1 blocks, bucketed primes)
/// and rayon for parallel chunk processing.
///
/// A near-miss is a window [n-k, n] containing exactly one non-governor,
/// where that non-governor fails at a barrier prime (not just a large prime).
pub fn scan_near_misses(
    start: u64,
    count: u64,
    k: u32,
    primes: &[u64],
    sieve: &PrimeSieve,
    progress: bool,
) -> Vec<NearMiss> {
    let window_size = k as u64 + 1;
    if count < window_size {
        return Vec::new();
    }

    // Split range into chunks with k overlap for cross-boundary windows.
    // Each chunk covers [chunk_start, chunk_start + chunk_len) and produces
    // near-misses with block top n in [chunk_start + k, chunk_start + chunk_len - 1].
    // The overlap ensures windows spanning chunk boundaries are seen by one chunk.
    let overlap = k as u64;
    let mut chunks: Vec<(u64, u64)> = Vec::new();

    let mut chunk_start = start;
    while chunk_start < start + count {
        let remaining = (start + count) - chunk_start;
        let chunk_len = remaining.min(CHUNK_SIZE + overlap);
        if chunk_len >= window_size {
            chunks.push((chunk_start, chunk_len));
        }
        // Advance by CHUNK_SIZE (not chunk_len) so chunks overlap by k
        chunk_start += CHUNK_SIZE;
    }

    let prime_data = build_prime_data(primes);

    let t_start = Instant::now();
    let chunks_done = AtomicU64::new(0);
    let total_chunks = chunks.len() as u64;

    let all_results: Vec<Vec<NearMiss>> = chunks
        .par_iter()
        .map(|&(c_start, c_len)| {
            let result = scan_chunk(c_start, c_len, k, &prime_data, sieve, start, count);

            if progress {
                let done = chunks_done.fetch_add(1, Ordering::Relaxed) + 1;
                if done % 100 == 0 || done == total_chunks {
                    let elapsed = t_start.elapsed().as_secs_f64();
                    let scanned = done * CHUNK_SIZE;
                    let rate = scanned as f64 / elapsed / 1e6;
                    let pct = done as f64 / total_chunks as f64 * 100.0;
                    eprintln!(
                        "[{:.1}%] {:.3}B scanned  {:.1}M/s  {:.1}s elapsed",
                        pct,
                        scanned as f64 / 1e9,
                        rate,
                        elapsed,
                    );
                }
            }

            result
        })
        .collect();

    // Flatten and deduplicate (overlapping chunks can find the same near-miss)
    let mut results: Vec<NearMiss> = all_results.into_iter().flatten().collect();
    results.sort_by_key(|nm| nm.n);
    results.dedup_by_key(|nm| nm.n);

    results
}

/// Scan a single chunk for near-misses using two-phase strip-then-classify.
///
/// Phase 1: Strip-only sieve (fast, no Kummer). rem[i] > 1 = definitely not governor.
/// Phase 2: Scan rem[] for promising windows (≤1 position with rem > 1).
///          Only re-factor the ~0.1% of positions in promising windows.
fn scan_chunk(
    chunk_start: u64,
    chunk_len: u64,
    k: u32,
    prime_data: &[PrimeData],
    sieve: &PrimeSieve,
    global_start: u64,
    global_count: u64,
) -> Vec<NearMiss> {
    let mut results = Vec::new();
    let window_size = k as u64 + 1;

    if chunk_len < window_size {
        return results;
    }

    // === PHASE 1: Strip-only sieve (hot path, no Kummer) ===
    let rem = strip_only_sieve(chunk_start, chunk_len as usize, prime_data);

    // === PHASE 2: Scan rem[] for promising near-miss windows ===
    // A "candidate non-governor" has rem > 1 (large prime factor).
    // A "candidate governor" has rem == 1 (all factors stripped, but v_p unknown).
    // We want windows with exactly 1 non-governor. Since rem > 1 is a SUBSET of
    // true non-governors, we look for windows with 0 or 1 positions having rem > 1.
    //
    // - 0 rem > 1: all positions are candidate governors. Might be a governor run,
    //   or might have a v_p failure. Need to classify each via is_governor_cold.
    //   If exactly 1 fails v_p → near-miss.
    // - 1 rem > 1: that position is definitely non-governor. The rest are candidate
    //   governors. If all pass v_p → near-miss (the rem > 1 position is the sole non-governor).
    // - 2+ rem > 1: definitely not a near-miss. Skip.

    let len = chunk_len as usize;
    let k_usize = k as usize;
    let global_end = global_start + global_count;

    // Sliding window tracking count of rem > 1 positions
    let mut rem_fail_count: u32 = 0;
    let mut rem_fail_positions: Vec<usize> = Vec::with_capacity(k_usize + 2);

    // Initialize first window [0, k]
    for i in 0..=k_usize {
        if i < len {
            let n = chunk_start + i as u64;
            if n <= 1 || rem[i] > 1 {
                rem_fail_count += 1;
                rem_fail_positions.push(i);
            }
        }
    }

    // Process windows
    let mut process_window = |right: usize, fail_count: u32, fail_positions: &[usize]| {
        let n_top = chunk_start + right as u64;
        if n_top < global_start + k as u64 || n_top >= global_end {
            return;
        }

        if fail_count > 1 {
            return; // 2+ definite non-governors → can't be a near-miss
        }

        let window_left = right - k_usize;

        if fail_count == 1 {
            // Exactly 1 position has rem > 1 (definite non-governor).
            // Check if all OTHER positions are true governors via v_p.
            let fail_idx = fail_positions[0];
            let m = chunk_start + fail_idx as u64;

            // Verify all candidate positions are actual governors
            for i in window_left..=right {
                if i == fail_idx { continue; }
                let val = chunk_start + i as u64;
                if !is_governor_cold(val, prime_data) {
                    return; // Second non-governor → not a near-miss
                }
            }

            // This is a near-miss. Check barrier-prime filter.
            if has_barrier_prime_failure(m, k, sieve) {
                results.push(NearMiss {
                    n: n_top,
                    m,
                    j: (n_top - m) as u32,
                });
            }
        } else {
            // fail_count == 0: all positions have rem == 1 (candidate governors).
            // Need v_p check on each to find if exactly 1 fails.
            let mut non_gov_idx: Option<usize> = None;
            let mut multi_fail = false;

            for i in window_left..=right {
                let val = chunk_start + i as u64;
                if !is_governor_cold(val, prime_data) {
                    if non_gov_idx.is_some() {
                        multi_fail = true;
                        break;
                    }
                    non_gov_idx = Some(i);
                }
            }

            if multi_fail || non_gov_idx.is_none() {
                return; // 0 or 2+ non-governors
            }

            let fail_idx = non_gov_idx.unwrap();
            let m = chunk_start + fail_idx as u64;

            if has_barrier_prime_failure(m, k, sieve) {
                results.push(NearMiss {
                    n: n_top,
                    m,
                    j: (n_top - m) as u32,
                });
            }
        }
    };

    // Check first window
    process_window(k_usize, rem_fail_count, &rem_fail_positions);

    // Slide window
    for right in (k_usize + 1)..len {
        let leaving = right - k_usize - 1;

        // Remove leaving element
        if chunk_start + leaving as u64 <= 1 || rem[leaving] > 1 {
            rem_fail_count -= 1;
            if !rem_fail_positions.is_empty() && rem_fail_positions[0] == leaving {
                rem_fail_positions.remove(0);
            }
        }

        // Add entering element
        let n_enter = chunk_start + right as u64;
        if n_enter <= 1 || rem[right] > 1 {
            rem_fail_count += 1;
            rem_fail_positions.push(right);
        }

        // Trim stale
        let window_left = right - k_usize;
        while !rem_fail_positions.is_empty() && rem_fail_positions[0] < window_left {
            rem_fail_positions.remove(0);
        }

        // Only call process_window when ≤ 1 definite non-governor in window
        if rem_fail_count <= 1 {
            process_window(right, rem_fail_count, &rem_fail_positions);
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scan_known_witness_no_near_misses() {
        let n = 339_949_252u64;
        let k = 8u32;
        let start = n - k as u64;
        let count = k as u64 + 1;
        let sieve = PrimeSieve::for_range(n + 100);

        let results = scan_near_misses(start, count, k, sieve.primes(), &sieve, false);
        assert!(
            results.is_empty(),
            "Known witness block should have no near-misses"
        );
    }

    #[test]
    fn test_scan_small_range() {
        let start = 100u64;
        let count = 1000u64;
        let k = 2u32;
        let sieve = PrimeSieve::for_range(start + count + 100);
        let checker = erdos396::GovernorChecker::with_sieve(sieve.clone());

        let results = scan_near_misses(start, count, k, sieve.primes(), &sieve, false);

        for nm in &results {
            assert!(nm.n >= start + k as u64);
            assert!(nm.n < start + count);
            assert_eq!(nm.m, nm.n - nm.j as u64);
            assert!(nm.j <= k);
            assert!(!checker.is_governor_fast(nm.m));
            for i in 0..=k {
                if i != nm.j {
                    assert!(
                        checker.is_governor_fast(nm.n - i as u64),
                        "Block term n-{} = {} should be a governor",
                        i,
                        nm.n - i as u64
                    );
                }
            }
        }
    }

    #[test]
    fn test_scan_matches_original_k2() {
        // Compare batch scanner against known results from the original scanner
        // at a range where we know pure candidates exist.
        let start = 40000u64;
        let count = 2000u64;
        let k = 2u32;
        let sieve = PrimeSieve::for_range(start + count + 100);

        let results = scan_near_misses(start, count, k, sieve.primes(), &sieve, false);

        // n=40664 is the smallest k=2 non-governor witness — it should appear
        let has_40664 = results.iter().any(|nm| nm.n == 40664);
        assert!(has_40664, "Should find the known k=2 non-governor witness at n=40664");
    }
}
