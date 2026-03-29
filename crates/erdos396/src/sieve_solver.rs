//! High-performance sieve solver — ported from search-lab.
//!
//! This module provides the same bucketed-sieve + sliding-window architecture
//! as the standalone search-lab binary, but as a library function callable from
//! the erdos396 CLI. It uses the shared [`PrimeData`] infrastructure from
//! [`crate::sieve`].
//!
//! Key difference from [`crate::search::parallel_search`]: this solver never
//! computes per-element `is_governor[i]`. Instead it strips all prime factors
//! into a `rem[]` array, scans for runs of `rem[j] <= 1`, and only calls
//! [`exact_check`] on the rare candidates. This makes the hot path ~100x
//! faster than the fused-sieve approach.

use crate::governor::vp_central_binom_dispatch;
use crate::sieve::{build_prime_data, PrimeData, PrimeSieve};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::time::Instant;

const CHUNK_SIZE: u64 = 1_048_576; // 1M — matches search-lab reference
const BLOCK_SIZE: u32 = 32768;
const BLOCK_SHIFT: u32 = 15;
const BLOCK_MASK: u32 = 0x7FFF;

// ---------------------------------------------------------------------------
// BucketItem / FastBucket — for bucketing large primes by target block
// ---------------------------------------------------------------------------
#[derive(Clone, Copy)]
struct BucketItem {
    p_idx: u32,
    offset: u32,
}

struct FastBucket {
    data: Vec<BucketItem>,
}

impl FastBucket {
    fn new() -> Self {
        FastBucket {
            data: Vec::with_capacity(16384),
        }
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.data.clear();
    }

    #[inline(always)]
    fn push(&mut self, p_idx: u32, offset: u32) {
        self.data.push(BucketItem { p_idx, offset });
    }
}

// ---------------------------------------------------------------------------
// Integer sqrt (correct for all u64)
// ---------------------------------------------------------------------------
#[inline(always)]
fn isqrt_u64(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut x = (n as f64).sqrt() as u64;
    while x < n / (x + 1).max(1) {
        x += 1;
    }
    while x > 0 && x > n / x {
        x -= 1;
    }
    x
}

// ---------------------------------------------------------------------------
// Const-generic strip with offset tracking (matches search-lab)
// ---------------------------------------------------------------------------
#[inline(always)]
fn process_prime<const P: u32>(start_j: &mut u32, w_block: u32, rem: *mut u64) {
    const fn inv_p_const(p: u64) -> u64 {
        let mut inv = p.wrapping_mul(3) ^ 2;
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        inv
    }

    const fn limit_const(p: u64) -> u64 {
        u64::MAX / p
    }

    let inv = inv_p_const(P as u64);
    let limit = limit_const(P as u64);

    let mut j = *start_j;
    while j < w_block {
        unsafe {
            let r = *rem.add(j as usize);
            let mut temp = r.wrapping_mul(inv);
            let mut q = temp.wrapping_mul(inv);
            if q <= limit {
                temp = q;
                loop {
                    q = temp.wrapping_mul(inv);
                    if q > limit {
                        break;
                    }
                    temp = q;
                }
            }
            *rem.add(j as usize) = temp;
        }
        j += P;
    }
    *start_j = j - w_block;
}

// ---------------------------------------------------------------------------
// Dynamic strip with precomputed inv_p/max_quot
// ---------------------------------------------------------------------------
#[inline(always)]
fn process_prime_dyn(
    p: u64,
    inv_p: u64,
    max_quot: u64,
    start_j: &mut u32,
    w_block: u32,
    rem: *mut u64,
) {
    let p32 = p as u32;
    let mut j = *start_j;
    while j < w_block {
        unsafe {
            let r = *rem.add(j as usize);
            let mut temp = r.wrapping_mul(inv_p);
            let mut q = temp.wrapping_mul(inv_p);
            if q <= max_quot {
                temp = q;
                loop {
                    q = temp.wrapping_mul(inv_p);
                    if q > max_quot {
                        break;
                    }
                    temp = q;
                }
            }
            *rem.add(j as usize) = temp;
        }
        j += p32;
    }
    *start_j = j - w_block;
}

// ---------------------------------------------------------------------------
// Dispatch macro
// ---------------------------------------------------------------------------
macro_rules! dispatch_strip {
    ($p:expr, $inv:expr, $limit:expr, $sj:expr, $w:expr, $rm:expr, [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => process_prime::<$prime>($sj, $w, $rm),)*
            _ => process_prime_dyn($p as u64, $inv, $limit, $sj, $w, $rm),
        }
    };
}

// ---------------------------------------------------------------------------
// Exact check — cold path, called only on rare candidate witnesses
// ---------------------------------------------------------------------------
#[cold]
#[inline(never)]
fn exact_check(n: u64, k: u64, prime_data: &[PrimeData]) -> bool {
    let two_k = 2 * k;
    let two_n = 2 * n;

    // Part 1: v_2 check (popcount-based Kummer's theorem for p=2)
    let n_ones = n.count_ones() as u64;
    let nu2_prod = ((n - k - 1).count_ones() as u64)
        .wrapping_sub(n_ones)
        .wrapping_add(k + 1);
    if n_ones < nu2_prod {
        return false;
    }

    // Part 2: Legendre formula for primes p in (2, 2k]
    let nk1 = n - k - 1;
    for pd in prime_data {
        let p = pd.p as u64;
        if p == 2 {
            continue;
        }
        if p > two_k {
            break;
        }
        let (mut nu_prod, mut nu_comb, mut pw) = (0u64, 0u64, p);
        loop {
            let vn = n / pw;
            nu_prod += vn - nk1 / pw;
            nu_comb += two_n / pw - (vn << 1);
            if pw > two_n / p {
                break;
            }
            pw *= p;
        }
        if nu_prod > nu_comb {
            return false;
        }
    }

    // Part 3: Per-element large prime check
    for i in 0..=k {
        let ni = n - i;
        let mut temp = ni;
        temp >>= temp.trailing_zeros();

        for pd in prime_data {
            let p = pd.p as u64;
            if p == 2 {
                continue;
            }

            if p <= two_k {
                let mut q = temp.wrapping_mul(pd.inv_p);
                while q <= pd.max_quot {
                    temp = q;
                    q = temp.wrapping_mul(pd.inv_p);
                }
                continue;
            }

            if p * p > temp {
                if temp > two_k {
                    let p_val = temp;
                    let mut nu_prod = 1u64;
                    let mut t2 = ni / p_val;
                    while t2 > 0 && t2 % p_val == 0 {
                        nu_prod += 1;
                        t2 /= p_val;
                    }

                    let mut nu_comb = 0u64;
                    let mut power = p_val;
                    loop {
                        let v_n = n / power;
                        nu_comb += two_n / power - 2 * v_n;
                        if power > two_n / p_val {
                            break;
                        }
                        power *= p_val;
                    }
                    if nu_prod > nu_comb {
                        return false;
                    }
                }
                break;
            }

            let q = temp.wrapping_mul(pd.inv_p);
            if q <= pd.max_quot {
                let mut nu_prod = 0u64;
                let mut q2 = q;
                loop {
                    nu_prod += 1;
                    temp = q2;
                    q2 = temp.wrapping_mul(pd.inv_p);
                    if q2 > pd.max_quot {
                        break;
                    }
                }

                let mut nu_comb = 0u64;
                let mut power = p;
                loop {
                    let v_n = n / power;
                    nu_comb += two_n / power - 2 * v_n;
                    if power > two_n / p {
                        break;
                    }
                    power *= p;
                }
                if nu_prod > nu_comb {
                    return false;
                }
            }
        }
    }

    true
}

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

/// Result from a sieve solver run.
pub struct SolveResult {
    /// Minimum witness found (u64::MAX if none)
    pub min_witness: u64,
    /// Total chunks processed
    pub chunks_processed: u64,
    /// Number of verified witnesses found
    pub witness_count: u64,
    /// Wall-clock duration
    pub duration: std::time::Duration,
}

// ---------------------------------------------------------------------------
// Solver — bucketed sieve architecture matching search-lab
// ---------------------------------------------------------------------------

/// Run the high-performance sieve solver for a single k value.
///
/// This is the search-lab solver integrated as a library function. It uses
/// the same bucketed-sieve + sliding-window architecture for maximum
/// throughput.
///
/// - `k`: target k (searches for runs of k+1 consecutive governors)
/// - `start_l`, `end_l`: search range
/// - `prime_data`: precomputed prime data (from [`build_prime_data`])
/// - `n_threads`: number of worker threads
/// - `bench_secs`: if > 0, stop after this many seconds (bench mode)
/// - `progress`: if true, emit progress lines to stderr
pub fn solve(
    k: u64,
    start_l: u64,
    end_l: u64,
    prime_data: &[PrimeData],
    n_threads: u32,
    bench_secs: f64,
    progress: bool,
) -> SolveResult {
    let k32 = k as u32;
    let two_k = 2 * k;
    let t_start = Instant::now();
    let bench_mode = bench_secs > 0.0;

    let global_min_n = AtomicU64::new(u64::MAX);
    let current_chunk = AtomicU64::new(0);
    let witness_count = AtomicU64::new(0);
    let time_up = AtomicBool::new(false);

    std::thread::scope(|s| {
        // Progress reporter thread
        if progress {
            let current_chunk_ref = &current_chunk;
            let global_min_ref = &global_min_n;
            let time_up_ref = &time_up;
            s.spawn(move || {
                let mut last_chunks: u64 = 0;
                loop {
                    std::thread::sleep(std::time::Duration::from_millis(500));
                    let chunks_now = current_chunk_ref.load(Relaxed);
                    let delta = chunks_now - last_chunks;
                    let candidates_delta = delta * CHUNK_SIZE;
                    let mcps = candidates_delta as f64 / 0.5 / 1e6;
                    let total_candidates = chunks_now * CHUNK_SIZE;
                    let elapsed = t_start.elapsed().as_secs_f64();
                    let cumulative_mcps = total_candidates as f64 / elapsed / 1e6;
                    eprintln!(
                        "P\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.4}",
                        k, mcps, cumulative_mcps, total_candidates as f64 / 1e9, elapsed
                    );
                    last_chunks = chunks_now;

                    if time_up_ref.load(Relaxed) {
                        break;
                    }
                    let l_chunk = start_l + chunks_now * CHUNK_SIZE;
                    if l_chunk >= end_l || l_chunk > global_min_ref.load(Relaxed) {
                        break;
                    }
                }
            });
        }

        for _ in 0..n_threads {
            s.spawn(|| {
                let buf_size = BLOCK_SIZE as usize + k as usize + 1;
                let mut rem = vec![0u64; buf_size];
                let mut prime_offsets: Vec<u32> = Vec::new();

                let num_buckets = (CHUNK_SIZE as u32).div_ceil(BLOCK_SIZE) + 1;
                let mut buckets: Vec<FastBucket> =
                    (0..num_buckets).map(|_| FastBucket::new()).collect();

                loop {
                    let chunk_id = current_chunk.fetch_add(1, Relaxed);
                    let l_chunk = start_l + chunk_id * CHUNK_SIZE;

                    if l_chunk >= end_l {
                        break;
                    }
                    if bench_mode {
                        if chunk_id % 100 == 0 && t_start.elapsed().as_secs_f64() >= bench_secs {
                            time_up.store(true, Relaxed);
                        }
                        if time_up.load(Relaxed) {
                            break;
                        }
                    } else if l_chunk > global_min_n.load(Relaxed) {
                        break;
                    }

                    let effective_chunk = CHUNK_SIZE.min(end_l - l_chunk);
                    let r_chunk = l_chunk + effective_chunk + k - 1;
                    let chunk_w = (r_chunk - l_chunk + 1) as u32;

                    let max_p = isqrt_u64(r_chunk.saturating_mul(2)).max(two_k) + 1;
                    let chunk_total_primes =
                        prime_data.partition_point(|pd| (pd.p as u64) <= max_p);

                    if prime_offsets.len() < chunk_total_primes {
                        prime_offsets.resize(chunk_total_primes, 0);
                    }

                    let first_large_prime_idx = prime_data[..chunk_total_primes]
                        .partition_point(|pd| pd.p <= BLOCK_SIZE);

                    // Compute initial offsets via Barrett reduction
                    for idx in 1..chunk_total_primes {
                        let pd = &prime_data[idx];
                        let p = pd.p as u64;
                        let num = l_chunk + p - 1;
                        let start_c =
                            (((num as u128).wrapping_mul(pd.magic as u128)) >> 64) as u64
                                >> pd.shift as u64;
                        prime_offsets[idx] = (start_c * p - l_chunk) as u32;
                    }

                    // Bucket large primes by target block
                    for bucket in buckets.iter_mut() {
                        bucket.clear();
                    }
                    for idx in first_large_prime_idx..chunk_total_primes {
                        let p = prime_data[idx].p;
                        let mut sj = prime_offsets[idx];
                        while sj < chunk_w {
                            let block_idx = sj >> BLOCK_SHIFT;
                            let offset = sj & BLOCK_MASK;
                            buckets[block_idx as usize].push(idx as u32, offset);
                            sj += p;
                        }
                    }

                    let mut overlap: u32 = 0;
                    let mut scan_j: u32 = 0;

                    let mut block_start: u32 = 0;
                    while block_start < chunk_w {
                        let block_idx = block_start >> BLOCK_SHIFT;
                        let block_l = l_chunk + block_start as u64;
                        let block_r = r_chunk.min(block_l + BLOCK_SIZE as u64 - 1);
                        let num_new = (block_r - block_l + 1) as u32;

                        // Initialize new elements: strip factor of 2
                        {
                            let mut x = block_l;
                            for idx_new in 0..num_new {
                                unsafe {
                                    *rem.as_mut_ptr().add((overlap + idx_new) as usize) =
                                        x >> x.trailing_zeros();
                                }
                                x += 1;
                            }
                        }

                        let rem_ptr = unsafe { rem.as_mut_ptr().add(overlap as usize) };

                        // Process small primes (skip p=2)
                        for idx in 1..first_large_prime_idx {
                            let p = prime_data[idx].p;
                            let mut sj = prime_offsets[idx];

                            if sj < num_new {
                                let inv = prime_data[idx].inv_p;
                                let limit = prime_data[idx].max_quot;

                                dispatch_strip!(
                                    p,
                                    inv,
                                    limit,
                                    &mut sj,
                                    num_new,
                                    rem_ptr,
                                    [
                                        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                                        59, 61, 67, 71, 73, 79, 83, 89, 97
                                    ]
                                );
                                prime_offsets[idx] = sj;
                            } else {
                                prime_offsets[idx] = sj - num_new;
                            }
                        }

                        // Process bucketed large primes
                        {
                            let b = &buckets[block_idx as usize];
                            let b_count = b.data.len();
                            let b_ptr = b.data.as_ptr();
                            let pd_ptr = prime_data.as_ptr();
                            for i in 0..b_count {
                                unsafe {
                                    if i + 8 < b_count {
                                        let future_idx = (*b_ptr.add(i + 8)).p_idx as usize;
                                        let ptr = pd_ptr.add(future_idx) as *const u8;
                                        #[cfg(target_arch = "aarch64")]
                                        std::arch::asm!(
                                            "prfm pldl2keep, [{ptr}]",
                                            ptr = in(reg) ptr,
                                            options(nostack, preserves_flags)
                                        );
                                        #[cfg(target_arch = "x86_64")]
                                        std::arch::x86_64::_mm_prefetch(
                                            ptr as *const i8,
                                            std::arch::x86_64::_MM_HINT_T1,
                                        );
                                    }

                                    let item = *b_ptr.add(i);
                                    let pd = &*pd_ptr.add(item.p_idx as usize);
                                    let r = *rem_ptr.add(item.offset as usize);
                                    let mut temp = r.wrapping_mul(pd.inv_p);
                                    let mut q = temp.wrapping_mul(pd.inv_p);
                                    if q <= pd.max_quot {
                                        temp = q;
                                        loop {
                                            q = temp.wrapping_mul(pd.inv_p);
                                            if q > pd.max_quot {
                                                break;
                                            }
                                            temp = q;
                                        }
                                    }
                                    *rem_ptr.add(item.offset as usize) = temp;
                                }
                            }
                        }

                        // Scan for consecutive governors (sliding window)
                        let w_search = overlap + num_new;
                        while scan_j + k32 < w_search {
                            let mut i = k32 as i32;
                            while i >= 0
                                && unsafe { *rem.as_ptr().add((scan_j + i as u32) as usize) } <= 1
                            {
                                i -= 1;
                            }

                            if i < 0 {
                                let cand_n =
                                    block_l - overlap as u64 + scan_j as u64 + k;
                                if bench_mode {
                                    if cand_n > k
                                        && exact_check(
                                            cand_n,
                                            k,
                                            &prime_data[..chunk_total_primes],
                                        )
                                    {
                                        witness_count.fetch_add(1, Relaxed);
                                        let mut current = global_min_n.load(Relaxed);
                                        while cand_n < current {
                                            match global_min_n.compare_exchange_weak(
                                                current, cand_n, Relaxed, Relaxed,
                                            ) {
                                                Ok(_) => break,
                                                Err(c) => current = c,
                                            }
                                        }
                                    }
                                } else if cand_n > k
                                    && cand_n < global_min_n.load(Relaxed)
                                    && exact_check(
                                        cand_n,
                                        k,
                                        &prime_data[..chunk_total_primes],
                                    )
                                {
                                    witness_count.fetch_add(1, Relaxed);
                                    let mut current = global_min_n.load(Relaxed);
                                    while cand_n < current {
                                        match global_min_n.compare_exchange_weak(
                                            current, cand_n, Relaxed, Relaxed,
                                        ) {
                                            Ok(_) => break,
                                            Err(c) => current = c,
                                        }
                                    }
                                }
                                scan_j += 1;
                            } else {
                                scan_j += i as u32 + 1;
                            }
                        }

                        // Carry overlap for cross-boundary detection
                        if w_search >= k32 {
                            let src = (w_search - k32) as usize;
                            unsafe {
                                let p = rem.as_mut_ptr();
                                std::ptr::copy(p.add(src), p, k32 as usize);
                            }
                            scan_j -= w_search - k32;
                            overlap = k32;
                        } else {
                            overlap = w_search;
                        }

                        block_start += BLOCK_SIZE;
                    }
                }
            });
        }
    });

    let duration = t_start.elapsed();
    SolveResult {
        min_witness: global_min_n.load(Relaxed),
        chunks_processed: current_chunk.load(Relaxed),
        witness_count: witness_count.load(Relaxed),
        duration,
    }
}

/// Convenience: sieve primes and build PrimeData in one call.
pub fn build_solver_primes(limit: u64) -> Vec<PrimeData> {
    let sieve = PrimeSieve::new(limit);
    build_prime_data(sieve.primes())
}
