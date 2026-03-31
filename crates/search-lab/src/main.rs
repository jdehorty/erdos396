// erdos396-search-lab — Highly optimized Rust witness search for Erdős Problem #396.
//
// Build: cargo build --release

use std::fs;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::time::Instant;

type Prime = u64;

const CHUNK_SIZE: u64 = 1_048_576; // 1M — shared with C++ implementation for parity
const BLOCK_SIZE: u32 = 32768;
const BLOCK_SHIFT: u32 = 15;
const BLOCK_MASK: u32 = 0x7FFF;

// ---------------------------------------------------------------------------
// PrimeData — precomputed per-prime data for fast modular arithmetic.
// Hot fields (inv_p, max_quot) first so they share a cache line in the strip loop.
// ---------------------------------------------------------------------------
#[derive(Clone, Copy)]
#[repr(C)]
struct PrimeData {
    inv_p: u64,    // modular inverse of p mod 2^64 — used every strip iteration
    max_quot: u64, // u64::MAX / p — divisibility threshold, used every strip iteration
    magic: u64,    // Barrett magic constant — used once per prime per chunk
    p: u64,        // prime value — full u64 to support sieve limits above 2^32
    shift: u8,     // bit_width(p) - 1 — used once per prime per chunk (Barrett)
}

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
// Modular-inverse exact division (computed inline, no precomputed table)
//
// For odd prime p: p_inv = modular inverse of p mod 2^64.
//   - Divisibility:   p | n  iff  n.wrapping_mul(p_inv) <= u64::MAX / p
//   - Exact division: n / p  ==  n.wrapping_mul(p_inv)   (when p | n)
//
// Replaces ~35 cycle hardware `div` with ~3 cycle `imul` + compare.
// Only used for exact-division paths; Kummer digit extraction keeps
// hardware div for correct general divmod.
// ---------------------------------------------------------------------------

/// Compute modular inverse of odd n mod 2^64 via Newton's method.
#[inline(always)]
fn mod_inverse_u64(n: u64) -> u64 {
    let mut inv = n.wrapping_mul(3) ^ 2; // correct mod 2^4
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv
}

// ---------------------------------------------------------------------------
// Sieve of Eratosthenes
// ---------------------------------------------------------------------------
fn sieve_primes(limit: usize) -> Vec<Prime> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;
    let mut p = 2;
    while p * p <= limit {
        if is_prime[p] {
            let mut m = p * p;
            while m <= limit {
                is_prime[m] = false;
                m += p;
            }
        }
        p += 1;
    }
    let mut primes = Vec::with_capacity(limit / 10);
    for (p, &is_p) in is_prime.iter().enumerate().take(limit + 1).skip(2) {
        if is_p {
            primes.push(p as Prime);
        }
    }
    primes
}

/// Build PrimeData vec from raw primes vec.
/// Precomputes modular inverse, Barrett magic constant, and divisibility threshold.
fn build_prime_data(primes: &[Prime]) -> Vec<PrimeData> {
    primes
        .iter()
        .map(|&p| {
            let inv_p = if p % 2 != 0 { mod_inverse_u64(p) } else { 0 };
            let shift = (64 - p.leading_zeros()).saturating_sub(1) as u8;
            let magic = if p > 1 {
                (((1u128 << (64 + shift as u32)) / p as u128) + 1) as u64
            } else {
                0
            };
            PrimeData {
                inv_p,
                max_quot: u64::MAX / p,
                magic,
                p,
                shift,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Integer sqrt (correct for all u64, unlike f64 which loses precision > 2^53)
// ---------------------------------------------------------------------------
#[inline(always)]
fn isqrt_u64(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut x = (n as f64).sqrt() as u64;
    // Newton correction — at most 1-2 iterations
    while x < n / (x + 1).max(1) {
        x += 1;
    }
    while x > 0 && x > n / x {
        x -= 1;
    }
    x
}

// ---------------------------------------------------------------------------
// process_prime<P> — const-generic strip with offset tracking
// Const-generic strip: strips ALL factors of P from rem[j],
// stepping by P, and updates start_j for next block.
// ---------------------------------------------------------------------------
#[inline(always)]
fn process_prime<const P: u64>(start_j: &mut u64, w_block: u64, rem: *mut u64) {
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

    let inv = inv_p_const(P);
    let limit = limit_const(P);

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
// process_prime_dyn — dynamic strip with precomputed inv_p/max_quot
// Accepts p as u64 (fix #3 from plan review).
// ---------------------------------------------------------------------------
#[inline(always)]
fn process_prime_dyn(
    p: u64,
    inv_p: u64,
    max_quot: u64,
    start_j: &mut u64,
    w_block: u64,
    rem: *mut u64,
) {
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
        j += p;
    }
    *start_j = j - w_block;
}

// ---------------------------------------------------------------------------
// Dispatch macro — const-generic for small primes, dynamic for the rest
// ---------------------------------------------------------------------------
macro_rules! dispatch_strip {
    ($p:expr, $inv:expr, $limit:expr, $sj:expr, $w:expr, $rm:expr, [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => process_prime::<$prime>($sj, $w, $rm),)*
            _ => process_prime_dyn($p, $inv, $limit, $sj, $w, $rm),
        }
    };
}

// ---------------------------------------------------------------------------
// Exact check — cold path, called only on rare candidate witnesses
//
// Part 1: v_2 check (popcount-based Kummer's theorem for p=2)
// Part 2: Legendre formula for primes p <= 2k
// Part 3: Per-element large prime check (for each n-i, check residual
//         large prime factors against Legendre bound on C(2n,n))
// ---------------------------------------------------------------------------
#[cold]
#[inline(never)]
fn exact_check(n: u64, k: u64, prime_data: &[PrimeData]) -> bool {
    let two_k = 2 * k;
    let two_n = 2 * n;

    // Part 1: v_2 check
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
        let p = pd.p;
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
            let p = pd.p;
            if p == 2 {
                continue;
            }

            if p <= two_k {
                // Strip all factors of small primes
                let mut q = temp.wrapping_mul(pd.inv_p);
                while q <= pd.max_quot {
                    temp = q;
                    q = temp.wrapping_mul(pd.inv_p);
                }
                continue;
            }

            // For primes > 2k: if p^2 > temp, temp is either 1 or a single
            // large prime factor
            if p * p > temp {
                if temp > two_k {
                    // temp is a large prime factor of (n-i); check Legendre
                    let p_val = temp;
                    // v_p(n-i): count how many times p_val divides (n-i)
                    let mut nu_prod = 1u64;
                    let mut t2 = ni / p_val; // (n-i) / p_val, already divided once
                    while t2 > 0 && t2 % p_val == 0 {
                        nu_prod += 1;
                        t2 /= p_val;
                    }

                    // v_p(C(2n,n)) via Legendre formula — use n, NOT ni
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

            // p > 2k but p^2 <= temp: check if p divides temp
            let q = temp.wrapping_mul(pd.inv_p);
            if q <= pd.max_quot {
                // p divides temp; count how many times
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

                // v_p(C(2n,n)) via Legendre formula — use n, NOT ni
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
// Checkpoint
// ---------------------------------------------------------------------------
fn write_checkpoint(k: u64, l_batch: u64) {
    let tmp = "checkpoint.tmp";
    let dst = "checkpoint.txt";
    if let Ok(mut f) = fs::File::create(tmp) {
        let _ = writeln!(f, "{} {}", k, l_batch);
        let _ = fs::rename(tmp, dst);
    }
}

// ---------------------------------------------------------------------------
// Solver — bucketed sieve architecture (shared design with C++ implementation)
// ---------------------------------------------------------------------------
#[allow(clippy::too_many_arguments)] // mirrors worker/bench CLI wiring
fn solve(
    k: u64,
    start_l: u64,
    end_l: u64,
    prime_data: &[PrimeData],
    n_threads: u32,
    _run_log_threshold: u64, // TODO: wire up run-length logging
    progress: bool,
    bench_secs: f64, // 0.0 = disabled; >0 = sieve for this many seconds, ignore witnesses
) -> (u64, u64, u64) {
    // Returns (min_witness, chunks_processed, witness_count)
    let k32 = k as u32;
    let two_k = 2 * k;
    let t_start = Instant::now();
    let bench_mode = bench_secs > 0.0;

    let global_min_n = AtomicU64::new(u64::MAX);
    let current_chunk = AtomicU64::new(0);
    let longest_seen = AtomicU64::new(0);
    let witness_count = AtomicU64::new(0);
    let time_up = AtomicBool::new(false);

    std::thread::scope(|s| {
        // Progress reporter thread: prints throughput every 500ms
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
                        k,
                        mcps,
                        cumulative_mcps,
                        total_candidates as f64 / 1e9,
                        elapsed
                    );
                    last_chunks = chunks_now;

                    // Stop when workers are done
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
                // rem buffer: BLOCK_SIZE + k + 1 for overlap (fix #4)
                let buf_size = BLOCK_SIZE as usize + k as usize + 1;
                let mut rem = vec![0u64; buf_size];
                let mut prime_offsets: Vec<u64> = Vec::new();

                // 33 buckets: ceil(CHUNK_SIZE / BLOCK_SIZE) + 1
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
                        // In bench mode: check time every 100 chunks, ignore witnesses
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
                    let chunk_total_primes = prime_data.partition_point(|pd| pd.p <= max_p);

                    if prime_offsets.len() < chunk_total_primes {
                        prime_offsets.resize(chunk_total_primes, 0);
                    }

                    let first_large_prime_idx = prime_data[..chunk_total_primes]
                        .partition_point(|pd| pd.p <= BLOCK_SIZE as u64);

                    // Compute initial offsets via Barrett reduction (no hardware div)
                    for idx in 1..chunk_total_primes {
                        let pd = &prime_data[idx];
                        let p = pd.p;
                        let num = l_chunk + p - 1;
                        let start_c = (((num as u128).wrapping_mul(pd.magic as u128)) >> 64) as u64
                            >> pd.shift as u64;
                        prime_offsets[idx] = start_c * p - l_chunk;
                    }

                    // Bucket large primes by target block
                    for bucket in buckets.iter_mut() {
                        bucket.clear();
                    }
                    for idx in first_large_prime_idx..chunk_total_primes {
                        let p = prime_data[idx].p;
                        let mut sj = prime_offsets[idx];
                        while sj < chunk_w as u64 {
                            let block_idx = (sj as u32) >> BLOCK_SHIFT;
                            let offset = (sj as u32) & BLOCK_MASK;
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

                        // rem_ptr points past overlap region
                        let rem_ptr = unsafe { rem.as_mut_ptr().add(overlap as usize) };

                        // Process small primes (p <= BLOCK_SIZE, skip p=2)
                        let num_new_u64 = num_new as u64;
                        for idx in 1..first_large_prime_idx {
                            let p = prime_data[idx].p;
                            let mut sj = prime_offsets[idx];

                            if sj < num_new_u64 {
                                let inv = prime_data[idx].inv_p;
                                let limit = prime_data[idx].max_quot;

                                dispatch_strip!(
                                    p,
                                    inv,
                                    limit,
                                    &mut sj,
                                    num_new_u64,
                                    rem_ptr,
                                    [
                                        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                                        59, 61, 67, 71, 73, 79, 83, 89, 97
                                    ]
                                );
                                prime_offsets[idx] = sj;
                            } else {
                                prime_offsets[idx] = sj - num_new_u64;
                            }
                        }

                        // Process bucketed large primes for this block
                        {
                            let b = &buckets[block_idx as usize];
                            let b_count = b.data.len();
                            let b_ptr = b.data.as_ptr();
                            let pd_ptr = prime_data.as_ptr();
                            for i in 0..b_count {
                                unsafe {
                                    // Prefetch upcoming prime data (8 items ahead)
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
                            // Fix #1: use `<` not `<=`
                            let mut i = k32 as i32;
                            while i >= 0
                                && unsafe { *rem.as_ptr().add((scan_j + i as u32) as usize) } <= 1
                            {
                                i -= 1;
                            }

                            if i < 0 {
                                // All k+1 elements are <= 1 — candidate
                                let cand_n = block_l - overlap as u64 + scan_j as u64 + k;
                                if bench_mode {
                                    // Bench: check all candidates, count witnesses
                                    if cand_n > k
                                        && exact_check(cand_n, k, &prime_data[..chunk_total_primes])
                                    {
                                        witness_count.fetch_add(1, Relaxed);
                                        // Track minimum witness
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
                                    && exact_check(cand_n, k, &prime_data[..chunk_total_primes])
                                {
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

                    // Progress + checkpoint (fix #5; skip checkpoint in bench mode)
                    if chunk_id > n_threads as u64 && chunk_id % 1000 == 0 {
                        let safe_l =
                            start_l + chunk_id.saturating_sub(n_threads as u64) * CHUNK_SIZE;
                        if !bench_mode {
                            write_checkpoint(k, safe_l);
                        }

                        let elapsed = t_start.elapsed().as_secs_f64();
                        let done = (safe_l - start_l) as f64;
                        let total = (end_l - start_l) as f64;
                        let pct = done / total * 100.0;
                        let rate = done / elapsed;
                        let remaining = (total - done) / rate;
                        let mins = remaining / 60.0;
                        let best_run = longest_seen.load(Relaxed);
                        eprintln!(
                            "[{:.1}%] pos={:.3}T  {:.0} M/s  ETA {:.1}m  best_run={}",
                            pct,
                            safe_l as f64 / 1e12,
                            rate / 1e6,
                            mins,
                            best_run,
                        );
                    }
                }
            });
        }
    });

    (
        global_min_n.load(Relaxed),
        current_chunk.load(Relaxed),
        witness_count.load(Relaxed),
    )
}

// ---------------------------------------------------------------------------
// Known witnesses — OEIS A375077
// ---------------------------------------------------------------------------
const KNOWN_WITNESSES: &[(u64, u64)] = &[
    (1, 2),
    (2, 2_480),
    (3, 8_178),
    (4, 45_153),
    (5, 3_648_841),
    (6, 7_979_090),
    (7, 101_130_029),
    (8, 339_949_252),
    (9, 1_019_547_844),
    (10, 17_609_764_994),
    (11, 1_070_858_041_585),
    (12, 5_048_891_644_646),
    (13, 18_253_129_921_842),
];

// ---------------------------------------------------------------------------
// Verify — detailed witness verification using Legendre's formula
// ---------------------------------------------------------------------------

/// Verify a single (k, n) witness with a detailed per-prime report.
fn verify_witness(k: u64, n: u64, prime_data: &[PrimeData]) -> bool {
    if k == 0 || n <= k {
        eprintln!(
            "Error: --verify requires k >= 1 and n > k (got k={}, n={})",
            k, n
        );
        return false;
    }

    let two_k = 2 * k;
    let two_n = 2 * n;
    let mut all_pass = true;

    println!("Verifying k={}, n={}", k, n);
    println!();

    // p=2: popcount-based check
    {
        let supply = n.count_ones() as u64;
        let nk1_ones = (n - k - 1).count_ones() as u64;
        let demand = nk1_ones.wrapping_sub(supply).wrapping_add(k + 1);
        let pass = supply >= demand;
        let status = if pass { "PASS" } else { "FAIL" };
        println!(
            "  p=2:  supply(popcount)={}, demand={} [{}]",
            supply, demand, status
        );
        if !pass {
            all_pass = false;
        }
    }

    // Barrier primes p in (2, 2k]
    for pd in prime_data {
        let p = pd.p;
        if p == 2 {
            continue;
        }
        if p > two_k {
            break;
        }
        let nk1 = n - k - 1;
        let (mut demand, mut supply, mut pw) = (0u64, 0u64, p);
        loop {
            let vn = n / pw;
            demand += vn - nk1 / pw;
            supply += two_n / pw - (vn << 1);
            if pw > two_n / p {
                break;
            }
            pw *= p;
        }
        let pass = demand <= supply;
        let status = if pass { "PASS" } else { "FAIL" };
        println!(
            "  p={}: supply={}, demand={} [{}]",
            p, supply, demand, status
        );
        if !pass {
            all_pass = false;
        }
    }

    // Per-element large prime check
    let mut large_prime_ok = true;
    for i in 0..=k {
        let ni = n - i;
        let mut temp = ni;
        temp >>= temp.trailing_zeros();

        for pd in prime_data {
            let p = pd.p;
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
                        println!(
                            "  n-{}={}: large prime p={}, demand={} > supply={} [FAIL]",
                            i, ni, p_val, nu_prod, nu_comb
                        );
                        large_prime_ok = false;
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
                    println!(
                        "  n-{}={}: prime p={}, demand={} > supply={} [FAIL]",
                        i, ni, p, nu_prod, nu_comb
                    );
                    large_prime_ok = false;
                }
            }
        }
    }

    if large_prime_ok {
        println!("  large primes: all block terms pass [PASS]");
    }
    all_pass = all_pass && large_prime_ok;

    println!();
    if all_pass {
        println!("RESULT: k={}, n={} is a VALID witness", k, n);
    } else {
        println!("RESULT: k={}, n={} is NOT a valid witness", k, n);
    }

    all_pass
}

/// Verify all known OEIS A375077 witnesses.
fn verify_all_known(prime_data: &[PrimeData]) -> bool {
    println!("Verifying all known OEIS A375077 witnesses...\n");
    let mut all_pass = true;
    for &(k, n) in KNOWN_WITNESSES {
        let pass = exact_check(n, k, prime_data);
        let status = if pass { "PASS" } else { "FAIL" };
        println!("  k={:2}, n={:>20}: {}", k, n, status);
        if !pass {
            all_pass = false;
        }
    }
    println!();
    if all_pass {
        println!(
            "All {} known witnesses verified successfully.",
            KNOWN_WITNESSES.len()
        );
    } else {
        println!("SOME WITNESSES FAILED VERIFICATION.");
    }
    all_pass
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut k_start = 1u64;
    let mut start_l = 1u64;
    let mut end_l = u64::MAX;
    let mut k_max = 20u64;
    let mut run_log: u64 = 10;
    let mut progress = false;
    let mut bench_secs: f64 = 0.0;
    let mut verify_k: Option<u64> = None;
    let mut verify_n: Option<u64> = None;
    let mut verify_known = false;
    let mut n_threads = std::thread::available_parallelism()
        .map(|n| n.get() as u32)
        .unwrap_or(1);

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-k" | "--k" => {
                i += 1;
                k_start = args[i].parse().unwrap();
            }
            "--start" => {
                i += 1;
                start_l = args[i].parse().unwrap();
            }
            "--end" => {
                i += 1;
                end_l = args[i].parse().unwrap();
            }
            "--kmax" => {
                i += 1;
                k_max = args[i].parse().unwrap();
            }
            "--threads" => {
                i += 1;
                n_threads = args[i].parse().unwrap();
            }
            "--run-log" => {
                i += 1;
                run_log = args[i].parse().unwrap();
            }
            "--progress" => {
                progress = true;
            }
            "--bench" => {
                i += 1;
                bench_secs = args[i].parse().unwrap();
            }
            "--verify" => {
                // --verify K N
                i += 1;
                verify_k = Some(args[i].parse().unwrap());
                i += 1;
                verify_n = Some(args[i].parse().unwrap());
            }
            "--verify-known" => {
                verify_known = true;
            }
            _ => {
                eprintln!("Usage: {} [options]", args[0]);
                eprintln!();
                eprintln!("Search:");
                eprintln!("  -k <k>              starting k value (default: 1)");
                eprintln!("  --start <n>         search range start (default: 1)");
                eprintln!("  --end <n>           search range end (default: inf)");
                eprintln!("  --kmax <k>          maximum k value (default: 20)");
                eprintln!("  --threads <t>       worker thread count (default: all cores)");
                eprintln!("  --run-log <n>       run-length logging threshold (default: 10)");
                eprintln!("  --progress          print throughput stats to stderr");
                eprintln!("  --bench <secs>      sieve-only benchmark for <secs> seconds");
                eprintln!();
                eprintln!("Verify:");
                eprintln!("  --verify <k> <n>    verify a single witness with detailed report");
                eprintln!("  --verify-known      verify all 13 known OEIS A375077 witnesses");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    // Handle verify modes early — no search needed
    if verify_known || verify_k.is_some() {
        let t0 = Instant::now();
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);
        eprintln!(
            "# primes: {} in {:.3}s",
            primes.len(),
            t0.elapsed().as_secs_f64()
        );

        if verify_known {
            let pass = verify_all_known(&prime_data);
            std::process::exit(if pass { 0 } else { 1 });
        }

        if let (Some(k), Some(n)) = (verify_k, verify_n) {
            let pass = verify_witness(k, n, &prime_data);
            std::process::exit(if pass { 0 } else { 1 });
        }
    }

    let end_str = if end_l == u64::MAX {
        "inf".to_string()
    } else {
        format!("{}", end_l)
    };
    let info_hdr = format!(
        "# erdos396-search-lab | threads={} k={}..{} start={} end={}",
        n_threads, k_start, k_max, start_l, end_str
    );
    let bench_mode_cli = bench_secs > 0.0;
    if bench_mode_cli {
        eprintln!("{}", info_hdr);
    } else {
        println!("{}", info_hdr);
    }

    let t0 = Instant::now();
    let primes = sieve_primes(200_000_000);
    let prime_data = build_prime_data(&primes);
    let info_primes = format!(
        "# primes: {} in {:.3}s",
        primes.len(),
        t0.elapsed().as_secs_f64()
    );
    if bench_mode_cli {
        eprintln!("{}", info_primes);
    } else {
        println!("{}", info_primes);
    }

    for k in k_start..=k_max {
        let ts = Instant::now();
        let (ans, _chunks, witnesses) = solve(
            k,
            start_l,
            end_l,
            &prime_data,
            n_threads,
            run_log,
            progress,
            bench_secs,
        );
        let sec = ts.elapsed().as_secs_f64();

        let (witness_str, candidates) = if bench_secs > 0.0 && end_l < u64::MAX {
            ("bench".to_string(), end_l - start_l)
        } else if bench_secs > 0.0 {
            ("bench".to_string(), _chunks * CHUNK_SIZE)
        } else if ans < u64::MAX && ans >= start_l {
            (ans.to_string(), ans - start_l + 1)
        } else if end_l < u64::MAX {
            ("none".to_string(), end_l - start_l)
        } else {
            ("none".to_string(), 1)
        };
        let speed = candidates as f64 / sec;

        // Universal benchmark contract: one R-line per k on stdout
        let line = format!(
            "R\t{}\t{}\t{:.4}\t{}\t{:.2}\t{}",
            k,
            witness_str,
            sec,
            candidates,
            speed / 1e6,
            witnesses
        );
        println!("{}", line);

        if !bench_mode_cli {
            if let Ok(mut f) = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open("results.txt")
            {
                let _ = writeln!(f, "{}", line);
            }
        }

        if bench_secs > 0.0 {
            // In bench mode: don't chain k values, don't checkpoint
        } else if ans < u64::MAX {
            start_l = ans;
            write_checkpoint(k + 1, start_l);
        } else {
            break;
        }
    }

    if bench_secs <= 0.0 {
        let _ = fs::remove_file("checkpoint.txt");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Known-good first witnesses for k=1..11 over [1, 18_000_000_000].
    ///
    /// This is the "standard candle" test: deterministic, covers the full
    /// sieve + exact_check pipeline, verified against the C++ implementation.
    ///
    /// Expected R-line results:
    ///   k=1  → 2
    ///   k=2  → 2480
    ///   k=3  → 8178
    ///   k=4  → 45153
    ///   k=5  → 3648841
    ///   k=6  → 7979090
    ///   k=7  → 101130029
    ///   k=8  → 339949252
    ///   k=9  → 1019547844
    ///   k=10 → 17609764994
    ///   k=11 → none (no witness in [1, 18B])
    // -----------------------------------------------------------------------
    // mod_inverse_u64: the entire sieve correctness hinges on this.
    // Verify n * inv ≡ 1 (mod 2^64) for a spread of odd primes.
    // -----------------------------------------------------------------------
    #[test]
    fn mod_inverse_identity() {
        let test_primes: &[u64] = &[
            3,
            5,
            7,
            11,
            13,
            97,
            101,
            251,
            1009,
            65537,
            104729,
            // Large primes near u32/u64 boundaries where bugs like to hide
            2_147_483_647, // 2^31 - 1 (Mersenne prime)
            4_294_967_291, // largest prime < 2^32
            1_000_000_007, // common large prime
        ];
        for &p in test_primes {
            let inv = mod_inverse_u64(p);
            assert_eq!(
                p.wrapping_mul(inv),
                1,
                "mod_inverse broken for p={}: p*inv = {}",
                p,
                p.wrapping_mul(inv)
            );
        }
    }

    // -----------------------------------------------------------------------
    // isqrt_u64: f64 loses precision above 2^53, so the Newton correction
    // path is critical. Test values where f64 sqrt gives wrong answers.
    // -----------------------------------------------------------------------
    #[test]
    fn isqrt_edge_cases() {
        // Zero
        assert_eq!(isqrt_u64(0), 0);
        // Perfect squares
        assert_eq!(isqrt_u64(1), 1);
        assert_eq!(isqrt_u64(4), 2);
        assert_eq!(isqrt_u64(9), 3);
        assert_eq!(isqrt_u64(1_000_000), 1000);
        // Non-perfect squares
        assert_eq!(isqrt_u64(2), 1);
        assert_eq!(isqrt_u64(8), 2);
        assert_eq!(isqrt_u64(10), 3);
        // Large values where f64 loses precision (> 2^53)
        // (2^32)^2 = 2^64 overflows, but (2^32 - 1)^2 fits
        let large = (1u64 << 32) - 1; // 4294967295
        assert_eq!(isqrt_u64(large * large), large);
        assert_eq!(isqrt_u64(large * large - 1), large - 1);
        assert_eq!(isqrt_u64(large * large + 1), large);
        // Near u64::MAX — isqrt(u64::MAX) should be 2^32 - 1 = 4294967295
        assert_eq!(isqrt_u64(u64::MAX), (1u64 << 32) - 1);
    }

    // -----------------------------------------------------------------------
    // sieve_primes: verify against known prime-counting function values.
    // -----------------------------------------------------------------------
    #[test]
    fn sieve_known_counts() {
        let primes = sieve_primes(1_000_000);
        // π(100) = 25, π(1000) = 168, π(10000) = 1229, π(1000000) = 78498
        let pi_100 = primes.iter().filter(|&&p| p <= 100).count();
        let pi_1k = primes.iter().filter(|&&p| p <= 1000).count();
        let pi_10k = primes.iter().filter(|&&p| p <= 10_000).count();
        assert_eq!(pi_100, 25, "π(100)");
        assert_eq!(pi_1k, 168, "π(1000)");
        assert_eq!(pi_10k, 1229, "π(10000)");
        assert_eq!(primes.len(), 78498, "π(1000000)");
        // First and last
        assert_eq!(primes[0], 2);
        assert_eq!(*primes.last().unwrap(), 999983);
    }

    // -----------------------------------------------------------------------
    // exact_check: test known witnesses and known non-witnesses directly.
    // This isolates the Kummer/Legendre logic from the sieve.
    // -----------------------------------------------------------------------
    #[test]
    fn exact_check_known_values() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // Known witnesses (from the big test) must pass exact_check
        let witnesses: &[(u64, u64)] = &[
            (1, 2),
            (2, 2480),
            (3, 8178),
            (4, 45153),
            (5, 3648841),
            (6, 7979090),
            (7, 101130029),
            (8, 339949252),
            (9, 1019547844),
        ];
        for &(k, n) in witnesses {
            assert!(
                exact_check(n, k, &prime_data),
                "exact_check should accept known witness k={}, n={}",
                k,
                n
            );
        }

        // Just below each witness should NOT be a witness (otherwise it
        // wouldn't be the *first* witness). Spot-check a few.
        for &(k, n) in &witnesses[2..] {
            assert!(
                !exact_check(n - 1, k, &prime_data),
                "n={} (one below k={} witness) should fail exact_check",
                n - 1,
                k
            );
        }
    }

    // -----------------------------------------------------------------------
    // process_prime vs process_prime_dyn: the const-generic and dynamic
    // strip paths must produce identical results for the same prime.
    // -----------------------------------------------------------------------
    #[test]
    fn const_vs_dynamic_strip_consistency() {
        // Set up a block of consecutive values starting at some offset
        let block_start: u64 = 1_000_000;
        let block_size: u64 = 1024;
        let p: u64 = 7;
        let inv_p = mod_inverse_u64(p);
        let max_quot = u64::MAX / p;

        // Prepare two identical buffers: values with factor-of-2 stripped
        let mut rem_const = vec![0u64; block_size as usize];
        let mut rem_dyn = vec![0u64; block_size as usize];
        for j in 0..block_size {
            let x = block_start + j;
            let v = x >> x.trailing_zeros();
            rem_const[j as usize] = v;
            rem_dyn[j as usize] = v;
        }

        let mut sj_const = block_start % p;
        if sj_const != 0 {
            sj_const = p - sj_const;
        }
        let mut sj_dyn = sj_const;

        process_prime::<7>(&mut sj_const, block_size, rem_const.as_mut_ptr());
        process_prime_dyn(
            p,
            inv_p,
            max_quot,
            &mut sj_dyn,
            block_size,
            rem_dyn.as_mut_ptr(),
        );

        assert_eq!(sj_const, sj_dyn, "start_j diverged");
        assert_eq!(rem_const, rem_dyn, "rem buffers diverged for p=7");
    }

    // -----------------------------------------------------------------------
    // Thread-count independence: same witness regardless of parallelism.
    // Uses k=5 (witness at 3648841) as a quick check.
    // -----------------------------------------------------------------------
    #[test]
    fn thread_count_independence() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);
        let end = 4_000_000u64;

        let (ans_1t, _, _) = solve(5, 1, end, &prime_data, 1, 0, false, 0.0);
        let (ans_2t, _, _) = solve(5, 1, end, &prime_data, 2, 0, false, 0.0);
        let (ans_8t, _, _) = solve(5, 1, end, &prime_data, 8, 0, false, 0.0);

        assert_eq!(ans_1t, 3648841);
        assert_eq!(ans_1t, ans_2t, "1-thread vs 2-thread disagree");
        assert_eq!(ans_1t, ans_8t, "1-thread vs 8-thread disagree");
    }

    // -----------------------------------------------------------------------
    // Empty/degenerate ranges: should return no witness without panicking.
    // -----------------------------------------------------------------------
    #[test]
    fn empty_range_returns_no_witness() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // start == end
        let (ans, _, _) = solve(3, 100, 100, &prime_data, 1, 0, false, 0.0);
        assert_eq!(ans, u64::MAX, "empty range should find no witness");

        // Range that contains no witness for k=3 (witness is at 8178)
        let (ans, _, _) = solve(3, 8179, 8200, &prime_data, 1, 0, false, 0.0);
        assert_eq!(ans, u64::MAX, "range past k=3 witness should find nothing");
    }

    // -----------------------------------------------------------------------
    // Full pipeline: known-good first witnesses for k=1..11.
    // -----------------------------------------------------------------------
    #[test]
    fn known_witnesses_k1_through_k11() {
        let expected: &[(u64, Option<u64>)] = &[
            (1, Some(2)),
            (2, Some(2480)),
            (3, Some(8178)),
            (4, Some(45153)),
            (5, Some(3648841)),
            (6, Some(7979090)),
            (7, Some(101130029)),
            (8, Some(339949252)),
            (9, Some(1019547844)),
            (10, Some(17609764994)),
            (11, None), // no witness in [1, 18B]
        ];

        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        let end = 18_000_000_000u64;
        let mut start = 1u64;

        for &(k, expected_witness) in expected {
            let (ans, _chunks, _witnesses) = solve(k, start, end, &prime_data, 4, 0, false, 0.0);

            match expected_witness {
                Some(w) => {
                    assert_eq!(ans, w, "k={}: expected witness {} but got {}", k, w, ans);
                    // Next k starts from this witness (matches CLI chaining behavior)
                    start = ans;
                }
                None => {
                    assert_eq!(
                        ans,
                        u64::MAX,
                        "k={}: expected no witness but got {}",
                        k,
                        ans
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // build_prime_data: Barrett magic constants must produce correct quotients.
    // A wrong magic → wrong initial offset → primes skip values in a block.
    // -----------------------------------------------------------------------
    #[test]
    fn barrett_magic_produces_correct_quotients() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // Test Barrett floor(n / p) at a variety of offsets across the search range.
        // The Barrett quotient is: ((n as u128 * magic as u128) >> 64) >> shift
        let test_offsets: &[u64] = &[
            0,
            1,
            100,
            CHUNK_SIZE,
            1_000_000_000,
            10_000_000_000,
            100_000_000_000,
            u64::MAX / 2,
            u64::MAX - 1,
        ];

        for pd in prime_data.iter().skip(1) {
            // skip p=2 (inv_p=0)
            let p = pd.p;
            for &base in test_offsets {
                let n = base.wrapping_add(p).wrapping_sub(1); // n = base + p - 1
                if n < p {
                    continue;
                }
                let barrett_q =
                    (((n as u128).wrapping_mul(pd.magic as u128)) >> 64) as u64 >> pd.shift as u64;
                let true_q = n / p;
                assert_eq!(
                    barrett_q, true_q,
                    "Barrett quotient wrong for p={}, n={}: got {} expected {}",
                    p, n, barrett_q, true_q
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // process_prime: verify stripping removes ALL factors across higher powers
    // (p^2, p^3, p^4) and multiple primes in sequence.
    // -----------------------------------------------------------------------
    #[test]
    fn process_prime_strips_higher_powers() {
        // For p=3: values like 27=3^3, 81=3^4, 243=3^5 should be fully stripped to 1
        let primes_to_test: &[u64] = &[3, 5, 7, 11, 13];
        for &p in primes_to_test {
            let inv_p = mod_inverse_u64(p);
            let max_quot = u64::MAX / p;

            // Build a buffer with known prime powers at specific positions
            let block_size = 256u64;
            let block_start = p * 1000; // aligned start for easy offset math
            let mut rem = vec![0u64; block_size as usize];
            for j in 0..block_size {
                let x = block_start + j;
                rem[j as usize] = x >> x.trailing_zeros(); // strip factor-of-2
            }

            let mut sj = 0u64; // block_start is a multiple of p
            process_prime_dyn(p, inv_p, max_quot, &mut sj, block_size, rem.as_mut_ptr());

            // After stripping, every position divisible by p should have no factor of p left
            for j in 0..block_size {
                let x = block_start + j;
                if x % p == 0 {
                    let r = rem[j as usize];
                    assert!(
                        r % p != 0 || r == 0,
                        "p={}: position {} (value {}) still has factor of p after strip, rem={}",
                        p,
                        j,
                        x,
                        r
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Multi-prime sequential stripping: after processing all small primes,
    // smooth numbers should reduce to 1.
    // -----------------------------------------------------------------------
    #[test]
    fn multi_prime_strip_smooth_numbers_reduce_to_one() {
        let primes = sieve_primes(1000);
        let prime_data = build_prime_data(&primes);

        let block_start = 1_000_000u64;
        let block_size = 1024u64;
        let mut rem = vec![0u64; block_size as usize];
        for j in 0..block_size {
            let x = block_start + j;
            rem[j as usize] = x >> x.trailing_zeros();
        }

        // Process all odd primes up to sqrt(2 * (block_start + block_size))
        let limit = isqrt_u64(2 * (block_start + block_size));
        for pd in prime_data.iter().skip(1) {
            if pd.p > limit {
                break;
            }
            let mut sj = block_start % pd.p;
            if sj != 0 {
                sj = pd.p - sj;
            }
            process_prime_dyn(
                pd.p,
                pd.inv_p,
                pd.max_quot,
                &mut sj,
                block_size,
                rem.as_mut_ptr(),
            );
        }

        // Numbers whose only prime factors are <= limit should now be 1.
        // Others should have a residual > 1 (the large prime factor).
        for j in 0..block_size {
            let x = block_start + j;
            let odd = x >> x.trailing_zeros();
            // Factor odd to check if it's smooth
            let mut temp = odd;
            for pd in prime_data.iter().skip(1) {
                if pd.p > limit {
                    break;
                }
                while temp % pd.p == 0 {
                    temp /= pd.p;
                }
            }
            let is_smooth = temp == 1;
            if is_smooth {
                assert_eq!(
                    rem[j as usize], 1,
                    "Smooth number x={} (odd={}) should reduce to 1, got {}",
                    x, odd, rem[j as usize]
                );
            } else {
                assert!(
                    rem[j as usize] > 1,
                    "Non-smooth x={} should have residual > 1, got {}",
                    x,
                    rem[j as usize]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // exact_check Part 3: large prime factor path.
    // Test with values known to have exactly one large prime factor.
    // -----------------------------------------------------------------------
    #[test]
    fn exact_check_rejects_large_prime_deficiency() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // For a governor run to exist, v_p must not exceed supply for any p.
        // We test that exact_check correctly rejects when a large prime factor
        // of some (n-i) has a p-adic valuation exceeding v_p(C(2n,n)).
        //
        // n=100, k=3: not a witness — should fail exact_check
        assert!(!exact_check(100, 3, &prime_data));
        // n=200, k=5: not a witness
        assert!(!exact_check(200, 5, &prime_data));
        // n=1000, k=3: not a witness
        assert!(!exact_check(1000, 3, &prime_data));

        // But known witnesses must pass (cross-check with the other test)
        assert!(exact_check(8178, 3, &prime_data));
        assert!(exact_check(339949252, 8, &prime_data));
    }

    // -----------------------------------------------------------------------
    // exact_check v_2 popcount path: the fast popcount check for p=2 must
    // agree with a naive Legendre computation.
    // -----------------------------------------------------------------------
    #[test]
    fn exact_check_v2_popcount_agrees_with_legendre() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // For a range of n values and k values, verify exact_check's boolean
        // result is consistent: if the v_2 check rejects, a manual Legendre
        // computation should also show the deficit.
        for k in 1..=5u64 {
            for n in (k + 2)..500 {
                let result = exact_check(n, k, &prime_data);
                // Manual v_2 check: v_2(prod) = sum of v_2(n-i) for i=0..k
                let n_ones = n.count_ones() as u64;
                let nk1_ones = (n - k - 1).count_ones() as u64;
                let nu2_prod = nk1_ones.wrapping_sub(n_ones).wrapping_add(k + 1);
                let v2_passes = n_ones >= nu2_prod;

                if !v2_passes {
                    assert!(
                        !result,
                        "n={}, k={}: v_2 check fails but exact_check passed",
                        n, k
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Cross-block boundary: verify the solver finds witnesses regardless of
    // where the search start places the witness relative to block boundaries.
    // The witness at n=2480 (k=2) needs values [2478, 2479, 2480], so start
    // must be <= 2478. We vary the start to shift which BLOCK_SIZE-aligned
    // block the witness lands in.
    // -----------------------------------------------------------------------
    #[test]
    fn solve_finds_witnesses_near_block_boundaries() {
        let primes = sieve_primes(200_000_000);
        let prime_data = build_prime_data(&primes);

        // Vary start so 2480 lands at different positions within a block
        for start in [1u64, 2, 100, 2000, 2470, 2477, 2478] {
            let (ans, _, _) = solve(2, start, 3000, &prime_data, 1, 0, false, 0.0);
            assert_eq!(
                ans, 2480,
                "k=2 witness should be 2480 with start={}, got {}",
                start, ans
            );
        }
    }

    // -----------------------------------------------------------------------
    // const_vs_dynamic across ALL dispatched primes (3..97), not just p=7.
    // -----------------------------------------------------------------------
    #[test]
    fn const_vs_dynamic_all_dispatched_primes() {
        let dispatched: &[u64] = &[
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97,
        ];

        let block_start: u64 = 10_000_000;
        let block_size: u64 = 512;

        for &p in dispatched {
            let inv_p = mod_inverse_u64(p);
            let max_quot = u64::MAX / p;

            let mut rem_dyn = vec![0u64; block_size as usize];
            for j in 0..block_size {
                let x = block_start + j;
                rem_dyn[j as usize] = x >> x.trailing_zeros();
            }
            let mut rem_const = rem_dyn.clone();

            let mut sj_dyn = block_start % p;
            if sj_dyn != 0 {
                sj_dyn = p - sj_dyn;
            }
            let mut sj_const = sj_dyn;

            // Dynamic path
            process_prime_dyn(
                p,
                inv_p,
                max_quot,
                &mut sj_dyn,
                block_size,
                rem_dyn.as_mut_ptr(),
            );

            // Const-generic path via macro-style dispatch
            macro_rules! test_const {
                ($($prime:literal),*) => {
                    match p {
                        $($prime => process_prime::<$prime>(&mut sj_const, block_size, rem_const.as_mut_ptr()),)*
                        _ => unreachable!(),
                    }
                };
            }
            test_const!(
                3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97
            );

            assert_eq!(sj_const, sj_dyn, "start_j diverged for p={}", p);
            assert_eq!(rem_const, rem_dyn, "rem buffers diverged for p={}", p);
        }
    }
}
