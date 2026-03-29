// erdos396-speed-rust — Highly optimized Rust port of the C++ witness search.
//
// Build: cargo build --release

use std::fs;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::time::Instant;

type Prime = u32;

const CHUNK_SIZE: u64 = 1_048_576; // 1M — matches C++ reference
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
    p: u32,        // prime value — used for stepping
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
// Sieve of Eratosthenes (u32 primes — halves cache footprint)
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
            let p64 = p as u64;
            let inv_p = if p % 2 != 0 { mod_inverse_u64(p64) } else { 0 };
            // Barrett reduction: magic = ceil(2^(64+shift) / p)
            // shift = bit_width(p) - 1
            let shift = (32 - p.leading_zeros()).saturating_sub(1) as u8;
            let magic = if p > 1 {
                (((1u128 << (64 + shift as u32)) / p64 as u128) + 1) as u64
            } else {
                0
            };
            PrimeData {
                inv_p,
                max_quot: u64::MAX / p64,
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
// Matches C++ process_prime_p<p>: strips ALL factors of P from rem[j],
// stepping by P, and updates start_j for next block.
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

    // These are computed at compile time for each monomorphization
    let inv = inv_p_const(P as u64);
    let limit = limit_const(P as u64);

    let mut j = *start_j;
    while j < w_block {
        unsafe {
            let r = *rem.add(j as usize);
            // First division (always — we know p | (n-i) at this offset)
            let mut temp = r.wrapping_mul(inv);
            // Strip additional factors (rare)
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
// Dispatch macro — const-generic for small primes, dynamic for the rest
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
    // For each element n-i (i=0..k), strip small primes then check
    // any remaining large prime factor against Legendre bound.
    for i in 0..=k {
        let ni = n - i;
        let mut temp = ni;
        // Strip factor of 2
        temp >>= temp.trailing_zeros();

        for pd in prime_data {
            let p = pd.p as u64;
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
// Solver — bucketed sieve architecture matching C++ reference
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
                let mut prime_offsets: Vec<u32> = Vec::new();

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
                    let chunk_total_primes =
                        prime_data.partition_point(|pd| (pd.p as u64) <= max_p);

                    if prime_offsets.len() < chunk_total_primes {
                        prime_offsets.resize(chunk_total_primes, 0);
                    }

                    // Find first large prime (p > BLOCK_SIZE)
                    let first_large_prime_idx =
                        prime_data[..chunk_total_primes].partition_point(|pd| pd.p <= BLOCK_SIZE);

                    // Compute initial offsets via Barrett reduction (no hardware div)
                    for idx in 1..chunk_total_primes {
                        let pd = &prime_data[idx];
                        let p = pd.p as u64;
                        let num = l_chunk + p - 1;
                        // Barrett: start_c = ceil(num / p) = floor((num * magic) >> (64 + shift))
                        let start_c = (((num as u128).wrapping_mul(pd.magic as u128)) >> 64) as u64
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

                        // rem_ptr points past overlap region
                        let rem_ptr = unsafe { rem.as_mut_ptr().add(overlap as usize) };

                        // Process small primes (p <= BLOCK_SIZE, skip p=2)
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
            _ => {
                eprintln!(
                    "Usage: {} [-k K] [--start N] [--end N] [--kmax KMAX] [--threads T] [--run-log N] [--progress] [--bench SECS]",
                    args[0]
                );
                std::process::exit(1);
            }
        }
        i += 1;
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
        } else {
            if ans < u64::MAX {
                start_l = ans;
            }
            write_checkpoint(k + 1, start_l);
        }
    }

    if bench_secs <= 0.0 {
        let _ = fs::remove_file("checkpoint.txt");
    }
}
