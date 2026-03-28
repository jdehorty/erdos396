// erdos396-speed-rust — Highly optimized Rust port of the C++ witness search.
//
// Build: cargo build --release

use std::fs;
use std::io::Write;
use std::sync::atomic::{AtomicI32, AtomicU64, AtomicU32, Ordering::Relaxed};
use std::time::Instant;

type Prime = u32;

const CHUNK_SIZE: u64 = 10_000_000;
const BLOCK_SIZE: usize = 32768; // L2-resident (~256 KiB for u64)
const NUM_SEQ: u64 = 2;

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
    for p in 2..=limit {
        if is_prime[p] {
            primes.push(p as Prime);
        }
    }
    primes
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
// Exact check — cold path, called only on rare candidate witnesses
// ---------------------------------------------------------------------------
#[cold]
#[inline(never)]
fn exact_check(n: u64, k: u64, small_primes: &[Prime]) -> bool {
    let n_ones = n.count_ones() as u64;
    let nu2_prod = ((n - k - 1).count_ones() as u64)
        .wrapping_sub(n_ones)
        .wrapping_add(k + 1);
    if n_ones < nu2_prod {
        return false;
    }
    let two_n = n << 1;
    let nk1 = n - k - 1;
    for &p32 in small_primes {
        let p = p32 as u64;
        if p == 2 {
            continue;
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
    true
}

// ---------------------------------------------------------------------------
// Strip kernel — for primes <= 2k (no governor check needed)
// ---------------------------------------------------------------------------
#[inline(always)]
fn strip_prime<const P: u64>(start_j: usize, w: usize, rem: *mut u64) {
    let ps = P as usize;
    let mut j = start_j;
    while j < w {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                let mut t = r / P;
                while t % P == 0 {
                    t /= P;
                }
                *rem.add(j) = t;
            }
        }
        j += ps;
    }
}

#[inline(always)]
fn strip_prime_dyn(p: u64, start_j: usize, w: usize, rem: *mut u64) {
    let ps = p as usize;
    let p_inv = mod_inverse_u64(p);
    let max_quot = u64::MAX / p;
    let mut j = start_j;
    while j < w {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                let mut t = r.wrapping_mul(p_inv);
                while t.wrapping_mul(p_inv) <= max_quot {
                    t = t.wrapping_mul(p_inv);
                }
                *rem.add(j) = t;
            }
        }
        j += ps;
    }
}

// ---------------------------------------------------------------------------
// Kummer helpers
// ---------------------------------------------------------------------------
#[inline(always)]
fn strip_factors_const<const P: u64>(mut r: u64) -> u64 {
    r /= P;
    while r % P == 0 {
        r /= P;
    }
    r
}

#[inline(always)]
fn strip_factors_dyn(mut r: u64, p: u64) -> u64 {
    let p_inv = mod_inverse_u64(p);
    let max_quot = u64::MAX / p;
    r = r.wrapping_mul(p_inv);
    while r.wrapping_mul(p_inv) <= max_quot {
        r = r.wrapping_mul(p_inv);
    }
    r
}

/// Kummer check on c = n/p: returns true if v_p(C(2n,n)) > v_p(c).
/// Since v_p(n) = v_p(c) + 1 and v_p(C(2n,n)) = carries(c+c in base p),
/// the governor condition v_p(C(2n,n)) >= v_p(n) becomes carries > v_p(c).
#[inline(always)]
fn kummer_ok_const<const P: u64>(mut c: u64) -> bool {
    let mut nu = 0u64;
    while c % P == 0 {
        nu += 1;
        c /= P;
    }
    let mut carries = 0u64;
    let mut carry = 0u64;
    while c != 0 {
        let d = c % P;
        let ge = (d + d + carry >= P) as u64;
        carries += ge;
        if carries > nu {
            return true;
        }
        carry = ge;
        c /= P;
    }
    false
}

#[inline(always)]
fn kummer_ok_dyn(mut c: u64, p: u64) -> bool {
    let p_inv = mod_inverse_u64(p);
    let max_quot = u64::MAX / p;
    // Part 1: strip factors of p (exact division via modular inverse)
    let mut nu = 0u64;
    while c.wrapping_mul(p_inv) <= max_quot {
        nu += 1;
        c = c.wrapping_mul(p_inv);
    }
    // Part 2: base-p digit extraction (hardware div for correct general divmod)
    let mut carries = 0u64;
    let mut carry = 0u64;
    while c != 0 {
        let d = c % p;
        let ge = (d + d + carry >= p) as u64;
        carries += ge;
        if carries > nu {
            return true;
        }
        carry = ge;
        c /= p;
    }
    false
}

// ---------------------------------------------------------------------------
// Three-phase Kummer kernel — for primes > 2k
//
// Phase 1: c < ceil(P/2) → 0 carries, but v_p(n) >= 1 → reject
// Phase 2: ceil(P/2) <= c < P → 1 carry, v_p(n) = 1 → accept, single strip
// Phase 3: c >= P → full Kummer check
// ---------------------------------------------------------------------------
#[inline(always)]
fn kummer_prime<const P: u64>(
    mut c: u64,
    mut j: usize,
    w: usize,
    rem: *mut u64,
) {
    let ps = P as usize;
    let half = (P + 1) >> 1;

    // Phase 1: unconditional reject
    while j < w && c < half {
        unsafe {
            if *rem.add(j) != 0 {
                *rem.add(j) = 0;
            }
        }
        j += ps;
        c += 1;
    }

    // Phase 2: unconditional accept + single division (v_p(n) = 1)
    while j < w && c < P {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                *rem.add(j) = r / P;
            }
        }
        j += ps;
        c += 1;
    }

    // Phase 3: full Kummer
    while j < w {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                if kummer_ok_const::<P>(c) {
                    *rem.add(j) = strip_factors_const::<P>(r);
                } else {
                    *rem.add(j) = 0;
                }
            }
        }
        j += ps;
        c += 1;
    }
}

fn kummer_prime_dyn(p: u64, mut c: u64, mut j: usize, w: usize, rem: *mut u64) {
    let ps = p as usize;
    let half = (p + 1) >> 1;
    let p_inv = mod_inverse_u64(p);

    while j < w && c < half {
        unsafe {
            if *rem.add(j) != 0 {
                *rem.add(j) = 0;
            }
        }
        j += ps;
        c += 1;
    }
    while j < w && c < p {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                // Exact division via modular inverse (r is divisible by p)
                *rem.add(j) = r.wrapping_mul(p_inv);
            }
        }
        j += ps;
        c += 1;
    }
    while j < w {
        unsafe {
            let r = *rem.add(j);
            if r != 0 {
                if kummer_ok_dyn(c, p) {
                    *rem.add(j) = strip_factors_dyn(r, p);
                } else {
                    *rem.add(j) = 0;
                }
            }
        }
        j += ps;
        c += 1;
    }
}

// ---------------------------------------------------------------------------
// Dispatch macros — strip (small) and kummer (large) kernels
// ---------------------------------------------------------------------------
macro_rules! dispatch_strip {
    ($p:expr, $sj:expr, $w:expr, $rm:expr, [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => strip_prime::<$prime>($sj, $w, $rm),)*
            _ => strip_prime_dyn($p as u64, $sj, $w, $rm),
        }
    };
}

macro_rules! dispatch_kummer {
    ($p:expr, $sc:expr, $sj:expr, $w:expr, $rm:expr, [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => kummer_prime::<$prime>($sc, $sj, $w, $rm),)*
            _ => kummer_prime_dyn($p as u64, $sc, $sj, $w, $rm),
        }
    };
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
// Solver
// ---------------------------------------------------------------------------
fn solve(
    k: u64,
    start_l: u64,
    end_l: u64,
    primes: &[Prime],
    n_threads: u32,
    run_log_threshold: u64,
) -> u64 {
    let cpb = n_threads as u64 * NUM_SEQ;
    let k2 = 2 * k;
    let mut l_batch = start_l;
    let t_start = Instant::now();

    // Pre-partition: strip-only primes (p <= 2k, skip p=2) and Kummer primes
    let strip_end = primes.partition_point(|&p| (p as u64) <= k2);
    let strip_primes = if strip_end > 1 { &primes[1..strip_end] } else { &[] as &[Prime] };
    let exact_primes = &primes[..strip_end]; // for exact_check (includes p=2)

    let best = AtomicU64::new(u64::MAX);
    let longest_seen = AtomicU32::new(0);

    loop {
        // Check end bound
        if l_batch >= end_l {
            return best.load(Relaxed);
        }

        let next_chunk = AtomicI32::new(0);

        std::thread::scope(|s| {
            for _ in 0..n_threads {
                s.spawn(|| {
                    let mut rem_buf = vec![0u64; BLOCK_SIZE];

                    loop {
                        let cid = next_chunk.fetch_add(1, Relaxed);
                        if cid >= cpb as i32 {
                            break;
                        }
                        let l = l_batch + cid as u64 * CHUNK_SIZE;
                        if l >= end_l {
                            continue;
                        }
                        let r = l + CHUNK_SIZE + k - 1;
                        if l > best.load(Relaxed) {
                            continue;
                        }

                        let max_p = isqrt_u64(r.saturating_mul(2)).max(k2);
                        let kummer_end =
                            primes.partition_point(|&p| (p as u64) <= max_p);
                        let kummer_primes = &primes[strip_end..kummer_end];

                        let mut consec = 0u64;
                        let mut local_best = best.load(Relaxed);

                        let mut bl = l;
                        while bl <= r {
                            let br = r.min(bl + BLOCK_SIZE as u64 - 1);
                            let w = (br - bl + 1) as usize;
                            let rm = rem_buf.as_mut_ptr();

                            // Init: fold p=2 strip into initialization
                            // rem[j] = odd part of n (strips all factors of 2)
                            {
                                let mut n = bl;
                                for j in 0..w {
                                    unsafe {
                                        *rm.add(j) = n >> n.trailing_zeros();
                                    }
                                    n += 1;
                                }
                            }

                            // Strip-only primes (p <= 2k, p > 2)
                            for &p32 in strip_primes {
                                let p = p32 as u64;
                                let sc = (bl + p - 1) / p;
                                let sj = (sc * p - bl) as usize;
                                if sj >= w {
                                    continue;
                                }
                                dispatch_strip!(
                                    p, sj, w, rm,
                                    [3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
                                     37, 41, 43, 47, 53, 59, 61, 67, 71,
                                     73, 79, 83, 89, 97, 101, 103, 107,
                                     109, 113, 127, 131, 137, 139, 149,
                                     151, 157, 163, 167, 173, 179, 181,
                                     191, 193, 197, 199]
                                );
                            }

                            // Kummer primes (p > 2k)
                            for &p32 in kummer_primes {
                                let p = p32 as u64;
                                let sc = (bl + p - 1) / p;
                                let sj = (sc * p - bl) as usize;
                                if sj >= w {
                                    continue;
                                }
                                dispatch_kummer!(
                                    p, sc, sj, w, rm,
                                    [3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
                                     37, 41, 43, 47, 53, 59, 61, 67, 71,
                                     73, 79, 83, 89, 97, 101, 103, 107,
                                     109, 113, 127, 131, 137, 139, 149,
                                     151, 157, 163, 167, 173, 179, 181,
                                     191, 193, 197, 199]
                                );
                            }

                            // Scan for consecutive governors
                            let mut n = bl;
                            for j in 0..w {
                                if n >= local_best {
                                    break;
                                }
                                unsafe {
                                    if *rm.add(j) == 1 {
                                        consec += 1;
                                        if consec >= k + 1
                                            && n > k
                                            && exact_check(n, k, exact_primes)
                                        {
                                            let prev =
                                                best.fetch_min(n, Relaxed);
                                            local_best = prev.min(n);
                                        }
                                    } else {
                                        if consec >= run_log_threshold {
                                            let run_end = n - 1;
                                            let run_start = run_end - consec + 1;
                                            let prev_longest = longest_seen.fetch_max(consec as u32, Relaxed);
                                            let marker = if consec as u32 > prev_longest { " *** NEW BEST ***" } else { "" };
                                            eprintln!("  run={} at [{}, {}]{}", consec, run_start, run_end, marker);
                                        }
                                        consec = 0;
                                    }
                                }
                                n += 1;
                            }

                            bl += BLOCK_SIZE as u64;
                        }
                    }
                });
            }
        });

        let result = best.load(Relaxed);
        if result != u64::MAX {
            return result;
        }
        l_batch += cpb * CHUNK_SIZE;
        write_checkpoint(k, l_batch);

        // Progress + ETA
        let elapsed = t_start.elapsed().as_secs_f64();
        let done = (l_batch - start_l) as f64;
        let total = (end_l - start_l) as f64;
        let pct = done / total * 100.0;
        let rate = done / elapsed;
        let remaining = (total - done) / rate;
        let mins = remaining / 60.0;
        eprintln!(
            "[{:.1}%] pos={:.3}T  {:.0} M/s  ETA {:.1}m  best_run={}",
            pct,
            l_batch as f64 / 1e12,
            rate / 1e6,
            mins,
            longest_seen.load(Relaxed),
        );
    }
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
            _ => {
                eprintln!(
                    "Usage: {} [-k K] [--start N] [--end N] [--kmax KMAX] [--threads T] [--run-log N]",
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
    println!(
        "erdos396-speed-rust | threads={} k={}..{} start={} end={} run_log>={}",
        n_threads, k_start, k_max, start_l, end_str, run_log
    );

    let t0 = Instant::now();
    let primes = sieve_primes(200_000_000);
    println!(
        "Primes: {} in {:.3}s\n",
        primes.len(),
        t0.elapsed().as_secs_f64()
    );

    for k in k_start..=k_max {
        let ts = Instant::now();
        let ans = solve(k, start_l, end_l, &primes, n_threads, run_log);
        let sec = ts.elapsed().as_secs_f64();
        let speed =
            if ans >= start_l { (ans - start_l + 1) as f64 } else { 1.0 } / sec;
        let line = format!(
            "k={:2}  n={:15}  {:.4}s  {:.2} M/s",
            k, ans, sec, speed / 1e6
        );
        println!("{}", line);

        if let Ok(mut f) = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open("results.txt")
        {
            let _ = writeln!(f, "{}", line);
        }

        start_l = ans;
        write_checkpoint(k + 1, start_l);
    }

    let _ = fs::remove_file("checkpoint.txt");
}
