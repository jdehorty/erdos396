//! High-performance governor sieve using search-lab architecture.
//!
//! Key optimizations:
//! 1. Modular-inverse exact division (~3 cycles vs ~35 for hardware div)
//! 2. Barrett fast-divmod for Kummer carry counting (no hardware div)
//! 3. L1-cache-sized block processing (32KB blocks)
//! 4. Bucketed large primes (avoid iterating all primes per block)
//! 5. Const-generic dispatch for first 24 primes
//! 6. Single fused pass: strip factors + check v_p inline (no redundant second pass)
//! 7. Barrett reduction for initial offset computation

const BLOCK_SIZE: usize = 32768;
const BLOCK_SHIFT: u32 = 15;
const BLOCK_MASK: u32 = 0x7FFF;

// ---------------------------------------------------------------------------
// PrimeData — precomputed per-prime data for fast modular arithmetic
// ---------------------------------------------------------------------------
#[derive(Clone, Copy)]
#[repr(C)]
pub struct PrimeData {
    // Hot fields first (used every strip iteration)
    inv_p: u64,    // modular inverse of p mod 2^64 — for exact division
    max_quot: u64, // u64::MAX / p — divisibility threshold
    // Barrett fields (used for Kummer and offset computation)
    barrett_magic: u64, // ceil(2^(64+shift) / p) — for fast divmod
    barrett_shift: u8,  // bit_width(p) - 1
    // Prime value
    pub p: u64,
    // Kummer thresholds (precomputed for branchless carry logic)
    half: u64,    // p / 2 (floor)
    half_up: u64, // ceil(p / 2)
}

/// Compute modular inverse of odd n mod 2^64 via Newton's method.
#[inline(always)]
fn mod_inverse_u64(n: u64) -> u64 {
    let mut inv = n.wrapping_mul(3) ^ 2;
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
    inv
}

/// Barrett fast divmod: compute (n / p, n % p) without hardware division.
/// Uses precomputed magic constant: q = ((n as u128 * magic) >> (64 + shift)) as u64
/// Then r = n - q * p. May need a single correction step.
#[inline(always)]
fn barrett_divmod(n: u64, p: u64, magic: u64, shift: u8) -> (u64, u64) {
    let q = (((n as u128).wrapping_mul(magic as u128)) >> 64) as u64 >> shift as u64;
    let r = n - q * p;
    // Correction: r might be >= p (off by 1)
    if r >= p {
        (q + 1, r - p)
    } else {
        (q, r)
    }
}

/// Kummer carry count using Barrett fast-divmod: v_p(C(2n, n)) without hardware div.
#[inline(always)]
fn kummer_barrett(n: u64, pd: &PrimeData) -> u64 {
    let mut carries = 0u64;
    let mut remaining = n;
    let mut carry = 0u64;
    while remaining > 0 {
        let (quot, digit) = barrett_divmod(remaining, pd.p, pd.barrett_magic, pd.barrett_shift);
        let threshold = if carry == 0 { pd.half_up } else { pd.half };
        carry = (digit >= threshold) as u64;
        carries += carry;
        remaining = quot;
    }
    carries
}

/// Dispatch Kummer to const-generic for small primes, Barrett for the rest.
/// p=2 uses popcount (1 cycle).
#[inline(always)]
fn kummer_dispatch(n: u64, pd: &PrimeData) -> u64 {
    use erdos396::governor::{vp_central_binom_kummer_const, vp_central_binom_p2};
    match pd.p {
        2 => vp_central_binom_p2(n),
        3 => vp_central_binom_kummer_const::<3>(n),
        5 => vp_central_binom_kummer_const::<5>(n),
        7 => vp_central_binom_kummer_const::<7>(n),
        11 => vp_central_binom_kummer_const::<11>(n),
        13 => vp_central_binom_kummer_const::<13>(n),
        17 => vp_central_binom_kummer_const::<17>(n),
        19 => vp_central_binom_kummer_const::<19>(n),
        23 => vp_central_binom_kummer_const::<23>(n),
        29 => vp_central_binom_kummer_const::<29>(n),
        31 => vp_central_binom_kummer_const::<31>(n),
        37 => vp_central_binom_kummer_const::<37>(n),
        41 => vp_central_binom_kummer_const::<41>(n),
        43 => vp_central_binom_kummer_const::<43>(n),
        47 => vp_central_binom_kummer_const::<47>(n),
        _ => kummer_barrett(n, pd),
    }
}

/// Build PrimeData vec from raw primes.
pub fn build_prime_data(primes: &[u64]) -> Vec<PrimeData> {
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
                barrett_magic: magic,
                barrett_shift: shift,
                p,
                half: p / 2,
                half_up: p / 2 + (p & 1),
            }
        })
        .collect()
}

/// Compute initial offset (first multiple of p >= lo) using Barrett reduction.
#[inline(always)]
fn barrett_offset(lo: u64, pd: &PrimeData) -> u64 {
    if lo == 0 {
        return pd.p;
    }
    let (_q, r) = barrett_divmod(lo, pd.p, pd.barrett_magic, pd.barrett_shift);
    if r == 0 {
        0
    } else {
        pd.p - r
    }
}

// ---------------------------------------------------------------------------
// Bucketed large primes
// ---------------------------------------------------------------------------
#[derive(Clone, Copy)]
struct BucketItem {
    pd_idx: u32,
    offset: u32,
}

struct Bucket {
    data: Vec<BucketItem>,
}

impl Bucket {
    fn new() -> Self {
        Bucket {
            data: Vec::with_capacity(256),
        }
    }
    #[inline(always)]
    fn push(&mut self, pd_idx: u32, offset: u32) {
        self.data.push(BucketItem { pd_idx, offset });
    }
}

// ---------------------------------------------------------------------------
// Const-generic fused strip + v_p check — unsafe raw pointers, no bounds checks
// ---------------------------------------------------------------------------
#[inline(always)]
fn strip_prime_const<const P: u64>(
    start_j: &mut u64,
    block_len: u64,
    rem: *mut u64,
    is_gov: *mut bool,
    lo: u64,
    block_offset: usize,
    pd: &PrimeData,
) {
    const fn inv_const(p: u64) -> u64 {
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

    let inv = inv_const(P);
    let limit = limit_const(P);
    let mut j = *start_j;
    while j < block_len {
        unsafe {
            let ji = j as usize;
            if !*is_gov.add(ji) {
                j += P;
                continue;
            }
            let r = *rem.add(ji);
            let q = r.wrapping_mul(inv);
            if q <= limit {
                let mut exp = 1u32;
                let mut temp = q;
                loop {
                    let q2 = temp.wrapping_mul(inv);
                    if q2 > limit {
                        break;
                    }
                    temp = q2;
                    exp += 1;
                }
                *rem.add(ji) = temp;
                let n = lo + block_offset as u64 + j;
                let supply = kummer_dispatch(n, pd);
                if (exp as u64) > supply {
                    *is_gov.add(ji) = false;
                }
            }
        }
        j += P;
    }
    *start_j = j - block_len;
}

/// Dynamic fused strip + v_p check — unsafe raw pointers, no bounds checks.
#[inline(always)]
fn strip_prime_dyn_fused(
    pd: &PrimeData,
    start_j: &mut u64,
    block_len: u64,
    rem: *mut u64,
    is_gov: *mut bool,
    lo: u64,
    block_offset: usize,
) {
    let mut j = *start_j;
    while j < block_len {
        unsafe {
            let ji = j as usize;
            if !*is_gov.add(ji) {
                j += pd.p;
                continue;
            }
            let r = *rem.add(ji);
            let q = r.wrapping_mul(pd.inv_p);
            if q <= pd.max_quot {
                let mut exp = 1u32;
                let mut temp = q;
                loop {
                    let q2 = temp.wrapping_mul(pd.inv_p);
                    if q2 > pd.max_quot {
                        break;
                    }
                    temp = q2;
                    exp += 1;
                }
                *rem.add(ji) = temp;
                let n = lo + block_offset as u64 + j;
                let supply = kummer_dispatch(n, pd);
                if (exp as u64) > supply {
                    *is_gov.add(ji) = false;
                }
            }
        }
        j += pd.p;
    }
    *start_j = j - block_len;
}

macro_rules! dispatch_strip_fused {
    ($p:expr, $pd:expr, $sj:expr, $blen:expr, $rem:expr, $is_gov:expr, $lo:expr, $bo:expr,
     [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => strip_prime_const::<$prime>($sj, $blen, $rem, $is_gov, $lo, $bo, $pd),)*
            _ => strip_prime_dyn_fused($pd, $sj, $blen, $rem, $is_gov, $lo, $bo),
        }
    };
}

// ---------------------------------------------------------------------------
// Fast governor sieve — single fused pass
// ---------------------------------------------------------------------------

/// Compute governor membership for [lo, lo+len) using:
/// - Modular-inverse exact division for factor stripping
/// - Barrett fast-divmod for Kummer carry counting
/// - L1-cache-sized 32KB block processing
/// - Bucketed large primes
/// - Single fused pass: strip + v_p check inline
pub fn fast_governor_sieve(lo: u64, len: usize, prime_data: &[PrimeData]) -> Vec<bool> {
    if len == 0 {
        return vec![];
    }

    let max_n = lo + len as u64 - 1;
    let sieve_limit = isqrt_2n(max_n);

    let total_primes = prime_data.partition_point(|pd| pd.p <= sieve_limit);
    let first_large_idx =
        prime_data[..total_primes].partition_point(|pd| pd.p <= BLOCK_SIZE as u64);

    // Per-prime offsets using Barrett reduction (no hardware div)
    let mut prime_offsets = vec![0u64; total_primes];
    // p=2 offset: special case
    if total_primes > 0 {
        prime_offsets[0] = if lo % 2 == 0 { 0 } else { 1 };
    }
    for idx in 1..total_primes {
        prime_offsets[idx] = barrett_offset(lo, &prime_data[idx]);
    }

    let mut rem = vec![0u64; len];
    let mut is_governor = vec![true; len];

    // Handle n <= 1
    for (i, gov) in is_governor.iter_mut().enumerate().take(len) {
        let n = lo + i as u64;
        if n <= 1 {
            *gov = n == 1;
        }
    }

    // Bucket large primes by target block
    let num_blocks = len.div_ceil(BLOCK_SIZE);
    let num_buckets = num_blocks + 1;
    let mut buckets: Vec<Bucket> = (0..num_buckets).map(|_| Bucket::new()).collect();

    for idx in first_large_idx..total_primes {
        let p = prime_data[idx].p;
        let mut sj = prime_offsets[idx];
        while (sj as usize) < len {
            let block_idx = (sj as u32) >> BLOCK_SHIFT;
            let offset = (sj as u32) & BLOCK_MASK;
            buckets[block_idx as usize].push(idx as u32, offset);
            sj += p;
        }
    }

    // === Single fused pass: strip factors + check v_p inline ===
    let mut block_start: usize = 0;
    while block_start < len {
        let block_end = (block_start + BLOCK_SIZE).min(len);
        let block_len = block_end - block_start;
        let block_lo = lo + block_start as u64;
        let block_idx = block_start >> BLOCK_SHIFT as usize;

        // Initialize rem[]: strip factor of 2 and check v_2 inline
        let pd2 = &prime_data[0]; // p=2
        for i in 0..block_len {
            let gi = block_start + i;
            if !is_governor[gi] {
                continue;
            }
            let n = block_lo + i as u64;
            if n == 0 {
                continue;
            }
            let tz = n.trailing_zeros();
            rem[gi] = n >> tz;
            if tz > 0 {
                let supply = kummer_dispatch(n, pd2);
                if (tz as u64) > supply {
                    is_governor[gi] = false;
                }
            }
        }

        // Raw pointers for this block (eliminates bounds checks in hot loops)
        let rem_ptr = unsafe { rem.as_mut_ptr().add(block_start) };
        let gov_ptr = unsafe { is_governor.as_mut_ptr().add(block_start) };

        // Fused strip + v_p check for small odd primes
        for idx in 1..first_large_idx {
            let pd = &prime_data[idx];
            let p = pd.p;
            let mut sj = prime_offsets[idx];

            if sj < block_len as u64 {
                dispatch_strip_fused!(
                    p,
                    pd,
                    &mut sj,
                    block_len as u64,
                    rem_ptr,
                    gov_ptr,
                    lo,
                    block_start,
                    [
                        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                        73, 79, 83, 89, 97
                    ]
                );
                prime_offsets[idx] = sj;
            } else {
                prime_offsets[idx] = sj - block_len as u64;
            }
        }

        // Fused strip + v_p check for bucketed large primes (with prefetching)
        {
            let bucket = &buckets[block_idx];
            let b_count = bucket.data.len();
            let b_ptr = bucket.data.as_ptr();
            let pd_ptr = prime_data.as_ptr();
            for i in 0..b_count {
                unsafe {
                    // Prefetch PrimeData 8 items ahead to hide L2/L3 latency
                    if i + 8 < b_count {
                        let future_idx = (*b_ptr.add(i + 8)).pd_idx as usize;
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
                    let ji = item.offset as usize;
                    if !*gov_ptr.add(ji) {
                        continue;
                    }
                    let pd = &*pd_ptr.add(item.pd_idx as usize);
                    let r = *rem_ptr.add(ji);
                    let q = r.wrapping_mul(pd.inv_p);
                    if q <= pd.max_quot {
                        let mut exp = 1u32;
                        let mut temp = q;
                        loop {
                            let q2 = temp.wrapping_mul(pd.inv_p);
                            if q2 > pd.max_quot {
                                break;
                            }
                            temp = q2;
                            exp += 1;
                        }
                        *rem_ptr.add(ji) = temp;
                        let n = lo + block_start as u64 + item.offset as u64;
                        let supply = kummer_dispatch(n, pd);
                        if (exp as u64) > supply {
                            *gov_ptr.add(ji) = false;
                        }
                    }
                }
            }
        }

        block_start = block_end;
    }

    // Post-sieve: remaining > 1 means large prime factor → not a governor
    for i in 0..len {
        if is_governor[i] && rem[i] > 1 {
            is_governor[i] = false;
        }
    }

    is_governor
}

// ---------------------------------------------------------------------------
// Strip-only sieve (Phase 1 of two-phase approach)
// No Kummer, no v_p check — just reduce each number to its residual.
// ---------------------------------------------------------------------------

/// Strip-only const-generic — no v_p check, no is_gov tracking.
#[inline(always)]
fn strip_only_const<const P: u64>(start_j: &mut u64, block_len: u64, rem: *mut u64) {
    const fn inv_const(p: u64) -> u64 {
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

    let inv = inv_const(P);
    let limit = limit_const(P);
    let mut j = *start_j;
    while j < block_len {
        unsafe {
            let r = *rem.add(j as usize);
            let mut temp = r.wrapping_mul(inv);
            if temp <= limit {
                loop {
                    let q = temp.wrapping_mul(inv);
                    if q > limit {
                        break;
                    }
                    temp = q;
                }
                *rem.add(j as usize) = temp;
            }
        }
        j += P;
    }
    *start_j = j - block_len;
}

#[inline(always)]
fn strip_only_dyn(
    inv_p: u64,
    max_quot: u64,
    p: u64,
    start_j: &mut u64,
    block_len: u64,
    rem: *mut u64,
) {
    let mut j = *start_j;
    while j < block_len {
        unsafe {
            let r = *rem.add(j as usize);
            let mut temp = r.wrapping_mul(inv_p);
            if temp <= max_quot {
                loop {
                    let q = temp.wrapping_mul(inv_p);
                    if q > max_quot {
                        break;
                    }
                    temp = q;
                }
                *rem.add(j as usize) = temp;
            }
        }
        j += p;
    }
    *start_j = j - block_len;
}

macro_rules! dispatch_strip_only {
    ($p:expr, $pd:expr, $sj:expr, $blen:expr, $rem:expr, [$($prime:literal),* $(,)?]) => {
        match $p {
            $($prime => strip_only_const::<$prime>($sj, $blen, $rem),)*
            _ => strip_only_dyn($pd.inv_p, $pd.max_quot, $p, $sj, $blen, $rem),
        }
    };
}

/// Strip-only sieve: reduce each number in [lo, lo+len) to its residual
/// after dividing out all primes up to sqrt(2*max_n).
///
/// Returns rem[] where rem[i] == 1 means "candidate governor" (all small prime
/// factors stripped, no large factor remaining). rem[i] > 1 means "definitely
/// not a governor" (has a large prime factor).
///
/// NO Kummer calls, NO v_p checks. This is the fast hot path.
pub fn strip_only_sieve(lo: u64, len: usize, prime_data: &[PrimeData]) -> Vec<u64> {
    if len == 0 {
        return vec![];
    }

    let max_n = lo + len as u64 - 1;
    let sieve_limit = isqrt_2n(max_n);

    let total_primes = prime_data.partition_point(|pd| pd.p <= sieve_limit);
    let first_large_idx =
        prime_data[..total_primes].partition_point(|pd| pd.p <= BLOCK_SIZE as u64);

    let mut prime_offsets = vec![0u64; total_primes];
    if total_primes > 0 {
        prime_offsets[0] = if lo % 2 == 0 { 0 } else { 1 };
    }
    for idx in 1..total_primes {
        prime_offsets[idx] = barrett_offset(lo, &prime_data[idx]);
    }

    let mut rem = vec![0u64; len];

    // Bucket large primes
    let num_blocks = len.div_ceil(BLOCK_SIZE);
    let mut buckets: Vec<Bucket> = (0..num_blocks + 1).map(|_| Bucket::new()).collect();
    for idx in first_large_idx..total_primes {
        let p = prime_data[idx].p;
        let mut sj = prime_offsets[idx];
        while (sj as usize) < len {
            let block_idx = (sj as u32) >> BLOCK_SHIFT;
            let offset = (sj as u32) & BLOCK_MASK;
            buckets[block_idx as usize].push(idx as u32, offset);
            sj += p;
        }
    }

    // Process in L1-cache-sized blocks
    let mut block_start: usize = 0;
    while block_start < len {
        let block_end = (block_start + BLOCK_SIZE).min(len);
        let block_len = block_end - block_start;
        let block_lo = lo + block_start as u64;
        let block_idx = block_start >> BLOCK_SHIFT as usize;

        // Initialize rem[]: strip factor of 2 via trailing_zeros
        for i in 0..block_len {
            let n = block_lo + i as u64;
            if n == 0 {
                unsafe {
                    *rem.as_mut_ptr().add(block_start + i) = 0;
                }
            } else {
                unsafe {
                    *rem.as_mut_ptr().add(block_start + i) = n >> n.trailing_zeros();
                }
            }
        }

        let rem_ptr = unsafe { rem.as_mut_ptr().add(block_start) };

        // Strip small odd primes
        for idx in 1..first_large_idx {
            let pd = &prime_data[idx];
            let p = pd.p;
            let mut sj = prime_offsets[idx];
            if sj < block_len as u64 {
                dispatch_strip_only!(
                    p,
                    pd,
                    &mut sj,
                    block_len as u64,
                    rem_ptr,
                    [
                        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                        73, 79, 83, 89, 97
                    ]
                );
                prime_offsets[idx] = sj;
            } else {
                prime_offsets[idx] = sj - block_len as u64;
            }
        }

        // Strip bucketed large primes with prefetching
        {
            let bucket = &buckets[block_idx];
            let b_count = bucket.data.len();
            let b_ptr = bucket.data.as_ptr();
            let pd_ptr = prime_data.as_ptr();
            for i in 0..b_count {
                unsafe {
                    if i + 8 < b_count {
                        let future_idx = (*b_ptr.add(i + 8)).pd_idx as usize;
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
                    let pd = &*pd_ptr.add(item.pd_idx as usize);
                    let r = *rem_ptr.add(item.offset as usize);
                    let mut temp = r.wrapping_mul(pd.inv_p);
                    if temp <= pd.max_quot {
                        loop {
                            let q = temp.wrapping_mul(pd.inv_p);
                            if q > pd.max_quot {
                                break;
                            }
                            temp = q;
                        }
                        *rem_ptr.add(item.offset as usize) = temp;
                    }
                }
            }
        }

        block_start = block_end;
    }

    rem
}

/// Full governor check for a single number. Used on the cold path for
/// the ~0.1% of positions that appear in promising near-miss windows.
/// Re-factors n from scratch using modular-inverse division.
pub fn is_governor_cold(n: u64, prime_data: &[PrimeData]) -> bool {
    if n <= 1 {
        return n == 1;
    }

    let sieve_limit = isqrt_2n(n);
    let mut remaining = n;

    // p=2
    let tz = remaining.trailing_zeros();
    if tz > 0 {
        remaining >>= tz;
        let supply = kummer_dispatch(n, &prime_data[0]);
        if (tz as u64) > supply {
            return false;
        }
    }

    // Odd primes
    for pd in &prime_data[1..] {
        if pd.p > sieve_limit {
            break;
        }
        if pd.p * pd.p > remaining {
            break;
        }

        let q = remaining.wrapping_mul(pd.inv_p);
        if q <= pd.max_quot {
            let mut exp = 1u64;
            let mut temp = q;
            loop {
                let q2 = temp.wrapping_mul(pd.inv_p);
                if q2 > pd.max_quot {
                    break;
                }
                temp = q2;
                exp += 1;
            }
            remaining = temp;
            let supply = kummer_dispatch(n, pd);
            if exp > supply {
                return false;
            }
            if remaining == 1 {
                return true;
            }
        }
    }

    // Remaining large prime
    if remaining > 1 {
        let sqrt_2n = isqrt_2n(n);
        if remaining > sqrt_2n {
            return false;
        }
        // Single prime <= sqrt(2n), exp=1. Check Kummer with hardware div (rare path).
        let half = remaining / 2;
        let half_up = half + (remaining & 1);
        let mut carries = 0u64;
        let mut rem_n = n;
        let mut carry = 0u64;
        while rem_n > 0 {
            let digit = rem_n % remaining;
            let threshold = if carry == 0 { half_up } else { half };
            carry = (digit >= threshold) as u64;
            carries += carry;
            rem_n /= remaining;
        }
        if 1 > carries {
            return false;
        }
    }

    true
}

/// Integer sqrt of 2n, exact.
#[inline]
fn isqrt_2n(n: u64) -> u64 {
    let two_n = 2u128 * n as u128;
    let mut x = (two_n as f64).sqrt() as u64;
    while (x as u128) * (x as u128) > two_n {
        x -= 1;
    }
    while ((x + 1) as u128) * ((x + 1) as u128) <= two_n {
        x += 1;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_inverse() {
        let inv3 = mod_inverse_u64(3);
        assert_eq!(3u64.wrapping_mul(inv3), 1);
        let inv7 = mod_inverse_u64(7);
        assert_eq!(7u64.wrapping_mul(inv7), 1);
        let inv97 = mod_inverse_u64(97);
        assert_eq!(97u64.wrapping_mul(inv97), 1);
    }

    #[test]
    fn test_barrett_divmod() {
        let primes = erdos396::PrimeSieve::for_range(1000);
        let pd = build_prime_data(primes.primes());
        for pdata in &pd[1..] {
            // skip p=2 (inv_p=0)
            let p = pdata.p;
            for &n in &[0u64, 1, p - 1, p, p + 1, 2 * p, 1000000007, 10000000000u64] {
                let (q, r) = barrett_divmod(n, p, pdata.barrett_magic, pdata.barrett_shift);
                assert_eq!(q, n / p, "quotient mismatch: n={}, p={}", n, p);
                assert_eq!(r, n % p, "remainder mismatch: n={}, p={}", n, p);
            }
        }
    }

    #[test]
    fn test_kummer_barrett_matches_dispatch() {
        let primes_sieve = erdos396::PrimeSieve::for_range(1000);
        let pd = build_prime_data(primes_sieve.primes());
        // Test primes > 47 (where Barrett kicks in instead of const-generic)
        let test_primes: Vec<&PrimeData> = pd.iter().filter(|p| p.p > 47 && p.p < 200).collect();
        for pdata in test_primes {
            for &n in &[100u64, 10000, 1000000, 10000000000u64] {
                let barrett_result = kummer_barrett(n, pdata);
                let dispatch_result = erdos396::governor::vp_central_binom_dispatch(n, pdata.p);
                assert_eq!(
                    barrett_result, dispatch_result,
                    "Kummer mismatch: n={}, p={}: barrett={}, dispatch={}",
                    n, pdata.p, barrett_result, dispatch_result
                );
            }
        }
    }

    #[test]
    fn test_fast_sieve_matches_fused() {
        use erdos396::prefilter::FusedBatchResult;

        let lo = 10_000_000_000u64;
        let len = 100_000usize;
        let primes = erdos396::PrimeSieve::for_range(lo + len as u64 + 100);
        let prime_data = build_prime_data(primes.primes());

        let erdos_prime_data = erdos396::build_prime_data(primes.primes());
        let fused = FusedBatchResult::compute(lo, len, &erdos_prime_data);
        let fast = fast_governor_sieve(lo, len, &prime_data);

        for (i, (f, g)) in fused.is_governor.iter().zip(fast.iter()).enumerate() {
            assert_eq!(
                f,
                g,
                "Mismatch at lo+{}={}: fused={}, fast={}",
                i,
                lo + i as u64,
                f,
                g
            );
        }
    }

    #[test]
    fn test_fast_sieve_small_range() {
        let lo = 2u64;
        let len = 10000usize;
        let primes = erdos396::PrimeSieve::for_range(lo + len as u64 + 100);
        let prime_data = build_prime_data(primes.primes());

        let erdos_prime_data = erdos396::build_prime_data(primes.primes());
        let fused = erdos396::prefilter::FusedBatchResult::compute(lo, len, &erdos_prime_data);
        let fast = fast_governor_sieve(lo, len, &prime_data);

        for (i, (f, g)) in fused.is_governor.iter().zip(fast.iter()).enumerate() {
            assert_eq!(
                f,
                g,
                "Mismatch at {}={}: fused={}, fast={}",
                i,
                lo + i as u64,
                f,
                g
            );
        }
    }
}
