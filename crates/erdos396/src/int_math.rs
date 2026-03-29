//! Integer-only math utilities.
//!
//! This module exists to avoid floating-point in correctness-critical code paths.
//! In particular, we rely on exact floor square roots for:
//! - trial division bounds (`p ≤ ⌊√n⌋`)
//! - Sanna's `√(2n)` large-prime-factor barrier

/// Floor square root: `isqrt(n) = ⌊√n⌋`.
///
/// Uses a branch-free bit-by-bit algorithm (Hacker's Delight style) that works
/// for all `u128` values without overflow.
#[inline]
pub fn isqrt_u128(n: u128) -> u128 {
    if n < 2 {
        return n;
    }

    // `bit` is the highest power of four ≤ n.
    let mut bit: u128 = 1u128 << 126;
    while bit > n {
        bit >>= 2;
    }

    let mut remainder = n;
    let mut root = 0u128;

    while bit != 0 {
        if remainder >= root + bit {
            remainder -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }

    root
}

/// Floor square root: `⌊√n⌋` for `u64`.
#[inline]
pub fn isqrt_u64(n: u64) -> u64 {
    isqrt_u128(n as u128) as u64
}

/// Floor square root: `⌊√(2n)⌋` computed exactly (no floating point).
#[inline]
pub fn isqrt_2n_u64(n: u64) -> u64 {
    isqrt_u128(2u128 * (n as u128)) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isqrt_u128_small() {
        for n in 0u128..=10_000 {
            let r = isqrt_u128(n);
            assert!(r * r <= n, "r^2 > n for n={}, r={}", n, r);
            assert!((r + 1) * (r + 1) > n, "(r+1)^2 <= n for n={}, r={}", n, r);
        }
    }

    #[test]
    fn test_isqrt_u64_edges() {
        // Perfect squares and neighbors
        for &x in &[0u64, 1, 2, 3, 4, 15, 16, 17, 24, 25, 26, 1_000_000] {
            let r = isqrt_u64(x);
            assert!(r * r <= x);
            assert!((r + 1) * (r + 1) > x);
        }

        // Large values: ensure no overflow and the inequality contract holds.
        for &x in &[u64::MAX, u64::MAX - 1, u64::MAX / 2, u64::MAX / 2 - 1] {
            let r = isqrt_u64(x);
            assert!((r as u128) * (r as u128) <= x as u128);
            assert!(((r + 1) as u128) * ((r + 1) as u128) > x as u128);
        }
    }

    #[test]
    fn test_isqrt_2n_u64_matches_definition() {
        for &n in &[0u64, 1, 2, 3, 4, 10, 123_456, 18_253_129_921_842] {
            let r = isqrt_2n_u64(n);
            let two_n = 2u128 * (n as u128);
            assert!((r as u128) * (r as u128) <= two_n);
            assert!(((r + 1) as u128) * ((r + 1) as u128) > two_n);
        }
    }
}
