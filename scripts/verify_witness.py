#!/usr/bin/env python3
"""
Independent (stdlib-only) verifier for Erdős Problem #396 witnesses.

Given (k, n), checks:
    n(n-1)...(n-k) | C(2n, n)

Method:
  - Factor all terms in the block {n-k, ..., n} (trial division using a prime sieve)
  - For each prime p appearing in the block, compare:
        demand_p = Σ v_p(n-i)         (i=0..k)
        supply_p = v_p(C(2n,n))       (Legendre: v_p((2n)!)-2v_p(n!))

This script is intended as an independent cross-check of Claim A (witness validity),
separate from the Rust implementation.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


KNOWN_WITNESSES: List[Tuple[int, int]] = [
    (1, 2),
    (2, 2480),
    (3, 8178),
    (4, 45153),
    (5, 3648841),
    (6, 7979090),
    (7, 101130029),
    (8, 339949252),
    (9, 17609764993),
    (10, 17609764994),
    (11, 1070858041585),
    (12, 5048891644646),
    (13, 18253129921842),
]


def primes_up_to(n: int) -> List[int]:
    """All primes p with 2 <= p <= n (stdlib-only sieve)."""
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start : n + 1 : step] = b"\x00" * (((n - start) // step) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def vp_factorial(n: int, p: int) -> int:
    """Legendre: v_p(n!)."""
    v = 0
    power = p
    while power <= n:
        v += n // power
        power *= p
    return v


def vp_central_binom(n: int, p: int) -> int:
    """v_p(C(2n,n)) via Legendre."""
    return vp_factorial(2 * n, p) - 2 * vp_factorial(n, p)


def factorize(m: int, primes: List[int]) -> Dict[int, int]:
    """Trial-division factorization using the provided primes."""
    if m <= 1:
        return {}
    out: Dict[int, int] = {}
    x = m
    for p in primes:
        if p * p > x:
            break
        if x % p != 0:
            continue
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        out[p] = e
    if x > 1:
        out[x] = out.get(x, 0) + 1
    return out


@dataclass(frozen=True)
class WitnessCheckResult:
    k: int
    n: int
    is_witness: bool
    failing_prime: int | None
    demand: Dict[int, int]
    supply: Dict[int, int]


def check_witness(k: int, n: int) -> WitnessCheckResult:
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if n <= k:
        return WitnessCheckResult(k=k, n=n, is_witness=False, failing_prime=None, demand={}, supply={})

    block_lo = n - k
    block_hi = n
    limit = math.isqrt(block_hi)
    primes = primes_up_to(limit)

    demand: Counter[int] = Counter()
    for t in range(block_lo, block_hi + 1):
        for p, e in factorize(t, primes).items():
            demand[p] += e

    supply: Dict[int, int] = {}
    failing: int | None = None
    for p in sorted(demand.keys()):
        s = vp_central_binom(n, p)
        supply[p] = s
        if demand[p] > s and failing is None:
            failing = p

    return WitnessCheckResult(
        k=k,
        n=n,
        is_witness=failing is None,
        failing_prime=failing,
        demand=dict(demand),
        supply=supply,
    )


def main(argv: Iterable[str]) -> int:
    parser = argparse.ArgumentParser(description="Verify Erdős 396 witnesses (stdlib-only).")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--known", action="store_true", help="Verify the known witness table.")
    mode.add_argument("-k", type=int, help="Witness parameter k.")
    parser.add_argument("-n", type=int, help="Candidate n.")
    parser.add_argument("--verbose", action="store_true", help="Print per-prime details.")
    args = parser.parse_args(list(argv))

    if args.known:
        all_ok = True
        for k, n in KNOWN_WITNESSES:
            r = check_witness(k, n)
            status = "PASS" if r.is_witness else "FAIL"
            print(f"k={k:2d}, n={n:14d}: {status}")
            if not r.is_witness:
                all_ok = False
                if r.failing_prime is not None:
                    p = r.failing_prime
                    print(f"  failing p={p}: demand={r.demand.get(p, 0)}, supply={r.supply.get(p, 0)}")
        return 0 if all_ok else 1

    if args.k is None or args.n is None:
        parser.error("When not using --known, both -k and -n are required.")

    r = check_witness(args.k, args.n)
    if r.is_witness:
        print(f"PASS: k={r.k}, n={r.n} is a witness")
    else:
        print(f"FAIL: k={r.k}, n={r.n} is NOT a witness")
        if r.failing_prime is not None:
            p = r.failing_prime
            print(f"  failing p={p}: demand={r.demand.get(p, 0)}, supply={r.supply.get(p, 0)}")

    if args.verbose:
        for p in sorted(r.demand.keys()):
            print(f"  p={p}: demand={r.demand[p]}, supply={r.supply[p]}")

    return 0 if r.is_witness else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

