#!/usr/bin/env python3
"""
Independent verifier for `validate_report_*.json`.

This script checks the *coverage invariants* recorded by the Rust `validate` binary:
count, sum, and XOR of all scanned n values (with the `n <= k` prefix excluded).

It is intentionally stdlib-only, so reviewers can run it without dependencies.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable, List


def _as_int(x: Any, *, field: str) -> int:
    if isinstance(x, int):
        return x
    if isinstance(x, str):
        try:
            return int(x, 10)
        except ValueError as e:
            raise ValueError(f"{field}: expected integer string, got {x!r}") from e
    raise TypeError(f"{field}: expected int or str, got {type(x).__name__}")


def sum_range(a: int, b: int) -> int:
    """Sum of integers in [a, b)."""
    if b <= a:
        return 0
    count = b - a
    first_plus_last = a + (b - 1)
    if count % 2 == 0:
        return (count // 2) * first_plus_last
    return count * (first_plus_last // 2)


def xor_upto(n: int) -> int:
    """XOR of integers in [0, n] inclusive."""
    r = n & 3
    if r == 0:
        return n
    if r == 1:
        return 1
    if r == 2:
        return n + 1
    return 0


def xor_range(a: int, b: int) -> int:
    """XOR of integers in [a, b)."""
    if b <= a:
        return 0
    hi = xor_upto(b - 1)
    lo = xor_upto(a - 1) if a > 0 else 0
    return hi ^ lo


def primes_up_to(n: int) -> List[int]:
    """All primes p with 2 <= p <= n."""
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


class ReportError(Exception):
    pass


def check_report(path: Path) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))

    for key in ("target_k", "start", "end", "scan_start", "barrier_primes"):
        if key not in data:
            raise ReportError(f"Missing required field {key!r}")

    k = _as_int(data["target_k"], field="target_k")
    start = _as_int(data["start"], field="start")
    end = _as_int(data["end"], field="end")
    scan_start = _as_int(data["scan_start"], field="scan_start")

    if k < 1:
        raise ReportError(f"target_k must be >= 1, got {k}")
    if end <= start:
        raise ReportError(f"end must be > start, got start={start}, end={end}")

    expected_scan_start = max(start, k + 1)
    if scan_start != expected_scan_start:
        raise ReportError(
            f"scan_start mismatch: report {scan_start}, expected {expected_scan_start}"
        )

    if scan_start >= end:
        expected_checked = 0
        expected_sum = 0
        expected_xor = 0
    else:
        expected_checked = end - scan_start
        expected_sum = sum_range(scan_start, end)
        expected_xor = xor_range(scan_start, end)

    checked = _as_int(data.get("checked"), field="checked")
    expected_checked_report = _as_int(data.get("expected_checked"), field="expected_checked")
    if checked != expected_checked or expected_checked_report != expected_checked:
        raise ReportError(
            f"checked mismatch: report checked={checked}, expected_checked={expected_checked_report}; expected both {expected_checked}"
        )

    sum_checked = _as_int(data.get("sum_checked"), field="sum_checked")
    expected_sum_report = _as_int(
        data.get("expected_sum_checked"), field="expected_sum_checked"
    )
    if sum_checked != expected_sum or expected_sum_report != expected_sum:
        raise ReportError(
            f"sum mismatch: report sum={sum_checked}, expected_sum={expected_sum_report}; expected both {expected_sum}"
        )

    xor_checked = _as_int(data.get("xor_checked"), field="xor_checked")
    expected_xor_report = _as_int(
        data.get("expected_xor_checked"), field="expected_xor_checked"
    )
    if xor_checked != expected_xor or expected_xor_report != expected_xor:
        raise ReportError(
            f"xor mismatch: report xor={xor_checked}, expected_xor={expected_xor_report}; expected both {expected_xor}"
        )

    barrier_primes = data["barrier_primes"]
    if not isinstance(barrier_primes, list) or not all(
        isinstance(p, int) for p in barrier_primes
    ):
        raise ReportError("barrier_primes must be a list of integers")
    expected_barrier = primes_up_to(2 * k)
    if barrier_primes != expected_barrier:
        raise ReportError(
            f"barrier_primes mismatch: report {barrier_primes}, expected {expected_barrier} (all primes p < 2k+1)"
        )

    witnesses = data.get("witnesses", [])
    if not isinstance(witnesses, list) or not all(isinstance(w, int) for w in witnesses):
        raise ReportError("witnesses must be a list of integers")
    if witnesses != sorted(witnesses):
        raise ReportError("witnesses list is not sorted")
    for a, b in zip(witnesses, witnesses[1:]):
        if a == b:
            raise ReportError("witnesses list contains duplicates")
    for w in witnesses:
        if not (start <= w < end):
            raise ReportError(f"witness {w} is outside range [{start}, {end})")
        if w <= k:
            raise ReportError(f"witness {w} violates necessary condition w >= k+1")


def _load_core_fields(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        "path": path,
        "k": _as_int(data["target_k"], field="target_k"),
        "start": _as_int(data["start"], field="start"),
        "end": _as_int(data["end"], field="end"),
        "scan_start": _as_int(data["scan_start"], field="scan_start"),
        "barrier_primes": data["barrier_primes"],
        "checked": _as_int(data["checked"], field="checked"),
        "sum_checked": _as_int(data["sum_checked"], field="sum_checked"),
        "xor_checked": _as_int(data["xor_checked"], field="xor_checked"),
        "witnesses": data.get("witnesses", []),
    }


def check_partition(paths: List[Path]) -> None:
    """
    Check that multiple reports (same k) form a contiguous, non-overlapping partition
    of the scanned range (using each report's [scan_start, end) interval).
    """
    if len(paths) < 2:
        return

    reports = [_load_core_fields(p) for p in paths]

    k0 = reports[0]["k"]
    barrier0 = reports[0]["barrier_primes"]
    for r in reports[1:]:
        if r["k"] != k0:
            raise ReportError(f"partition check requires same k: got {k0} and {r['k']}")
        if r["barrier_primes"] != barrier0:
            raise ReportError("partition check requires identical barrier_primes in all reports")

    nonempty = [r for r in reports if r["scan_start"] < r["end"]]
    if not nonempty:
        return

    nonempty.sort(key=lambda r: (r["scan_start"], r["end"]))
    for a, b in zip(nonempty, nonempty[1:]):
        if a["end"] != b["scan_start"]:
            raise ReportError(
                f"scanned intervals not contiguous: {a['path']} ends at {a['end']}, "
                f"but {b['path']} scan_start is {b['scan_start']}"
            )

    global_scan_start = nonempty[0]["scan_start"]
    global_end = nonempty[-1]["end"]

    agg_checked = sum(r["checked"] for r in nonempty)
    agg_sum = sum(r["sum_checked"] for r in nonempty)
    agg_xor = 0
    for r in nonempty:
        agg_xor ^= r["xor_checked"]

    expected_checked = global_end - global_scan_start
    expected_sum = sum_range(global_scan_start, global_end)
    expected_xor = xor_range(global_scan_start, global_end)

    if agg_checked != expected_checked:
        raise ReportError(
            f"partition checked mismatch: aggregated {agg_checked}, expected {expected_checked} "
            f"(scan [{global_scan_start}, {global_end}), k={k0})"
        )
    if agg_sum != expected_sum or agg_xor != expected_xor:
        raise ReportError(
            f"partition invariant mismatch: sum {agg_sum} (expected {expected_sum}), "
            f"xor {agg_xor} (expected {expected_xor})"
        )


def main(argv: Iterable[str]) -> int:
    parser = argparse.ArgumentParser(description="Check validate_report_*.json invariants")
    parser.add_argument("reports", nargs="+", type=Path, help="Path(s) to report JSON files")
    parser.add_argument(
        "--check-partition",
        action="store_true",
        help="Also check that the reports form a contiguous, non-overlapping partition "
        "of the scanned range (by [scan_start, end) intervals) and that invariants "
        "aggregate correctly.",
    )
    args = parser.parse_args(list(argv))

    ok = True
    for report in args.reports:
        try:
            check_report(report)
        except (OSError, json.JSONDecodeError, ReportError, TypeError, ValueError) as e:
            ok = False
            print(f"{report}: FAIL: {e}", file=sys.stderr)
        else:
            print(f"{report}: OK")

    if ok and args.check_partition:
        try:
            check_partition(list(args.reports))
        except (OSError, json.JSONDecodeError, ReportError, TypeError, ValueError) as e:
            ok = False
            print(f"partition: FAIL: {e}", file=sys.stderr)
        else:
            print("partition: OK")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
