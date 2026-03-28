#!/usr/bin/env python3
"""
Independent verifier for `search_report_*.json`.

This script checks the *coverage invariants* recorded by the Rust `erdos396` binary:
count, sum, and XOR of all scanned n values in [start, end).

It is intentionally stdlib-only, so reviewers can run it without dependencies.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable, List, Tuple


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


class ReportError(Exception):
    pass


def _load_worker_intervals(workers: list[dict[str, Any]]) -> List[Tuple[int, int, str]]:
    out: List[Tuple[int, int, str]] = []
    for idx, w in enumerate(workers):
        w_start = _as_int(w.get("start"), field=f"worker_stats[{idx}].start")
        w_end = _as_int(w.get("end"), field=f"worker_stats[{idx}].end")
        out.append((w_start, w_end, f"worker_stats[{idx}]"))
    return out


def check_report(path: Path) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))

    for key in ("target_k", "start", "end", "checked", "sum_checked", "xor_checked", "worker_stats"):
        if key not in data:
            raise ReportError(f"Missing required field {key!r}")

    k = _as_int(data["target_k"], field="target_k")
    start = _as_int(data["start"], field="start")
    end = _as_int(data["end"], field="end")
    if k < 1:
        raise ReportError(f"target_k must be >= 1, got {k}")
    if end <= start:
        raise ReportError(f"end must be > start, got start={start}, end={end}")

    expected_checked = end - start
    expected_sum = sum_range(start, end)
    expected_xor = xor_range(start, end)

    checked = _as_int(data.get("checked"), field="checked")
    expected_checked_report = _as_int(data.get("expected_checked"), field="expected_checked")
    if checked != expected_checked or expected_checked_report != expected_checked:
        raise ReportError(
            f"checked mismatch: report checked={checked}, expected_checked={expected_checked_report}; expected both {expected_checked}"
        )

    sum_checked = _as_int(data.get("sum_checked"), field="sum_checked")
    expected_sum_report = _as_int(data.get("expected_sum_checked"), field="expected_sum_checked")
    if sum_checked != expected_sum or expected_sum_report != expected_sum:
        raise ReportError(
            f"sum mismatch: report sum={sum_checked}, expected_sum={expected_sum_report}; expected both {expected_sum}"
        )

    xor_checked = _as_int(data.get("xor_checked"), field="xor_checked")
    expected_xor_report = _as_int(data.get("expected_xor_checked"), field="expected_xor_checked")
    if xor_checked != expected_xor or expected_xor_report != expected_xor:
        raise ReportError(
            f"xor mismatch: report xor={xor_checked}, expected_xor={expected_xor_report}; expected both {expected_xor}"
        )

    worker_stats = data.get("worker_stats")
    if not isinstance(worker_stats, list) or not all(isinstance(w, dict) for w in worker_stats):
        raise ReportError("worker_stats must be a list of objects")

    # Check worker partition of [start, end).
    intervals = _load_worker_intervals(worker_stats)
    intervals.sort(key=lambda t: (t[0], t[1]))

    if intervals[0][0] != start:
        raise ReportError(f"worker partition starts at {intervals[0][0]}, expected {start}")
    if intervals[-1][1] != end:
        raise ReportError(f"worker partition ends at {intervals[-1][1]}, expected {end}")
    for (a0, a1, a_name), (b0, b1, b_name) in zip(intervals, intervals[1:]):
        if a1 != b0:
            raise ReportError(
                f"worker intervals not contiguous: {a_name} ends at {a1}, but {b_name} starts at {b0}"
            )
        if not (a0 < a1 and b0 < b1):
            raise ReportError(f"invalid worker interval: {a_name}=[{a0},{a1}), {b_name}=[{b0},{b1})")

    # Check per-worker invariants aggregate correctly.
    agg_checked = 0
    agg_sum = 0
    agg_xor = 0
    for idx, w in enumerate(worker_stats):
        w_start = _as_int(w.get("start"), field=f"worker_stats[{idx}].start")
        w_end = _as_int(w.get("end"), field=f"worker_stats[{idx}].end")
        w_checked = _as_int(w.get("checked"), field=f"worker_stats[{idx}].checked")
        w_sum = _as_int(w.get("sum_checked"), field=f"worker_stats[{idx}].sum_checked")
        w_xor = _as_int(w.get("xor_checked"), field=f"worker_stats[{idx}].xor_checked")

        w_expected_checked = w_end - w_start
        w_expected_sum = sum_range(w_start, w_end)
        w_expected_xor = xor_range(w_start, w_end)
        if w_checked != w_expected_checked:
            raise ReportError(
                f"worker_stats[{idx}] checked mismatch: got {w_checked}, expected {w_expected_checked}"
            )
        if w_sum != w_expected_sum or w_xor != w_expected_xor:
            raise ReportError(
                f"worker_stats[{idx}] invariant mismatch: sum {w_sum} (expected {w_expected_sum}), "
                f"xor {w_xor} (expected {w_expected_xor})"
            )

        agg_checked += w_checked
        agg_sum += w_sum
        agg_xor ^= w_xor

    if agg_checked != expected_checked:
        raise ReportError(f"worker aggregate checked mismatch: got {agg_checked}, expected {expected_checked}")
    if agg_sum != expected_sum or agg_xor != expected_xor:
        raise ReportError(
            f"worker aggregate invariant mismatch: sum {agg_sum} (expected {expected_sum}), "
            f"xor {agg_xor} (expected {expected_xor})"
        )

    # Optional: check candidate/witness lists if present.
    def _check_sorted_unique_ints(field: str) -> list[int]:
        vals = data.get(field, [])
        if not isinstance(vals, list) or not all(isinstance(x, int) for x in vals):
            raise ReportError(f"{field} must be a list of integers")
        if vals != sorted(vals):
            raise ReportError(f"{field} list is not sorted")
        for a, b in zip(vals, vals[1:]):
            if a == b:
                raise ReportError(f"{field} list contains duplicates")
        for x in vals:
            if not (start <= x < end):
                raise ReportError(f"{field} value {x} is outside range [{start}, {end})")
        return vals

    candidates = _check_sorted_unique_ints("candidates")
    witnesses = _check_sorted_unique_ints("witnesses")
    false_positives = _check_sorted_unique_ints("false_positives")

    if data.get("verify_candidates") is True:
        cand_set = set(candidates)
        for w in witnesses:
            if w not in cand_set:
                raise ReportError(f"witness {w} not present in candidates list")
        for fp in false_positives:
            if fp not in cand_set:
                raise ReportError(f"false positive {fp} not present in candidates list")


def main(argv: Iterable[str]) -> int:
    parser = argparse.ArgumentParser(description="Check search_report_*.json invariants")
    parser.add_argument("report", type=Path, help="Path to report JSON file")
    args = parser.parse_args(list(argv))

    try:
        check_report(args.report)
    except (OSError, json.JSONDecodeError, ReportError, TypeError, ValueError) as e:
        print(f"{args.report}: FAIL: {e}", file=sys.stderr)
        return 1
    else:
        print(f"{args.report}: OK")
        return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
