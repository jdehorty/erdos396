#!/usr/bin/env python3
"""
Theorem-facing parquet exploration for the affine / prefix-core strategy.

This script is intentionally focused on:

- audited `run_window_classifications` for normalized bad-density statistics
- factorized small-prime family comparisons
- dynamic prefix-height proxies that mimic the first B=0 prefix-core face
- lightweight forward-looking checks on fresh `run_windows` data

It is not a proof tool. It is a reproducible hypothesis-generation tool for
finding theorem-friendly candidate class families.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import pyarrow as pa
import pyarrow.dataset as ds


CLASSIFICATION_SCHEMA = pa.schema(
    [
        ("schema_version", pa.string()),
        ("classification_id", pa.string()),
        ("window_id", pa.string()),
        ("window_length", pa.uint32()),
        ("window_k", pa.uint32()),
        ("window_start", pa.uint64()),
        ("window_end", pa.uint64()),
        ("window_n", pa.uint64()),
        ("audit_status", pa.string()),
        ("failing_prime", pa.uint64()),
        ("demand", pa.uint64()),
        ("supply", pa.uint64()),
        ("verifier_mode", pa.string()),
        ("verified_at_utc", pa.string()),
        ("build_git_hash", pa.string()),
    ]
)


RUN_WINDOW_SCHEMA = pa.schema(
    [
        ("schema_version", pa.string()),
        ("window_id", pa.string()),
        ("parent_run_id", pa.string()),
        ("checkpoint_id", pa.string()),
        ("source_label", pa.string()),
        ("archive_role", pa.string()),
        ("server", pa.string()),
        ("source_path", pa.string()),
        ("campaign_target_k", pa.uint32()),
        ("worker_id", pa.uint32()),
        ("max_run_start", pa.uint64()),
        ("max_run_end", pa.uint64()),
        ("max_run_length", pa.uint32()),
        ("window_length", pa.uint32()),
        ("window_k", pa.uint32()),
        ("window_offset", pa.uint32()),
        ("window_start", pa.uint64()),
        ("window_end", pa.uint64()),
        ("window_n", pa.uint64()),
        ("range_bucket_t", pa.uint64()),
        ("recorded_status", pa.string()),
        ("recorded_failing_prime", pa.uint64()),
        ("recorded_demand", pa.uint64()),
        ("recorded_supply", pa.uint64()),
        ("is_unique_coverage", pa.bool_()),
    ]
)


@dataclass(frozen=True)
class Family:
    name: str
    pred: Callable[[int], bool]


FAMILIES: tuple[Family, ...] = (
    Family("(3,2)", lambda n: (n % 4, n % 9) == (3, 2)),
    Family("(3,8)", lambda n: (n % 4, n % 9) == (3, 8)),
    Family("top2={(3,2),(3,8)}", lambda n: (n % 4, n % 9) in {(3, 2), (3, 8)}),
    Family(
        "top4={(3,2),(3,5),(3,7),(3,8)}",
        lambda n: (n % 4, n % 9) in {(3, 2), (3, 5), (3, 7), (3, 8)},
    ),
    Family("mod4=3 & mod9 in {2,8}", lambda n: n % 4 == 3 and n % 9 in {2, 8}),
    Family(
        "mod4=3 & mod9 in {2,5,7,8}",
        lambda n: n % 4 == 3 and n % 9 in {2, 5, 7, 8},
    ),
    Family("mod9 in {2,8}", lambda n: n % 9 in {2, 8}),
    Family("mod9 in {2,5,7,8}", lambda n: n % 9 in {2, 5, 7, 8}),
)


def classifications_dataset(root: Path) -> ds.Dataset:
    return ds.dataset(
        str(root / "data/run_corpus/v1/run_window_classifications"),
        format="parquet",
        schema=CLASSIFICATION_SCHEMA,
        partitioning="hive",
    )


def run_windows_dataset(root: Path) -> ds.Dataset:
    return ds.dataset(
        str(root / "data/run_corpus/v1/run_windows"),
        format="parquet",
        schema=RUN_WINDOW_SCHEMA,
        partitioning="hive",
    )


def audited_slice(dataset: ds.Dataset, length: int) -> tuple[list[int], list[str], list[int | None]]:
    table = dataset.scanner(
        columns=["window_n", "audit_status", "failing_prime"],
        filter=ds.field("window_length") == length,
    ).to_table()
    return (
        table["window_n"].to_pylist(),
        table["audit_status"].to_pylist(),
        table["failing_prime"].to_pylist(),
    )


def family_report(ns: list[int], status: list[str], primes: list[int | None], min_support: float = 0.01) -> list[str]:
    total = len(ns)
    base = {
        "all_fp": sum(1 for s in status if s == "false_positive") / total,
        "p2": sum(1 for s, p in zip(status, primes) if s == "false_positive" and p == 2) / total,
        "p3": sum(1 for s, p in zip(status, primes) if s == "false_positive" and p == 3) / total,
        "p5": sum(1 for s, p in zip(status, primes) if s == "false_positive" and p == 5) / total,
    }
    lines = [f"base={{{', '.join(f'{k}={base[k]:.8f}' for k in ('all_fp','p2','p3','p5'))}}}"]
    for family in FAMILIES:
        idx = [i for i, n in enumerate(ns) if family.pred(n)]
        support = len(idx) / total
        if support < min_support:
            continue
        rates = {
            "all_fp": sum(1 for i in idx if status[i] == "false_positive") / len(idx),
            "p2": sum(1 for i in idx if status[i] == "false_positive" and primes[i] == 2) / len(idx),
            "p3": sum(1 for i in idx if status[i] == "false_positive" and primes[i] == 3) / len(idx),
            "p5": sum(1 for i in idx if status[i] == "false_positive" and primes[i] == 5) / len(idx),
        }
        lifts = {
            key: (rates[key] / base[key] if base[key] else float("nan"))
            for key in rates
        }
        lines.append(
            f"{family.name}: support={support:.6f} "
            f"lift_all={lifts['all_fp']:.3f} lift_p2={lifts['p2']:.3f} "
            f"lift_p3={lifts['p3']:.3f} lift_p5={lifts['p5']:.3f}"
        )
    return lines


def dynamic_proxy_report(
    ns: list[int],
    status: list[str],
    primes: list[int | None],
    length: int,
    support_threshold: float,
) -> list[str]:
    d2 = int(math.log(length, 2))
    d3 = int(math.log(length, 3))
    m2 = 2 ** (d2 + 1)
    m3 = 3 ** (d3 + 1)
    total = len(ns)
    global_fp = sum(1 for s in status if s == "false_positive")
    global_23 = sum(
        1 for s, p in zip(status, primes) if s == "false_positive" and p in (2, 3)
    )
    counts: dict[tuple[int, int], int] = {}
    fp: dict[tuple[int, int], int] = {}
    fp23: dict[tuple[int, int], int] = {}
    fp5: dict[tuple[int, int], int] = {}

    for n, s, p in zip(ns, status, primes):
        key = (n % m2, n % m3)
        counts[key] = counts.get(key, 0) + 1
        if s == "false_positive":
            fp[key] = fp.get(key, 0) + 1
            if p in (2, 3):
                fp23[key] = fp23.get(key, 0) + 1
            elif p == 5:
                fp5[key] = fp5.get(key, 0) + 1

    ranked = []
    for key, c in counts.items():
        support = c / total
        if support < support_threshold:
            continue
        rate = fp.get(key, 0) / c
        rate23 = fp23.get(key, 0) / c
        ranked.append(
            (
                rate23 / (global_23 / total) if global_23 else float("nan"),
                rate / (global_fp / total) if global_fp else float("nan"),
                -support,
                key,
                support,
                fp.get(key, 0),
                fp23.get(key, 0),
                fp5.get(key, 0),
            )
        )
    ranked.sort()

    lines = [
        f"dynamic_proxy: m2={m2} m3={m3} support_threshold={support_threshold:.4f} "
        f"global23={global_23/total:.6f} globalfp={global_fp/total:.6f}"
    ]
    for row in ranked[:10]:
        lines.append(
            f"  key={row[3]} support={row[4]:.6f} lift23={row[0]:.3f} "
            f"lift={row[1]:.3f} fp={row[5]} fp23={row[6]} fp5={row[7]}"
        )
    return lines


def fresh_run_windows_report(root: Path) -> list[str]:
    dataset = run_windows_dataset(root)
    queries = [(12, 13), (13, 14)]
    family_names = [
        "mod4=3 & mod9 in {2,8}",
        "mod4=3 & mod9 in {2,5,7,8}",
        "mod9 in {2,8}",
        "mod9 in {2,5,7,8}",
    ]
    family_map = {f.name: f.pred for f in FAMILIES if f.name in family_names}
    lines: list[str] = []
    for campaign_target_k, window_length in queries:
        table = dataset.scanner(
            columns=["window_n"],
            filter=(ds.field("campaign_target_k") == campaign_target_k)
            & (ds.field("window_length") == window_length),
        ).to_table()
        ns = table["window_n"].to_pylist()
        total = len(ns)
        lines.append(
            f"fresh_run_windows: campaign_target_k={campaign_target_k} "
            f"window_length={window_length} total={total}"
        )
        for name in family_names:
            count = sum(1 for n in ns if family_map[name](n))
            support = count / total if total else 0.0
            lines.append(f"  {name}: support={support:.6f} count={count}")
    return lines


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/Users/jdehorty/src/witness-search-erdos396"),
        help="Repository root containing data/run_corpus/v1",
    )
    parser.add_argument(
        "--lengths",
        type=int,
        nargs="*",
        default=[7, 8, 9, 10, 11, 12, 13],
        help="Audited window lengths to inspect",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    dataset = classifications_dataset(args.root)
    thresholds = {7: 0.01, 8: 0.003, 9: 0.003, 10: 0.002, 11: 0.002, 12: 0.0, 13: 0.0}

    for length in args.lengths:
        ns, status, primes = audited_slice(dataset, length)
        print(f"=== audited_length={length} rows={len(ns)} ===")
        for line in family_report(ns, status, primes):
            print(line)
        print()
        threshold = thresholds.get(length, 0.01)
        for line in dynamic_proxy_report(ns, status, primes, length, threshold):
            print(line)
        print()

    for line in fresh_run_windows_report(args.root):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
