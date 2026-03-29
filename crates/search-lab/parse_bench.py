#!/usr/bin/env python3
"""Parse benchmark TSV output into markdown table rows.

Reads lines starting with 'R\t' (the universal benchmark contract) and
formats throughput results into a markdown table.

R-line fields: R, k, witness, elapsed_secs, candidates, speed_mcps[, witnesses_found]

Usage: parse_bench.py FILE   or   ... | parse_bench.py -
"""
import sys

# Sweep start positions (must match justfile STARTS array)
STARTS = {
    8:  169_974_626,
    9:  509_773_922,
    10: 8_804_882_497,
    11: 1_062_858_041_585,
    12: 5_042_391_644_646,
    13: 18_247_629_921_842,
}

_input = sys.stdin if len(sys.argv) < 2 or sys.argv[1] == "-" else open(sys.argv[1])

for line in _input:
    line = line.strip()

    # R-line TSV format (universal benchmark contract)
    if not line.startswith("R\t"):
        continue
    fields = line.split("\t")
    # Fields: R, k, witness, elapsed_secs, candidates, speed_mcps[, witnesses_found]
    k = int(fields[1])
    elapsed = float(fields[3])
    candidates = int(fields[4])
    speed = float(fields[5])
    witnesses = int(fields[6]) if len(fields) > 6 else 0

    if k in STARTS:
        s = STARTS[k]
        e = s + candidates
        print("| %d | `[%s, %s)` | %s | %d | %.4f | %.1f |" % (
            k,
            format(s, ","),
            format(e, ","),
            format(candidates, ","),
            witnesses,
            elapsed,
            speed,
        ))
        sys.stderr.write("  k=%-4d  %.4fs  %.1f M/s  witnesses=%d\n" % (k, elapsed, speed, witnesses))
