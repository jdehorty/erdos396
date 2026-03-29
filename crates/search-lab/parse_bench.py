#!/usr/bin/env python3
"""Parse benchmark TSV output into markdown table rows.

Reads lines starting with 'R\t' (the universal benchmark contract) and
formats throughput results into a markdown table.

R-line fields: R, k, witness, elapsed_secs, candidates, speed_mcps[, witnesses_found]
"""
import sys

# Known witness positions (start of each k's search range for benchmarking)
WITNESS = {
    8:  339_949_252,
    9:  1_019_547_844,
    10: 17_609_764_994,
    11: 1_070_858_041_585,
    12: 5_048_891_644_646,
    13: 18_253_129_921_842,
}

for line in open(sys.argv[1]):
    line = line.strip()

    # --bench mode outputs markdown rows directly: pass them through
    if line.startswith("| ") and not line.startswith("| k"):
        print(line)
        # Extract k and speed for stderr progress
        parts = [p.strip().strip('`') for p in line.split("|") if p.strip()]
        if len(parts) >= 6:
            sys.stderr.write("  k=%-4s  %ss  %s M/s\n" % (parts[0], parts[4], parts[5]))
        continue

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

    if k in WITNESS:
        s = WITNESS[k]
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
