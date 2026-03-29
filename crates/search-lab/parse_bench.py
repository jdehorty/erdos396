#!/usr/bin/env python3
"""Parse benchmark TSV output into markdown table rows.

Reads lines starting with 'R\t' (the universal benchmark contract) and
computes throughput from elapsed time and known range sizes.
"""
import sys

# k=9,10: full witness ranges (chains from k=8 witness)
# k=11:   capped at 50B candidates (witness at 1.07T is too far)
# k=13:   same total candidates as k=9+10+11 combined (~67.3B)
RANGES = {
    9:  (339949252,        1019547845),
    10: (1019547844,       17609764995),
    11: (17609764994,      67609764994),
    13: (18185829921842,   18253129921842),
}

for line in open(sys.argv[1]):
    if not line.startswith("R\t"):
        continue
    fields = line.strip().split("\t")
    # Fields: R, k, witness, elapsed_secs, candidates_checked, speed_mcps
    k = int(fields[1])
    elapsed = float(fields[3])
    if k in RANGES:
        s, e = RANGES[k]
        r = e - s
        speed = r / elapsed / 1e6
        print("| %d | `[%d, %d)` | %s | %.4f | %.1f |" % (k, s, e, format(r, ","), elapsed, speed))
        sys.stderr.write("  k=%-4d  %.4fs  %.1f M/s\n" % (k, elapsed, speed))
