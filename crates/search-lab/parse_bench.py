#!/usr/bin/env python3
"""Parse Rust benchmark output into markdown table rows.

Computes throughput from elapsed time and the known range (not the binary's
self-reported speed, which is wrong when --end caps the search before finding
a witness).
"""
import re, sys

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
    m = re.match(r'k=\s*(\d+)\s+n=\s*\d+\s+([\d.]+)s', line)
    if m:
        k, elapsed = int(m.group(1)), float(m.group(2))
        if k in RANGES:
            s, e = RANGES[k]
            r = e - s
            speed = r / elapsed / 1e6
            print("| %d | `[%d, %d)` | %s | %.3f | %.1f |" % (k, s, e, format(r, ","), elapsed, speed))
            sys.stderr.write("  k=%-4d  %.3fs  %.1f M/s\n" % (k, elapsed, speed))
