#!/usr/bin/env python3
"""Parse Rust/C++ benchmark output into markdown table rows."""
import re, sys

RANGES = {
    8:  (1,              339949253),
    9:  (339949252,      1019547845),
    10: (1019547844,     17609764995),
    13: (18243129921842, 18253129921843),
}

for line in open(sys.argv[1]):
    m = re.match(r'k=\s*(\d+)\s+n=\s*\d+\s+([\d.]+)s\s+([\d.]+)\s+M/s', line)
    if m:
        k, elapsed, speed = int(m.group(1)), m.group(2), m.group(3)
        if k in RANGES:
            s, e = RANGES[k]
            r = e - s
            print("| %d | `[%d, %d)` | %s | %s | %s |" % (k, s, e, format(r, ","), elapsed, speed))
            sys.stderr.write("  k=%-4d  %ss  %s M/s\n" % (k, elapsed, speed))
