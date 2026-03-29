#!/usr/bin/env python3
"""Parse C++ reference benchmark output into markdown table rows.

Computes throughput from the binary's self-reported elapsed time and the known
range. For capped runs (k=11, k=13) the binary's reported speed is based on
witness distance, which is wrong when no witness is found.
"""
import re, sys

RANGES = {
    9:  (339949252,        1019547845),
    10: (1019547844,       17609764995),
    11: (17609764994,      67609764994),
    13: (18185829921842,   18253129921842),
}

# C++ output: k =  9 | min n =  1019547844 | Time:  0.2040 s | Speed:  3330.88 M candidates/s
for line in open(sys.argv[1]):
    m = re.match(r'k\s*=\s*(\d+)\s*\|.*Time:\s*([\d.]+)\s*s', line)
    if m:
        k, elapsed = int(m.group(1)), float(m.group(2))
        if k in RANGES:
            s, e = RANGES[k]
            r = e - s
            speed = r / elapsed / 1e6
            print("| %d | `[%d, %d)` | %s | %.4f | %.1f |" % (k, s, e, format(r, ","), elapsed, speed))
            sys.stderr.write("  k=%-4d  %.4fs  %.1f M/s\n" % (k, elapsed, speed))
