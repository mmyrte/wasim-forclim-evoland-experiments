#!/usr/bin/env python3
"""Create temperature_merged.dat and precipitation_merged.dat for the Rhone golden case.

The Rhone sample reads gridded RhiresD/TabsD inputs that are not part of the sample.
The golden case uses station regression instead, from the two regional station files
in the sample (<var>_reg1.dat, <var>_reg2.dat), merged column-wise into one file per
variable. Both files have the same time stamps.

usage: merge_station_files.py <unpacked sample-run-rhone dir> [output dir, default: this script's dir]
"""
import os
import sys


def merge(ts_dir, var, out_dir):
    def read(name):
        with open(os.path.join(ts_dir, name), encoding="latin-1") as f:
            return f.read().replace("\r", "").rstrip("\n").split("\n")

    a, b = read(f"{var}_reg1.dat"), read(f"{var}_reg2.dat")
    if len(a) != len(b):
        sys.exit(f"{var}: line counts differ ({len(a)} vs {len(b)})")
    out = [var.upper()]
    for la, lb in zip(a[1:], b[1:]):
        fa, fb = la.rstrip("\t").split("\t"), lb.rstrip("\t").split("\t")
        if fa[0] != "YY" and fa[:4] != fb[:4]:
            sys.exit(f"{var}: time stamps differ: {fa[:4]} vs {fb[:4]}")
        out.append("\t".join(fa + fb[4:]))
    with open(os.path.join(out_dir, f"{var}_merged.dat"), "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    case = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))
    for var in ("temperature", "precipitation"):
        merge(os.path.join(case, "Input", "timeseries"), var, out_dir)
