#!/usr/bin/env python3
"""Aggregate WaSiM hourly station/time series files to daily values.

Each output row is stamped YYYY MM DD 24 (end of the day). Values are the mean
of the valid hourly values of that calendar day (as stamped in the input), or
the sum for files given with --sum (precipitation). Values <= -999 or >= 998
(WaSiM's usual nodata markers) are ignored; a day without valid values gets -9999.

usage: hourly_to_daily.py OUT_DIR [--sum FILE ...] [--mean FILE ...]
"""
import argparse
import os


def aggregate(path, out_dir, how):
    with open(path, encoding="latin-1") as f:
        lines = f.read().replace("\r", "").rstrip("\n").split("\n")
    header, rows = [], []
    for line in lines:
        fields = line.split()
        if len(fields) >= 4 and all(x.lstrip("-").isdigit() for x in fields[:4]):
            rows.append(fields)
        elif not rows:
            header.append(line)
    out = list(header)
    day, acc = None, None
    def flush():
        vals = []
        for s, n in acc:
            if n == 0:
                vals.append("-9999")
            else:
                vals.append(f"{(s if how == 'sum' else s / n):.4g}")
        out.append("\t".join([day[0], day[1], day[2], "24"] + vals))
    for fields in rows:
        key = (fields[0], fields[1].zfill(2), fields[2].zfill(2))
        if key != day:
            if day is not None:
                flush()
            day, acc = key, [[0.0, 0] for _ in fields[4:]]
        for k, v in enumerate(fields[4:]):
            x = float(v)
            if -999 < x < 998:
                acc[k][0] += x
                acc[k][1] += 1
    if day is not None:
        flush()
    with open(os.path.join(out_dir, os.path.basename(path)), "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"{os.path.basename(path)}: {len(rows)} hourly rows -> {len(out) - len(header)} days ({how})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--sum", nargs="*", default=[])
    ap.add_argument("--mean", nargs="*", default=[])
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    for p in a.sum:
        aggregate(p, a.out_dir, "sum")
    for p in a.mean:
        aggregate(p, a.out_dir, "mean")
