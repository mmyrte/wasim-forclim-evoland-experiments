#!/usr/bin/env python3
"""Which lines of a WaSiM function did a reference run execute?

A port can only be validated on the code paths its reference run exercises. This tool reads the
coverage data of a WaSiM build with --coverage and reports, per function, how many executable
lines ran, and with --missed which did not. A reference case (real or synthetic) is good enough
for a port when every line of the ported branch has run.

Build without OpenMP: gcov does not attribute the outlined OpenMP loop bodies to source lines.

  cd wasim-eth/src
  make MPI=0 DEBUG=1 ARCH=--coverage OMPFLAGS= OMPLIBS= CXXTAG=cov TARGET=wasim-cov
  # run the reference case with ./wasim-cov (each run adds to the counters in build/*-cov/*.gcda)
  python3 validation/wasim_coverage.py wasim-eth/src wasim-eth/src/build/omp-bundled-cov \\
      interzep.cc:interzeptionsmodellklasse::run etpetr.cc:verdunstungsmodellklasse::snow --missed

Delete the .gcda files to start counting afresh.

usage: wasim_coverage.py SRC_DIR OBJECT_DIR FILE:CLASS::FUNCTION [...] [--missed] [--range FILE:FIRST-LAST ...]
"""
import argparse
import os
import re
import subprocess


def gcov_lines(src_dir, obj_dir, source):
    """{line number: count or None (not executable)} for one source file."""
    out = subprocess.run(["gcov", "-t", "-o", obj_dir, source], cwd=src_dir, capture_output=True, text=True).stdout
    lines, keep = {}, False
    for line in out.splitlines():
        m = re.match(r"\s*([^:]+):\s*(\d+):(.*)", line)
        if not m:
            continue
        count, number, text = m.group(1).strip(), int(m.group(2)), m.group(3)
        if number == 0:
            if text.startswith("Source:"):
                keep = text[len("Source:"):] == source
            continue
        if keep:
            lines[number] = (None if count == "-" else 0 if count in ("#####", "=====") else int(count.rstrip("*")), text)
    return lines


def function_range(src_dir, source, name):
    """First and last line of a function definition 'Class::function' (closing brace in column 0)."""
    text = open(os.path.join(src_dir, source), encoding="utf-8", errors="replace").read().splitlines()
    pattern = re.compile(r"^\S.*\b" + re.escape(name) + r"\s*\(")
    for i, line in enumerate(text):
        if pattern.match(line):
            for j in range(i + 1, len(text)):
                if text[j].startswith("}"):
                    return i + 1, j + 1
    raise SystemExit(f"{name} not found in {source}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src_dir")
    ap.add_argument("obj_dir")
    ap.add_argument("functions", nargs="*", help="FILE:Class::function")
    ap.add_argument("--range", action="append", default=[], help="FILE:FIRST-LAST")
    ap.add_argument("--missed", action="store_true", help="list the executable lines that did not run")
    args = ap.parse_args()
    targets = []
    for spec in args.functions:
        source, name = spec.split(":", 1)
        targets.append((source, name, *function_range(args.src_dir, source, name)))
    for spec in args.range:
        source, span = spec.split(":", 1)
        first, last = map(int, span.split("-"))
        targets.append((source, f"lines {first}-{last}", first, last))
    cache = {}
    for source, name, first, last in targets:
        if source not in cache:
            cache[source] = gcov_lines(args.src_dir, args.obj_dir, source)
        lines = {n: v for n, v in cache[source].items() if first <= n <= last and v[0] is not None}
        ran = sum(1 for c, _ in lines.values() if c > 0)
        print(f"{source}:{first}-{last} {name}: {ran} of {len(lines)} executable lines ran"
              f" ({100 * ran / max(len(lines), 1):.1f}%)")
        if args.missed:
            for n, (c, text) in sorted(lines.items()):
                if c == 0:
                    print(f"    {n:6d}  {text.strip()[:110]}")


if __name__ == "__main__":
    main()
