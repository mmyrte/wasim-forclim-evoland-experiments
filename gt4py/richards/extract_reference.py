"""Extract unchanged WaSiM helper definitions into a minimal, isolated C ABI.

Usage: python extract_reference.py /path/to/wasim-eth/src/unsatzon.cc /scratch/reference.cc
Compile the output with clang++ -O3 -ffp-contract=off -shared -fPIC.
The adapter supplies new_thetas' only class state: levels and GWLEVEL[i][j].
No sample data or WaSiM source is stored alongside this extraction script.
"""

import argparse
import hashlib
import json
import pathlib


def extract(text, signature):
    lines = text.splitlines(keepends=True)
    start = next(n for n, line in enumerate(lines) if line.startswith(signature))
    end = next(n for n in range(start + 1, len(lines)) if lines[n].startswith("}"))
    return "".join(lines[start : end + 1]), start + 1, end + 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    args = parser.parse_args()
    text = args.source.read_text()
    definitions = [
        extract(text, sig)
        for sig in ("double unsatzonmodellklasse::new_thetas(", "double unsatzonmodellklasse::quflow(")
    ]
    header = """#include <cmath>
using std::fabs;
class unsatzonmodellklasse {
public:
    int levels;
    double **GWLEVEL;
    double new_thetas(int,int,int,int,int&,double*&,double*&,double*&,double*&,
        double*&,double*&,double*&,double*&,long&,double);
    double quflow(double,double,double,double,double);
};
"""
    adapter = """
extern "C" double reference(int levels, int offset, int iteration, int maxlevel,
    double gwlevel, long *weight, double dt, double *gwin, double *gwex,
    double *theta, double *average, double *old, double *flux,
    double *interflow, double *thickness) {
    double *row = &gwlevel;
    unsatzonmodellklasse soil;
    soil.levels = levels;
    soil.GWLEVEL = &row;
    return soil.new_thetas(offset, 0, 0, iteration, maxlevel, gwin, gwex,
        theta, average, old, flux, interflow, thickness, *weight, dt);
}
extern "C" double darcy(double k, double h1, double h2, double dz1, double dz2) {
    unsatzonmodellklasse soil;
    return soil.quflow(k,h1,h2,dz1,dz2);
}
"""
    args.output.write_text(header + "\n".join(d[0] for d in definitions) + adapter)
    metadata = {
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "functions": [
            {"first": first, "last": last, "sha256": hashlib.sha256(code.encode()).hexdigest()}
            for code, first, last in definitions
        ],
    }
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
