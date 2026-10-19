#!/usr/bin/env python3
"""Compare the outputs of two Richards proxy runs (e.g. richards_ref.cpp and another formulation).

Prints, per output, the largest absolute difference over all cells (and steps), and for the
fluxes also the difference of the catchment totals, the measure that matters for discharge.

usage: compare.py CASE_DIR REF_OUT OTHER_OUT [--steps T]
"""
import argparse
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("ref")
    ap.add_argument("other")
    ap.add_argument("--steps", type=int, default=0)
    args = ap.parse_args()
    case = np.load(os.path.join(args.case_dir, "case.npz"))
    ncell, _, nlev = case["params"].shape
    load = lambda d, k, shape: np.fromfile(os.path.join(d, k + ".f64")).reshape(shape)  # noqa: E731
    th_r, th_o = load(args.ref, "theta", (ncell, nlev)), load(args.other, "theta", (ncell, nlev))
    nsteps = args.steps or np.fromfile(os.path.join(args.other, "drain.f64")).size // ncell
    print(f"theta after {nsteps} steps: max |diff| {np.abs(th_r - th_o).max():.3g}, "
          f"bit-identical cells {np.mean(np.all(th_r == th_o, axis=1)) * 100:.1f} %")
    for k in ("drain", "excess", "eta"):
        r = load(args.ref, k, (-1, ncell))[:nsteps]
        o = load(args.other, k, (-1, ncell))[:nsteps]
        tot_r, tot_o = r.sum(), o.sum()
        print(f"{k:7s} max |diff| per cell-step {np.abs(r - o).max():.3g} mm; total {tot_r:.6g} vs {tot_o:.6g} mm "
              f"({(tot_o - tot_r) / tot_r * 100 if tot_r else 0.0:+.3g} %)")
    for k in ("substeps", "iters"):
        r = load(args.ref, k, (-1, ncell))[:nsteps]
        o = load(args.other, k, (-1, ncell))[:nsteps]
        print(f"{k:8s} mean {r.mean():.3f} vs {o.mean():.3f}; cell-steps that differ {np.mean(r != o) * 100:.2f} %")


if __name__ == "__main__":
    main()
