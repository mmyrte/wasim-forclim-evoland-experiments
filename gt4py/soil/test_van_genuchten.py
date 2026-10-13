#!/usr/bin/env python3
"""Compare wasim-gt's soil hydraulic functions (wasim_gt.soil.van_genuchten) with WaSiM's C++ originals on random soils and inputs.

Builds the C++ reference from the WaSiM sources (extract_reference.sh), evaluates h_calc, k_calc
and get_theta for N random soil horizons with one random query each, and reports the largest
differences. Includes edge cases: table suctions, h <= 0, theta at or beyond thetar and thetas,
frozen soil (SE < 1), converted thetar/thetas.

usage: test_van_genuchten.py WASIM_SRC_DIR [--n 20000] [--backend gt:cpu_ifirst]
"""
import argparse
import os
import subprocess
import tempfile
import time

import numpy as np
import gt4py.cartesian.gtscript as gs
import gt4py.storage as gts
from gt4py.cartesian.gtscript import PARALLEL, computation, interval

from wasim_gt.soil.van_genuchten import H_TAB, get_theta, h_calc, k_calc

HERE = os.path.dirname(os.path.abspath(__file__))


def build_stencil(backend):
    @gs.stencil(backend=backend)
    def vg(htab: gs.Field[gs.IJ, (np.float64, (31,))],
           ts: gs.Field[np.float64], tr: gs.Field[np.float64], alpha: gs.Field[np.float64],
           n: gs.Field[np.float64], m: gs.Field[np.float64], tau: gs.Field[np.float64],
           theta: gs.Field[np.float64], se: gs.Field[np.float64], h: gs.Field[np.float64],
           trx: gs.Field[np.float64], tsx: gs.Field[np.float64],
           out_h: gs.Field[np.float64], out_k: gs.Field[np.float64], out_t: gs.Field[np.float64]):
        with computation(PARALLEL), interval(...):
            out_h = h_calc(htab, theta, ts, tr, alpha, n, m, trx, tsx)
            out_k = k_calc(htab, theta, se, ts, tr, alpha, n, m, tau, trx, tsx)
            out_t = get_theta(htab, h, ts, tr, alpha, n, m, trx, tsx)
    return vg


def random_cases(count, rng):
    ts = rng.uniform(0.3, 0.6, count)
    tr = rng.uniform(0.01, 0.15, count)
    alpha = 10 ** rng.uniform(-0.5, 1.3, count)   # 1/m
    n = rng.uniform(1.05, 3.0, count)
    tau = np.where(rng.random(count) < 0.5, 0.5, rng.uniform(-1.0, 1.0, count))
    theta = rng.uniform(tr - 0.02, ts + 0.02, count)
    edge = rng.random(count)
    theta = np.where(edge < 0.05, tr, np.where(edge < 0.10, ts, theta))
    h = 10 ** rng.uniform(-5, 2.5, count)
    pick = rng.random(count)
    h = np.where(pick < 0.10, H_TAB[rng.integers(1, 31, count)], h)
    h = np.where((pick >= 0.10) & (pick < 0.13), 0.0, h)
    h = np.where((pick >= 0.13) & (pick < 0.15), -0.1, h)
    se = np.where(rng.random(count) < 0.6, 1.0, rng.uniform(0.0, 1.0, count))
    conv = rng.random(count)
    trx = np.where(conv < 0.5, -1.0, np.where(conv < 0.7, tr, tr * rng.uniform(0.5, 1.5, count)))
    tsx = np.where(conv < 0.5, -1.0, np.where(conv < 0.7, ts, ts * rng.uniform(0.9, 1.1, count)))
    return dict(ts=ts, tr=tr, alpha=alpha, n=n, m=1.0 - 1.0 / n, tau=tau, theta=theta, h=h, se=se, trx=trx, tsx=tsx)


def reference(src, c, workdir):
    cc = os.path.join(workdir, "reference.cc")
    exe = os.path.join(workdir, "reference")
    with open(cc, "w") as f:
        subprocess.run([os.path.join(HERE, "extract_reference.sh"), src], stdout=f, check=True)
    # WaSiM's Linux builds do not contract to FMA (no -march), so neither does the reference
    subprocess.run(["g++", "-O2", "-ffp-contract=off", "-o", exe, cc], check=True)
    x = lambda v: float(v).hex()  # exact transfer of the inputs
    lines = [f"{x(c['ts'][i])} {x(c['tr'][i])} {x(c['alpha'][i])} {x(c['n'][i])} {x(c['tau'][i])}" for i in range(len(c["ts"]))]
    lines.append("END")
    for i in range(len(c["ts"])):
        extra = f"{x(c['se'][i])} {x(c['trx'][i])} {x(c['tsx'][i])}"
        lines += [f"{i + 1} h {x(c['theta'][i])} 1.0 {x(c['trx'][i])} {x(c['tsx'][i])}",
                  f"{i + 1} k {x(c['theta'][i])} {extra}",
                  f"{i + 1} t {x(c['h'][i])} 1.0 {x(c['trx'][i])} {x(c['tsx'][i])}"]
    out = subprocess.run([exe], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True).stdout.split()
    r = np.array(out, dtype="f8").reshape(-1, 3)
    return r[:, 0], r[:, 1], r[:, 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--backend", default="gt:cpu_ifirst")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    c = random_cases(args.n, rng)
    with tempfile.TemporaryDirectory() as tmp:
        ref_h, ref_k, ref_t = reference(args.src, c, tmp)

    shape = (args.n, 1, 1)
    st = lambda a: gts.from_array(np.asarray(a, dtype="f8").reshape(shape), backend=args.backend, aligned_index=(0, 0, 0))
    htab = gts.zeros((args.n, 1), (np.float64, (31,)), backend=args.backend, aligned_index=(0, 0))
    htab[...] = H_TAB[None, None, :]
    outs = [st(np.zeros(args.n)) for _ in range(3)]
    t0 = time.perf_counter()
    vg = build_stencil(args.backend)
    t1 = time.perf_counter()
    vg(htab, *(st(c[k]) for k in ("ts", "tr", "alpha", "n", "m", "tau", "theta", "se", "h", "trx", "tsx")), *outs)
    t2 = time.perf_counter()
    mine = [np.asarray(o).ravel() for o in outs]
    print(f"{args.n} random soil horizons and queries, backend {args.backend} (build {t1 - t0:.1f} s, run {1000 * (t2 - t1):.1f} ms)")
    worst = 0.0
    for name, a, b in zip(("h_calc", "k_calc", "get_theta"), mine, (ref_h, ref_k, ref_t)):
        same = np.sum(a == b)
        rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
        i = int(np.argmax(rel))
        worst = max(worst, rel[i])
        print(f"  {name:9s}: {same} of {args.n} bit-identical, largest relative difference {rel[i]:.3g}"
              + (f" (case {i}: {a[i]!r} vs {b[i]!r})" if rel[i] > 0 else ""))
    return 0 if worst == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
