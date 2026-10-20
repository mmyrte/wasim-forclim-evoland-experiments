#!/usr/bin/env python3
"""Richards proxy written directly in DaCe: a parallel map over cells, WaSiM's loops inside.

Same algorithm and loop structure as richards_ref.cpp (adaptive sub-steps per cell, iterate to
convergence), expressed as a DaCe program: `dace.map` over the cells, ordinary Python loops
(`for`, `while`) over sub-steps, iterations and levels inside it. DaCe turns the loops into
control flow inside each map iteration, so the per-cell adaptivity is kept and nothing returns
to Python during a step.

usage: richards_dace.py CASE_DIR OUT_DIR [--steps T] [--loop]   (threads: OMP_NUM_THREADS)
   --loop: also run the time loop inside the DaCe program (one call for all steps)
"""
import argparse
import math
import os
import time

import dace
import numpy as np

N = dace.symbol("N")
L = dace.symbol("L")
T = dace.symbol("T")
KS, THS, THR, ALPHA, INV_N, M, INV_M, DZ, ZC = range(9)


@dace.program
def column_step(p: dace.float64[N, 9, L], th: dace.float64[N, L], rain: dace.float64[N], et: dace.float64[N],
                drain: dace.float64[N], excess: dace.float64[N], eta: dace.float64[N], nsub: dace.float64[N],
                iters: dace.float64[N], DT: dace.float64, max_err: dace.float64, max_iter: dace.int64,
                min_step: dace.float64, root_levels: dace.int64):
    for c in dace.map[0:N]:
        old = np.empty(L, dtype=np.float64)
        avg = np.empty(L, dtype=np.float64)
        h1 = np.empty(L, dtype=np.float64)
        hh = np.empty(L, dtype=np.float64)
        kk = np.empty(L, dtype=np.float64)
        q = np.empty(L + 1, dtype=np.float64)
        dz = np.empty(L, dtype=np.float64)
        for k in range(L):
            dz[k] = p[c, DZ, k]  # local copy: DaCe 2.0.0a10 drops the +1 in p[c, DZ, k + 1] (dace/known-issues/)
        P = rain[c]
        E = et[c]
        dtmin = DT * 1.0
        for k in range(L - 1):
            kr = 1.0
            if P <= 0:
                s = (th[c, k] - p[c, THR, k]) / (p[c, THS, k] - p[c, THR, k])
                if s < 1e-6:
                    s = 1e-6
                if s > 1.0:
                    s = 1.0
                t = 1.0 - ((1.0 - ((s) ** (p[c, INV_M, k]))) ** (p[c, M, k]))
                kr = math.sqrt(s) * t * t
            if kr > 0:
                d = p[c, DZ, k] / (kr * p[c, KS, k])
                if d < dtmin:
                    dtmin = d
        if dtmin < min_step:
            dtmin = min_step
        n = int(math.ceil(DT / dtmin))
        dt = DT / n
        Psub = P / n
        Esub = E / n
        sum_drain = 0.0
        sum_excess = 0.0
        sum_eta = 0.0
        sum_iter = 0
        for sub in range(n):
            for k in range(root_levels):
                take = Esub / root_levels / 1000.0 / p[c, DZ, k]
                avail = th[c, k] - p[c, THR, k]
                if avail < 0:
                    avail = 0.0
                if take > avail:
                    take = avail
                th[c, k] = th[c, k] - take
                sum_eta = sum_eta + take * p[c, DZ, k] * 1000.0
            for k in range(L):
                old[k] = th[c, k]
                avg[k] = th[c, k]
                s = (old[k] - p[c, THR, k]) / (p[c, THS, k] - p[c, THR, k])
                if s < 1e-6:
                    s = 1e-6
                if s > 1.0:
                    s = 1.0
                h1[k] = ((((s) ** (-p[c, INV_M, k])) - 1.0) ** (p[c, INV_N, k])) / p[c, ALPHA, k]
            z = 0
            maxe = max_err * 1.0
            err = 1.0
            while err > maxe and z < max_iter:
                for k in range(L):
                    s = (avg[k] - p[c, THR, k]) / (p[c, THS, k] - p[c, THR, k])
                    if s < 1e-6:
                        s = 1e-6
                    if s > 1.0:
                        s = 1.0
                    h2 = ((((s) ** (-p[c, INV_M, k])) - 1.0) ** (p[c, INV_N, k])) / p[c, ALPHA, k]
                    hh[k] = p[c, ZC, k] - 0.5 * (h1[k] + h2)
                    t = 1.0 - ((1.0 - ((s) ** (p[c, INV_M, k]))) ** (p[c, M, k]))
                    kk[k] = p[c, KS, k] * (math.sqrt(s) * t * t)
                q[0] = Psub / (1000.0 * dt)
                for k in range(L - 1):
                    q[k + 1] = math.sqrt(kk[k] * kk[k + 1]) * (hh[k] - hh[k + 1]) / (0.5 * (dz[k] + dz[k + 1]))
                q[L] = kk[L - 1]
                for k in range(L):
                    allowed = q[k] + (old[k] - p[c, THR, k]) * p[c, DZ, k] / dt
                    if q[k + 1] > allowed:
                        q[k + 1] = allowed
                for kb in range(L):
                    k = L - 1 - kb
                    w = old[k] + (q[k] - q[k + 1]) * dt / p[c, DZ, k]
                    if w > p[c, THS, k]:
                        q[k] = q[k + 1] + (p[c, THS, k] - old[k]) * p[c, DZ, k] / dt
                z = z + 1
                if z == 21:
                    maxe = max_err * 5.0
                weight = z - 3
                err = 0.0
                for k in range(L):
                    last = th[c, k]
                    if z <= 3:
                        avg[k] = th[c, k]
                    v = old[k] + (q[k] - q[k + 1]) * dt / p[c, DZ, k]
                    if v < 0.0:
                        v = 0.0
                    elif v > 1.0:
                        v = 1.0
                    th[c, k] = v
                    e = abs(last - v)
                    if e > err:
                        err = e
                    if z > 3:
                        avg[k] = (avg[k] * weight + v) / (weight + 1)
                    else:
                        avg[k] = (avg[k] + v) / 2.0
            sum_excess = sum_excess + (Psub - q[0] * dt * 1000.0)
            sum_drain = sum_drain + q[L] * dt * 1000.0
            sum_iter = sum_iter + z
        drain[c] = sum_drain
        excess[c] = sum_excess
        eta[c] = sum_eta
        nsub[c] = n
        iters[c] = sum_iter


@dace.program
def all_steps(p: dace.float64[N, 9, L], th: dace.float64[N, L], rain: dace.float64[T, N], et: dace.float64[T, N],
              drain: dace.float64[T, N], excess: dace.float64[T, N], eta: dace.float64[T, N], nsub: dace.float64[T, N],
              iters: dace.float64[T, N], DT: dace.float64, max_err: dace.float64, max_iter: dace.int64,
              min_step: dace.float64, root_levels: dace.int64):
    for t in range(T):
        column_step(p, th, rain[t], et[t], drain[t], excess[t], eta[t], nsub[t], iters[t], DT, max_err, max_iter,
                    min_step, root_levels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--steps", type=int, default=0)
    ap.add_argument("--loop", action="store_true")
    args = ap.parse_args()
    case = np.load(os.path.join(args.case_dir, "case.npz"))
    p = np.array(case["params"], order="C")
    th = np.ascontiguousarray(case["theta"]).copy()
    nsteps = args.steps or case["rain"].shape[0]
    ncell, _, nlev = p.shape
    rain, et = np.array(case["rain"][:nsteps], order="C"), np.array(case["et"][:nsteps], order="C")
    out = {k: np.zeros((nsteps, ncell)) for k in ("drain", "excess", "eta", "substeps", "iters")}
    scal = dict(DT=float(case["dt_step"]), max_err=float(case["max_err"]), max_iter=int(case["max_iter"]),
                min_step=float(case["min_step"]), root_levels=int(case["root_levels"]))
    t0 = time.perf_counter()
    if args.loop:
        sdfg = all_steps.to_sdfg()
        sdfg.specialize(dict(L=nlev))
        prog = sdfg.compile()
        t_compile = time.perf_counter() - t0
        t0 = time.perf_counter()
        prog(p=p, th=th, rain=rain, et=et, drain=out["drain"], excess=out["excess"], eta=out["eta"],
             nsub=out["substeps"], iters=out["iters"], N=ncell, T=nsteps, **scal)
        ms = 1000 * (time.perf_counter() - t0) / nsteps
    else:
        sdfg = column_step.to_sdfg()
        sdfg.specialize(dict(L=nlev))  # fixed-size local arrays (registers/stack, not heap)
        prog = sdfg.compile()
        t_compile = time.perf_counter() - t0
        res = {k: np.zeros(ncell) for k in out}
        rain_t, et_t = np.zeros(ncell), np.zeros(ncell)
        t0 = time.perf_counter()
        for t in range(nsteps):
            rain_t[:], et_t[:] = rain[t], et[t]
            prog(p=p, th=th, rain=rain_t, et=et_t, drain=res["drain"], excess=res["excess"], eta=res["eta"],
                 nsub=res["substeps"], iters=res["iters"], N=ncell, **scal)
            for k in out:
                out[k][t] = res[k]
        ms = 1000 * (time.perf_counter() - t0) / nsteps
    os.makedirs(args.out_dir, exist_ok=True)
    th.tofile(os.path.join(args.out_dir, "theta.f64"))
    for k, a in out.items():
        a.tofile(os.path.join(args.out_dir, k + ".f64"))
    print(f"{ncell} cells, {nlev} levels, {nsteps} steps, dace{' (time loop inside)' if args.loop else ''}: "
          f"{ms:.2f} ms per step (compiled in {t_compile:.0f} s)")


if __name__ == "__main__":
    main()
