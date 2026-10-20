#!/usr/bin/env python3
"""Richards proxy in gt4py.cartesian: the same lockstep ticks as richards_next.py, as one stencil.

The levels are K; the column sweeps are FORWARD/BACKWARD computations that read the level
above or below. Cells are I (J has size 1). Per-cell values are kept on all levels, as in the
next version, so the two do the same work. A gtscript stencil cannot loop over K inside a
`while`, so the per-cell loops become a loop of stencil calls in Python, as in next.

usage: richards_cartesian.py CASE_DIR OUT_DIR [--mode adaptive|fixed] [--substeps N] [--iters I]
                             [--check C] [--backend gt:cpu_ifirst|gt:cpu_kfirst|dace:cpu] [--steps T]
"""
import argparse
import os
import time

import numpy as np
import gt4py.cartesian.gtscript as gs
from gt4py.cartesian.gtscript import BACKWARD, FORWARD, PARALLEL, ceil, computation, interval, sqrt

FT = gs.Field[np.float64]
PARAMS = ("ks", "ths", "thr", "alpha", "inv_n", "m", "inv_m", "dz", "zc")
STATE = ("th", "old", "avg", "h1", "z", "s", "done", "drain", "excess", "eta", "iters")


def substeps_def(th: FT, rain: FT, ks: FT, ths: FT, thr: FT, m: FT, inv_m: FT, dz: FT, n: FT,
                 dt_step: float, min_step: float):
    with computation(PARALLEL), interval(...):
        s = (th - thr) / (ths - thr)
        if s < 1e-6:
            s = 1e-6
        if s > 1.0:
            s = 1.0
        t = 1.0 - (1.0 - s ** inv_m) ** m
        kr = sqrt(s) * t * t
        if rain > 0.0:
            kr = 1.0
        d = dt_step
        if kr > 0.0:
            d = dz / (kr * ks)
    with computation(PARALLEL), interval(-1, None):
        d = dt_step  # the lowest level does not count (as in WaSiM)
    with computation(FORWARD):
        with interval(0, 1):
            dmin = d
        with interval(1, None):
            dmin = d if d < dmin[0, 0, -1] else dmin[0, 0, -1]
    with computation(BACKWARD):
        with interval(-1, None):
            dtot = dmin
        with interval(0, -1):
            dtot = dtot[0, 0, 1]
    with computation(PARALLEL), interval(...):
        if dtot > dt_step:
            dtot = dt_step
        if dtot < min_step:
            dtot = min_step
        n = ceil(dt_step / dtot)


def tick_def(th: FT, old: FT, avg: FT, h1: FT, z: FT, s: FT, done: FT, drain: FT, excess: FT, eta: FT, iters: FT,
             n: FT, rain: FT, et: FT, root: FT,
             ks: FT, ths: FT, thr: FT, alpha: FT, inv_n: FT, m: FT, inv_m: FT, dz: FT, zc: FT,
             th_o: FT, old_o: FT, avg_o: FT, h1_o: FT, z_o: FT, s_o: FT, done_o: FT, drain_o: FT, excess_o: FT,
             eta_o: FT, iters_o: FT,
             dt_step: float, max_err: float, max_iter: float, root_levels: float):
    with computation(PARALLEL), interval(...):
        dt = dt_step / n
        psub = rain / n
        esub = et / n
        live = done == 0.0
        fresh = live and z == 0.0
        take = 0.0
        if root > 0.0:
            take = esub / root_levels / 1000.0 / dz
        avail = th - thr
        if avail < 0.0:
            avail = 0.0
        if take > avail:
            take = avail
        th0 = th
        eta1 = eta
        old1 = old
        avg1 = avg
        h11 = h1
        if fresh:
            th0 = th - take
            eta1 = eta + take * dz * 1000.0
            old1 = th0
            avg1 = th0
            sa = (th0 - thr) / (ths - thr)
            if sa < 1e-6:
                sa = 1e-6
            if sa > 1.0:
                sa = 1.0
            h11 = (sa ** (-inv_m) - 1.0) ** inv_n / alpha
        sb = (avg1 - thr) / (ths - thr)
        if sb < 1e-6:
            sb = 1e-6
        if sb > 1.0:
            sb = 1.0
        h2 = (sb ** (-inv_m) - 1.0) ** inv_n / alpha
        hh = zc - 0.5 * (h11 + h2)
        tb = 1.0 - (1.0 - sb ** inv_m) ** m
        kk = ks * (sqrt(sb) * tb * tb)
        q0 = psub / (1000.0 * dt)
    with computation(BACKWARD):
        with interval(-1, None):
            qraw = kk
        with interval(0, -1):
            qraw = sqrt(kk * kk[0, 0, 1]) * (hh - hh[0, 0, 1]) / (0.5 * (dz + dz[0, 0, 1]))
    with computation(FORWARD):
        with interval(0, 1):
            qtop = q0
            allowed = qtop + (old1 - thr) * dz / dt
            qbot = allowed if qraw > allowed else qraw
        with interval(1, None):
            qtop = qbot[0, 0, -1]
            allowed = qtop + (old1 - thr) * dz / dt
            qbot = allowed if qraw > allowed else qraw
    with computation(BACKWARD):
        with interval(-1, None):
            qb2 = qbot
            w = old1 + (qtop - qb2) * dt / dz
            qt2 = qbot + (ths - old1) * dz / dt if w > ths else qtop
        with interval(0, -1):
            qb2 = qt2[0, 0, 1]
            w = old1 + (qtop - qb2) * dt / dz
            qt2 = qb2 + (ths - old1) * dz / dt if w > ths else qtop
    with computation(PARALLEL), interval(...):
        z1 = z + 1.0
        weight = z1 - 3.0
        avg2 = avg1
        if z1 <= 3.0:
            avg2 = th0
        thn = old1 + (qt2 - qb2) * dt / dz
        if thn < 0.0:
            thn = 0.0
        elif thn > 1.0:
            thn = 1.0
        e = abs(th0 - thn)
        avg3 = (avg2 * weight + thn) / (weight + 1.0)
        if z1 <= 3.0:
            avg3 = (avg2 + thn) / 2.0
    with computation(FORWARD):
        with interval(0, 1):
            emax = e
            qsurf = qt2
        with interval(1, None):
            emax = e if e > emax[0, 0, -1] else emax[0, 0, -1]
            qsurf = qsurf[0, 0, -1]
    with computation(BACKWARD):
        with interval(-1, None):
            err = emax
            qdeep = qb2
        with interval(0, -1):
            err = err[0, 0, 1]
            qdeep = qdeep[0, 0, 1]
    with computation(PARALLEL), interval(...):
        maxe = max_err
        if z1 >= 21.0:
            maxe = max_err * 5.0
        finished = err <= maxe or z1 >= max_iter
        end = live and finished
        drain_o = drain
        excess_o = excess
        iters_o = iters
        s_o = s
        done_o = done
        th_o = th
        old_o = old
        avg_o = avg
        h1_o = h1
        z_o = z
        eta_o = eta1
        if end:
            drain_o = drain + qdeep * dt * 1000.0
            excess_o = excess + (psub - qsurf * dt * 1000.0)
            iters_o = iters + z1
            s_o = s + 1.0
            if s_o >= n:
                done_o = 1.0
        if live:
            th_o = thn
            old_o = old1
            avg_o = avg3
            h1_o = h11
            z_o = z1
            if end:
                z_o = 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--mode", default="adaptive", choices=["adaptive", "fixed"])
    ap.add_argument("--substeps", type=int, default=1)
    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--check", type=int, default=1)
    ap.add_argument("--backend", default="gt:cpu_ifirst")
    ap.add_argument("--steps", type=int, default=0)
    args = ap.parse_args()

    case = np.load(os.path.join(args.case_dir, "case.npz"))
    params = case["params"]
    ncell, _, nlev = params.shape
    nsteps = args.steps or case["rain"].shape[0]
    shape = (ncell, 1, nlev)
    as3 = lambda a: np.array(np.broadcast_to(np.asarray(a, dtype=np.float64).reshape(ncell, 1, -1), shape), order="C")  # noqa: E731
    p = {k: as3(params[:, i, :]) for i, k in enumerate(PARAMS)}
    root = as3(np.broadcast_to(np.arange(nlev) < int(case["root_levels"]), (ncell, nlev)))
    dt_step, max_err = float(case["dt_step"]), float(case["max_err"])
    max_iter = float(case["max_iter"]) if args.mode == "adaptive" else float(args.iters)
    rl = float(case["root_levels"])
    t0 = time.perf_counter()
    sub_st = gs.stencil(backend=args.backend, definition=substeps_def, name="proxy_substeps")
    tick_st = gs.stencil(backend=args.backend, definition=tick_def, name="proxy_tick")
    t_compile = time.perf_counter() - t0

    zeros = lambda: np.zeros(shape)  # noqa: E731
    th = as3(case["theta"])
    bufs = [{k: zeros() for k in STATE}, {k: zeros() for k in STATE}]
    start = {k: zeros() for k in STATE}
    n = zeros()
    out = {k: np.zeros((nsteps, ncell)) for k in ("drain", "excess", "eta", "substeps", "iters")}
    budget = args.substeps * args.iters if args.mode == "fixed" else None
    ticks_total, t_run, t_first = 0, 0.0, 0.0
    for t in range(nsteps):
        rain, et = as3(np.repeat(case["rain"][t][:, None], nlev, 1)), as3(np.repeat(case["et"][t][:, None], nlev, 1))
        t1 = time.perf_counter()
        for k in STATE:
            start[k][...] = 0.0
        start["th"][...] = th
        if budget is None:
            sub_st(th, rain, p["ks"], p["ths"], p["thr"], p["m"], p["inv_m"], p["dz"], n, dt_step=dt_step,
                   min_step=float(case["min_step"]))
        else:
            n[...] = float(args.substeps)
        cur, ticks, i = start, 0, 0
        while True:
            nxt = bufs[i % 2]
            tick_st(*(cur[k] for k in STATE), n, rain, et, root, *(p[k] for k in PARAMS), *(nxt[k] for k in STATE),
                    dt_step=dt_step, max_err=max_err, max_iter=max_iter, root_levels=rl)
            cur, i, ticks = nxt, i + 1, ticks + 1
            if budget is not None:
                if ticks >= budget:
                    break
            elif ticks % args.check == 0 and cur["done"].min() > 0.0:
                break
        th = cur["th"].copy()
        t_run += time.perf_counter() - t1
        if t == 0:
            t_first = t_run
        ticks_total += ticks
        out["drain"][t], out["excess"][t] = cur["drain"][:, 0, 0], cur["excess"][:, 0, 0]
        out["eta"][t] = cur["eta"][:, 0, :int(case["root_levels"])].sum(axis=1)
        out["substeps"][t], out["iters"][t] = n[:, 0, 0], cur["iters"][:, 0, 0]
    os.makedirs(args.out_dir, exist_ok=True)
    np.ascontiguousarray(th[:, 0, :]).tofile(os.path.join(args.out_dir, "theta.f64"))
    for k, a in out.items():
        a.tofile(os.path.join(args.out_dir, k + ".f64"))
    warm = (t_run - t_first) / max(nsteps - 1, 1)
    print(f"{ncell} cells, {nlev} levels, {nsteps} steps, {args.backend}, mode {args.mode}"
          f"{'' if args.mode == 'adaptive' else f' N={args.substeps} I={args.iters}'}: "
          f"{1000 * warm:.2f} ms per step (warm; compile {t_compile:.0f} s), "
          f"{ticks_total / nsteps:.1f} ticks per step, {1000 * warm / (ticks_total / nsteps):.3f} ms per tick")


if __name__ == "__main__":
    main()
