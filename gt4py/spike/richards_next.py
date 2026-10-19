#!/usr/bin/env python3
"""Richards proxy in gt4py.next: lockstep iterations ("ticks") with masks, the column sweeps as scans.

Field operators have no loops, so the per-cell loops of WaSiM (sub-steps, iterations until
convergence) become a loop of ticks around one compiled program. In each tick every cell that
is not done does one iteration of its current sub-step; a cell that converges finishes the
sub-step and starts its next one in the next tick. The levels are a vertical dimension K; the
sweeps over a column are scan operators. Per-cell values (sub-step count, error, fluxes at the
top and bottom) are kept on all levels of the column (broadcast), since next has no reduction
over K.

Modes:
  adaptive       Courant sub-steps per cell, iterate to convergence (max_iter); after every
                 `--check` ticks the host asks whether all cells are done. Same algorithm and
                 same results as richards_ref.cpp.
  fixed          N sub-steps for all cells (--substeps) and at most I iterations (--iters):
                 exactly N*I ticks per step, no host check. Cells that converge earlier wait.

usage: richards_next.py CASE_DIR OUT_DIR [--mode adaptive|fixed] [--substeps N] [--iters I]
                        [--check C] [--backend gtfn_cpu|embedded|dace_cpu] [--steps T]
"""
import argparse
import os
import time

import numpy as np
import gt4py.next as gtx
from gt4py.next import abs, ceil, float64, maximum, minimum, power, sqrt, where  # noqa: F401

C = gtx.Dimension("C")
K = gtx.Dimension("K", kind=gtx.DimensionKind.VERTICAL)
F = gtx.Field[gtx.Dims[C, K], float64]
PARAMS = ("ks", "ths", "thr", "alpha", "inv_n", "m", "inv_m", "dz", "zc")


@gtx.field_operator
def saturation(th: F, thr: F, ths: F) -> F:
    s = (th - thr) / (ths - thr)
    s = where(s < 1e-6, 1e-6, s)
    return where(s > 1.0, 1.0, s)


@gtx.field_operator
def suction(th: F, thr: F, ths: F, alpha: F, inv_n: F, inv_m: F) -> F:
    s = saturation(th, thr, ths)
    return power(power(s, -inv_m) - 1.0, inv_n) / alpha


@gtx.field_operator
def krel(th: F, thr: F, ths: F, m: F, inv_m: F) -> F:
    s = saturation(th, thr, ths)
    t = 1.0 - power(1.0 - power(s, inv_m), m)
    return sqrt(s) * t * t


# --- column sweeps -------------------------------------------------------------------------------

@gtx.scan_operator(axis=K, forward=False, init=(0.0, 0.0, 0.0, 0.0))
def darcy(state: tuple[float64, float64, float64, float64], hh: float64, kk: float64, dz: float64
          ) -> tuple[float64, float64, float64, float64]:
    """Flux at the bottom of each level (down positive); free drainage below the lowest level."""
    _, hh_below, k_below, dz_below = state
    q = kk if dz_below == 0.0 else sqrt(kk * k_below) * (hh - hh_below) / (0.5 * (dz + dz_below))
    return q, hh, kk, dz


@gtx.scan_operator(axis=K, forward=True, init=(0.0, 0.0, 0.0))
def drying(state: tuple[float64, float64, float64], q0: float64, qbot: float64, old: float64, thr: float64,
           dz: float64, dt: float64) -> tuple[float64, float64, float64]:
    """Top-down: no level below its residual water content. Returns (q at top, q at bottom)."""
    _, qbot_above, started = state
    qtop = qbot_above if started > 0.0 else q0
    allowed = qtop + (old - thr) * dz / dt
    return qtop, (allowed if qbot > allowed else qbot), 1.0


@gtx.scan_operator(axis=K, forward=False, init=(0.0, 0.0, 0.0))
def saturating(state: tuple[float64, float64, float64], qtop: float64, qbot_in: float64, old: float64,
               ths: float64, dz: float64, dt: float64) -> tuple[float64, float64, float64]:
    """Bottom-up: no level above saturation."""
    qtop_below, _, started = state
    qbot = qtop_below if started > 0.0 else qbot_in
    w = old + (qtop - qbot) * dt / dz
    return (qbot + (ths - old) * dz / dt if w > ths else qtop), qbot, 1.0


@gtx.scan_operator(axis=K, forward=True, init=(0.0, 0.0, 0.0))
def down_max(state: tuple[float64, float64, float64], e: float64, qtop: float64
             ) -> tuple[float64, float64, float64]:
    """Running maximum of e downwards, and the flux at the top of the column on every level."""
    emax, q0, started = state
    return (e if e > emax else emax), (q0 if started > 0.0 else qtop), 1.0


@gtx.scan_operator(axis=K, forward=False, init=(0.0, 0.0, 0.0))
def up_bottom(state: tuple[float64, float64, float64], emax: float64, qbot: float64
              ) -> tuple[float64, float64, float64]:
    """The column's maximum and the flux at its bottom, on every level."""
    e, qn, started = state
    return (e if started > 0.0 else emax), (qn if started > 0.0 else qbot), 1.0


@gtx.scan_operator(axis=K, forward=True, init=0.0)
def down_min(state: float64, x: float64) -> float64:
    return x if ((state == 0.0) | (x < state)) else state


@gtx.scan_operator(axis=K, forward=False, init=(0.0, 0.0))
def up_last(state: tuple[float64, float64], x: float64) -> tuple[float64, float64]:
    v, started = state
    return (v if started > 0.0 else x), 1.0


# --- one step's start: Courant sub-steps ---------------------------------------------------------
# (gtfn needs each scan as its own stencil: a scan cannot take another scan's output inside one
# field operator, so the steps below are separate field operators called from a program)

@gtx.field_operator
def courant_dt(th: F, rain: F, ks: F, ths: F, thr: F, m: F, inv_m: F, dz: F, below: F, dt_step: float64) -> F:
    """Per level the Courant time step dz / K (the lowest level does not count, as in WaSiM)."""
    kr = where(rain > 0.0, 1.0, krel(th, thr, ths, m, inv_m))
    return where((kr > 0.0) & (below > 0.0), dz / (kr * ks), dt_step)


@gtx.field_operator
def substep_count(dmin: F, dt_step: float64, min_step: float64) -> F:
    d = minimum(dmin, dt_step)
    d = where(d < min_step, min_step, d)
    return ceil(dt_step / d)


@gtx.program
def substeps(th: F, rain: F, ks: F, ths: F, thr: F, m: F, inv_m: F, dz: F, below: F, dt_step: float64,
             min_step: float64, t_d: F, t_run: F, t_min: F, t_flag: F, n: F):
    """Per cell (on all levels) the number of sub-steps by WaSiM's Courant criterion."""
    courant_dt(th, rain, ks, ths, thr, m, inv_m, dz, below, dt_step, out=t_d)
    down_min(t_d, out=t_run)
    up_last(t_run, out=(t_min, t_flag))
    substep_count(t_min, dt_step, min_step, out=n)


# --- one tick: three field operators in one program ---------------------------------------------
# (one large field operator for the whole tick did not finish compiling in 15 minutes: gt4py's
# type inference on the inlined tree; three stages compile in about a minute)

@gtx.field_operator
def prepare(th: F, old: F, avg: F, h1: F, z: F, done: F, eta: F, n: F, rain: F, et: F, root: F,
            ths: F, thr: F, alpha: F, inv_n: F, inv_m: F, dz: F, dt_step: float64, root_levels: float64,
            ) -> tuple[F, F, F, F, F, F]:
    """Start of a sub-step where z == 0: transpiration from the root levels, old = avg = theta."""
    esub = et / n
    fresh = (done == 0.0) & (z == 0.0)
    take = where(root > 0.0, esub / root_levels / 1000.0 / dz, 0.0)
    avail = maximum(th - thr, 0.0)
    take = where(take > avail, avail, take)
    th0 = where(fresh, th - take, th)
    return (th0, where(fresh, eta + take * dz * 1000.0, eta), where(fresh, th0, old), where(fresh, th0, avg),
            where(fresh, suction(th0, thr, ths, alpha, inv_n, inv_m), h1), dt_step / n)


@gtx.field_operator
def heads(old: F, avg: F, h1: F, dt: F, n: F, rain: F, ks: F, ths: F, thr: F, alpha: F, inv_n: F, m: F,
          inv_m: F, zc: F) -> tuple[F, F, F]:
    """Hydraulic heads and conductivities of the iteration (WaSiM: hydraulic_heads_iterative), top inflow."""
    h2 = suction(avg, thr, ths, alpha, inv_n, inv_m)
    return zc - 0.5 * (h1 + h2), ks * krel(avg, thr, ths, m, inv_m), rain / n / (1000.0 * dt)


@gtx.field_operator
def new_thetas(th0: F, old: F, avg: F, z: F, qtop: F, qbot: F, dt: F, dz: F) -> tuple[F, F, F]:
    """New water contents and their moving average (WaSiM new_thetas); per level change e."""
    z1 = z + 1.0
    weight = z1 - 3.0
    early = z1 <= 3.0
    avg2 = where(early, th0, avg)
    thn = old + (qtop - qbot) * dt / dz
    thn = where(thn < 0.0, 0.0, where(thn > 1.0, 1.0, thn))
    return thn, abs(th0 - thn), where(early, (avg2 + thn) / 2.0, (avg2 * weight + thn) / (weight + 1.0))


@gtx.field_operator
def finish(thn: F, avg3: F, old: F, h1: F, z: F, s: F, done: F, drain: F, excess: F, eta1: F, iters: F,
           err: F, qsurf: F, qdeep: F, dt: F, n: F, rain: F, th_in: F, old_in: F, avg_in: F, h1_in: F,
           max_err: float64, max_iter: float64) -> tuple[F, F, F, F, F, F, F, F, F, F, F]:
    """Convergence; at the end of a sub-step the fluxes go into the sums."""
    live = done == 0.0
    z1 = z + 1.0
    maxe = where(z1 >= 21.0, max_err * 5.0, max_err)
    end = live & ((err <= maxe) | (z1 >= max_iter))
    s1 = where(end, s + 1.0, s)
    return (where(live, thn, th_in), where(live, old, old_in), where(live, avg3, avg_in), where(live, h1, h1_in),
            where(live, where(end, 0.0, z1), z), s1, where(end & (s1 >= n), 1.0, done),
            where(end, drain + qdeep * dt * 1000.0, drain),
            where(end, excess + (rain / n - qsurf * dt * 1000.0), excess),
            eta1, where(end, iters + z1, iters))


@gtx.program
def tick(th: F, old: F, avg: F, h1: F, z: F, s: F, done: F, drain: F, excess: F, eta: F, iters: F,
         n: F, rain: F, et: F, root: F,
         ks: F, ths: F, thr: F, alpha: F, inv_n: F, m: F, inv_m: F, dz: F, zc: F,
         dt_step: float64, max_err: float64, max_iter: float64, root_levels: float64,
         t_th0: F, t_eta1: F, t_old1: F, t_avg1: F, t_h11: F, t_dt: F, t_hh: F, t_kk: F, t_q0: F,
         t_qraw: F, t_c1: F, t_c2: F, t_c3: F, t_qtop: F, t_qbot: F, t_c4: F, t_qtop2: F, t_qbot2: F, t_c5: F,
         t_thn: F, t_e: F, t_avg3: F, t_emax: F, t_qsurf: F, t_c6: F, t_err: F, t_qdeep: F, t_c7: F,
         th_o: F, old_o: F, avg_o: F, h1_o: F, z_o: F, s_o: F, done_o: F, drain_o: F, excess_o: F, eta_o: F,
         iters_o: F):
    prepare(th, old, avg, h1, z, done, eta, n, rain, et, root, ths, thr, alpha, inv_n, inv_m, dz, dt_step,
            root_levels, out=(t_th0, t_eta1, t_old1, t_avg1, t_h11, t_dt))
    heads(t_old1, t_avg1, t_h11, t_dt, n, rain, ks, ths, thr, alpha, inv_n, m, inv_m, zc, out=(t_hh, t_kk, t_q0))
    darcy(t_hh, t_kk, dz, out=(t_qraw, t_c1, t_c2, t_c3))
    drying(t_q0, t_qraw, t_old1, thr, dz, t_dt, out=(t_qtop, t_qbot, t_c4))
    saturating(t_qtop, t_qbot, t_old1, ths, dz, t_dt, out=(t_qtop2, t_qbot2, t_c5))
    new_thetas(t_th0, t_old1, t_avg1, z, t_qtop2, t_qbot2, t_dt, dz, out=(t_thn, t_e, t_avg3))
    down_max(t_e, t_qtop2, out=(t_emax, t_qsurf, t_c6))
    up_bottom(t_emax, t_qbot2, out=(t_err, t_qdeep, t_c7))
    finish(t_thn, t_avg3, t_old1, t_h11, z, s, done, drain, excess, t_eta1, iters, t_err, t_qsurf, t_qdeep, t_dt, n,
           rain, th, old, avg, h1, max_err, max_iter,
           out=(th_o, old_o, avg_o, h1_o, z_o, s_o, done_o, drain_o, excess_o, eta_o, iters_o))


STAGES = (prepare, heads, darcy, drying, saturating, new_thetas, down_max, up_bottom, finish)


def tick_split(ops, th, old, avg, h1, z, s, done, drain, excess, eta, iters, n, rain, et, root,
               ks, ths, thr, alpha, inv_n, m, inv_m, dz, zc, dt_step, max_err, max_iter, root_levels,
               t_th0, t_eta1, t_old1, t_avg1, t_h11, t_dt, t_hh, t_kk, t_q0, t_qraw, t_c1, t_c2, t_c3, t_qtop, t_qbot,
               t_c4, t_qtop2, t_qbot2, t_c5, t_thn, t_e, t_avg3, t_emax, t_qsurf, t_c6, t_err, t_qdeep, t_c7,
               th_o, old_o, avg_o, h1_o, z_o, s_o, done_o, drain_o, excess_o, eta_o, iters_o, offset_provider):
    """The program `tick`, as one call per stage from Python (--split)."""
    o = offset_provider
    ops[0](th, old, avg, h1, z, done, eta, n, rain, et, root, ths, thr, alpha, inv_n, inv_m, dz, dt_step, root_levels,
           out=(t_th0, t_eta1, t_old1, t_avg1, t_h11, t_dt), offset_provider=o)
    ops[1](t_old1, t_avg1, t_h11, t_dt, n, rain, ks, ths, thr, alpha, inv_n, m, inv_m, zc, out=(t_hh, t_kk, t_q0),
           offset_provider=o)
    ops[2](t_hh, t_kk, dz, out=(t_qraw, t_c1, t_c2, t_c3), offset_provider=o)
    ops[3](t_q0, t_qraw, t_old1, thr, dz, t_dt, out=(t_qtop, t_qbot, t_c4), offset_provider=o)
    ops[4](t_qtop, t_qbot, t_old1, ths, dz, t_dt, out=(t_qtop2, t_qbot2, t_c5), offset_provider=o)
    ops[5](t_th0, t_old1, t_avg1, z, t_qtop2, t_qbot2, t_dt, dz, out=(t_thn, t_e, t_avg3), offset_provider=o)
    ops[6](t_e, t_qtop2, out=(t_emax, t_qsurf, t_c6), offset_provider=o)
    ops[7](t_emax, t_qbot2, out=(t_err, t_qdeep, t_c7), offset_provider=o)
    ops[8](t_thn, t_avg3, t_old1, t_h11, z, s, done, drain, excess, t_eta1, iters, t_err, t_qsurf, t_qdeep, t_dt, n,
           rain, th, old, avg, h1, max_err, max_iter,
           out=(th_o, old_o, avg_o, h1_o, z_o, s_o, done_o, drain_o, excess_o, eta_o, iters_o), offset_provider=o)


def backend_of(name):
    if name == "embedded":
        return None
    if name == "gtfn_cpu":
        return gtx.gtfn_cpu
    if name == "dace_cpu":
        from gt4py.next.program_processors.runners.dace import run_dace_cpu
        return run_dace_cpu
    raise SystemExit(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--mode", default="adaptive", choices=["adaptive", "fixed"])
    ap.add_argument("--substeps", type=int, default=1)
    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--check", type=int, default=1, help="adaptive: ticks between host checks")
    ap.add_argument("--backend", default="gtfn_cpu")
    ap.add_argument("--steps", type=int, default=0)
    ap.add_argument("--split", action="store_true", help="call the stages of a tick one by one instead of the program")
    args = ap.parse_args()

    case = np.load(os.path.join(args.case_dir, "case.npz"))
    params, theta = case["params"], case["theta"]
    ncell, _, nlev = params.shape
    nsteps = args.steps or case["rain"].shape[0]
    dom = gtx.domain({C: ncell, K: nlev})
    fld = lambda a: gtx.as_field(dom, np.ascontiguousarray(a, dtype=np.float64))  # noqa: E731
    col = lambda v: fld(np.broadcast_to(np.asarray(v, dtype=np.float64)[:, None], (ncell, nlev)))  # noqa: E731
    p = {name: fld(params[:, i, :]) for i, name in enumerate(PARAMS)}
    root = fld(np.broadcast_to(np.arange(nlev) < int(case["root_levels"]), (ncell, nlev)))
    below = fld(np.broadcast_to(np.arange(nlev) < nlev - 1, (ncell, nlev)))
    dt_step, max_err = float(case["dt_step"]), float(case["max_err"])
    max_iter = float(case["max_iter"]) if args.mode == "adaptive" else float(args.iters)
    n_fixed = 0.0 if args.mode == "adaptive" else float(args.substeps)

    b = backend_of(args.backend)
    step_op = tick if b is None else tick.with_backend(b)
    if args.split:
        ops = [op if b is None else op.with_backend(b) for op in STAGES]
        step_op = lambda *a, offset_provider: tick_split(ops, *a, offset_provider=offset_provider)  # noqa: E731
    sub_op = substeps if b is None else substeps.with_backend(b)

    zeros = lambda: gtx.zeros(dom, dtype=float64)  # noqa: E731
    temps = [zeros() for _ in range(28)]
    sub_temps = [zeros() for _ in range(4)]
    th = fld(theta)
    names = ("th", "old", "avg", "h1", "z", "s", "done", "drain", "excess", "eta", "iters")
    bufs = [{k: zeros() for k in names}, {k: zeros() for k in names}]
    start = {k: zeros() for k in names}
    n = zeros()
    n_const = gtx.full(dom, n_fixed, dtype=float64)
    out_drain, out_excess, out_eta, out_sub, out_iters = (np.zeros((nsteps, ncell)) for _ in range(5))
    budget = args.substeps * args.iters if args.mode == "fixed" else None
    ticks_total, t_run, t_first = 0, 0.0, 0.0
    for t in range(nsteps):
        rain, et = col(case["rain"][t]), col(case["et"][t])
        t0 = time.perf_counter()
        start["th"] = th
        if n_fixed > 0.0:
            n = n_const
        else:
            sub_op(th, rain, p["ks"], p["ths"], p["thr"], p["m"], p["inv_m"], p["dz"], below, dt_step,
                   float(case["min_step"]), *sub_temps, n, offset_provider={})
        cur, ticks, i = start, 0, 0
        while True:
            nxt = bufs[i % 2]
            step_op(*(cur[k] for k in names), n, rain, et, root, *(p[k] for k in PARAMS),
                    dt_step, max_err, max_iter, float(case["root_levels"]), *temps, *(nxt[k] for k in names),
                    offset_provider={})
            cur, i, ticks = nxt, i + 1, ticks + 1
            if budget is not None:
                if ticks >= budget:
                    break
            elif ticks % args.check == 0 and cur["done"].asnumpy().min() > 0.0:
                break
        th = gtx.as_field(dom, cur["th"].asnumpy().copy())
        t_run += time.perf_counter() - t0
        if t == 0:
            t_first = t_run
        ticks_total += ticks
        out_drain[t], out_excess[t] = cur["drain"].asnumpy()[:, 0], cur["excess"].asnumpy()[:, 0]
        out_eta[t] = cur["eta"].asnumpy()[:, :int(case["root_levels"])].sum(axis=1)
        out_sub[t], out_iters[t] = n.asnumpy()[:, 0], cur["iters"].asnumpy()[:, 0]
    os.makedirs(args.out_dir, exist_ok=True)
    np.ascontiguousarray(th.asnumpy()).tofile(os.path.join(args.out_dir, "theta.f64"))
    for name, a in (("drain", out_drain), ("excess", out_excess), ("eta", out_eta), ("substeps", out_sub),
                    ("iters", out_iters)):
        a.tofile(os.path.join(args.out_dir, name + ".f64"))
    warm = (t_run - t_first) / max(nsteps - 1, 1)
    print(f"{ncell} cells, {nlev} levels, {nsteps} steps, {args.backend}, mode {args.mode}"
          f"{'' if args.mode == 'adaptive' else f' N={args.substeps} I={args.iters}'}{', split' if args.split else ''}: "
          f"{1000 * warm:.2f} ms per step (warm; first step {t_first:.1f} s), "
          f"{ticks_total / nsteps:.1f} ticks per step, {1000 * warm / (ticks_total / nsteps):.3f} ms per tick")


if __name__ == "__main__":
    main()
