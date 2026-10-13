#!/usr/bin/env python3
"""GT4Py port of WaSiM's snow evaporation (verdunstungsmodellklasse::snow), checked against grid dumps.

The evapotranspiration module takes snow evaporation from the snow model's solid storage
(src/etpetr.cc, verdunstungsmodellklasse::snow), after interception and the snow model's first
part have run. It writes three grids it does not own: SSNOW (snow model), and through ETP and
EVAP its own potentials that the soil model uses later (gt4py/coupling.md).

Conventional snow methods (T-index etc.), cells with snow:
  for each vegetation layer: EVRS += min(SSNOW, ETP_layer), taken from both SSNOW and ETP_layer
  then the same for the soil evaporation EVAP
Energy-balance snow methods (5, 6, 7, 11): EVRS comes from the snow model; where there is snow,
ETP of every vegetation layer and EVAP are set to 0.
Not ported: the glacierized fraction of cells (SSNOWGL, EVRSGL), spin-up HRU copies.
Note: the layer-0 total ETP is not reduced here (WaSiM does the same; interception does reduce it).

Inputs: dump point after_interception (WaSiM order for interception method 1: snow run1,
interception, snow evaporation; for method 11: interception, snow run1, snow evaporation, so use
after_snow1 there). Reference: after_snowevap.

usage: snow_evaporation.py DUMP_DIR [--energy-balance] [--pre POINT] [--backend numpy|gtfn_cpu]
"""
import argparse
import time

import numpy as np
import gt4py.next as gtx
from gt4py.next import where, float64

from wasim_dump import Dump

I = gtx.Dimension("I")
J = gtx.Dimension("J")
F = gtx.Field[gtx.Dims[I, J], float64]
M = gtx.Field[gtx.Dims[I, J], bool]


@gtx.field_operator
def take(store: F, demand: F, use: M, evrs: F) -> tuple[F, F, F]:
    """Evaporate min(store, demand) from the snow store where use is set."""
    e = where(store < demand, store, demand)
    e = where(use, e, 0.0)
    return store - e, demand - e, evrs + e


@gtx.field_operator
def zero_where(x: F, cond: M) -> F:
    return where(cond, 0.0, x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--pre", default="after_interception")
    ap.add_argument("--energy-balance", action="store_true", help="snow method 5, 6, 7 or 11")
    ap.add_argument("--backend", default="numpy", choices=["numpy", "gtfn_cpu"])
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    pre, post = args.pre, "after_snowevap"
    nl = args.layers
    ssnow_in = d(pre, "snowstorage_solid")
    nsteps, nrows, ncols = ssnow_in.shape
    active = ssnow_in[0] > -9999
    backend = gtx.gtfn_cpu if args.backend == "gtfn_cpu" else None
    take_op = take.with_backend(backend) if backend else take
    zero_op = zero_where.with_backend(backend) if backend else zero_where
    dom = gtx.domain({I: nrows, J: ncols})
    fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))
    msk = lambda a: gtx.as_field(dom, a & active)

    def call(op, *a, nout):
        out = tuple(gtx.zeros(dom, dtype=float64) for _ in range(nout))
        op(*a, out=out if nout > 1 else out[0], offset_provider={})
        return out if nout > 1 else out[0]

    names = ["snowstorage_solid", "potential_evaporation"] + [f"pot_transpiration{l}" for l in range(1, nl + 1)]
    if not args.energy_balance:
        names.append("snow_evaporation")
    worst = dict.fromkeys(names, 0.0)
    t_run = 0.0
    for t in range(nsteps):
        start = time.perf_counter()
        ssnow = fld(ssnow_in[t])
        snow = ssnow_in[t] > 0.0
        evap = fld(d(pre, "potential_evaporation")[t])
        res = {}
        if args.energy_balance:
            for l in range(1, nl + 1):
                use = msk(snow & (d(pre, f"landuse{l}")[t] > -9999))
                res[f"pot_transpiration{l}"] = call(zero_op, fld(d(pre, f"pot_transpiration{l}")[t]), use, nout=1)
            res["potential_evaporation"] = call(zero_op, evap, msk(snow), nout=1)
            res["snowstorage_solid"] = ssnow
        else:
            evrs = gtx.zeros(dom, dtype=float64)
            for l in range(1, nl + 1):
                use = msk(snow & (d(pre, f"landuse{l}")[t] > -9999))
                ssnow, etp, evrs = call(take_op, ssnow, fld(d(pre, f"pot_transpiration{l}")[t]), use, evrs, nout=3)
                res[f"pot_transpiration{l}"] = etp
            ssnow, evap, evrs = call(take_op, ssnow, evap, msk(snow), evrs, nout=3)
            res.update(snowstorage_solid=ssnow, potential_evaporation=evap, snow_evaporation=evrs)
        t_run += time.perf_counter() - start
        for name in names:
            diff = np.abs(res[name].asnumpy() - d(post, name)[t])[active]
            worst[name] = max(worst[name], float(diff.max()))

    print(f"grid {nrows} x {ncols}, {active.sum()} active cells, {nsteps} steps, "
          f"{'energy-balance' if args.energy_balance else 'conventional'} snow method, backend {args.backend}")
    print(f"{1000 * t_run / nsteps:.1f} ms per step including host-device copies (first step includes compilation)")
    print("largest absolute difference to WaSiM over all steps and active cells [mm]:")
    for name, v in worst.items():
        print(f"  {name}: {v:.3g}")
    return 0 if max(worst.values()) <= 1e-12 else 1


if __name__ == "__main__":
    raise SystemExit(main())
