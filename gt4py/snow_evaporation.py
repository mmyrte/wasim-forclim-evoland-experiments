#!/usr/bin/env python3
"""Validate wasim-gt's snow evaporation (wasim_gt.processes.snow_evaporation) against WaSiM grid dumps.

Inputs: dump point after_interception (WaSiM order for interception method 1: snow run1,
interception, snow evaporation; for method 11: interception, snow run1, snow evaporation, so use
--pre after_snow1 there). Reference: after_snowevap, every step and active cell.

usage: snow_evaporation.py DUMP_DIR [--energy-balance] [--pre POINT] [--backend embedded|gtfn_cpu]
"""
import argparse
import time

import numpy as np
import gt4py.next as gtx

from wasim_gt.fields import I, J
from wasim_gt.processes.snow_evaporation import conventional, energy_balance

from wasim_dump import Dump


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--pre", default="after_interception")
    ap.add_argument("--energy-balance", action="store_true", help="snow method 5, 6, 7 or 11")
    ap.add_argument("--backend", default="embedded", choices=["embedded", "gtfn_cpu"])
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    pre, post = args.pre, "after_snowevap"
    nl = args.layers
    ssnow_in = d(pre, "snowstorage_solid")
    nsteps, nrows, ncols = ssnow_in.shape
    active = ssnow_in[0] > -9999
    dom = gtx.domain({I: nrows, J: ncols})
    fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))
    msk = lambda a: gtx.as_field(dom, a & active)

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
        etps = [fld(d(pre, f"pot_transpiration{l}")[t]) for l in range(1, nl + 1)]
        veg = [msk(d(pre, f"landuse{l}")[t] > -9999) for l in range(1, nl + 1)]
        snow_f = msk(snow)
        if args.energy_balance:
            etps, evap = energy_balance(etps, veg, evap, snow_f, dom, backend=args.backend)
            res = dict(snowstorage_solid=ssnow, potential_evaporation=evap)
        else:
            ssnow, etps, evap, evrs = conventional(ssnow, etps, veg, evap, snow_f, dom, backend=args.backend)
            res = dict(snowstorage_solid=ssnow, potential_evaporation=evap, snow_evaporation=evrs)
        res.update({f"pot_transpiration{l}": e for l, e in enumerate(etps, start=1)})
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
