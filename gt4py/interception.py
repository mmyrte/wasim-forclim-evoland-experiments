#!/usr/bin/env python3
"""GT4Py port of WaSiM's interception model (methods 1 and 11), checked against grid dumps.

What it reimplements (src/interzep.cc, interzeptionsmodellklasse::run, multi-layer branch):
for each vegetation layer from the top,
  1. capacity SImax = LAI * IntercepCap (land use table; SIdicke from the control file if missing)
  2. the uncovered fraction (1 - VBG) of the input falls through, the covered fraction fills SI
  3. SI above SImax overflows into the throughfall Q of the layer
  4. interception evaporation EI = min(SI, ETP of the layer); EI is taken from ETP, which
     belongs to the evapotranspiration module: this is one of the cross-module writes
     (gt4py/coupling.md). Methods 1 and 11 use ETP, 2 and 12 use EIP (not ported).
  5. the next layer's input is this layer's throughfall
Layer 0 holds the totals: EI0 = sum EI, SI0 = sum SI, ETP0 -= sum EI, Q0 = throughfall of the
lowest layer with vegetation. Cells with a glacier fraction get Q0 = Q0 (1 - ice) + P ice;
fully glacierized cells, ArcVeg, the snow canopy model (methods > 20) and tracers are not ported.

The inputs are the states WaSiM had just before the interception call (dump point after_etp,
because irrigation is off), and the results are compared with the dump after_interception.

usage: interception.py DUMP_DIR CONTROL_FILE [--backend numpy|gtfn_cpu]
"""
import argparse
import time

import numpy as np
import gt4py.next as gtx
from gt4py.next import where, minimum, float64

from wasim_dump import Dump, landuse_table

I = gtx.Dimension("I")
J = gtx.Dimension("J")
F = gtx.Field[gtx.Dims[I, J], float64]
M = gtx.Field[gtx.Dims[I, J], bool]


@gtx.field_operator
def layer_step(inp: F, lai: F, vbg: F, cap: F, etp: F, si: F, has_veg: M,
               etp0: F, ei0: F, si0: F, q0: F,
               ) -> tuple[F, F, F, F, F, F, F, F]:
    simax = lai * cap
    q = (1.0 - vbg) * inp
    si1 = si + vbg * inp
    over = si1 > simax
    q = where(over, q + (si1 - simax), q)
    si1 = where(over, simax, si1)
    water_limited = si1 <= etp
    ei = where(water_limited, si1, etp)
    etp_new = etp - ei
    si2 = where(water_limited, 0.0, si1 - ei)
    # cells without vegetation in this layer: nothing is intercepted and nothing passes on
    ei = where(has_veg, ei, 0.0)
    q = where(has_veg, q, 0.0)
    si2 = where(has_veg, si2, 0.0)
    etp_new = where(has_veg, etp_new, etp)
    return (ei, q, si2, etp_new,
            where(has_veg, etp0 - ei, etp0), ei0 + ei, si0 + si2, where(has_veg, q, q0))


@gtx.field_operator
def glacier_mix(q0: F, prec: F, ice: F) -> F:
    return where(ice > 0.0, q0 * (1.0 - ice) + prec * ice, q0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_file", help="the filled control file of the run (run_case.py writes it to the output directory)")
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--sidicke", type=float, default=0.35)
    ap.add_argument("--backend", default="numpy", choices=["numpy", "gtfn_cpu"])
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    pre, post = "after_etp", "after_interception"
    nl = args.layers
    cap_table = landuse_table(args.control_file, "IntercepCap")
    prec = d(pre, "precipitation")
    nsteps, nrows, ncols = prec.shape
    luse = [d(pre, f"landuse{l}") for l in range(1, nl + 1)]
    has_veg = [u > -9999 for u in luse]
    active = d(post, "interception_storage0")[0] > -9999

    def cap_grid(lu):
        c = np.full(lu.shape, args.sidicke)
        for code, v in cap_table.items():
            if v >= 0:
                c[lu == code] = v
        return c

    backend = gtx.gtfn_cpu if args.backend == "gtfn_cpu" else None
    step = layer_step.with_backend(backend) if backend else layer_step
    mix = glacier_mix.with_backend(backend) if backend else glacier_mix
    dom = gtx.domain({I: nrows, J: ncols})
    fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))
    msk = lambda a: gtx.as_field(dom, a & active)
    zero = gtx.zeros(dom, dtype=float64)

    def call(op, *a, nout):
        out = tuple(gtx.zeros(dom, dtype=float64) for _ in range(nout))
        op(*a, out=out if nout > 1 else out[0], offset_provider={})
        return out if nout > 1 else out[0]
    ice = zero  # no glaciers in the Doubs setup

    names = {"ei": "interception_evaporation", "q": "throughfall", "si": "interception_storage", "etp": "pot_transpiration"}
    worst = {(k, l): 0.0 for k in names for l in range(nl + 1)}
    t_run = 0.0
    for t in range(nsteps):
        start = time.perf_counter()
        etp0, ei0, si0, q0 = fld(d(pre, "pot_transpiration0")[t]), zero, zero, zero
        inp = fld(prec[t])
        res = {}
        for l in range(1, nl + 1):
            ei, q, si, etp, etp0, ei0, si0, q0 = call(step,
                inp, fld(d(pre, f"leaf_area_index{l}")[t]), fld(d(pre, f"vegetation_coverage_degree{l}")[t]),
                fld(cap_grid(luse[l - 1][t])), fld(d(pre, f"pot_transpiration{l}")[t]),
                fld(d(pre, f"interception_storage{l}")[t]), msk(has_veg[l - 1][t]),
                etp0, ei0, si0, q0, nout=8)
            res[l] = dict(ei=ei, q=q, si=si, etp=etp)
            inp = q
        res[0] = dict(ei=ei0, q=call(mix, q0, fld(prec[t]), ice, nout=1), si=si0, etp=etp0)
        t_run += time.perf_counter() - start
        for (k, l) in worst:
            ref = d(post, f"{names[k]}{l}")[t]
            diff = np.abs(res[l][k].asnumpy() - ref)[active]
            worst[(k, l)] = max(worst[(k, l)], float(diff.max()))

    print(f"grid {nrows} x {ncols}, {active.sum()} active cells, {nl} vegetation layers, {nsteps} steps, backend {args.backend}")
    print(f"{1000 * t_run / nsteps:.1f} ms per step including host-device copies (first step includes compilation)")
    print("largest absolute difference to WaSiM over all steps and active cells [mm]:")
    for (k, l), v in worst.items():
        print(f"  {names[k]}{l}: {v:.3g}")
    return 0 if max(worst.values()) <= 1e-12 else 1


if __name__ == "__main__":
    raise SystemExit(main())
