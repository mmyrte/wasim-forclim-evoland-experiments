#!/usr/bin/env python3
"""Validate wasim-gt's glacier mass balance and yearly glacier dynamics (wasim_gt.processes.glacier_dynamics)
against WaSiM grid dumps.

  - every step: the glacier mass balance GMB, from the last sub step's after_gwflow plus this step's
    snow, melt, slide and evaporation grids, compared with after_snow2
  - the evaluation step: the slides (snow_slides, validated in snow_energy_balance.py), the mass balance,
    then the yearly update with the glacier bodies set up from the start dump, compared with after_snow2:
    glacier fraction and codes (and their old values), mass balance and its reference, firn stack, and
    the snow storages off and on the glacier fraction
Every compared grid must match exactly, including NaN, infinity and nodata patterns.

The evaluation runs at 24:00 of the evaluation day of [ice_firn]. WaSiM tests the hour with an offset of 1
(60 * hour + 60 * hour_offset == 1440), so it is the step that ends at 23:00. The golden Rhone run (January)
does not reach it; the Rhone summer variant (1-3 July 2021) moves the date to 1 July.

Grids to dump at start, after_snow1, after_gwflow and after_snow2: [zone_grid], zonegrid_soilmodel,
[elevation_model], GlacierizedCells, GlacierizedCellsOld, GlacierCodes, GlacierCodesOld, GlacierMassBalance,
OldGlacierMassBalance, firn_WE_stack, snowstorage_solid, snowstorage_liquid, snow_storage_total,
SolidSnowStorage_onGlacier, LiquidSnowStorage_onGlacier, TotalSnowStorage_onGlacier, snow_fraction,
rain_fraction, SnowMelt_onGlacier, melt_from_firn, melt_from_ice, SnowEvaporation_onGlacier, input_mass,
mobile_mass, deposition, slope_angle, slidefraction1..4, max_deposition, elev_order_row, elev_order_col.

usage: glacier_dynamics.py DUMP_DIR CONTROL_FILE --start YYYY-MM-DD [--first-hour 1] [--cellsize 500]
"""
import argparse
import datetime as dt
import re

import numpy as np

from wasim_gt.processes.glacier_dynamics import DynamicsParameters, GlacierModel, mass_balance_step
from wasim_gt.processes.snow_slides import snow_slides

from snow_energy_balance import control_section, difference, snow_parameters
from wasim_dump import Dump

COMPARED = ("GlacierizedCells", "GlacierizedCellsOld", "GlacierCodes", "GlacierCodesOld", "GlacierMassBalance",
            "OldGlacierMassBalance", "snow_storage_total", "snowstorage_solid", "snowstorage_liquid",
            "TotalSnowStorage_onGlacier", "SolidSnowStorage_onGlacier", "LiquidSnowStorage_onGlacier")


def dynamics_parameters(control_file):
    """Evaluation date, the volume-area line and the new-glacier line of [ice_firn]."""
    eval_date = va = new = None
    for line in control_section(control_file, "ice_firn").splitlines():
        value = line.split("#", 1)[0].split()
        comment = line.partition("#")[2]
        if "month and day" in comment:
            eval_date = (int(value[0]), int(value[1]))
        elif "VAscaling" in comment:
            va = value
        elif len(value) == 2 and re.fullmatch(r"[0-9.]+", value[0]) and "SWE" in comment:
            new = value
    params = DynamicsParameters(va_scaling=float(va[0]), va_exponent=float(va[1]), elevation_bands=int(va[2]),
                                extra_weight_band0=int(va[3]), max_snow_new_glacier=float(new[0]),
                                rel_size_new_glacier=float(new[1]))
    return eval_date, params


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_file")
    ap.add_argument("--start", required=True, help="date of the first step, YYYY-MM-DD")
    ap.add_argument("--first-hour", type=int, default=1)
    ap.add_argument("--cellsize", type=float, default=500.0)
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    (month, day), params = dynamics_parameters(args.control_file)
    _, slide = snow_parameters(args.control_file)
    s0 = lambda name: d("start", name)[0]
    zone, zone_routing, dhm = s0("[zone_grid]"), s0("zonegrid_soilmodel"), s0("[elevation_model]")
    active = zone > -9999
    model = GlacierModel.from_grids(params, args.cellsize, dhm, zone_routing, zone, s0("GlacierizedCells"),
                                    s0("GlacierCodes"), d("start", "firn_WE_stack")[0][..., 0])
    nsteps = d("after_snow2", "GlacierMassBalance").shape[0]
    nsub = d("after_gwflow", "GlacierMassBalance").shape[0] // nsteps
    start = dt.datetime.fromisoformat(args.start)
    worst = {}

    def check(name, value, ref):
        worst[name] = max(worst.get(name, 0.0), difference(value, ref, active))

    evaluated = []
    for t in range(nsteps):
        gw = lambda name: d("after_gwflow", name)[(t + 1) * nsub - 1].copy()
        s1 = lambda name: d("after_snow1", name)[t]
        s2 = lambda name: d("after_snow2", name)[t]
        ice = gw("GlacierizedCells")
        gmb = mass_balance_step(gw("GlacierMassBalance"), ice, s1("snow_fraction"), s1("rain_fraction"),
                                s1("SnowMelt_onGlacier"), s1("melt_from_firn"), s1("melt_from_ice"),
                                gw("SnowEvaporation_onGlacier"), active, s2("deposition"), s2("mobile_mass"),
                                s2("input_mass"))
        end = start + dt.timedelta(hours=t + args.first_hour)
        if not (end.month == month and end.day == day and end.hour == 23):
            check("GlacierMassBalance (every step)", gmb, s2("GlacierMassBalance"))
            continue
        evaluated.append(t)
        glacier = tuple(gw(k) for k in ("TotalSnowStorage_onGlacier", "LiquidSnowStorage_onGlacier",
                                        "SolidSnowStorage_onGlacier"))
        (ssto, sliq, ssnow), _, (ssto_gl, sliq_gl, ssnow_gl) = snow_slides(
            gw("snow_storage_total"), gw("snowstorage_liquid"), gw("snowstorage_solid"), active, zone,
            s0("slope_angle"), [s0(f"slidefraction{k}") for k in range(1, 5)], s0("max_deposition"),
            s0("elev_order_row"), s0("elev_order_col"), ice=ice, glacier=glacier, **slide)
        ice_old, gid, gid_old = gw("GlacierizedCellsOld"), gw("GlacierCodes"), gw("GlacierCodesOld")
        gmbo, firn = gw("OldGlacierMassBalance"), gw("firn_WE_stack")
        model.yearly_update(ice, ice_old, gid, gid_old, gmb, gmbo, firn, ssto, ssnow, sliq, ssto_gl, ssnow_gl,
                            sliq_gl)
        result = dict(GlacierizedCells=ice, GlacierizedCellsOld=ice_old, GlacierCodes=gid, GlacierCodesOld=gid_old,
                      GlacierMassBalance=gmb, OldGlacierMassBalance=gmbo, snow_storage_total=ssto,
                      snowstorage_solid=ssnow, snowstorage_liquid=sliq, TotalSnowStorage_onGlacier=ssto_gl,
                      SolidSnowStorage_onGlacier=ssnow_gl, LiquidSnowStorage_onGlacier=sliq_gl)
        for k in COMPARED:
            check(k, result[k], s2(k))
        ref_firn = s2("firn_WE_stack")
        for k in range(firn.shape[-1]):
            check(f"firn_WE_stack[{k}]", firn[..., k], ref_firn[..., k])
        changed = int((ice != gw("GlacierizedCells"))[active].sum())
        print(f"evaluation at step {t}: {len(model.glaciers)} glaciers, glacier fraction changed in {changed} cells")

    print(f"{active.sum()} active cells, {nsteps} steps, evaluation steps {evaluated}")
    print("largest absolute difference to WaSiM (inf: NaN or nodata patterns differ):")
    for k, v in worst.items():
        print(f"  {k}: {v:.3g}")
    return 0 if evaluated and max(worst.values()) == 0.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
