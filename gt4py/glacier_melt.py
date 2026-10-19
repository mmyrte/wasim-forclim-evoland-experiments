#!/usr/bin/env python3
"""Validate wasim-gt's glacier melt (snow on the glacier fraction and ice/firn melt, glacier method 13)
against WaSiM grid dumps.

Per step, from the dump point before the snow model (after_irrigation in Rhone) to after_snow1:
  1. snow_energy_balance on all cells (the fraction without ice), as in snow_energy_balance.py, now
     compared on glacierized cells too
  2. snow_energy_balance(on_glacier=True) on cells with a glacier fraction, with the first call's
     precipitation, snow and rain fractions and albedo
  3. glacier_melt: ice and firn melt and the firn stack
Every compared grid must match exactly, including NaN, infinity and nodata patterns.

The golden Rhone run (January) never melts ice: every glacier carries snow. A variant runs 1-3 July
2021 from a copy of the Rhone states with the snow on glaciers set to 0 (ssno_onGlac_, sliq_onGlac_,
ssto_onGlac_), so that ice and firn melt. A second variant adds a constant debris grid (internal standard
grid debris_on_glaciers, readcode = 0, value 0.7) to reach the debris factor.

Grids to dump at start, after_meteo, after_irrigation and after_snow1: those of snow_energy_balance.py,
plus GlacierizedCells, global_radiation, firn_WE_stack, melt_from_ice, melt_from_firn,
SolidSnowStorage_onGlacier, LiquidSnowStorage_onGlacier, TotalSnowStorage_onGlacier, snowtemp_onGlacier,
surftemp_onGlacier, SnowEvaporation_onGlacier, SnowMelt_onGlacier.

usage: glacier_melt.py DUMP_DIR CONTROL_FILE --start YYYY-MM-DD [--first-hour 1] [--pre after_irrigation]
                       [--input precipitation] [--backend gt:cpu_ifirst]
"""
import argparse
import datetime as dt
import re

import gt4py.next as gtx
import numpy as np

from wasim_gt.fields import I, J
from wasim_gt.processes.glacier_melt import GlacierMeltParameters, glacier_melt, radiation_coefficient
from wasim_gt.processes.snow_energy_balance import snow_energy_balance

from snow_energy_balance import control_section, difference, landuse_albedo, snow_parameters
from wasim_dump import Dump

GLACIER_SNOW = {"snowstorage_solid": "SolidSnowStorage_onGlacier", "snowstorage_liquid": "LiquidSnowStorage_onGlacier",
                "snow_storage_total": "TotalSnowStorage_onGlacier", "snowtemp": "snowtemp_onGlacier",
                "surftemp": "surftemp_onGlacier", "snow_evaporation": "SnowEvaporation_onGlacier",
                "snowcover_outflow": "SnowMelt_onGlacier"}
FIRST_CALL = ("snowstorage_solid", "snowstorage_liquid", "snow_storage_total", "snowtemp", "surftemp",
              "snow_age", "albedo", "snowcover_outflow", "melt_from_snowcover", "days_snow", "snow_evaporation")


def ice_firn_parameters(control_file):
    """[ice_firn] numbers in reading order: method, t-index factors (ice, firn, snow), melt factor MF,
    radiation coefficients (ice min, ice max, snow min, snow max), ... One value per line is supported."""
    numbers = []
    for line in control_section(control_file, "ice_firn").splitlines():
        line = line.split("#", 1)[0].strip()
        if line and re.fullmatch(r"[-+0-9.eE\s]+", line):
            numbers.append([float(v) for v in line.split()])
    if int(numbers[0][0]) != 13:
        raise SystemExit(f"glacier method {numbers[0][0]}: this check covers method 13 only")
    if any(len(n) != 1 for n in numbers[1:9]):
        raise SystemExit("per-subbasin glacier parameters are not supported by this check")
    v = [n[0] for n in numbers[1:9]]
    return dict(mf=v[3], ice_min=v[4], ice_max=v[5], snow_min=v[6], snow_max=v[7])


def debris_scaling(control_file):
    """The global debris scaling factor of [ice_firn] (the numeric line whose comment mentions the debris grid)."""
    for line in control_section(control_file, "ice_firn").splitlines():
        value, _, comment = line.partition("#")
        if "debris" in comment and re.fullmatch(r"\s*[-+0-9.eE]+\s*", value):
            return float(value)
    raise SystemExit("no debris scaling factor in [ice_firn]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_file")
    ap.add_argument("--start", required=True, help="date of the first step, YYYY-MM-DD")
    ap.add_argument("--first-hour", type=int, default=1)
    ap.add_argument("--pre", default="after_irrigation")
    ap.add_argument("--input", default="precipitation")
    ap.add_argument("--backend", default="gt:cpu_ifirst")
    ap.add_argument("--step-minutes", type=float, default=60.0)
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    snow, _ = snow_parameters(args.control_file)
    glac = ice_firn_parameters(args.control_file)
    debris = np.ones(d("start", "[zone_grid]")[0].shape)
    scaling = 1.0
    if d.has("start", "debris_on_glaciers"):  # WaSiM reads the scaling factor only when the grid exists
        grid = d("start", "debris_on_glaciers")[0]
        debris = np.where(grid > -9999 + 1, grid, 1.0)
        scaling = debris_scaling(args.control_file)
    params = GlacierMeltParameters(t0=snow.t0, debris_scaling=scaling)
    active = d("start", "[zone_grid]")[0] > -9999
    lu_alb, unobserved = landuse_albedo(d, active, args.pre, args.first_hour)
    elevation = d("start", "[elevation_model]")[0]
    nsteps, nrows, ncols = d(args.pre, "snowstorage_solid").shape
    dom = gtx.domain({I: nrows, J: ncols})
    fld = lambda a: gtx.as_field(dom, np.ascontiguousarray(np.broadcast_to(a, (nrows, ncols)), dtype=np.float64))
    start = dt.datetime.fromisoformat(args.start)

    worst = {}
    def check(name, value, ref, cells):
        worst[name] = max(worst.get(name, 0.0), difference(value, ref, cells))

    counts = dict(ice_melt=0, firn_melt=0)
    for t in range(nsteps):
        pre = lambda name: d(args.pre, name)[t]
        post = lambda name: d("after_snow1", name)[t]
        ice = pre("GlacierizedCells")
        glacier = active & (ice > 0)
        no_snow = ((pre("lake_codes") > 0) & (pre("ponding_storage") > 0)).astype(float)
        base = dict(temperature=pre("corr_temperature"), deposition_index=pre("deposition_index"),
                    air_humidity=pre("air_humidity"), cloudiness=pre("cloudiness"),
                    radiation=pre("AspectCorrectedGlobalRadiation"), wind_speed=pre("wind_speed"), elevation=elevation,
                    landuse_albedo=lu_alb[t], no_snow=no_snow)
        first = snow_energy_balance(dict(base, precipitation=pre(args.input), active=active.astype(float),
                                         snowstorage_solid=pre("snowstorage_solid"),
                                         snowstorage_liquid=pre("snowstorage_liquid"), snowtemp=pre("snowtemp"),
                                         surftemp=pre("surftemp"), snow_age=pre("snow_age")),
                                    snow, args.step_minutes, backend=args.backend)
        for k in FIRST_CALL:
            check(k, first[k], post(k), active)
        gl = snow_energy_balance(dict(base, precipitation=first["precipitation"], active=glacier.astype(float),
                                      snowstorage_solid=pre(GLACIER_SNOW["snowstorage_solid"]),
                                      snowstorage_liquid=pre(GLACIER_SNOW["snowstorage_liquid"]),
                                      snowtemp=pre(GLACIER_SNOW["snowtemp"]), surftemp=pre(GLACIER_SNOW["surftemp"]),
                                      snow_age=pre("snow_age"), snow_fraction=first["snow_fraction"],
                                      rain_fraction=first["rain_fraction"], albedo_used=first["albedo_used"]),
                                 snow, args.step_minutes, backend=args.backend, on_glacier=True)
        for k, name in GLACIER_SNOW.items():
            check(name, gl[k], post(name), glacier)

        # the date of the step's end time: at 24:00 WaSiM already uses the next day (checked against the dumps)
        day = start + dt.timedelta(hours=t + args.first_hour)
        a_ice = radiation_coefficient(glac["ice_min"], glac["ice_max"], day.day, day.month)
        a_snow = radiation_coefficient(glac["snow_min"], glac["snow_max"], day.day, day.month)
        stack = pre("firn_WE_stack")
        layers = [fld(stack[..., k]) for k in range(1, stack.shape[-1])]
        qice, qfirn, new_layers, total = glacier_melt(
            fld(pre("corr_temperature")), fld(gl["snowstorage_solid"]), fld(gl["melt"]), layers, fld(stack[..., 0]),
            fld(pre("global_radiation")), fld(glac["mf"]), fld(a_ice), fld(a_snow), fld(debris),
            gtx.as_field(dom, glacier), params, args.step_minutes, dom,
            backend=None if args.backend in ("embedded", "gt:cpu_ifirst") else args.backend)
        check("melt_from_ice", qice.asnumpy(), post("melt_from_ice"), glacier)
        check("melt_from_firn", qfirn.asnumpy(), post("melt_from_firn"), glacier)
        ref_stack = post("firn_WE_stack")
        check("firn_WE_stack[0]", total.asnumpy(), ref_stack[..., 0], glacier)
        for k, layer in enumerate(new_layers, start=1):
            check(f"firn_WE_stack[{k}]", layer.asnumpy(), ref_stack[..., k], glacier)
        counts["ice_melt"] += int((post("melt_from_ice")[glacier] > 0).sum())
        counts["firn_melt"] += int((post("melt_from_firn")[glacier] > 0).sum())

    print(f"{active.sum()} active cells, {nsteps} steps, backend {args.backend}, {counts}, "
          f"{unobserved} cell-steps with an unobserved land-use albedo")
    print("largest absolute difference to WaSiM (inf: NaN or nodata patterns differ):")
    for k, v in worst.items():
        print(f"  {k}: {v:.3g}")
    return 0 if max(worst.values()) == 0.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
