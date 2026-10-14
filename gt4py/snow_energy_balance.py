#!/usr/bin/env python3
"""Validate wasim-gt's snow method 7 (wasim_gt.processes.snow_energy_balance and snow_slides) against
WaSiM grid dumps.

Part 1 (run1: accumulation, energy balance, melt, albedo) is driven by the state at the dump point
before the snow call and compared with after_snow1. With interception before snow (methods 11/12, as
in Doubs) that point is after_interception and the snow input is throughfall0; with methods 1/2 it is
after_irrigation (or after_etp) and the input is precipitation. Part 2 (run2: gravitational slides) runs
after the sub-step loop of the soil and groundwater models and before their post-runs, so it is driven by
the last sub step's after_gwflow and compared with after_snow2. Every step, every active cell, exact.

Grids to dump (WASIM_DUMP_GRIDS), at points start,after_meteo,after_interception,after_snow1,after_gwflow,after_snow2:
  corr_temperature,throughfall0,precipitation,wind_speed,air_humidity,AspectCorrectedGlobalRadiation,
  cloudiness,deposition_index,lake_codes,ponding_storage,[elevation_model],[zone_grid],
  snowstorage_solid,snowstorage_liquid,snow_storage_total,snowtemp,surftemp,snow_age,days_snow,albedo,
  snow_fraction,rain_fraction,snowcover_outflow,melt_from_snowcover,snow_evaporation,
  slope_angle,slidefraction1,slidefraction2,slidefraction3,slidefraction4,max_deposition,
  elev_order_row,elev_order_col,input_mass,mobile_mass,deposition

The golden Doubs setup never reaches the slide threshold (35 degrees; its steepest cell has 34).
Lowering it in the control file ("minimum slope for creating slides", e.g. to 20) exercises the slides.

usage: snow_energy_balance.py DUMP_DIR CONTROL_FILE [--pre after_interception] [--input throughfall0]
                              [--ice GlacierizedCells] [--no-lakes] [--backend gt:cpu_ifirst]

  Doubs (interception method 11, before snow):  defaults
  Rhone (interception method 1, after snow; glaciers):  --pre after_irrigation --input precipitation --ice GlacierizedCells
  With --ice, also dump GlacierizedCells and the snow on glaciers (TotalSnowStorage_onGlacier,
  LiquidSnowStorage_onGlacier, SolidSnowStorage_onGlacier) at after_gwflow and after_snow2.
"""
import argparse
import re
import time

import numpy as np

from wasim_gt.processes.snow_energy_balance import SnowParameters, snow_energy_balance
from wasim_gt.processes.snow_slides import snow_slides

from wasim_dump import Dump


GLACIER_SNOW = ("TotalSnowStorage_onGlacier", "LiquidSnowStorage_onGlacier", "SolidSnowStorage_onGlacier")


def control_section(control_file, name):
    text = open(control_file, encoding="utf-8", errors="replace").read()
    return text.split(f"[{name}]", 1)[1].split("\n[", 1)[0]


def snow_parameters(control_file):
    """Numbers of [snow_model] in the order of schneemodellklasse's reader (method 7 and its relatives)."""
    lines = []
    for line in control_section(control_file, "snow_model").splitlines():
        line = line.split("#", 1)[0].strip()
        if line and re.fullmatch(r"[-+0-9.eE\s]+", line):
            lines.append([float(v) for v in line.split()])
    method = int(lines[0][0])
    if method != 7:
        raise SystemExit(f"snow method {method}: this check covers method 7 only")
    # method, transzone, T0R, T0, CWH, CRFR, C0[12], C1[12], C2[12], z0, RMFMIN, RMFMAX, ALBmin, ALBmax,
    # then (after the file names) max deposition slope, deposition scaling, min slide slope,
    # slide mass fraction, LWINcorr, LWOUTcorr, energy balance scaling
    tail = [v[0] for v in lines[14:]]
    if len(tail) != 7:
        raise SystemExit(f"unexpected [snow_model] layout: {lines}")
    p = SnowParameters(transition=lines[1][0], t0r=lines[2][0], t0=lines[3][0], cwh=lines[4][0],
                       albedo_min=lines[12][0], albedo_max=lines[13][0],
                       lwin_corr=tail[4], lwout_corr=tail[5], enbal_scaling=tail[6])
    return p, dict(min_slide_slope=tail[2], slide_mass_fraction=tail[3])


def landuse_albedo(d, active, pre, first_hour):
    """Land-use albedo per step, exactly as the snow model gets it, and the number of cells where it
    could not be observed.

    run1 asks the land-use table for the albedo of each cell (LUSEtabletype::GetAlbedo with the day
    and the elevation) twice with the same arguments: before the step if SWE <= 2.5 mm and after it if
    SWE <= 5 mm. Its interpolation between the table's fix points rounds (a constant 0.15 becomes
    0.15 +- 1 ulp, depending on day and elevation), and porting it belongs to the vegetation
    parameters (int_tvar_params; issue #2 takes them as inputs too). So this check reads the value
    back: on every day, from the after_snow1 albedo of a step where the cell ended with SWE <= 5 mm.
    Cells never snow-free on a day fall back to the value int_tvar_params wrote to the albedo grid at
    the day's first after_meteo, which can differ in the last bit; they only matter if their SWE
    before a step is <= 2.5 mm.
    """
    meteo, snow1 = d("after_meteo", "albedo"), d("after_snow1", "albedo")
    swe1 = d("after_snow1", "snow_storage_total")
    nsteps = meteo.shape[0]
    day = (np.arange(nsteps) + first_hour) // 24
    out, unobserved = [None] * nsteps, 0
    for k in np.unique(day):
        steps = np.flatnonzero(day == k)
        value = meteo[steps[0]].copy()
        seen = np.zeros(active.shape, bool)
        for t in steps:
            free = active & (swe1[t] <= 5.0)
            if (seen & free & (value != snow1[t])).any():
                raise SystemExit(f"land-use albedo not constant within day {k}")
            value[free], seen = snow1[t][free], seen | free
        for t in steps:
            unobserved += int((active & ~seen & (d(pre, "snow_storage_total")[t] <= 2.5)).sum())
            out[t] = value
    return out, unobserved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_file", help="the filled control file of the run (run_case.py writes it to the output directory)")
    ap.add_argument("--pre", default="after_interception", help="dump point just before the snow model, part 1")
    ap.add_argument("--input", default="throughfall0", help="the snow model's input grid (throughfall0 or precipitation)")
    ap.add_argument("--backend", default="gt:cpu_ifirst")
    ap.add_argument("--no-lakes", action="store_true", help="the lake model is off")
    ap.add_argument("--ice", default=None, help="glacier grid (GlacierizedCells or ice_firn): cells with ice are not "
                    "ported (issue #5) and are left out of the comparison")
    ap.add_argument("--step-minutes", type=float, default=60.0)
    ap.add_argument("--first-hour", type=int, default=1, help="hour of day of the first step (to group steps by day)")
    args = ap.parse_args()

    d = Dump(args.dump_dir)
    params, slide = snow_parameters(args.control_file)
    print(params, slide)
    pre, post = args.pre, "after_snow1"
    ssnow0 = d(pre, "snowstorage_solid")
    nsteps = ssnow0.shape[0]
    active = d("start", "[zone_grid]")[0] > -9999
    lu_alb, unobserved = landuse_albedo(d, active, pre, args.first_hour)
    elevation = d("start", "[elevation_model]")[0]
    check = active.copy()  # cells compared with WaSiM
    if args.ice:
        ice = d("start", args.ice)[0]
        check &= ~(ice > 0)
        print(f"{(active & (ice > 0)).sum()} glacierized cells left out")

    compare1 = {"precipitation": args.input, "snow_fraction": "snow_fraction", "rain_fraction": "rain_fraction",
                "snowstorage_solid": "snowstorage_solid", "snowstorage_liquid": "snowstorage_liquid",
                "snow_storage_total": "snow_storage_total", "snowtemp": "snowtemp", "surftemp": "surftemp",
                "snow_age": "snow_age", "albedo": "albedo", "snowcover_outflow": "snowcover_outflow",
                "melt_from_snowcover": "melt_from_snowcover", "days_snow": "days_snow",
                "snow_evaporation": "snow_evaporation"}
    worst1 = dict.fromkeys(compare1, 0.0)
    t_run = 0.0
    for t in range(nsteps):
        g = lambda name, point=pre: d(point, name)[t]
        no_snow = np.zeros(active.shape)
        if not args.no_lakes:
            no_snow = ((g("lake_codes") > 0) & (g("ponding_storage") > 0)).astype(float)
        inputs = dict(temperature=g("corr_temperature"), precipitation=g(args.input),
                      deposition_index=g("deposition_index"), air_humidity=g("air_humidity"),
                      cloudiness=g("cloudiness"), radiation=g("AspectCorrectedGlobalRadiation"),
                      wind_speed=g("wind_speed"), elevation=elevation, landuse_albedo=lu_alb[t], no_snow=no_snow,
                      active=active.astype(float), snowstorage_solid=g("snowstorage_solid"),
                      snowstorage_liquid=g("snowstorage_liquid"), snowtemp=g("snowtemp"),
                      surftemp=g("surftemp"), snow_age=g("snow_age"))
        start = time.perf_counter()
        res = snow_energy_balance(inputs, params, args.step_minutes, backend=args.backend)
        t_run += time.perf_counter() - start
        for k, name in compare1.items():
            diff = np.abs(res[k] - g(name, post))[check]
            worst1[k] = max(worst1[k], float(diff.max()))

    worst2 = dict.fromkeys(["snow_storage_total", "snowstorage_liquid", "snowstorage_solid",
                            "input_mass", "mobile_mass", "deposition"], 0.0)
    zone = d("start", "[zone_grid]")[0]
    fractions = [d("start", f"slidefraction{k}")[0] for k in range(1, 5)]
    static = dict(active=active, zone=zone, slope=d("start", "slope_angle")[0], slide_fractions=fractions,
                  deposition_max=d("start", "max_deposition")[0], elev_order_row=d("start", "elev_order_row")[0],
                  elev_order_col=d("start", "elev_order_col")[0], **slide)
    if args.ice:  # the slides are ported for glacierized cells too: compare all cells, and the snow on glaciers
        static["ice"] = d("start", args.ice)[0]
        for k in GLACIER_SNOW:
            worst2[k] = 0.0
    n_slides = 0
    nsub = d("after_gwflow", "snow_storage_total").shape[0] // nsteps  # records per step: one per sub step
    for t in range(nsteps):
        g = lambda name, point: d(point, name)[t]
        last = lambda name: d("after_gwflow", name)[(t + 1) * nsub - 1]
        glacier = tuple(last(k) for k in GLACIER_SNOW) if args.ice else None
        states, slide_grids, glacier_states = snow_slides(
            last("snow_storage_total"), last("snowstorage_liquid"), last("snowstorage_solid"), glacier=glacier, **static)
        n_slides += int((slide_grids[0][active] > 0).sum())
        results = (*states, *slide_grids, *(glacier_states or ()))
        for k, v in zip(worst2, results):
            ref = g(k, "after_snow2")
            cells = active & (ref > -9999) if k in GLACIER_SNOW else (active if args.ice else check)
            if args.ice and k in ("snow_storage_total", "snowstorage_liquid", "snowstorage_solid"):
                cells = active & ~(static["ice"] >= 1.0)  # fully glacierized cells have no snow off the glacier
            worst2[k] = max(worst2[k], float(np.abs(v - ref)[cells].max()))

    print(f"{active.sum()} active cells, {nsteps} steps, backend {args.backend}, "
          f"{1000 * t_run / nsteps:.1f} ms per step for part 1 (first step includes compilation)")
    print(f"part 1 ({pre} -> after_snow1), largest absolute difference to WaSiM over all steps and active cells"
          f" ({unobserved} cell-steps used an unobserved land-use albedo):")
    for k, v in worst1.items():
        print(f"  {compare1[k]}: {v:.3g}")
    print(f"part 2, slides (after_gwflow, last sub step -> after_snow2), {n_slides} cell-steps with a slide:")
    for k, v in worst2.items():
        print(f"  {k}: {v:.3g}")
    return 0 if max(*worst1.values(), *worst2.values()) == 0.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
