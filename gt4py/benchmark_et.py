#!/usr/bin/env python3
"""Warm per-call timing of wasim-gt's radiation correction, potential ET and vegetation parameters.

Takes the inputs of one daylight step from a WaSiM grid dump (as potential_et.py does) and calls
each process repeatedly, after warm-up calls that include any compilation, so the numbers are
compute plus Python and host-device overhead, without compilation and without dump I/O.
Compare with WaSiM's per-module times from WASIM_TIMING (wasim-eth) on the same case and a
single thread: see ../validation/performance.md.

usage: benchmark_et.py DUMP_DIR CONTROL_INTERN START [--step T] [--repeat N] [--backend embedded|gtfn_cpu]
       [--cellsize M]
"""
import argparse
import datetime
import re
import time

import numpy as np
import gt4py.next as gtx

from wasim_gt.fields import I, J
from wasim_gt.processes import potential_et as pet
from wasim_gt.processes import vegetation as veg
from wasim_gt.processes.radiation_correction import radiation_correction, shadow, sub_intervals
from wasim_gt.processes.sun import Site

from potential_et import multilayer_table, numbers, section
from vegetation import landuse_table
from wasim_dump import Dump


def timed(fn, repeat, warmup=2):
    for _ in range(warmup):
        fn()
    t = time.perf_counter()
    for _ in range(repeat):
        fn()
    return 1000.0 * (time.perf_counter() - t) / repeat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_intern")
    ap.add_argument("start", help="first step of the dump, e.g. 2021-06-13T01")
    ap.add_argument("--step", type=int, default=11, help="dump step to use (default 11: midday of the first day)")
    ap.add_argument("--repeat", type=int, default=20)
    ap.add_argument("--backend", default="embedded", choices=["embedded", "gtfn_cpu"])
    ap.add_argument("--cellsize", type=float)
    args = ap.parse_args()

    text = open(args.control_intern, encoding="latin-1").read()
    coords, rc, et, lake = (section(text, s) for s in ("coordinates", "radiation_correction", "evapotranspiration",
                                                         "lake_model"))
    step = float(et[1])
    site = Site(float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3]),
                float(section(text, "output_interval")[3].split()[0]) / 60.0)
    amp = numbers(et[30])
    lakes_on = int(lake[0]) != 0
    params = pet.Parameters(cr=tuple(numbers(et[25])), dt_koeff=numbers(et[28])[0], dtmax=tuple(numbers(et[29])[:12]),
                            ampl_day=amp[0], ampl_annual=amp[1], albedo_open_water=float(lake[2]) if lakes_on else 0.0,
                            z0_open_water=float(lake[3]) if lakes_on else 0.0, lakes=lakes_on)
    d = Dump(args.dump_dir)
    p0, p1, t = "after_meteo", "after_expokorr", args.step
    active = d(p0, "[zone_grid]")[0] > -9999
    dem = np.array(d(p0, "[elevation_model]")[0])
    nrows, ncols = dem.shape
    cellsize = args.cellsize
    if cellsize is None:
        cellsize = float(re.search(r"cellsize\s+([0-9.]+)", open(section(text, "elevation_model")[0], "rb")
                                   .read(400).decode("latin-1")).group(1))
    dom = gtx.domain({I: nrows, J: ncols})
    raw = lambda a: gtx.as_field(dom, np.array(a))  # noqa: E731
    fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))  # noqa: E731
    date, hour = args.start.split("T")
    now = datetime.datetime.strptime(date, "%Y-%m-%d") + datetime.timedelta(hours=int(hour), minutes=step * t)

    slope, aspect = fld(d(p0, "slope_angle")[t]), fld(d(p0, "slope_aspect")[t])
    temp0, ssd0 = fld(d(p0, "temperature")[t]), fld(d(p0, "sunshine_duration")[t])
    rc_call = lambda: radiation_correction(  # noqa: E731
        site, now.day, now.month, now.hour, now.minute, step, int(rc[2]), float(rc[5]), int(rc[11]), dem, slope, aspect,
        cellsize, temp0, ssd0, active, dom, backend=args.backend)
    subs, _ = sub_intervals(site, now.day, now.month, now.hour, now.minute, step, int(rc[11]))
    shade_call = lambda: shadow(dem, active, subs[0].azimuth, subs[0].altitude, cellsize)  # noqa: E731

    ml = multilayer_table(text)
    luse0 = d(p1, "landuse")[0]
    kl, ls = np.full(luse0.shape, 0.3), np.full(luse0.shape, 10.0)
    for code, (k, s) in ml.items():
        kl[luse0 == code], ls[luse0 == code] = k, s
    layers = [pet.Layer(landuse=raw(d(p1, f"landuse{n}")[t]), lai=raw(d(p1, f"leaf_area_index{n}")[t]),
                        coverage=fld(d(p1, f"vegetation_coverage_degree{n}")[t]),
                        rsc=fld(d(p1, f"SurfaceCanopyResistance{n}")[t]), rsi=fld(d(p1, f"SurfaceIntercepResistance{n}")[t]),
                        z0=fld(d(p1, f"RoughnessLength{n}")[t])) for n in (1, 2)]
    tk = fld(d(p1, "corr_temperature")[t])
    ins = [fld(d(p1, g)[t]) for g in ("global_radiation", "sunshine_duration", "radiation_correction", "air_humidity",
                                       "wind_speed", "albedo")]
    rest = dict(dhm=fld(dem), rse=fld(d(p1, "SurfaceEvaporationResistance")[t]), layers=layers, k_lambert=fld(kl),
                lai_scale=fld(ls), lake=raw(d(p1, "lake_codes")[t]), pond=raw(d(p1, "ponding_storage")[t]))
    ssd_mean = np.array(d(p1, "SSDTemp")[t])

    def et_call():
        state = pet.State(ssd_mean=raw(ssd_mean))
        return pet.potential_et(site, now.day, now.month, now.hour, step, params, tk, tk, *ins,
                                state=state, active=active, domain=dom, backend=args.backend, **rest)

    table = landuse_table(text)
    luse = [np.array(d(p0, f"landuse{n}")[0]) for n in (1, 2)]
    veg_call = lambda: veg.vegetation_parameters(table, luse, dem, active, now.day, now.month, now.year)  # noqa: E731

    print(f"{active.sum()} active cells ({nrows} x {ncols}), step {t} = {now:%Y-%m-%d %H:%M}, backend {args.backend}, "
          f"{args.repeat} calls after warm-up")
    print(f"  radiation correction:   {timed(rc_call, args.repeat):8.2f} ms per call")
    print(f"    of which shading walk: {timed(shade_call, args.repeat):7.2f} ms (host Python)")
    print(f"  potential ET:            {timed(et_call, args.repeat):8.2f} ms per call")
    print(f"  vegetation parameters:   {timed(veg_call, args.repeat):8.2f} ms per call (host numpy, once a day)")


if __name__ == "__main__":
    raise SystemExit(main())
