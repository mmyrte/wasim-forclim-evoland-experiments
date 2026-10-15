#!/usr/bin/env python3
"""Validate wasim-gt's table-based phenology (wasim_gt.processes.vegetation) against WaSiM grid dumps.

WaSiM sets the vegetation parameters at the start date and recomputes them on the last step of
each day (defaultini::int_tvar_params), before the dump point after_meteo. For every step this
script compares the dumped parameters with wasim-gt's: at the first step with the initialisation,
on day-end steps with the update for that day (the initialisation uses the elevation model as read
from its file, before the lake model adjusts it), and on other steps it checks that WaSiM left them
unchanged. The albedo is compared on day-end steps only, because the snow model changes it in
between.

The land-use table is read from the expanded control file (comments removed) the way WaSiM
reads it (LUSEtabletype::initialize): a parameter is the first occurrence of its name in the
record, missing rs_interception and rs_evaporation get WaSiM's defaults (20 and 200), and the
code at the start of a record is not used for the lookup: WaSiM replaces the codes in the
land-use grids by the record's position (defltuzr.cc:1825), so the dumped grids hold positions.

usage: vegetation.py DUMP_DIR CONTROL_INTERN START [--layers N]
   START: first step as in WaSiM's input data, e.g. 2021-01-01T01 (hourly) or 2021-01-01T24 (daily)
"""
import argparse
import datetime
import re

import numpy as np

from wasim_gt.processes import vegetation as veg

from wasim_dump import Dump

GRIDS = {"lai": "leaf_area_index", "vcf": "vegetation_coverage_degree", "rsc": "SurfaceCanopyResistance",
         "rsi": "SurfaceIntercepResistance", "z0": "RoughnessLength", "root_depth": "root_depth"}
KEYS = {"albedo": "Albedo", "rsc": "rsc", "lai": "LAI", "z0": "Z0", "vcf": "VCF", "root_depth": "RootDepth",
        "alt_dep": "AltDep"}


def section(text, name):
    body = text.split(f"[{name}]", 1)[1].split("\n[", 1)[0]
    lines = [ln.split("#", 1)[0].strip() for ln in body.splitlines()]
    return [ln for ln in lines if ln]


def initial_dem(text, shape):
    """The elevation model as WaSiM reads it at the start (before the lake model lowers lake
    cells): an ASCII grid, or WaSiM's binary grid (12 float32 header values, rows north first)."""
    path = section(text, "elevation_model")[0]
    try:
        return np.loadtxt(path, skiprows=6)
    except ValueError:
        raw = np.fromfile(path, dtype="<f4")
        return raw[12:12 + shape[0] * shape[1]].reshape(shape).astype(np.float64)


def wasim_param(record, key, n, cast=float):
    """LUSEtabletype::stripParameters: the n values between '=' and ';' after the first `key`."""
    buf = record[record.find(key):]
    text = buf[buf.find("=") + 1:buf.find(";")]
    values = [cast(x) for x in re.split(r"[\t ,;]+", text.strip()) if x]
    if len(values) < n:
        raise SystemExit(f"{key}: {len(values)} values, {n} expected")
    return tuple(values[:n])


def landuse_table(text):
    """{position: LandUse} for the VariableDayCount records of [landuse_table], in file order."""
    body = text.split("[landuse_table]", 1)[1]
    count = int(body.split()[0])
    # as fscanf("%d") + fgets: the records start on the line after the count
    lines = body[body.index(body.split()[0]) + len(body.split()[0]):].splitlines(keepends=True)[1:]
    table, i = {}, 0
    for pos in range(1, count + 1):
        record = lines[i]
        i += 1
        if not ("VariableDayCount" in record and "method" in record and "{" in record):
            raise SystemExit(f"land use {pos}: only VariableDayCount records are ported")
        while "}" not in record:
            record += lines[i]
            i += 1
        buf = record[record.find("JulDays"):]
        days = tuple(int(x) for x in re.split(r"[\t ,;]+", buf[buf.find("=") + 1:buf.find(";")].strip()) if x)
        n = len(days)
        kw = {k: wasim_param(record, key, n) for k, key in KEYS.items()}
        if "rs_interception" in record:
            kw["rsi"] = wasim_param(record, "rs_interception", n)
        if "rs_evaporation" in record:
            kw["rse"] = wasim_param(record, "rs_evaporation", n)
        table[pos] = veg.LandUse(jul_days=days, **kw)
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_intern")
    ap.add_argument("start", help="first step as in WaSiM's input data, e.g. 2021-01-01T01")
    ap.add_argument("--layers", type=int, default=2)
    args = ap.parse_args()

    text = open(args.control_intern, encoding="latin-1").read()
    table = landuse_table(text)
    step = float(section(text, "evapotranspiration")[1])
    hour_offset = float(section(text, "output_interval")[3].split()[0]) / 60.0

    d = Dump(args.dump_dir)
    p = "after_meteo"
    active = d(p, "[zone_grid]")[0] > -9999
    dhm = d(p, "[elevation_model]")[0]
    luse = [d(p, f"landuse{n}")[0] for n in range(1, args.layers + 1)]
    # the initialisation uses the elevation model as read, before the lake model lowers lake cells
    # (lake initialisation method 2); the daily updates use the lowered one
    dem_file = initial_dem(text, dhm.shape)
    nsteps = d(p, "leaf_area_index1").shape[0]
    grids = {(k, n): f"{g}{n}" for k, g in GRIDS.items() for n in range(1, args.layers + 1)}
    grids[("rse", 0)] = "SurfaceEvaporationResistance"
    grids[("albedo", 0)] = "albedo"
    data = {g: d(p, g) for g in grids.values()}

    date, hour = args.start.split("T")
    start = datetime.datetime.strptime(date, "%Y-%m-%d") + datetime.timedelta(hours=int(hour))
    hour24 = int(hour) == 24
    worst = {g: 0.0 for g in grids.values()}
    unchanged = {g: 0 for g in grids.values()}
    updates = 0
    first = None
    for t in range(nsteps):
        now = start + datetime.timedelta(minutes=step * t)
        day, month, year, hr = now.day, now.month, now.year, now.hour
        if hour24:
            prev = now - datetime.timedelta(days=1)
            day, month, year, hr = prev.day, prev.month, prev.year, 24
        if first is None:
            first = (day, month, year)
        if veg.is_day_end(hr, 0, hour_offset):
            res = veg.vegetation_parameters(table, luse, dhm, active, day, month, year)
            updates += 1
            check_albedo = True
        elif t == 0:
            res = veg.vegetation_parameters(table, luse, dem_file, active, *first, initial=True)
            check_albedo = False
        else:
            for (key, n), g in grids.items():
                if key != "albedo" and not np.array_equal(data[g][t][active], data[g][t - 1][active]):
                    unchanged[g] += 1
            continue
        for (key, n), g in grids.items():
            if key == "albedo" and not check_albedo:
                continue
            v = res[key] if n == 0 else res["layers"][n - 1][key]
            worst[g] = max(worst[g], float(np.abs(v - data[g][t])[active].max()))

    print(f"{active.sum()} active cells, {nsteps} steps of {step:g} min from {args.start}, "
          f"{updates} day-end updates, {len(table)} land uses")
    print("largest absolute difference to WaSiM (initialisation and day-end steps); steps where WaSiM changed "
          "the grid between updates:")
    for g in grids.values():
        print(f"  {g}: {worst[g]:.3g}; {unchanged[g]}")
    return 0 if max(worst.values()) == 0.0 and not any(unchanged.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
