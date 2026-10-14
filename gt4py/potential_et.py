#!/usr/bin/env python3
"""Validate wasim-gt's radiation correction and potential ET against WaSiM grid dumps.

wasim_gt.processes.radiation_correction (expokorrklasse::run) is driven by the states at dump
point after_meteo and compared with after_expokorr; wasim_gt.processes.potential_et (Penman-
Monteith, verdunstungsmodellklasse::potentiell) is driven by after_expokorr and compared with
after_etp. Every step and every active cell. Parameters come from the expanded control file of
the run (control.wasim.intern, written by WaSiM next to the control file).

The scalar state of ET (ssdFlag, ssdCnt: counting the bright steps of a day) is carried by this
script from the first step, so the run must start at the start of the WaSiM run. Grids are taken
from WaSiM before each call, so each step is checked on its own.

Bit-identity needs the C library's exp, log, pow, sin and cos:
- gtfn_cpu calls them directly;
- embedded (numpy) uses its own vectorised exp, log and pow on AVX-512 machines, which differ in
  the last bit for 2-5% of arguments. NPY_DISABLE_CPU_FEATURES=X86_V4 (numpy 2.5; the AVX-512
  group) gives the C library's results; this script warns if it is not set.

usage: potential_et.py DUMP_DIR CONTROL_INTERN START [--steps N] [--cellsize M] [--backend embedded|gtfn_cpu]
   START: date and hour of the first step as in WaSiM's input data, e.g. 2021-01-01T01 (hourly;
          midnight is hour 0 of the next day) or 2021-01-01T24 (daily)
"""
import argparse
import datetime
import os
import re
import time

import numpy as np
import gt4py.next as gtx

from wasim_gt.fields import I, J
from wasim_gt.processes.sun import Site
from wasim_gt.processes.radiation_correction import radiation_correction
from wasim_gt.processes import potential_et as pet

from wasim_dump import Dump

GRIDS = [
    ("etp", 0, "pot_transpiration0"), ("etp", 1, "pot_transpiration1"), ("etp", 2, "pot_transpiration2"),
    ("eip", 0, "potential_InterceptionEvap0"), ("eip", 1, "potential_InterceptionEvap1"),
    ("eip", 2, "potential_InterceptionEvap2"), ("evap", None, "potential_evaporation"),
    ("rgex", None, "AspectCorrectedGlobalRadiation"), ("cloudiness", None, "cloudiness"),
    ("global_radiation", None, "global_radiation"), ("ssd_mean", None, "SSDTemp"),
    ("temperature", None, "corr_temperature"),
]


class WaSiMTime:
    """A step labelled with hour 24 of its day, as in WaSiM's daily inputs."""
    def __init__(self, day, hour):
        self.day, self.month, self.hour, self.minute = day.day, day.month, hour, 0


def section(text, name):
    """Non-empty lines of a control-file section, comments removed."""
    body = text.split(f"[{name}]", 1)[1].split("\n[", 1)[0]
    lines = [ln.split("#", 1)[0].strip() for ln in body.splitlines()]
    return [ln for ln in lines if ln]


def numbers(line):
    return [float(x) for x in line.split()]


def multilayer_table(text):
    """{position: (k_extinct, LAI_scale)} from [multilayer_landuse]. WaSiM replaces the codes in the
    land-use grid by the record's position (defltuzr.cc:1878), so the dumped grid holds positions."""
    body = text.split("[multilayer_landuse]", 1)[1].split("\n[", 1)[0]
    out = {}
    for pos, m in enumerate(re.finditer(r"^\s*(\d+)\s+.*?\{(.*?)\}", body, re.M | re.S), start=1):
        k = re.search(r"k_extinct\s*=\s*([-+0-9.eE]+)", m.group(2))
        s = re.search(r"LAI_scale\s*=\s*([-+0-9.eE]+)", m.group(2))
        out[pos] = (float(k.group(1)) if k else 0.3, float(s.group(1)) if s else 10.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_dir")
    ap.add_argument("control_intern", help="control.wasim.intern of the run (WaSiM writes it next to the control file)")
    ap.add_argument("start", help="first step as WaSiM counts it, e.g. 2021-01-01T01")
    ap.add_argument("--steps", type=int, default=0, help="number of steps (default: all in the dump)")
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--cellsize", type=float, help="cell size in m (default: from the elevation model's ASCII header)")
    ap.add_argument("--backend", default="embedded", choices=["embedded", "gtfn_cpu"])
    args = ap.parse_args()
    if args.backend == "embedded" and "NPY_DISABLE_CPU_FEATURES" not in os.environ:
        print("warning: NPY_DISABLE_CPU_FEATURES not set; numpy's vectorised exp/log/pow may differ from the C library")

    text = open(args.control_intern, encoding="latin-1").read()
    coords = section(text, "coordinates")
    rc = section(text, "radiation_correction")
    et = section(text, "evapotranspiration")
    lake = section(text, "lake_model")
    step = float(et[1])
    if int(et[2]) != 1:
        raise SystemExit("only ET method 1 (Penman-Monteith) is ported")
    out_iv = section(text, "output_interval")
    site = Site(float(coords[0]), float(coords[1]), float(coords[2]), float(coords[3]),
                float(out_iv[3].split()[0]) / 60.0)
    rc_method, spreiz, count, split = int(rc[2]), float(rc[5]), int(rc[10]), int(rc[11])
    if count != 1:
        raise SystemExit("the interval counter of the radiation correction must be 1")
    lakes_on = int(lake[0]) != 0
    amp = numbers(et[30])
    params = pet.Parameters(cr=tuple(numbers(et[25])), dt_koeff=numbers(et[28])[0],
                            dtmax=tuple(numbers(et[29])[:12]), ampl_day=amp[0], ampl_annual=amp[1],
                            albedo_open_water=float(lake[2]) if lakes_on else 0.0,
                            z0_open_water=float(lake[3]) if lakes_on else 0.0, lakes=lakes_on)
    ml = multilayer_table(text) if args.layers > 1 else {}

    d = Dump(args.dump_dir)
    p0, p1, p2 = "after_meteo", "after_expokorr", "after_etp"
    temp0 = d(p0, "temperature")
    nsteps, nrows, ncols = temp0.shape
    if args.steps:
        nsteps = min(nsteps, args.steps)
    active = d(p0, "[zone_grid]")[0] > -9999
    dem = d(p0, "[elevation_model]")[0]
    cellsize = args.cellsize
    if cellsize is None:  # from an ASCII grid header
        m = re.search(r"cellsize\s+([0-9.]+)", open(section(text, "elevation_model")[0], "rb").read(400).decode("latin-1"))
        if not m:
            raise SystemExit("the elevation model is not an ASCII grid: give --cellsize")
        cellsize = float(m.group(1))
    dom = gtx.domain({I: nrows, J: ncols})
    raw = lambda a: gtx.as_field(dom, a)  # noqa: E731
    fld = lambda a: gtx.as_field(dom, np.where(active, a, 0.0))  # noqa: E731

    luse0 = d(p1, "landuse")[0]
    kl = np.full(luse0.shape, 0.3)
    ls = np.full(luse0.shape, 10.0)
    for code, (k, s) in ml.items():
        kl[luse0 == code] = k
        ls[luse0 == code] = s
    state = pet.State(ssd_mean=None)

    worst_rc = {g: 0.0 for g in ("radiation_correction", "Ieff_max", "corr_temperature", "Shadow")}
    worst = {}
    t_rc = t_et = 0.0
    date, hour = args.start.split("T")
    start = datetime.datetime.strptime(date, "%Y-%m-%d") + datetime.timedelta(hours=int(hour))
    hour24 = int(hour) == 24  # daily data: the step of a day is labelled hour 24 of that day
    for t in range(nsteps):
        now = start + datetime.timedelta(minutes=step * t)
        if hour24:
            now = WaSiMTime(now - datetime.timedelta(days=1), 24)
        # --- radiation correction
        t1 = time.perf_counter()
        r = radiation_correction(site, now.day, now.month, now.hour, now.minute, step, rc_method, spreiz, split,
                                 dem, fld(d(p0, "slope_angle")[t]), fld(d(p0, "slope_aspect")[t]), cellsize,
                                 fld(temp0[t]), fld(d(p0, "sunshine_duration")[t]), active, dom,
                                 backend=args.backend)
        t_rc += time.perf_counter() - t1
        for k, g in (("exkorr", "radiation_correction"), ("ieffmax", "Ieff_max"), ("tkorr", "corr_temperature"),
                     ("shadow", "Shadow")):
            v = r[k] if isinstance(r[k], np.ndarray) else r[k].asnumpy()
            worst_rc[g] = max(worst_rc[g], float(np.abs(v - d(p1, g)[t])[active].max()))
        # --- potential ET, from WaSiM's state after the radiation correction
        if state.ssd_mean is None:
            state.ssd_mean = raw(d(p1, "SSDTemp")[t])
        layers = [pet.Layer(landuse=raw(d(p1, f"landuse{n}")[t]), lai=raw(d(p1, f"leaf_area_index{n}")[t]),
                            coverage=fld(d(p1, f"vegetation_coverage_degree{n}")[t]),
                            rsc=fld(d(p1, f"SurfaceCanopyResistance{n}")[t]),
                            rsi=fld(d(p1, f"SurfaceIntercepResistance{n}")[t]),
                            z0=fld(d(p1, f"RoughnessLength{n}")[t]))
                  for n in range(1, args.layers + 1)]
        tk = d(p1, "corr_temperature")[t]
        tn = d(p1, "temperature_night")[t] if step > 720.0 else tk
        t1 = time.perf_counter()
        res = pet.potential_et(
            site, now.day, now.month, now.hour, step, params, fld(tk), fld(tn),
            fld(d(p1, "global_radiation")[t]), fld(d(p1, "sunshine_duration")[t]),
            fld(d(p1, "radiation_correction")[t]), fld(d(p1, "air_humidity")[t]), fld(d(p1, "wind_speed")[t]),
            fld(d(p1, "albedo")[t]), fld(dem), fld(d(p1, "SurfaceEvaporationResistance")[t]),
            layers, fld(kl), fld(ls), raw(d(p1, "lake_codes")[t]), raw(d(p1, "ponding_storage")[t]),
            state, active, dom, backend=args.backend)
        t_et += time.perf_counter() - t1
        res["ssd_mean"] = state.ssd_mean
        grids = GRIDS + ([("temperature_night", None, "temperature_night")] if step > 720.0 else [])
        for key, lay, g in grids:
            v = res[key] if lay is None else res[key][lay]
            diff = np.abs(v.asnumpy() - d(p2, g)[t])[active]
            worst[g] = max(worst.get(g, 0.0), float(diff.max()))

    print(f"grid {nrows} x {ncols}, {active.sum()} active cells, {args.layers} vegetation layers, "
          f"{nsteps} steps of {step:g} min from {args.start}, backend {args.backend}")
    print(f"ms per step including host-device copies (first step includes compilation): "
          f"radiation correction {1000 * t_rc / nsteps:.1f}, potential ET {1000 * t_et / nsteps:.1f}")
    print("largest absolute difference to WaSiM over all steps and active cells:")
    for g, v in {**worst_rc, **worst}.items():
        print(f"  {g}: {v:.3g}")
    return 0 if max({**worst_rc, **worst}.values()) == 0.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
