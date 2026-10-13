#!/usr/bin/env python3
"""GT4Py prototype of WaSiM's T-index snow model (snow method 1), checked against WaSiM output.

What it reimplements (src/snowmod.cc, calc_Accu_and_Melt and run2, hourly time steps):
  1. split precipitation into snow and rain with a linear transition around T0R
     (no snow on lake cells with water in the pond)
  2. accumulate: SSNOW += snow, SLIQ += rain
  3. melt if T > T0 and SSNOW > 0: M = C0 * (T - T0) / steps_per_day, at most SSNOW
  4. move melt to the liquid store; outflow = SLIQ above the holding capacity CWH * SSNOW
  5. subtract snow evaporation (EVRS). WaSiM computes it in the evapotranspiration module
     (etpetr.cc, verdunstungsmodellklasse::snow) after the snow model, from the potential
     evapotranspiration of each vegetation layer and the soil evaporation. The prototype takes
     WaSiM's EVRS grid as an input instead of recomputing it.
Method 1 has no refreezing (WaSiM only refreezes for methods 3 and 4).

Inputs come from a WaSiM run that wrote, at every time step, the corrected temperature
(tcor), precipitation (prec), snow evaporation (etrs) and the snow storages and outflow (ssno, sliq, qsno) as
netCDF; see README.md. The prototype starts from the same initial states and
uses WaSiM's own temperature and precipitation fields, so any difference comes from the
snow step itself.

usage: snow_tindex.py WASIM_OUTPUT_DIR INITIAL_STATE_DIR LAKE_GRID [--backend numpy|gtfn_cpu]
"""
import argparse
import time

import netCDF4
import numpy as np

import gt4py.next as gtx
from gt4py.next import float64

from wasim_gt.fields import I, J
from wasim_gt.processes.snow_tindex import snow_program


def read_ascii_grid(path):
    with open(path) as f:
        header = {}
        for _ in range(6):
            key, value = f.readline().split()[:2]
            header[key.lower()] = float(value)
        data = np.loadtxt(f)
    return data, header["nodata_value"]


def read_nc_series(path, name):
    with netCDF4.Dataset(path) as d:
        # WaSiM writes (time, x, y) with y from north to south; transpose to (time, row, col)
        return np.asarray(d.variables[name][:], dtype="f8").transpose(0, 2, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wasim_out")
    ap.add_argument("state_dir")
    ap.add_argument("lake_grid")
    ap.add_argument("--backend", default="numpy", choices=["numpy", "gtfn_cpu"])
    args = ap.parse_args()

    # parameters of the Doubs setup ([snow_model], method 1); C0 is the same for all months there
    params = dict(t0r=1.0, transition=1.0, t0=1.0, c0=1.6, cwh=0.05, steps_per_day=24.0)

    temp = read_nc_series(f"{args.wasim_out}/tcor.grd.nc", "tcor")
    prec = read_nc_series(f"{args.wasim_out}/prec.grd.nc", "prec")
    evrs = read_nc_series(f"{args.wasim_out}/etrs.grd.nc", "etrs")
    ref_ssno = read_nc_series(f"{args.wasim_out}/ssno.grd.nc", "ssno")
    ref_sliq = read_nc_series(f"{args.wasim_out}/sliq.grd.nc", "sliq")
    ref_qsno = read_nc_series(f"{args.wasim_out}/qsno.grd.nc", "qsno")
    ssno0, nodata = read_ascii_grid(f"{args.state_dir}/ssno.grd")
    sliq0, _ = read_ascii_grid(f"{args.state_dir}/sliq.grd")
    lakes, lake_nodata = read_ascii_grid(args.lake_grid)

    valid = (ssno0 != nodata) & (ref_ssno[0] > -9998)
    lake = (lakes > 0) & (lakes != lake_nodata)
    nsteps, nrows, ncols = temp.shape
    print(f"grid {nrows} x {ncols}, {valid.sum()} active cells, {lake.sum()} lake cells, {nsteps} steps")

    backend = None if args.backend == "numpy" else gtx.gtfn_cpu
    program = snow_program.with_backend(backend) if backend else snow_program
    domain = {I: nrows, J: ncols}
    as_field = lambda a: gtx.as_field(gtx.domain(domain), np.where(valid, a, 0.0))
    ssnow = as_field(ssno0)
    sliq = as_field(sliq0)
    # WaSiM sets the snow fraction to 0 on lake cells only while their pond holds water;
    # the pond state is not compared here, so lake cells are excluded from the comparison
    no_snow = gtx.as_field(gtx.domain(domain), lake)
    out = [gtx.zeros(gtx.domain(domain), dtype=float64) for _ in range(3)]

    compare = valid & ~lake
    worst = {"ssno": 0.0, "sliq": 0.0, "qsno": 0.0}
    t_kernel = 0.0
    for t in range(nsteps):
        start = time.perf_counter()
        program(as_field(temp[t]), as_field(prec[t]), as_field(evrs[t]), ssnow, sliq, no_snow, **params,
                ssnow_out=out[0], sliq_out=out[1], outflow=out[2], offset_provider={})
        if t == 0:
            t_first = time.perf_counter() - start
        else:
            t_kernel += time.perf_counter() - start
        result = [o.asnumpy() for o in out]
        for name, mine, ref in zip(worst, result, (ref_ssno[t], ref_sliq[t], ref_qsno[t])):
            worst[name] = max(worst[name], float(np.abs(mine - ref)[compare].max()))
        ssnow, sliq = as_field(result[0]), as_field(result[1])

    swe_end = (result[0] + result[1])[compare]
    print(f"backend {args.backend}: first step (includes compilation) {t_first:.2f} s, "
          f"then {1000 * t_kernel / (nsteps - 1):.2f} ms per step")
    print(f"snow water equivalent at the end: mean {swe_end.mean():.3f} mm, max {swe_end.max():.3f} mm")
    print("largest absolute difference to WaSiM over all steps and non-lake cells [mm]:")
    for name, value in worst.items():
        print(f"  {name}: {value:.3g}")


if __name__ == "__main__":
    main()
