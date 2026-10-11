# Reimplementing WaSiM's physics in GT4Py: feasibility

**Short answer: not trivial, but feasible for a chosen subset of processes.** A one-to-one port of all of WaSiM is a multi-year effort. A "WaSiM-lite" covering only the processes our Swiss water-balance runs use is a realistic project of a few person-months, if every process is validated against WaSiM as it is ported. The prototype here shows that such validation works to machine precision. The main cost is not the equations but the implicit coupling between modules.

As of 2026-10-11, based on WaSiM 11.01.01 (wasim-eth) and GT4Py 1.2.3.

## What was tried

`snow_tindex.py` reimplements WaSiM's T-index snow model (snow method 1) as a GT4Py field operator and steps it through the same 240 hourly steps (Doubs, 2021-01-01 to 01-11) as a WaSiM reference run. WaSiM's own temperature, precipitation and snow evaporation fields drive it, from the same initial states.

| Check | Result |
| --- | --- |
| Largest difference in snow storage, liquid storage and snow outflow, all steps, 5,990 non-lake cells | 1e-13 mm (rounding) |
| Kernel size | about 15 lines of GT4Py |
| Time per step, Doubs grid (44,421 cells), `numpy` backend | 5 ms |
| Time per step, compiled `gtfn_cpu` backend (after 45 s compilation) | 17 ms |

The first attempt differed by up to 5 mm. The cause was that WaSiM's evapotranspiration module subtracts snow evaporation from the snow store in a separate pass after the snow model (`etpetr.cc`, `verdunstungsmodellklasse::snow`). That coupling is invisible from the snow code. It only showed up because every step was compared against WaSiM.

The compiled backend is slower here because a kernel this small is dominated by per-call overhead from Python. A port would have to run the time loop, or at least many fused processes per call, inside compiled code.

## Why a full port is not trivial

| Aspect | What the code looks like | Consequence for a port |
| --- | --- | --- |
| Size | about 44,000 lines of physics (Richards zone 13,600, routing/lakes 6,700, snow/glacier 5,800, meteo interpolation 4,500, surface routing 3,600, groundwater 2,500, ET 2,400), plus 21,000 lines of infrastructure | port selectively, not wholesale |
| Method variants | 14 snow methods, 12 interpolation methods, 4 ET methods, many switches in the soil model (ponding, drainage, clay layer, macropores, irrigation, heat transfer, layered snow, tracers) | choose one variant per process for our setup and drop the rest |
| Implicit coupling | about 140 shared grids, looked up by name (`gridinfo[getindex(...)]`); modules modify each other's state (e.g. ET writes the snow store) | write down the data flow per process before porting it; validate per step |
| Richards zone (60% of run time) | per-cell column solver with adaptive sub-steps (down to 120 s), table lookups for van Genuchten functions, many branches | fits GT4Py's vertical `scan_operator`, but adaptive sub-stepping per column has to become masked fixed iterations; this is the hardest and most valuable part |
| Groundwater | 2D Gauss-Seidel iteration | a natural stencil, but Gauss-Seidel is sequential; a parallel version (Jacobi, red-black or CG) changes the numerics slightly |
| Routing, lakes, reservoir rules | tree of subbasins, channel routing per link, rule tables | not grid-shaped; keep in ordinary code (it takes about 1% of run time) |
| Meteo interpolation (11% of run time) | IDW and regression over all stations for every cell and step | a dense operation that ports easily; precomputing station weights would help regardless of a port |
| Zonal statistics, I/O | per-zone reductions, many file formats | outside the stencil part; use xarray/netCDF |

## Where a port would pay off

At 500 m, all of hydrological Switzerland is about 0.2 million active cells, too small to keep a GPU busy on one run. The case for GT4Py is **ensembles**: calibration runs on Euler could be vectorized across many parameter sets as an extra field dimension, so one GPU or CPU node runs hundreds of members at once. WaSiM cannot do that without major changes. Daily time steps help too: they reduce the step count 24-fold, although the Richards sub-steps keep the saving smaller (Rhone, one thread: about 30 min per year hourly, 6.7 min daily).

## Suggested path

1. **Keep WaSiM for production** of the Swiss setup and build the config generator anyway. The validation method used here depends on WaSiM runs as the reference.
2. **Port process by process**, in order of value for a daily water balance: snow (done for T-index), interception, Penman-Monteith ET, Richards column (one soil option set), groundwater, then simple routing. Each process gets a per-step comparison against a WaSiM run that writes its inputs and outputs, as in `doubs_snow_tindex.wasim.in`.
3. **Decision point after the Richards column**: if one soil configuration reproduces WaSiM within an agreed tolerance, the rest is mostly routine. If not, the port stays a research prototype.

Rough effort, for orientation only: interception and ET a few days each, Richards column several weeks, groundwater one to two weeks, routing and lakes in plain Python or C++ another few weeks.

## Running the prototype

```sh
# 1. WaSiM reference run with per-step output (wasim-eth, tests/golden/run_case.py)
python3 tests/golden/run_case.py --binary <wasim> --template gt4py/doubs_snow_tindex.wasim.in \
    --case-dir path/to/Doubs --out runs/doubs_snow_tindex
# 2. prototype against it (pip install gt4py netCDF4)
python3 gt4py/snow_tindex.py runs/doubs_snow_tindex path/to/Doubs/StateIni_2021 \
    path/to/Doubs/Input/grids/doubs500.lake_IDs2 --backend numpy
```

The Doubs setup is `sample-run-doubs.zip` (upload pending, see `golden/README.md`).
