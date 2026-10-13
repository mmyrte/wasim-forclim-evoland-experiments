# Reimplementing WaSiM's physics in GT4Py: feasibility

**Short answer: not trivial, but feasible for a chosen subset of processes.** A one-to-one port of all of WaSiM would be very large. A "WaSiM-lite" covering only the processes our Swiss water-balance runs use is realistic for an LLM agent: about 2.5–4 million tokens of agent work (see [`coupling.md`](coupling.md#effort-in-tokens-and-compute)), if every process is validated against WaSiM as it is ported. The ports so far reproduce WaSiM bit for bit. The main cost is not the equations but the implicit coupling between modules, which is now documented in [`coupling.md`](coupling.md).

**Update 2026-10-13:** module coupling documented, no blockers found; interception, snow evaporation and the soil hydraulic functions ported bit-identically; WaSiM has a grid dump hook for per-process validation. Details in `coupling.md`.

First written 2026-10-11, based on WaSiM 11.01.01 (wasim-eth) and GT4Py 1.2.3.

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
| Richards zone (60% of run time) | per-cell column solver with adaptive sub-steps (down to 120 s), table lookups for van Genuchten functions, many branches | fits `gt4py.cartesian` column stencils (K intervals, `while` loops); per-cell sub steps cost GPU efficiency, not correctness; the hardest and most valuable part. The lookup functions are ported (`soil/`) |
| Groundwater | 2-D implicit, line relaxation by rows (parallel over rows) | neighbour stencils plus tridiagonal solves along rows |
| Routing, lakes, reservoir rules | tree of subbasins, channel routing per link, rule tables | not grid-shaped; keep in ordinary code (it takes about 1% of run time) |
| Meteo interpolation (11% of run time) | IDW and regression over all stations for every cell and step | a dense operation that ports easily; precomputing station weights would help regardless of a port |
| Zonal statistics, I/O | per-zone reductions, many file formats | outside the stencil part; use xarray/netCDF |

## Where a port would pay off

At 500 m, all of hydrological Switzerland is about 0.2 million active cells, too small to keep a GPU busy on one run. The case for GT4Py is **ensembles**: calibration runs on Euler could be vectorized across many parameter sets as an extra field dimension, so one GPU or CPU node runs hundreds of members at once. WaSiM cannot do that without major changes. Daily time steps help too: they reduce the step count 24-fold, although the Richards sub-steps keep the saving smaller (Rhone, one thread: about 30 min per year hourly, 6.7 min daily).

## Suggested path

1. **Keep WaSiM for production** of the Swiss setup and build the config generator anyway. The validation method used here depends on WaSiM runs as the reference.
2. **Port process by process**, in order of value for a daily water balance: snow (T-index done), interception (done), snow evaporation (done), Penman-Monteith ET, Richards column (one soil option set; hydraulic functions done), groundwater, then simple routing. Each process is compared at every step with the state WaSiM had just before and after the call (grid dump hook, `coupling.md`).
3. **Decision point after the Richards column**: if one soil configuration reproduces WaSiM within an agreed tolerance, the rest is mostly routine. If not, the port stays a research prototype.

Effort is given in tokens of agent work, measured on the ports so far, in [`coupling.md`](coupling.md#effort-in-tokens-and-compute). GPU validation needs a CUDA node.

## Running the prototype

```sh
# 1. WaSiM reference run with per-step output (wasim-eth, tests/golden/run_case.py)
python3 tests/golden/run_case.py --binary <wasim> --template gt4py/doubs_snow_tindex.wasim.in \
    --case-dir path/to/Doubs --out runs/doubs_snow_tindex
# 2. prototype against it (pip install gt4py netCDF4)
python3 gt4py/snow_tindex.py runs/doubs_snow_tindex path/to/Doubs/StateIni_2021 \
    path/to/Doubs/Input/grids/doubs500.lake_IDs2 --backend numpy
```

The Doubs setup is `sample-run-doubs.zip`.

The newer ports read grid dumps instead of netCDF output; see `coupling.md` ("How the ports are validated"):

```sh
python3 gt4py/interception.py <dump dir> <run dir>/control.wasim          # golden Doubs setup, interception method 11
python3 gt4py/snow_evaporation.py <dump dir>                              # T-index setup (doubs_snow_tindex.wasim.in)
python3 gt4py/snow_evaporation.py <dump dir> --energy-balance --pre after_snow1   # golden setup, snow method 7
python3 gt4py/soil/test_van_genuchten.py path/to/wasim-eth/src           # needs g++
```
