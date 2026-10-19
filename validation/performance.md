# Performance per module: WaSiM and wasim-gt

Status 2026-10-18. The aim is a per-module comparison of legacy WaSiM and the wasim-gt ports on the same case, so that the speed-up of the port (or its cost) can be judged module by module, and later for whole time steps, ensembles and GPUs.

## How it is measured

**WaSiM.** wasim-eth (branch of mmyrte/wasim-eth#2) has a timing hook: `WASIM_TIMING=<file>` sums the wall time between consecutive grid-dump points of the time loop and writes one row per segment at exit. Each segment is one module call, e.g. `after_expokorr -> after_etp` is the potential ET. Without the variable nothing changes (the golden outputs are identical). Run single-threaded:

```sh
OMP_NUM_THREADS=1 WASIM_TIMING=timing.tsv python3 tests/golden/run_case.py --binary src/wasim \
    --template golden/doubs/control.wasim.in --case-dir cases/Doubs --out runs/timing_doubs
```

**wasim-gt.** Each process is called repeatedly on the inputs of one step taken from a grid dump, after warm-up calls that include compilation, single-threaded, on the same machine: `gt4py/benchmark_et.py` for the radiation correction, potential ET and the vegetation parameters. The times include the Python glue and host-device copies of the process function, but not compilation and not dump I/O. The other ports report their times in their PRs (below); they should get the same kind of benchmark script.

Both measure the same work per call, with one difference to keep in mind: WaSiM's module times also contain its statistics and output writing for that module (zonal statistics, grids written per step), which the ports do not do.

## WaSiM, per module

Single thread, clang 18 `-O3`, x86-64 (4 cores, this container), golden cases (10 days, hourly). ms per step:

| Segment (module) | Doubs (6,022 cells) | Rhone (21,584 cells) |
| --- | --- | --- |
| precipitation correction | 0.02 | 0.05 |
| meteo interpolation | 1.92 | 24.65 |
| vegetation parameters (int_tvar_params) | 0.24 | 0.57 |
| data assimilation, output of meteo grids | 0.20 | 0.67 |
| radiation correction | 0.43 | 0.87 |
| potential ET | 1.87 | 4.95 |
| interception | 0.66 | 1.96 |
| snow part 1 | 1.68 | 6.38 |
| snow evaporation | 0.14 | 0.39 |
| unsaturated zone (Richards) | 63.21 | 113.44 |
| groundwater | 1.27 | 3.32 |
| snow part 2 (slides, glacier) | 1.29 | 5.52 |
| post-run (statistics, output) | 6.44 | 11.74 |
| routing | 0.39 | 2.67 |
| rest of the loop (output, next step) | 0.82 | 1.81 |
| **total per step** | **80.6** | **179.0** |

(Interception runs before snow in Doubs, method 11, and after it in Rhone, method 1.) The Richards zone is 63–78 % of the run time.

## wasim-gt, per module

Warm, single thread, same machine unless noted. ms per call (= per step):

| Module | Backend | Doubs | Rhone | WaSiM Doubs / Rhone | Source |
| --- | --- | --- | --- | --- | --- |
| radiation correction (with shading) | embedded | 24.7 | 26.6 | 0.43 / 0.87 | `benchmark_et.py` |
| | gtfn_cpu | 16.5 | 23.8 | | |
| of which the shading walk (host Python) | | 13–19 | 17 | | |
| potential ET (2 layers) | embedded | 67.4 | 66.3 | 1.87 / 4.95 | `benchmark_et.py` |
| | gtfn_cpu | 26.6 | 32.9 | | |
| vegetation parameters (host numpy, once a day) | — | 7.6 | 16.6 | 0.24 / 0.57 per step, i.e. ~5.8 / 13.7 per day | `benchmark_et.py` |
| meteo interpolation (incl. host fits and validation driver) | embedded | 35 | 138 | 1.92 / 24.65 | wasim-gt#20, macOS arm64 |
| | gtfn_cpu | 25 | 112 | | |
| groundwater (per sub step) | embedded | 29 | 38 | 1.27 / 3.32 | wasim-gt#22, macOS arm64 |
| | gtfn_cpu | 22 | 31 | | |
| of which the Gauss-Seidel sweep (host Python) | | 6 | 14 | | |
| snow part 1, method 7 | gt:cpu_ifirst | about 3 | | 1.68 / 6.38 | wasim-gt#19, macOS arm64 |
| glacier dynamics (yearly, host) | — | | about 4 s per run | (not separable: inside snow part 2 on one day) | wasim-gt#24 |

## What the numbers say

- **On one CPU thread and at these grid sizes, the ports are 5–40 times slower than WaSiM per module.** The time hardly grows from Doubs to Rhone (3.6 times the cells): it is per-call overhead, not arithmetic. Potential ET, for example, is about ten small field-operator calls per layer and step plus their Python glue and host-device copies.
- **Host Python loops** cost 13–19 ms (shading walk) and 6–14 ms (Gauss-Seidel sweep) per call on their own.
- **Fixed already:** binding a field operator to a backend on every call cost about 1 s per potential-ET call on `gtfn_cpu`; wasim-gt now caches the bound operators (`fields.with_backend`), which brought it from 1,105 to 27 ms.
- These numbers are not yet the comparison that matters. The case for wasim-gt is many parameter sets at once (an ensemble dimension) on a GPU, where the per-call overhead is shared by all members, and the Richards column, which dominates WaSiM's run time and is not ported yet (#7).

## Next steps

- One benchmark script per ported process, all printing the same table, and one run that fills this file (also for the other modules on this machine).
- Fuse the field operators of a process into one program, then a whole time step (wasim-gt #14); measure again.
- Measure with an ensemble dimension (e.g. 1, 10, 100 members) on CPU and on a GPU (wasim-gt #17), against WaSiM runs in parallel processes, one per member.
- Replace the host loops (shading, slides, glacier walk, Gauss-Seidel) or precompute what depends only on the DEM and the time (shading).
