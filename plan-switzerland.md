# WaSiM for hydrological Switzerland: working plan

Status 2026-10-11, adapted from the project plan doc and its review comments. **Entry point for the second half of November: read "Where things stand", then pick from "Next steps".**

## Where things stand

| Item | State | Where |
| --- | --- | --- |
| macOS and unified Linux/macOS build | done; g++ and clang builds give identical results; legacy makefiles kept | wasim-eth `src/Makefile`, `README.md` |
| Golden tests (run a case, compare all outputs) | done; Doubs and Rhone, 10 days hourly | wasim-eth `tests/golden/`, here `golden/` |
| Reference outputs and Doubs setup | uploaded (git-lfs) | `golden/<case>/reference.zip`, `sample-run-doubs.zip` |
| Memory and runtime measurements | done for Doubs and Rhone; hourly vs daily; OpenMP and MPI | below, and wasim-eth `docs/domain-size-limits.md` |
| Performance profile and quick wins | done: profile below; LTO 3%; `OMP_SCHEDULE=dynamic,4` helps short winter runs but slows a full year, so no general recommendation yet | below |
| GT4Py feasibility | assessment plus a validated T-index snow prototype | `gt4py/README.md` |
| Code fixes | done: g++ builds could not write netCDF (overlapping `sprintf`), memory and format errors, netCDF buffer sizes, zone code bounds, UTF-8 sources, `--help` | wasim-eth `TODO.md`, commit history |
| Other problems found | logged, not fixed (they change results or need a decision) | below and wasim-eth `TODO.md` |
| National input data, routing, config generator | not started | next steps |

## Goal and scope

One harmonized WaSiM setup for **hydrological Switzerland at 500 m in LV95 (EPSG:2056)**, replacing the 17 regional setups that today run one after another.

- **Domain:** includes the foreign headwaters (Rhine with Vorarlberg and Baden-Württemberg, Inn, Ticino with Italian headwaters, Doubs, Rhône/Arve), because part of the research is on international cooperation in adaptation. Input data must cover the whole domain; calibration gauges may be restricted to Switzerland, since that is the domain the model is optimised for.
- **Domain definition:** take the catchment outlines from CAMELS-CH (331 catchments in hydrological Switzerland) and derive the 500 m grid from them, instead of reconstructing the extents of the 17 setups.
- **In scope:** low-hanging performance work, small adjustments to the WaSiM entry point, an external config generator (robust interpolation and path generation, legible syntax).
- **Out of scope for now:** new process modules.

### Two uses, two run modes

| Use | Drivers | Time step | Notes |
| --- | --- | --- | --- |
| Water balance, low flows and drought; coupling with ForClim and evoland | coarse RCM output including radiation, bias-corrected (possibly by weather type) | daily (monthly exchange with ForClim, decadal with evoland) | cheap enough for large calibration ensembles |
| Extreme precipitation and floods | convection-permitting models | hourly or shorter | event runs only, started from states of the daily run ("hydro boosting", analogous to ensemble boosting, Fischer et al.) |

WaSiM supports starting a run from states saved at another time step: it converts the routing and runoff-concentration storages between time steps (manual 4.5.3.1, p. 319). The daily Rhone run below started from hourly states this way. **Still to test:** an hourly event run started from daily states, compared with a continuous hourly run.

## Compute: measured, then extrapolated

All runs on 4 cores of this machine, single-node, WaSiM built with clang at `-O3`.

| Run | Active cells | Steps | 1 thread | 4 threads | Peak memory |
| --- | --- | --- | --- | --- | --- |
| Doubs, 10 days hourly | 6,022 | 240 | 20 s | 10 s | 187 MB |
| Rhone, 10 days hourly | 21,584 | 240 | 49 s | 27 s (static) / 20 s (`dynamic,4`) | 324 MB |
| Rhone, 2021 hourly | 21,584 | 8,760 | about 30 min (extrapolated from 10 days) | 17 min (static, partly on a shared machine); with `dynamic,4` over 40 min | 331 MB |
| Rhone, 2021 daily | 21,584 | 365 | 401 s (6.7 min) | 119 s (`dynamic,4`) | 324 MB |

- **Daily vs hourly:** a daily year costs about 4.5 times less than an hourly year (single thread: 6.7 vs about 30 min), not 24 times less. The Richards zone sub-steps adaptively (down to 120 s), so most of its work remains. Your recollection of "about 5 minutes per Rhone year" fits a daily run.
- **All of Switzerland at 500 m:** about 10 times the Rhone's active cells (about 0.22 million). Memory stays below 1 GB per process. A daily year should take roughly 20 minutes on 4 threads, or about 1.2 core-hours. Hourly costs about 4.5 times more.
- **Ensembles on Euler:** at about 1.2 core-hours per member and simulated year (daily), a 30-year, 100-member ensemble is about 3,600 core-hours. Many single-node runs in parallel use the cores better than one MPI run: MPI gave only 2.2× on 4 ranks for Doubs, and rank 0 always needs the full domain in memory.
- **Where the time goes** (Rhone, hourly): Richards zone 60% (van Genuchten table lookups alone 36%), meteo interpolation 11% (IDW over all stations for every cell and step), snow 6%, reading and writing about 5%. I/O is synchronous in the main loop but small here; batching it would gain at most a few percent unless many grids are written every step.

## Input data

| Input | Default | Alternatives and notes |
| --- | --- | --- |
| Elevation | Copernicus DEM GLO-30, resampled to 500 m | EEA-10 for deriving stream density; swissALTI3D inside CH |
| Soil | SoilGrids 250 m, clustered into a WaSiM soil table (van Genuchten via pedotransfer functions) | Swiss soil property maps ([Zenodo 18428152](https://zenodo.org/records/18428152)) as refinement, though harmonizing a national product with the foreign parts will be hard |
| Land cover | CORINE 2018, 44 classes remapped to WaSiM classes with `landuse/corine_to_wasim.rmp` (44 to 17) | watch the CORINE 2025 release; the class set depends on evoland: SSP-CH targets were elicited for 10 classes, and "static" classes need splitting (deglaciated areas, built-up areas abandoned under extreme scenarios) |
| Routing | rebuilt from scratch: flow directions from the DEM, subbasins at gauges and river mouths | at 500 m probably a few hundred to a few thousand subbasins make sense (the FOEN mosaic has over 22,000); represent small channels through drainage density; mind the code's limit of 16,000 zones and codes up to 32,767 |
| Meteorology | historical: MeteoSwiss grids; scenarios: RCM (daily) and convection-permitting models (events) | bias correction by weather type for the RCM |
| Streamflow | compiled dataset from Haritha (Uni Bern), if available | FOEN hydrometric stations on [opendata.swiss](https://opendata.swiss/) (rivers and lakes); CAMELS-CH; foreign agencies for the headwaters |
| Glaciers | inventories 1973, 2010, 2016, 2023; [GLAMOS](https://doi.glamos.ch/index.html) mass balances | calibrate on recent years: the rapid melt of the last few years avoids a conservative bias |
| Hydropower | timing derived from electricity market data (ENTSO-E, SFOE) instead of per-dam rules | the 17 setups contain hand-written reservoir and abstraction rules; some research may want to adapt them to today's market or to other retention rationales (e.g. irrigation, water seeding) |

## Control files and the config generator

The inherited control files have drifted apart (1,257 different `$set` names across the 17 files, 143 soil types, 17 land use tables). Every file starts with a block marked "protected part to be changed by the module adapter": run IDs, start and end dates, paths and the `readgrids` switch. A separate tool, the module adapter, rewrites that block when runs are chained (it is not in this repository). A generator replacing it would also replace the regex edits used so far.

First step, as you suggested: **a general generator that writes correct control files from JSON**, before deciding how to fill individual sections.

1. A schema for the run-level settings (period, time step, paths relative to one root, module switches, output codes). Paths are built in one place; this would have avoided several problems found this week (below).
2. Templates per section, filled from JSON. Positional sections keep their order; the template documents each line, as in wasim-eth `docs/sample_control_file.wasim`.
3. Golden tests as the check: a generated Doubs and Rhone control file must give the same outputs as the originals.
4. Later: parameter tables (soil, land use) and regionalized parameters as separate files; a parser for the 17 legacy files only if their parameter values are still needed.

## Next steps

1. **Find the Doubs reference outputs** from the legacy runs, to check the golden references against them.
2. **Domain:** build the 500 m LV95 grid and mask from CAMELS-CH catchments.
3. **Config generator:** JSON to control file for the run-level settings; golden round trip on Doubs and Rhone.
4. **Hydro boosting test:** hourly Rhone event from daily states vs a continuous hourly run.
5. **Routing from scratch** on the national grid; decide the subbasin density.
6. **Decide on the GT4Py track** after reading `gt4py/README.md`; if yes, the next process to port is the Richards column.

## Problems found while setting up the runs

Logged in wasim-eth `TODO.md`, not fixed yet. (Fixed this week: g++ builds failed to write any netCDF output because of overlapping `sprintf` calls; clang builds were unaffected.)

- Per-layer statistics files get the output directory prepended twice when `DefaultOutputDirectory` is relative (multi-layer land use); absolute paths work around it.
- A symbol followed by `/` instead of `//` (e.g. `$outroot/`) makes the control file parser loop forever instead of stopping with an error.
- `$sym//` at the end of a value drops the trailing `/`.
- Results depend on the OpenMP thread count and schedule for a few dozen cells (up to about 1 in soil moisture for Rhone), beyond rounding: likely a data race or row dependency in the soil loop.
- Building with `-march=native` (or `-march=haswell`, as in the legacy Linux makefiles) changes results through FMA contraction.
- MPI splits the domain by a timing benchmark at start-up, so the decomposition, and thus the results, can change between runs; a `MpiSpeedFactor.dat` file fixes it.
- The Rhone sample references `r500.s01`; the archive contains `r500.s01.asc`.
