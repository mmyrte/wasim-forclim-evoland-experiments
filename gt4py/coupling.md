# How WaSiM's modules interfere with each other, and how well they map to stencils

Status 2026-10-13, WaSiM 11.01.01 (wasim-eth), GT4Py 1.2.3. This is the first step of the GT4Py track: write down the hidden data flow between the processes before porting them. **Result: there are no blockers.** The coupling is dense but explicit once listed. Every process except a few sequential ones (snow slides, river routing) maps onto pointwise or column stencils, and three processes are already ported bit-identically (below).

## How modules exchange data

WaSiM has no state object that is passed from process to process. Modules share data in three ways:

1. **The grid table.** Every grid is registered under a name (`initialize_grid(..., identification, ...)`, about 170 names in the Doubs setup), and any module can fetch any grid by name: `X = gridinfo[getindex(_Name_)].pointer`.
2. **Address setters** in `modapuzr.cc`, outside the table. For example, `SetExpoKorrAdress` hands the radiation correction factors to ET and snow, `SetSnowStorAdress` hands the snow store to ET, and `SetRiverflow_tableAdress` hands the river flow table from the soil model to groundwater, routing and irrigation.
3. **In-place modification.** A module often changes another module's output grid instead of writing its own. Then the meaning of a grid depends on the point in the time step. For example, `pot_transpiration1` holds three different values within one step: the potential, the potential after interception, and the potential after snow evaporation.

The third point is what makes a port error-prone, and it explains why the first snow prototype was off by 5 mm. A port has to make each of these hand-overs an explicit input and output.

`analysis/coupling.py` finds the first kind statically. It lists for each module the grids it creates, the grids it fetches, and which fetched grids it writes. Pointer reuse and helper functions make its raw output noisy, so the table below was checked line by line in the code.

## Order within a time step

From `ModellApplication` in `modapuzr.cc`. Arrows mark the grids handed over in place.

| # | Call | Reads (from others) | Writes in place (owned by others) |
| --- | --- | --- | --- |
| 1 | precipitation correction (`precicorr->run`) | station precipitation; wind and temperature at the stations | station values (scaled by `b + a·wind`, rain and snow separately); with gridded or superposed input, the precipitation grid instead (`runOnSuperposedRegions`, after step 2) |
| 2 | meteo interpolation (`meteo[i]->run`) | station data | (creates the meteo grids) |
| 3 | data assimilation | observed grids | meteo grids |
| 4 | radiation/temperature correction (`expokorrektur->run`) | DEM, slope, aspect, sun position | creates correction factors and the shading grid; temperature |
| 5 | potential ET (`verdunstungsmodell->potentiell`) | radiation, temperature, humidity, wind, vegetation parameters, correction factors (setter), albedo (from snow) | radiation (`RG` corrected in place), night temperature |
| 6 | irrigation | soil moisture, river flow table (setter) | precipitation (`+= irrigation`), irrigation amounts (soil model's grid) |
| 7a | interception (methods 11, 12, 21, …: before snow) | precipitation, LAI, coverage, `ETP` per layer | `pot_transpiration` per layer and total (reduced by interception evaporation) |
| 7b | snow model, part 1 (`schneemodell->run1`) | throughfall (7a) or precipitation, temperature, radiation, albedo | input precipitation grid (`+= DEPOINDEX·P_SNOW`, wind redistribution), albedo; with layered snow: the soil model's snow layer stacks |
| 8 | snow evaporation (`verdunstungsmodell->snow`) | snow store (setter), `ETP`, `EVAP` | snow store `SSNOW` (and on glaciers `SSNOWGL`, `EVRSGL`); `ETP` per layer, `EVAP` |
| 9 | silting up | snow outflow, direct runoff | snow outflow (`+= SSSU`, `-= QDSU`) |
| 10 | sub-step loop, `n` times: | | |
| 10a | unsaturated zone (`unsatzonmodell->run`) | snow outflow or rain fraction, `ETP`, `EVAP`, groundwater, river flow table, ponding | `EVAR` (accumulated real evaporation, owned by ET), snow outflow, total snow storage, melt and energy balance (layered snow), lake ponding |
| 10b | surface routing | direct runoff, DEM | direct runoff (`QD`), silting-up storage |
| 10c | real ET (`real_uzr`) | `ETR`, `EVAR` from 10a | (statistics only) |
| 10d | groundwater | recharge, river flow table | baseflow (soil model's grid) |
| 11 | snow model, part 2 (`run2`) | soil model results | snow stores, glacier |
| 12 | post-runs of 10a–d | | |
| 13 | routing, lakes, reservoirs (`abflussmodell->run`) | total runoff per zone (setter), river flow table | pond storage and real ET on lake cells |

Interception methods 1 and 2 swap 7a and 7b. Real evapotranspiration is computed inside the soil model, and the ET module only does statistics for it.

## Cross-module writes, verified

| Writer | Grid (owner) | What happens | Where |
| --- | --- | --- | --- |
| interception | `pot_transpiration1..n`, `pot_transpiration0` (ET) | `ETP -= EI` per layer and in the total | `interzep.cc:1414ff` |
| snow evaporation (ET) | `snowstorage_solid` (snow), `ETP`, `EVAP` | conventional snow methods: `EVRS += min(SSNOW, ETP)` per layer, then from `EVAP`; energy-balance methods (5, 6, 7, 11): `ETP = EVAP = 0` under snow. The total `pot_transpiration0` is *not* reduced here, unlike in interception | `etpetr.cc:2299–2512` |
| snow (wind redistribution) | throughfall or precipitation (interception/meteo) | `P += DEPOINDEX · P_SNOW` | `snowmod.cc:2863` |
| snow | albedo (input grid, read by ET) | snow albedo by snow age, otherwise land-use albedo | `snowmod.cc:3863ff` |
| snow (layered) ↔ soil | snow layer stacks SWES, SWEL, density, thickness, temperature (soil) | snow adds snowfall and melt to the soil model's snow layers; the soil model computes melt and writes back outflow, total storage and effective melt | `snowmod.cc:2739ff`, `unsatzon.cc:8487`, `11369–11377` |
| soil | `real_evaporation` (ET) | accumulated per sub step | `unsatzon.cc:10796` |
| soil | `EnergyBalance` (snow) | capped to the energy needed for melt | `unsatzon.cc:8629` |
| silting up | snow outflow (snow) | `P += SSSU`, `P -= QDSU` | `siltingup.cc:981, 1003` |
| surface routing | `directflow` (soil) | reset and refilled for river cells | `surfrout.cc:3409ff` |
| groundwater | `baseflow` (soil) | exchange with the soil from head changes | `gwflow.cc:2327` |
| routing | `ponding_storage` (soil) | lake water level → pond depth | `abfluss.cc:6698` |
| ET | `ponding_storage` (soil) | evaporation from ponds | `etpetr.cc:2594` |
| ET | `global_radiation` (meteo) | corrected for exposure in place | `etpetr.cc:1956` |
| irrigation | precipitation (meteo), `irrigation_amounts` (soil) | `P += IRRIG` | `irrigate.cc:338ff, 623` |
| precipitation correction (gridded or superposed input only) | precipitation (meteo) | in place | `meteoint.cc:4694ff` |

Not in the grid table, but equally binding: the soil model sorts river-link cells into an extra row, and their exfiltration reaches other cells (`MPISyncRivLnkCells1`, `unsatzon.cc:8079`). This is a cell-to-cell exchange through a graph inside an otherwise column-local module.

## Access patterns and fit to GT4Py

| Process | Share of run time (Rhone, hourly) | Pattern | Fit |
| --- | --- | --- | --- |
| precipitation correction (gridded), snow evaporation, interception, silting up, irrigation (per cell) | small | pointwise (station-based correction: per station); interception and snow evaporation loop over vegetation layers | direct. **Interception and snow evaporation are ported** |
| meteo interpolation (IDW, regression) | 11% | every cell reduces over all stations | good: a dense reduction or matrix product; precomputed station weights help on any platform |
| radiation correction, shading | small | the correction is pointwise; the shading search walks along the sun direction over the DEM | correction: direct; shading depends only on DEM and time, so it can be precomputed or kept on the CPU |
| potential ET (Penman-Monteith) | small | pointwise per vegetation layer | direct |
| snow, T-index and energy balance | 6% | pointwise | direct (T-index ported earlier) |
| snow slides (methods 6–9, 12) | | cascade sorted by elevation | sequential: keep on the CPU, or rewrite as an iterative stencil |
| unsaturated zone (Richards) | 60% | column per cell; adaptive sub steps and `while` loops per cell; variable active depth (groundwater level, snow layers); tridiagonal heat transport per column | good for `gt4py.cartesian` (K intervals, `FORWARD`/`BACKWARD`, `while`). Divergent sub steps cost GPU efficiency but do not block. **Soil hydraulic functions ported** |
| river-link exfiltration in the soil model | | gather/scatter over a list of cell pairs | outside the stencil, as an explicit graph step |
| surface routing | (inactive in the samples) | neighbour stencil with gradients to the four neighbours, iterated until converged | good |
| groundwater | | 2-D implicit, line relaxation by rows, parallel over rows | medium: neighbour stencils plus tridiagonal solves along rows; the rows can be laid out as the vertical axis to use GT4Py's solvers |
| routing, lakes, reservoirs | about 1% | tree of subbasins | not grid-shaped; plain Python or C++ |
| zonal statistics | | reductions per zone | segmented reductions after the stencils |

## Ported so far

All are checked bit for bit against WaSiM: the processes over every step and active cell, the soil functions against WaSiM's own C++ code.

| Port | Kernel (wasim-gt) / check (here) | Reference | Result |
| --- | --- | --- | --- |
| interception, methods 1/11, multiple vegetation layers | `wasim_gt.processes.interception` / `interception.py` | golden Doubs setup, 48 hourly steps, grid dumps before and after the call | all 12 output grids (per-layer and total storage, evaporation, throughfall, remaining ETP) **bit-identical**, numpy and `gtfn_cpu` backends |
| snow evaporation, conventional and energy-balance branches | `wasim_gt.processes.snow_evaporation` / `snow_evaporation.py` | T-index Doubs setup and golden setup (snow method 7), 48 steps | snow store, snow evaporation, ETP per layer, EVAP: **bit-identical** |
| radiation and temperature correction (method 2, with shading) | `wasim_gt.processes.radiation_correction` / `potential_et.py` | Doubs: golden setup 2021-01-01, 48 hourly steps, and `doubs_et_june.wasim.in` (48 sunny hours in June); Rhone daily (`golden/rhone/daily_jan.wasim.in`, `daily_june.wasim.in`), 61 days; dumps `after_meteo` → `after_expokorr` | correction factor, `Ieff_max`, shadow, corrected temperature: **bit-identical**, numpy and `gtfn_cpu` |
| potential ET, Penman-Monteith, several vegetation layers, open water on ponded lake cells, hourly and daily | `wasim_gt.processes.potential_et` / `potential_et.py` | same runs, dumps `after_expokorr` → `after_etp` | ETP and EIP per layer and total, EVAP, corrected radiation, `AspectCorrectedGlobalRadiation`, cloudiness, `SSDTemp`, day and night temperature (daily): **bit-identical**, numpy and `gtfn_cpu` |
| vegetation parameters from the land-use table (`VariableDayCount`), initialisation and daily update | `wasim_gt.processes.vegetation` / `vegetation.py` | golden Doubs setup, 10 days hourly (initialisation and 10 updates); Rhone daily, all of 2021 (365 updates); dump `after_meteo` | LAI, VCF, rsc, rsi, z0 and root depth per layer, albedo, rse: **bit-identical** |
| dynamic phenology DP1–DP4: forcing, chilling, start of growth, multiple start times, leaf fall, and the vegetation parameters that follow | `wasim_gt.processes.dynamic_phenology` / `dynamic_phenology.py` | synthetic Rhone case `rhone_dynphen.wasim.in` (one land use per model), daily, 2021-01-01 to 2022-07-31; dump `after_meteo` | all phenology grids and vegetation parameters: **bit-identical** (577 days, state carried from the start) |
| soil hydraulic functions `h_calc`, `k_calc`, `get_theta` | `wasim_gt.soil.van_genuchten` / `soil/test_van_genuchten.py` | WaSiM's C++ functions, extracted verbatim (`soil/extract_reference.sh`), 220,000 random horizons and queries including edge cases | **bit-identical** (`gt:cpu_ifirst`) |

In WaSiM the soil functions interpolate in lookup tables of 30 entries per soil horizon. Storing these tables per cell and level would take several GB for Switzerland, so the port recomputes only the entries that each lookup touches. It uses the same floating-point expressions as WaSiM's table generation, so the results are identical. The search is a fixed five-step bisection, so every GPU thread does the same work. This replaces memory traffic by about a dozen `pow` calls per lookup, a good trade on a GPU.

Together with the earlier T-index prototype, the chain phenology, radiation correction, potential ET, interception, snow evaporation can now run from the meteo fields and the albedo under snow (from snow method 7, not yet ported).

The coverage tool (`../validation/wasim_coverage.py`) over the four reference runs shows every line of `getETP` run except net radiation as input (not ported), and in `potentiell` and `expokorrklasse::run` all lines of the ported branches except: the cap of the radiation at 1260 Wh/m² per bright hour; the albedo without a lake model; the sun below 1° altitude in a bright sub interval; the temperature formulas for `split` > 1; and branches that cannot run at Swiss latitudes (sunset after midnight, polar day and night) or with these settings (SSDTemp without a bright step).

## How the ports are validated

The wasim-eth debug hook `griddump` in `src/modapuzr.cc` writes named grids at fixed points of the time step. It is enabled by environment variables and does nothing without them; the golden test is unchanged.

```sh
WASIM_DUMP_GRIDS=precipitation,pot_transpiration1,interception_storage1 \
WASIM_DUMP_DIR=dumps/doubs WASIM_DUMP_POINTS=after_etp,after_interception \
python3 tests/golden/run_case.py --binary src/wasim --template ... --case-dir ... --out ...
# dumps/doubs/identifications.txt lists every grid name
python3 gt4py/interception.py dumps/doubs runs/.../control.wasim
```

Each port reads the state WaSiM had just before the call and compares its results with the state just after it. The comparison therefore isolates one process, even when other modules change the same grids. `wasim_dump.py` reads the dumps.

## Findings along the way (not fixed)

- Snow evaporation adds each cell's EVRS to the zonal statistics twice (`etpetr.cc:2464` and `2507`).
- Snow evaporation reduces ETP per layer but not the layer total `pot_transpiration0`, while interception reduces both.
- `get_theta` reads `h_tab[0]` for suctions ≤ 0, but the table generation never sets it (uninitialized memory; in practice the result is θs).
- `k_calc` searches its table with `anz_h` instead of `anz_k`. This is harmless only as long as both are 30.
- Potential ET: the LAI accumulated above a layer (`LAIkum`, `etpetr.cc:2042`) adds the LAI of the layer above even where that layer has no land use, i.e. adds -9999, which makes the aerodynamic resistance of the lower layer negative. Only cells whose top layer is empty but a lower one is not are affected; the Doubs and Rhone grids have none.
- Potential ET corrects the global radiation in place (`etpetr.cc:1956`) but caps only its copy `RGEX` at 1260 Wh/m² per bright hour, so later modules (snow) can see a corrected radiation above the cap.
- The constant `c64_0` is 25.0 (`etpetr.cc:1349`), the ratio of measurement height and roughness length for z0 > 2 m.
- Phenology: the initialisation and the daily update clamp differently and use different elevations on lake cells (before and after the lake model lowers them); two leap-year rules, both wrong for 2000; before the first fix point the values are interpolated from the last one without its altitude shift. Details in wasim-eth `TODO.md`.
- WaSiM replaces the land-use codes in the grids by the position of the record in the table (`defltuzr.cc:1825`), so dumped land-use grids hold positions, not codes.
- Dynamic phenology: the state grids are created with zeros, so the initialisation rules for unset (negative) values never apply to a fresh run: the start of growth is 0 and F* is 0 until the first `JDReset_TStart`, and DP2/DP4 skip their DP1 fallback, unlike the manual says. The day number counts February with 29 days in every year of a run that ends in a leap year. Details in wasim-eth `TODO.md`.

## GT4Py limits met

- `gt4py.next`: a field operator called outside a program needs `out=` even on the embedded backend.
- `gt4py.cartesian`: data dimensions can be indexed at run time only on `gt:*` and `dace:*` backends, not on `numpy`. No function calls in `if` conditions, no `//`.
- Bit-identical results on a CPU rely on no FMA contraction and the same `pow`, `exp` and `log`. numpy's vectorised versions (used by the embedded backend) differ from the C library in the last bit for 2–5% of arguments on AVX-512 machines; `NPY_DISABLE_CPU_FEATURES=X86_V4` (numpy 2.5) makes them agree. `gtfn_cpu` calls the C library directly. CUDA's `pow` and `exp` differ by up to about 2 ulp, so GPU runs have to be validated with a tolerance.

## Effort, in tokens and compute

The work is done by an LLM agent, so effort is measured in tokens and machine time rather than person-months. The measure used here is the session token counter: new tokens entering the context, meaning code and output read plus code written. Each agent step re-reads the context, mostly from cache, so billed input tokens are much higher, typically 20–100 times this count.

**Measured in this session:**

| Work | Tokens (counter) | C++ lines read closely |
| --- | --- | --- |
| coupling analysis (script, verification, this document) | about 170,000 | across all modules |
| dump hook in WaSiM | about 15,000 | |
| interception port and validation | about 25,000 | about 400 |
| snow evaporation port, both branches | about 12,000 | about 200 |
| soil hydraulic functions, C++ reference extraction, test | about 45,000 | about 250 |
| radiation correction and Penman-Monteith potential ET, hourly and daily; Rhone daily setup, synthetic case, coverage (later session) | about 330,000 | about 1,100 |
| table-based phenology, same session | about 100,000 | about 500 |
| dynamic phenology DP1–DP4, synthetic case, same session | about 120,000 | about 900 |

For a port, the rate works out to about 100–150 tokens per C++ line that has to be read, rewritten and validated. That covers debugging and the GT4Py frontend's quirks, and comes on top of a fixed cost of about 10,000 tokens per process for a dump run and a comparison script.

**Estimate for the rest:**

| Item | Relevant C++ lines (estimate) | Tokens (counter) |
| --- | --- | --- |
| potential ET (Penman-Monteith) and radiation correction | 1,000 | 0.1–0.2 M (**done**, with both phenologies: about 0.55 M; see above) |
| precipitation correction, meteo interpolation (IDW, regression) | 1,500 | 0.15–0.25 M |
| snow energy balance (method 7) without slides; slides on the CPU | 2,000 | 0.2–0.35 M |
| dynamic glacier model (Rhone) | 1,500–2,000 | 0.2–0.3 M |
| simple permafrost model | 100–200 | 0.03–0.05 M |
| Richards column, one option set (no tracers, no macropores), plus river links | 6,000–8,000 | 0.8–1.5 M |
| heat transport and layered snow in the column, only if the Swiss setup needs them (not in wasim-gt's scope) | 2,000–3,000 | 0.3–0.5 M |
| groundwater | 2,000 | 0.2–0.35 M |
| routing and lakes in Python | 3,000–4,000 | 0.3–0.6 M |
| driver: state object, time loop in compiled code, I/O, ensemble dimension, GPU tuning | | 0.3–0.5 M |
| **total** | | **about 2.5–4 M** |

That is roughly 10–20 sessions of this one's size, and the uncertainty is a factor of 2. The Richards column carries most of the risk: its many switches may cost more than the line count suggests.

**Machine time:** a dump run takes seconds (Doubs) to a minute (Rhone, 10 days). GT4Py compiles take 1 minute per stencil on `gtfn_cpu`/`gt:cpu` and longer on `dace`. The CPU work fits in this container. GPU validation and benchmarks need a CUDA node (Euler) with `gt:gpu`, `dace:gpu` or `gtfn_gpu`. The payoff is in calibration: many parameter sets run as an extra field dimension in one GPU job, where WaSiM needs one process per member.

## Next steps

The remaining work is split into issues, one per process, in the wasim-gt repository, under its tracking issue [mmyrte/wasim-gt#1](https://github.com/mmyrte/wasim-gt/issues/1). Validation work stays here: harness #16, synthetic cases #17, legacy tools #18.

1. Richards column for the Doubs/Rhone option set, in `gt4py.cartesian`, built on `soil/van_genuchten.py`; validate against `after_unsatzon` dumps per sub step.
2. On a GPU node: run the soil-function test and interception on `gt:gpu`/`dace:gpu`, measure throughput with an ensemble dimension, and set the validation tolerance.
