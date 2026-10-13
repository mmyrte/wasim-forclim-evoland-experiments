# What the GT4Py port covers, compared with all of WaSiM

Status 2026-10-13, WaSiM 11.01.01; section numbers refer to the manual ([wasim_2025_en.pdf](https://www.wasim.ch/downloads/doku/wasim/wasim_2025_en.pdf), chapter 2).

**Terms.** *Full WaSiM* is everything WaSiM 11.01.01 can compute: every process and every method variant in chapter 2 of the manual, plus its infrastructure. The *port scope* is the subset the GT4Py port aims at: **the processes and method variants switched on in our two reference setups, Doubs and Rhone** (`golden/*/control.wasim.in`). These are the variants the Swiss setup is expected to use.

The port scope contains no simplified physics. A process that is in scope is ported as WaSiM computes it, and checked bit for bit (see `coupling.md`). What is left out are whole processes and method variants that these setups do not use. Earlier documents called the port scope "WaSiM-lite"; that term is dropped because it suggests simplified physics.

The Swiss setup is still being defined (`plan-switzerland.md`). If it switches on a process that is out of scope here, for example heat transport for permafrost, that process moves into scope.

## Process by process

| Process (manual §) | Variants in WaSiM | Doubs | Rhone | Port scope | Issue / status |
| --- | --- | --- | --- | --- | --- |
| precipitation correction (2.2) | station-based; gridded for superposition | on | on | yes | #4 |
| meteo interpolation (2.3) | 12 methods: IDW, regression, IDW + regression, Thiessen, bilinear (2), splines (2), external grids, regression with internal preprocessing (10), IDW + method 10 (11), nearest neighbour with lapse rate | 10 and 11 | 10 and 11 | methods 10 and 11 (IDW is part of 11); reading external grids (9) for gridded MeteoSwiss input | #4 |
| regional superposition (2.3.13), scenarios (2.4) | | off | off | no | |
| radiation and temperature correction by topography (2.5) | | on | on | yes | #3 |
| potential ET (2.6) | Penman-Monteith; Hamon, Wendling, Haude (daily only); open water | Penman-Monteith | Penman-Monteith | Penman-Monteith, and open-water evaporation on lake cells (the lake model is on) | #3 |
| real ET (2.8) | Topmodel version; Richards version; one or several vegetation layers; wet reduction | Richards, 2 layers | Richards | Richards version, several layers | inside #6 |
| snow (2.9) | 14 methods: T-index, T-u-index, energy balance (Anderson, Braun, enhanced), each with or without slides and wind redistribution | 7 (enhanced energy balance + slides + wind) | 7 | method 7; method 1 as a prototype | #5; T-index done |
| snow evaporation | conventional; energy-balance branch | energy balance | energy balance | both | **done** |
| interception (2.10) | bucket methods 1/2 (after snow) and 11/12 (before snow); snow canopy with energy balance (21, 22, 31, 32) | 11 | 1 | 1 and 11 | **done** |
| infiltration model (2.11) | Green-Ampt (Topmodel version) | off | off | no | |
| Topmodel soil model (2.12) | | off | off | no | |
| silting up, surface routing (2.13) | | off | off | no (#10, #12 kept for event runs) | #10, #12 |
| unsaturated zone, Richards (2.14) | FDM with dynamic time step; surface interaction; ponding; drainage; clay layer; macropores; irrigation; groundwater coupling; baseflow | on: method 3, river infiltration and exfiltration (2), ponds | on: method 3, exfiltration only (1) | yes, the options the two setups use | #6, #7; soil functions **done** |
| groundwater flow (2.15) | one or more aquifers; with or without lateral transport | on | on | yes, flow only | #9 |
| discharge routing (2.16) | kinematic wave; observed inflows; shorter routing step | on | on | yes | #11 |
| routing-only mode (2.17) | | off | off | no | |
| irrigation (2.18) | | off | off | no (#12 if needed) | #12 |
| lake model and abstraction rules (2.19) | lakes; reservoirs; time-variant rules; hydraulic links; externally given runoff | on | on (47 rule sections) | yes | #11 |
| solute transport (2.21) | tracers in all modules | off | off | no | |
| glacier model (2.22) | static glaciers; dynamic glaciers (retreat and advance, volume-area relation); per-subbasin parameters | no glaciers | dynamic (method 13) | yes, as the Rhone setup uses it | #15 |
| dynamic phenology (2.23) | DP1–DP4 | table-based (`VariableDayCount`) | table-based | table-based vegetation only | inside #3 and #13 |
| simple permafrost model (2.24) | | section present, inactive | section present, inactive | no | |
| heat transport (2.25), spin-up with HRUs (2.26) | | off | off | no, unless the Swiss setup needs permafrost | #8 |
| layered snow model (2.27) | | off | off | no | #8 |
| data assimilation (2.28) | gridded; reservoirs | off | off | no | |
| ArcVeg (2.29) | | off | off | no | |
| MPI, external coupling | | | | no: the port replaces MPI by GPU and ensemble parallelism | |

## Infrastructure outside the processes

| Item | Port scope |
| --- | --- |
| control file parser | no: the port reads its own configuration; the planned control file generator (`plan-switzerland.md`) can write both formats |
| input formats | netCDF and ASCII grids, WaSiM station tables; no binary grids or GIF |
| output | netCDF grids and zonal statistics, comparable with WaSiM's files; not all of WaSiM's file types |
| states | reading WaSiM's initial states (`StateIni_*`); writing states and converting them between time steps (manual 4.5.3.1) is needed for "hydro boosting" |
| performance criteria (2.20) | no: computed in the calibration framework |

## What a complete replication would take

A complete replication would port every row above, including all method variants, plus the infrastructure: control file parser, all file formats, all output types, state handling, MPI or an equivalent, and the routing-only mode.

- **Size:** about 44,000 lines of physics and 21,000 of infrastructure (`README.md`). The port scope needs roughly 15,000–20,000 of the physics lines.
- **Tokens:** at the rate measured so far, 100–150 tokens per line read, rewritten and validated, the physics alone comes to 4.5–6.5 M tokens and the infrastructure to another 2–3 M. Port scope: 2.5–4 M (`coupling.md`).
- **The larger cost is test setups, not code.** Every method variant needs a WaSiM setup that switches it on, as a reference to validate against. We have two, Doubs and Rhone, and they cover only the variants listed as "on" above. For the rest, someone has to build or obtain a setup per variant: the 14 snow methods, 12 interpolation methods, heat transport, tracers, ArcVeg and others. The test cases in the manual and the WaSiM distribution may help.
- **Without such setups, a port of an unused variant cannot be validated.** It would be a translation without a check, and the T-index prototype showed what hidden couplings do to that.

Our recommendation: port the scope listed here, and add a variant only when a setup that uses it exists.
