# Richards column: phase-0 design and validation plan

Status 2026-10-19 (phase 0 of 2026-10-15, plus the `h_temp1` follow-up and a first ensemble benchmark). Reference: WaSiM 11.01.01, wasim-eth commit `50c9ad8ca8eaced47a5a421248f8dafb43ceb736`; GT4Py 1.2.3, arm64 CPU, one WaSiM thread and no FMA contraction. This document is analysis and scratch-prototype evidence. **No production Richards process has been ported and no wasim-gt file has been created.**

The helper stencil is feasible and exact in the tested cases. Two risks prevented treating the full routine as a proven independent column: conditional reuse of row scratch `h_temp1`, and the ordered river-link exchange in Doubs. The first is now settled: the reuse is reachable in a real year run and depends on the column order, and a fix is proposed in mmyrte/wasim-eth#3 (§4, "Real runs"). The installed frontend also rejects a `while` enclosing an entire K computation. Phase 1 can proceed as isolated helper packages; phase 2 needs an explicit decision on history/order and staged iteration orchestration.

## 1. Resolved options and ownership

Values come from the resolved controls of fresh 48-hour runs. The actual section is `[unsatzon_model]`, not `[unsatzon]`; `[soil_model]` is disabled (the older Topmodel process).

| Option | Doubs | Rhone | Scope / owner |
| --- | --- | --- | --- |
| Soil algorithm | Method 3, dynamic FDM | Method 3, dynamic FDM | #7 |
| Soil levels | 20 plus surface slot 0 | 15 plus surface slot 0 | #7, layout agreement with #14 |
| Module interval / substeps | 60 min / 1 | 60 min / 1 | #7 + driver #14 |
| Minimum cell step | 120 s | 300 s | #7; actual splitfactor follows Courant test |
| Nonlinear iteration | max error 1e-5, max 100 iterations | Same | #7; after iteration 20 threshold becomes 5e-5 |
| Vegetation | 2 layers, ordered extraction | 2 layers, ordered extraction | #7 real ET; potential ET elsewhere |
| Surface interaction | 2: exfiltration + river reinfiltration | 1: exfiltration only | Local surface exchange contract in #7, river-link graph/table availability in #8 |
| Ponds / lakes | Enabled / enabled | Enabled / enabled | Soil pond budget #7; reservoir/lake routing #12 |
| Explicit groundwater | On, one aquifer | On, one aquifer | #7 lower-boundary handover; lateral solver #10 |
| Macropores | Helper called, little wet-path coverage | Helper called, more wet-path coverage | **In #7**; do not drop because an earlier estimate said “no macropores” |
| Drainage / clay | Both disabled | Both disabled | Out for current options; retain separate future option request |
| Simple thaw-depth input | Not bound in these references | Not bound in these references | #6 produces optional input, #7 consumes conductivity reduction |
| Heat / layered snow | Off / off | Off / off | #9, excluded; no negative snow K offsets in #7 |
| Glacier fraction | None | Dynamic fractional mask | #7 consumes mask and clears fully ice-covered columns; glacier evolution elsewhere |
| Irrigation | Off | Off | #13, excluded |
| Tracers / ArcVeg / HRU spin-up | Off | Off | Excluded |
| Travel-time zones + linear reservoirs | Enabled | Enabled | #7 zonal phase; #12 receives runoff |

Pond, river, macropore and groundwater options are separately configured switches in the legacy routine. The primary column interface must expose each, rather than inferring all of them from a single soil-method flag.

## 2. Exact call structure and time scales

`modapuzr.cc:1291–1335` establishes the order:

1. Silting-up call and `after_siltingup` dump occur **once per outer step** (:1291–1292).
2. Soil `pre_run` (:1304), surface/ET/groundwater pre-runs, then obtain the module `substepcnt` (:1309).
3. For each module substep: soil `run` → `after_unsatzon` → surface routing → real-ET statistics → groundwater (`:1310–1323`).
4. Snow `run2` follows the last groundwater substep. `get_max_pumping`, soil `post_run`, the other post-runs, and `after_post_run` follow (:1326–1335).

The issue's proposed `after_siltingup → after_unsatzon` pair is therefore valid for the present **one-substep** setups only, and contains soil pre-run effects as well. It is not a before dump per module substep. A future multi-substep comparator must snapshot directly before each soil run, or reconstruct from the preceding full handover including groundwater. It must not pair a once-hourly before record with every after record.

`pre_run` (`unsatzon.cc:7732–7913`) binds MPI/head pointers lazily, sets date/time, shifts direct/interflow travel-time bins (:7850–7869), clears statistics and broadcasts river-flow availability. With heat off it does not initialize temperatures. Persistent travel-time storage already changes here.

Inside `run`:

- With surface interaction >1, `MPISyncRivLnkCells1` builds an ordered temporary river-link row before ordinary rows (:8077–8092). It is not merely MPI housekeeping: it runs in serial Doubs.
- Allocate row scratch once (:8111–8246), then visit increasing columns from FirstLastCol (:8248). Skip invalid/fully glacierized cells, except for the explicit clearing branch (:11383–11409).
- Reset per-cell totals, copy soil parameters and all persistent column states (:8281–8347).
- Compute the Courant candidate `dt_min` by depth, using relative conductivity when precipitation is zero; floor it at `dParameterMinStep`; choose `splitfactor = ceil(60*InV/substepcnt/dt_min)` (:8718–8756). This chooses a split count once for the cell call, not a continuously retried time step.
- `DT = 60*InV/substepcnt/splitfactor` is seconds; `timestep = 60000*InV/substepcnt/splitfactor` converts m/s to mm (:8792–8794). `number_of_substeps` records this **cell splitfactor**, not the module substep count.
- For each `splitcnt` (:8816), reset column scratch (:8990–9043), process snow bypass/pond input and pond evaporation (:9117–9231), macropores (:9237–9248), then vegetation layers in input order (:9251–9451).
- Apply the previous groundwater handover: distribute lateral inflow upward or drain layers; update the discrete GWLEVEL/start_level (:9538–9650).
- Reset theta iterates (:9858–9864). Each nonlinear iteration calls, in this order: `hydraulic_heads_iterative` → `specific_fluxes` (which calls `calc_surface_interaction`, `quflow`, then ordered `gradient_restriction`) → `flux_restrictions_opt` → increment `ziter` → `new_thetas` (:9873–9922). The first hydraulic call sees z=0; the first state-update call sees z=1. At z=21 the convergence threshold is multiplied by 5. The alternate unoptimized restriction helper is used only with >=200 levels or layered snow (:9897–9906), so neither reference calls it.
- Update continuous groundwater position and `GWtheta`, call `groundwater_recharge` (:10434–10576), compute infiltration excess and ponds, and accumulate all per-subinterval output fluxes (:10586–10799).
- After all cell splits, update groundwater storage coefficient S0 (:10890–10902), convert DELTAGW to interval mm (:11161), then write back levels 1..L of theta/heads/qu (:11321–11325).
- Scatter river-link results (:11700–11740), then reduce per-cell runoff into travel-time-zone and statistics arrays (:11898–12202). The final block is buffered control-cell output, not column physics.

`post_run` (`:13105–13402`) gathers river availability and zonal quantities, applies two persistent linear reservoirs, forms Qg for routing, then writes statistics/state outputs. Routing (#12) consumes Qg and the ordered river-flow table. Groundwater (#10) consumes per-cell groundwater distance and S0, not the zonal reservoir outputs or Qg. `get_max_pumping` is an irrigation-related call outside this physical scope.

## 3. State layout and explicit handovers

| Quantity | Storage in WaSiM | Persists? | Proposed ownership |
| --- | --- | --- | --- |
| theta, hh, hgeo, qu | Grid-table stacks, per cell and levels 0..L | Yes, copied in/out | ColumnState; define slot 0 and padded L+1 explicitly |
| GWLEVEL, GWSTAND, GWtheta | 2-D registered grids | Yes | ColumnState, discrete layer index + continuous distance + moisture must remain distinct |
| POND | Registered soil grid shared with ET/routing | Yes | PondState passed explicitly between owners |
| DELTAGW/baseflow | Registered soil grid overwritten by groundwater | Between calls, with changing units | Named flux handovers, no shared mutable grid |
| ETR per vegetation layer, EVAR | Registered grids owned by ET | Accumulate across module substeps; reset by ET pre-run | Explicit accumulated real-ET outputs |
| QDgrid, QIgrid, GWIN/GWEX, GWN, UPRISE, PERCOL, MACRO, INFEX, INFILT | Registered per-cell grids | Totals for the current module substep; reset before cell computation | ColumnFluxes, with units and accumulation interval |
| S0 | Groundwater parameter grid | Updated after soil cell computation | Explicit soil-to-groundwater coefficient output |
| QDd, QIi, zoneFZT_QD/QI | Member arrays indexed by zone, travel-time class and tracer | Yes; storage file reads/writes preserve them | ZonalRunoffState, outside column stencil |
| riverflow_table | Member structure / pointer setter to routing | Yes; availability and consumed amounts change in flow order | #8 graph state, explicit table interface |
| ST, LUtable, zone coefficients | Soil/vegetation/zone parameter tables | Read-only during column physics | Driver-owned immutable parameter tables |
| Per-row working arrays | Local pointers allocated at :8143–8230 | Across cells in the same row/call, not across run calls | Scratch; **h_temp1 exception below** |
| Per-iteration locals / splitfactor | Scalars in row/cell scopes | No physical persistence | Kernel/control scratch |

The primary interface should be `richards_column(ColumnState, ColumnParameters, SurfaceForcing, GroundwaterExchange, RiverExchange, dt) -> (ColumnState, ColumnFluxes, Handover, Diagnostics)`. A staged CPU reference may additionally need an explicit LegacyRowScratch for the history behavior until its reachability is settled. An ensemble dimension belongs to state/parameter fields, never to shared mutable scratch or one shared river-flow availability table. Distinguish the last accepted interface velocity (q_down) from accumulated transported water: multiplying the final velocity by the outer interval does not recover the sum over adaptive substeps. Return both final fields and explicitly accumulated transfers where transport/conservation needs them.

Foreign reads are P (selected rain/throughfall/snow outflow pointer), SNOW storage for bypass, residual ETP per vegetation layer, EVAP after snow evaporation, root depth and coverage, LAKE/MAXPOND/POND, prior baseflow exchange, groundwater head/distance/index, optional THAWDEPTH, river geometry and availability, ice fraction and cell area/zone codes. Binding and pointer selection are at :1384–1420,1580–1706 and `setInterceptionBeforeSnow` :7916–7935. P depends on upstream module choices; it must not be hard-coded as precipitation.

Foreign writes and units:

- ETR/EVAR accumulate in mm per module substep (:10789–10796); potential ETP/EVAP are read into local subinterval copies, not overwritten by ordinary Richards extraction.
- S0 is recomputed as drainable/fillable porosity with a 0.01 floor (:10890–10902). The groundwater port (wasim-gt#22) showed it cannot be read only at start.
- Groundwater's incoming `baseflow`/DELTAGW is m/s. Soil adjusts/distributes it and multiplies by `timestep*splitfactor` at :11161, so the same registered grid holds **mm after soil**. Groundwater later overwrites it with m/s again (`gwflow.cc:2327`). The port must use different named typed handovers.
- Soil changes its own GWLEVEL/GWSTAND/GWtheta, which groundwater reads as the initial head for its next lateral solve. `GWN` is diagnostic recharge in mm; groundwater binds but does not numerically consume that grid directly (wasim-gt#22).
- Surface-local exchange computes per-level gwin/gwex and totals. Doubs also increments `riverflow_table[link].infilt_amount` in flow order (:10407). #8 owns availability/ordering, while #7 consumes/returns the per-level exchange through a settled interface.
- Layered-snow writes to P, snow stores and EnergyBalance exist but are #9 and inactive here. Heat-derived frozen depth/colmation updates are also #9, not the simple #6 thaw-depth input.

`coupling.md` needs these additions/corrections in a later reviewed update: groundwater is Gauss-Seidel in both samples, not line relaxation; it does not read the river-flow table or GWN numerically; S0 and DELTAGW's phase-dependent units are missing; pre_run/FZT state changes and the two different substep counts are missing; `after_siltingup` is once-step; local macropore paths are reached and belong in scope; `qu[0]` is not written back by :11321–11325, so q_down level 0 is not a snapshot of the most recent surface flux; and the conditional row-history finding below invalidates an unconditional independent-column claim.

## 4. Audit of every row scratch allocation

All m_* storage has length `2*(levels+1+levelsSnow)` and is zero-initialized once per row (:8142 onward). Its alias is offset by levelsSnow. The current option set has levelsSnow=0, but the extra tail is still present. Most arrays are safely reloaded or assigned over every subsequently read physical index. Slots beyond L are not ordinary soil levels; zero sentinels must not become arbitrary host/GPU memory.

The table covers every allocation at :8143–8230 plus ETPsub/ETRsub. “No carry” is a source audit under the stated options, not a whole-model reordering experiment.

| Array / allocation line | First assignment before numerical use | Cross-cell finding |
| thetasat / :8143 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| thetares / :8144 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| dicken / :8145 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| theta_sl / :8146 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| hgeo_sl / :8147 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| hh_sl / :8148 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| qu_sl / :8149 | Full levels 0..L reloaded from tables/stacks at :8339–8347. | No carry in current options; extra tail stays zero. |
| temper_sl / :8154 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| a_etr / :8157 | Cleared through ST[0].getLevels() at :8303–8308; a_etr/a_evar cleared per vegetation layer at :9279–9283. | No carry. ST[0] holds maximum L (:1803–1804). |
| a_evar / :8158 | Cleared through ST[0].getLevels() at :8303–8308; a_etr/a_evar cleared per vegetation layer at :9279–9283. | No carry. ST[0] holds maximum L (:1803–1804). |
| etrdetail / :8159 | Cleared through ST[0].getLevels() at :8303–8308; a_etr/a_evar cleared per vegetation layer at :9279–9283. | No carry. ST[0] holds maximum L (:1803–1804). |
| a_dicke / :8160 | Cleared through ST[0].getLevels() at :8303–8308; a_etr/a_evar cleared per vegetation layer at :9279–9283. | No carry. ST[0] holds maximum L (:1803–1804). |
| theta_temp / :8161 | All 0..L set from theta_sl at :9037–9043 and :9858–9861. | No carry; iteration-local values intentionally persist within one solve. |
| h_temp1 / :8162 | Only k<=min(max(maxlevel,GWLEVEL)+2,L), k<L, at first z=0; :5220–5241. | Conditional retained history: a later extension via gwin/gwex can read stale row values at :5247. |
| theta_average / :8163 | All 0..L set from theta_sl at :9037–9043 and :9858–9861. | No carry; iteration-local values intentionally persist within one solve. |
| theta_sl_old / :8164 | All 0..L set from theta_sl at :9037–9043 and :9858–9861. | No carry; iteration-local values intentionally persist within one solve. |
| q_temp / :8165 | Zero allocation only. | No numerical reads found; dead scratch in this routine. |
| extractionweight / :8166 | Assigned and normalized for 1..maxlayer at :5073–5088 before :5143. | No carry over the range used by ET. |
| qlat / :8167 | All 0..L cleared at :9037–9040; gwin/gwex also cleared at :5297. | No stale carry; river-table side effects are separate explicit shared state. |
| gwin / :8168 | All 0..L cleared at :9037–9040; gwin/gwex also cleared at :5297. | No stale carry; river-table side effects are separate explicit shared state. |
| gwex / :8169 | All 0..L cleared at :9037–9040; gwin/gwex also cleared at :5297. | No stale carry; river-table side effects are separate explicit shared state. |
| interflow / :8170 | All 0..L cleared at :9037–9040; gwin/gwex also cleared at :5297. | No stale carry; river-table side effects are separate explicit shared state. |
| lateralinflow / :8171 | All 0..L cleared at :9037–9040; gwin/gwex also cleared at :5297. | No stale carry; river-table side effects are separate explicit shared state. |
| maxfromlr / :8173 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| qufromlr / :8174 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| iflfromlr / :8175 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| qlatfromlr / :8176 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| gwexfromlr / :8177 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| ratio / :8180 | Assigned by new_concentrations at :6917–6933. | Tracer preparation still executes with zero tracers; subsequent concentration use is absent. No water-state carry. |
| macrocapacity / :8181 | Assigned 1..k at :7344–7349 before reverse use at :7356–7394. | No carry on the used range; tail is not a capacity input. |
| macroinfiltration / :8182 | All 0..L cleared at :9042 before optional macropore writes. | No carry. |
| SE_sl / :8218 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| T_SE08 / :8219 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| qHeat_sl / :8220 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| lambda_sl / :8221 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| heatcap_sl / :8222 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| deltaT_sl / :8223 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| dEOverdT_sl / :8224 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| temper_sl_new / :8225 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| coeff_a / :8226 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| coeff_b / :8227 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| coeff_c / :8228 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| right_v / :8229 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| result_x / :8230 | Allocated only when do_TMod>0; heat/snow paths initialize/use it. | Not allocated/read by either reference option set; separate #9 audit required. |
| conc_lat, conc_gwin, conc_gwex / :8184–8186 | Concentration scratch passed to tracer helpers; index 0 not initialized here | No concentration reads affect water with tracer_count=0; do not generalize to a tracer port. |
| conc_sl_old / :8187–8190 | Each column array allocated/zeroed per row; tracer-specific copies precede tracer use | Tracer-only, excluded. |
| conc_sl / :8193 | Alias for each tracer reset at :8476–8479 | Tracer-only, excluded. |
| balance, maxbalance / :8195–8196 | Debug balance initialization/use is guarded by _WASIMDEBUG_ | Not physical state; excluded in these builds. |
| ETPsub, ETRsub / :8244–8246 | Every vegetation entry set at :8990–8994 before a split | No carry. |

Other row scalars were checked as well: startlayer is set to 1 at :8237 and not modified; splitfactor, splitfactor_heat and heatmultiple are cell locals (:8319–8322); extraSnowRunoff is reset at :8281. The unused constructor member scratch arrays near :3861 are shadowed by these row-local pointers. `specific_fluxes` separately allocates and copies ksat per call (:5474–5481). The local *_loc[200] arrays in `flux_restrictions_opt` copy the accessed prefix before use (:6199–6222) and copy outputs back (:6547–6556).

### Conditional history witness

**Yes, there is a conditional read of previous-cell scratch.** The important exception is h_temp1. Hydraulic initialization uses zero gwin/gwex, but `specific_fluxes` subsequently fills river exchange deeper in the column. On the next nonlinear iteration, the active-depth extension at :5222–5233 can reach levels whose h_temp1 was not written by that cell at z=0. :5247 then reads the row buffer. The groundwater-head cap at :5251–5257 often hides the effect, so a stale read alone is not sufficient evidence of a changed output.

A scratch standalone witness extracts `hydraulic_heads_iterative` **and the real van Genuchten lookup functions verbatim**. It first processes a deeper predecessor, then an otherwise identical current cell with 8 thin 0.1 m horizons, GWLEVEL=2, GWSTAND=-0.3 m, current theta=0.3, maximum ET level 2. Between z=0 and z=1 it activates gwex at levels 3–5 with 1e-8 m/s. Soil parameters are theta_s=.45, theta_r=.05, alpha=2, n=1.6, tau=.5. Only predecessor theta changes:

| Predecessor theta | Retained h_temp1[5] | Current hh[5] after extension |
| --- | --- | --- |
| 0.2 | 2.4841136363046044 | -2.1430799264643166 |
| 0.4 | 0.29404573335432127 | -1.0480459749891748 |

This proves a conditional helper-level dependency, **not that either 48-hour golden case reaches it**. The witness intentionally contains a transient mismatch between discrete GWLEVEL and continuous GWSTAND. Such a mismatch must be investigated at the real call sites: GWLEVEL is recomputed at :9590–9650, while GWSTAND is updated later at :10559. Reachability is settled by the real-run follow-up below.

If it is reachable, strict legacy identity needs explicit row-history/order or a sequential fallback for affected cells. Resetting or filling the buffer is a numerical change and must not be silently introduced as a cleanup. New arrays are allocated each run, so this scratch does not itself persist across module calls.

### Real runs (follow-up, 2026-10-19)

Instrumented local builds of 50c9ad8, one thread, no FMA: A unmodified, B a detector that only records (outputs identical to A), C the same as A but using the cell's own old-state suction where A reads a value it did not write.

| Run | Reads of an unwritten level | Reads that change a head | C against A |
| --- | ---: | ---: | --- |
| Doubs 48 h / 2021 | 42 / 264,875 | 0 / 0 | identical |
| Rhone 48 h | 2,361 | 0 | identical |
| Rhone 2021 (8,759 hourly steps) | 1,165,252 | 103, in 11 cells on 14 days | 195 of 1,213 output files differ |

- In all 103 effective reads the value comes from another column (or the row's initial zero), and the continuous groundwater elevation lies outside the discrete groundwater layer, so the cap does not apply. Rhone uses surface interaction 1: river reinfiltration is not needed for the effect.
- Size of the difference C − A in Rhone 2021: routed discharge up to 0.040 m³/s (0.11 %), annual volume at gauge 1 +938 m³ of 4.7·10⁹ m³; final theta up to 0.041 m³/m³ and groundwater head up to 0.53 m in single cells; annual actual ET unchanged in every zone.
- **Order dependence:** processing the columns of every ordinary row in reverse order (river row, reductions and all other order unchanged) changes 214 of 1,213 Rhone 2021 output files with A and none with C. Doubs and Rhone 48 h are unchanged by the reversal in both.
- Per-cell traces of the affected cells over the year: no negative or supersaturated theta, the same iteration counts in the event windows, the same sign changes of the fluxes in A and C. In the clearest case (72,27, 2021-11-02) A keeps a hydraulic head 8.1 m lower at a saturated level for three days. A per-call soil water budget leaves the same remainder (equal to −GWEX, up to 0.115 mm) in A and C; it is not caused by this behaviour.
- The comments at `unsatzon.cc:3681`, `8162` and `5238` describe `h_temp1` as the same cell's values of the last time step, reused only to save time; the manual's Eq. 2.14.8 averages a column's own water contents. No public WaSiM source newer than 11.01.01 was found.

**Consequence for the port.** With the fix (mmyrte/wasim-eth#3: levels a cell has not initialized in the current split are computed from its own `theta_sl_old`), WaSiM no longer depends on the column order, its outputs equal C in all four cases, and the 10-day goldens are unchanged. If the fix is accepted, wasim-gt validates against the fixed WaSiM and the soil columns are independent apart from the river graph (#8). If not, the port has to process the columns of a row in order.

Even if that path is unreachable in the current cases, full Doubs run is not a set of independent columns: `calc_surface_interaction` reads flow/speed/consumed amounts (:5393–5411), `run` updates consumed reinfiltration (:10407), and MPISyncRivLnkCells1/2 preserve flow order. The #8 graph phase remains outside parallel column launches. Zonal sums are also ordered shared reductions, not per-column scratch.

## 5. In-scope line set and coverage

The explicit set is reproducible using [richards/coverage.py](richards/coverage.py). [richards/line-set.json](richards/line-set.json) lists every retained executable line and all inclusive exclusion ranges with reasons. [richards/coverage.md](richards/coverage.md) lists **every unexecuted in-scope line**, grouped by a named regime and a reaching proposal; [richards/coverage.json](richards/coverage.json) is the machine-readable summary. Common control-flow guards remain in the conservative set even if one alternative is excluded. It is not a branch-level proof and does not estimate production code size directly.

Separate no-OpenMP coverage builds and counters were used. Whole run coverage reproduces the prior measurement: Doubs 687/2421 (28.38%), Rhone 840/2421 (34.70%). Narrowing to the stated physical/control option set gives:

| Function | Source range | In-scope executable | Doubs | Rhone | Union |
| --- | --- | --- | --- | --- | --- |
| pre_run | 7732–7913 | 62 | 60 | 60 | 60 |
| run | 7939–12296 | 654 | 476 | 555 | 555 |
| post_run | 13105–13402 | 41 | 34 | 40 | 41 |
| quflow | 4850–4853 | 2 | 2 | 2 | 2 |
| gradient_restriction | 4859–4965 | 66 | 65 | 66 | 66 |
| GetMeanTheta | 4969–5034 | 38 | 25 | 25 | 25 |
| evapotranspiration | 5038–5206 | 74 | 71 | 71 | 71 |
| hydraulic_heads_iterative | 5210–5264 | 30 | 30 | 30 | 30 |
| calc_surface_interaction | 5269–5457 | 111 | 104 | 71 | 104 |
| specific_fluxes | 5461–5832 | 100 | 78 | 78 | 78 |
| flux_restrictions_opt | 6174–6558 | 225 | 62 | 91 | 91 |
| new_thetas | 6563–6619 | 33 | 29 | 29 | 29 |
| groundwater_recharge | 6623–6703 | 32 | 21 | 21 | 21 |
| CalculateMacroporeInfiltration | 7320–7400 | 41 | 7 | 16 | 16 |

Total: 1189/1509 in-scope lines executed in the union. The difference is 320 documented in-scope gaps.

The numerical helper ranges above are retained except explicit negative-snow/heat, clay and drainage alternatives. `infiltration` :6707–6822 has no active call site in the stripped source: this older Green-Ampt routine is excluded, while actual method-3 surface infiltration, restrictions and macropores remain in scope. `flux_restrictions` :5836–6170 is excluded by the <200-layer, non-layered-snow guard. Heat helpers :7405–7729, updateSnowStorageLayers :12299–12507, ArcVeg :12511–13101, tracers :6869–7275, irrigation and HRU expansion are excluded. Graph gather/scatter :15893–17000 belongs to #8. Allocation, legacy files and diagnostic dumps are not numerical port work.

Large run exclusions include :8357–8716 (ArcVeg/layered snow/tracers/heat), :8818–8988 (legacy dump/replay), :9674–9855 (heat/snow), :9955–10430 (diagnostics plus river-table update), :10801–10888 (snow/ArcVeg), :10904–11158 (heat-derived grids), :11226–11318 (debug balance), :11411–11554 (glacier temperature search), :11557–11890 (diagnostics, deletion, graph scatter and HRU), and :12206–12296 (buffer output). Smaller exclusions are explicit in the line-set file. These exclusions do not remove the actual ET call at :9286–9295; comment markers embedded in // comments must not be treated as real block-comment delimiters.

| Gap regime | First, cheap copied-control variant | Next reference if needed |
| --- | --- | --- |
| Drought/wilting, root depth and ET redistribution | Deeper/shallow roots, dry/wet reduction thresholds and vegetation fractions | Sustained dry forcing and shallow-layer wilting column |
| Saturation, upward/downward flow, limiter cases | More intense precipitation, lower min step, tighter error and altered iteration limit | Synthetic water contents at residual/saturation with each sign combination of interface fluxes |
| Pond filling/spill and lake bypass | Vary max pond depth and lake/pond flags consistently | Saturated column receiving a known pulse; do not disable pond while leaving a lake model that requires it |
| Macropores | Lower PMacroThresh, vary MacroCapacity/Depth and reduction per metre in copied soil parameters | Positive rain pulse with only part of the pore capacity available |
| Rising/falling GW and new active levels | Vary initial groundwater state and module substeps | Upward/downward exchange crossing interfaces; extend to the h_temp1 witness with #8 |
| Simple frost input | Enable valid #6 thaw-depth/coefficient grids and reduced ksat setting | Synthetic thaw boundary through a horizon; heat/SE evolution remains #9 |
| River availability / deep contact | Vary colmation, river depth and inflow with #8 | Exhausted reinfiltration availability and deep-contact cell chain, preserving graph order |
| Daily steps and zonal routing | Vary travel-time bins/reservoir coefficients; prepare a 1440-minute control | Explicit daily forcing aggregation and two daily steps; never relabel hourly forcing as daily |

These variants are proposals for later coverage tasks, not runs claimed here. Phase 0 measured the two original 48-hour references and helper/control-flow experiments. No whole-column exactness or completed dry/daily/frost coverage is claimed.

## 6. Intermediate validation: evaluate both strategies

### (a) Selected-cell column trace

A trace is strongest for coupled helpers and the hidden-history question because it records the real intermediate inputs, active range and order. It must identify `(outer time, module substep, physical cell id, temporary river-row id/order, splitcnt, ziter, helper, before/after)` and preserve float64 bytes and explicit NaN/inf/nodata patterns. Log both computed bounds and initialized ranges. Include levels 0..L+1, h_temp1 on cell entry, theta iterates, qu[0], gwin/gwex, per-level ET/interflow, DT/timestep, maxlevel and current GWLEVEL/GWSTAND. Static soil/LU parameters need source/version hashes and one exact snapshot.

Proposed environment contract: WASIM_COLUMN_TRACE_DIR and WASIM_COLUMN_TRACE_CELLS=`row,col;row,col`. Resolve physical river cells through their mapping before matching the temporary row; a naive i,j filter would trace the wrong Doubs river cell. Limit records by step/split/iteration and fail visibly on output errors. Default-off must perform no extra numerical arithmetic. Binary arrays plus a versioned metadata index avoid precision loss and compiler-dependent struct padding. Capture existing local bounds without recomputing helper predicates.

Cost/risk: patch touches private row-local arrays and many call sites; the graph remaps addresses; tracing can perturb timing/order in threaded builds, so use only one thread. A future implementation requires a maintainer review and unchanged-source-versus-hook Doubs full golden check at zero tolerance, both with the trace disabled and enabled. No such hook exists yet, so no equivalence test is claimed.

### (b) Verbatim helper extraction — recommended immediate reference

[richards/extract_reference.py](richards/extract_reference.py) copies `new_thetas` :6563–6619 and `quflow` :4850–4853 unchanged into a tiny C ABI. It records source/function SHA-256 hashes and exact ranges. The adapter supplies new_thetas' only required members: levels and GWLEVEL[i][j]. Build with `-O3 -ffp-contract=off`; allocate a defined zero tail L+1, matching the row allocation. No hydrological data or extracted legacy source is committed with the script.

This extraction is the recommended oracle for phase-1 isolated helper contracts. Working evidence: 1024 random 20-level columns and 16384 random 15-level columns, varying iteration number, weight, active depth, wet/dry clipping and exchange tails, including explicit NaN/infinity patterns. Three recorded-derived cells per setup also match. Those six cells replay a hypothetical one-second helper update from final recorded theta/flux/geometry, with zero per-level river exchange; they are **not a trace of the actual new_thetas inputs** and do not validate the assembled solver. Locations: Doubs (9,157), (117,109), (205,25); Rhone (4,252), (101,143), (180,55). The sample values are not part of this repository (no hydrological data).

Extraction dependencies differ markedly:

| Helper package | Non-argument context to make explicit | Extraction suitability |
| --- | --- | --- |
| new_thetas / quflow | L, GWLEVEL / none | Excellent; working prototype |
| hydraulic_heads_iterative | GWLEVEL/GWSTAND, ST lookup, dSE effectiveness, h_temp1 initialized range | Good once history is explicit; witness uses real extracted lookup functions |
| GetMeanTheta / ET | Soil/LU parameters, L, soil moisture/root geometry and ordered vegetation extraction | Good; reuse exact van Genuchten reference |
| gradient_restriction | Soil lookup, dSE constants, neighbouring layer state and mutable interflow | Good, but preserve 51-step sign/delta search |
| flux_restrictions_opt | Groundwater/surface flags, bounds, 200-entry scratch prefix, sequential mutations | Medium; must preserve top/bottom order and all side effects |
| macropores | Soil threshold/capacity/depth table and mutable p/theta | Good after table expansion |
| calc_surface_interaction | River geometry, lake/pond state, availability table, zone mapping | Medium/poor as isolated physics; #8 context is essential |
| assembled run / zonal phase | Many persistent grids, table setters and ordered reductions | Prefer whole-process dumps plus selected-cell trace |

Recommendation: extraction first for phase 1, combined with the retained process dumps. Require a selected-cell trace or an equally complete call-state capture before phase-2 acceptance, especially to settle the h_temp1 path and local/graph river interface. Phase 0 prototyped extraction, not the optional hook. A future trace is a separate reviewed change, not permission to rewrite legacy arithmetic.

## 7. GT4Py feasibility evidence and limits

The scratch prototype is not part of this repository. `new_thetas` uses a K interval over levels 1..L, K-1 interface flux, and per-column active-depth masks. Host preparation reproduces the legacy active-bound extension and weight update. The stencil reproduces theta clipping, ordered moving-average operations and error differences. The host max uses positive comparisons so NaNs behave like the C++ `if(temp>th_error)` check. This is a helper prototype, not a full Richards state/API implementation.

| Batch | Backend | First call s | Warm kernel ms | Exact against extracted new_thetas |
| --- | --- | --- | --- | --- |
| 1024 × 20 | gt-cpu_ifirst | 0.008722 | 0.349548 | Yes, theta, average, max error and patterns |
| 1024 × 20 | numpy | 0.004475 | 0.515956 | Yes, theta, average, max error and patterns |
| 16384 × 15 | gt-cpu_ifirst | 0.010379 | 4.908408 | Yes, theta, average, max error and patterns |
| 16384 × 15 | numpy | 0.008607 | 4.393885 | Yes, theta, average, max error and patterns |

Warm timing averages 20 calls, excluding field construction and comparison, with theta/average restored before each call. First calls can include cached compilation/loading; large-batch first calls reuse the stencil cache and are not cold compiler benchmarks. In the final runs, compiled CPU is faster for the 1024-column batch, while NumPy is faster for the 16384-column batch. These helper timings do not establish an overall Richards speedup or GPU throughput.

Additional tests pass on both backends: per-column split counts 1..32 inside a masked while, variable active depth, ordered FORWARD prefix and BACKWARD suffix sweeps. The toy NumPy masked-while failure seen in earlier work was **not reproduced** with these kernels. This does not establish NumPy as a reference for a larger masked solver or runtime-indexed van Genuchten tables.

Composition tests on the installed frontend:

- `while` around `with computation(FORWARD)` is rejected with AssertionError. Computation blocks must remain the outer DSL structure. A passing per-point while alone does not prove a full adaptive column fits in one stencil.
- A write to the preceding K level (`value[0,0,-1] = value` within a forward interval) compiles and its one-step shift is exact. Alias/halo constraints for the more complex flux limiter still need package-specific tests.
- A complete nonlinear iteration needs full-column maximum/error and mixed forward/backward operations. The initial safe design is a host/program loop calling compiled helper stages, maintaining column-active masks, with GT4Py prefix/reverse passes for vertical reductions and ordered restrictions. Measure dispatch cost before fusing. A cell-dependent while must not accidentally become independent level-dependent convergence.
- Groundwater depth should mask a fixed physical K range. Keep surface slot 0 and the legacy zero sentinel L+1. Do not equate a soil horizon code with the physical K index; expand table parameters explicitly or use the already ported lookup representation.
- Runtime data-dimension indexing remains unavailable on Cartesian numpy (coupling.md); soil lookup calls in conditions must first be assigned to temporaries, and gtscript does not support // integer floor division. Preserve arithmetic grouping and explicit no-FMA flags. The effective Cartesian C++ extra args must include `-ffp-contract=off`.

Divergent split/iteration counts cost work, and preserving the row-history dependency (§4) would additionally require a sequential fallback; the proposed WaSiM fix removes that need. GPU measurements of a column proxy follow below.

**Ensemble throughput on a CPU (2026-10-19).** A first measurement with the ported snow model (method 7, part 1: one `gt4py.cartesian` stencil with per-cell `while` loops, the closest existing analogue of a Richards column), Rhone grid, hourly step, Apple M-series CPU with 5 performance and 10 efficiency cores, warm calls, ms per ensemble member and step:

| | 1 thread | 15 threads |
| --- | ---: | ---: |
| WaSiM (`WASIM_TIMING`, one run / 15 runs in parallel) | 5.4 | 0.36 |
| port, 1 member (`gt:cpu_ifirst`) | 4.9 | 0.66 |
| port, 16–64 members in one call (`gt:cpu_ifirst` / `gt:cpu_kfirst`) | 5.5–8.9 / 6.3–6.9 | 1.45 / 0.73–0.87 |

Per cell the port is as fast as WaSiM. On a CPU an ensemble in one process does not beat independent WaSiM processes, which are already embarrassingly parallel over runs. Any speed-up of the port therefore has to come from GPUs (or from avoiding WaSiM's per-step I/O and statistics); this needs a measurement on GPU hardware (#17) before the go/no-go decision. On the same machine the Richards zone takes 45 ms of 87 ms per Rhone step (52 %).

**GPU measurements of a column proxy (2026-10-20).** Three sessions on a desktop PC (NVIDIA RTX 3080 10 GB, Intel i7-10700, WSL Ubuntu 24.04, CUDA 12.9, gt4py 1.2.3, dace 2.0.0a10). The proxy is synthetic: it runs the stages of one nonlinear theta update with the split and iteration counts recorded from WaSiM (§4 workload), not adaptive convergence, and it omits the repeated limiter passes, the gradient search, ET, macropores and the couplings. It is not a Richards port and gives no full-model speed-up.

- *Reference on the same PC.* Eight parallel WaSiM runs need 15.91 ms per member and hourly step for the `after_siltingup → after_unsatzon` segment over January; the ordinary-column share gives a **column budget of 15.04 ms**. About 18 ms for a full year is a Mac-to-PC scaling assumption, not a measurement.
- *Snow as the pointwise reference.* With fields resident on the device, the snow port runs at 0.14–0.17 ms per member against 0.86 ms for eight parallel WaSiM runs (5–6×). Calls that allocate and copy fields every step are slower than WaSiM.
- *Lockstep.* Over 2021, a Rhone column needs 2.1 iterations per call on average, but a step needs 117 global rounds at the median and 733 at the worst, because every round runs over all columns until the slowest has finished. Only 1.3 % of the column work in a round is useful. A masked column still costs about 21 % of an active one.
- *Hydraulics.* The port's `h_calc`/`k_calc`/`get_theta` recompute the table entries a query needs. WaSiM's tables depend only on soil type and layer (50 types × 15 layers for Rhone, 0.74 MB; 7 × 20 for Doubs, 0.14 MB), so they can be built once. GT4Py reads them as K-only fields with a flattened data dimension and a runtime soil index. The lookup is bit-identical to the CPU port on `gt:cpu_kfirst`, `gt:gpu` and `dace:gpu` when the GPU build disables fused multiply-add (`--fmad=false`); with the default build the difference reached 1.4e-12. Pointwise, the lookup costs 0.095 ms per member against 5.8 ms recomputing and 0.82 ms for the analytic formulas (`gt:gpu`, 15 levels, 64 members).
- *Proxy step times* (ms per member and step, `gt:gpu`, 16 members, fused driver):

| | Rhone median | Rhone worst | Doubs median | Doubs worst |
| --- | ---: | ---: | ---: | ---: |
| recomputing tables | 61.6 | 122.0 | 60.8 | 132.3 |
| precomputed tables | 7.2 | 39.0 | 3.3 | 16.1 |
| precomputed tables, compacting driver | 17.6 | 32.5 | 6.7 | 12.7 |

- *Annual estimate.* Interpolated over the global rounds of all 8,759 hourly steps of 2021, the fused proxy with precomputed tables averages about 11 ms (16 members) and 9 ms (64 members) for Rhone, and 4–5 ms for Doubs; `dace:gpu` is 20–45 % slower. The estimate rests on four measured steps; a regression over them did not pass its own quality check for Rhone and gave a range of 7.5–18.8 ms.
- *Compaction.* Gathering the unfinished columns into dense buffers and scattering them back costs more than it saves (Rhone, 16 members: 223 of 282 ms are data movement) and wins only in the worst steps. This rules out that implementation, not compaction in general.

So precomputed tables remove the hydraulics cost, and on this card the proxy column reaches 1.4–1.6× (Rhone) and 3.2–3.6× (Doubs) the throughput of eight parallel WaSiM runs, measured against the column budget. The rest of the column's work is not included, so the margin is not a decision. The remaining limit is the global lockstep. GT4Py does not allow a `while` around a computation over K, and it synchronises after every frontend call (126 calls per median Rhone step). A column kernel that iterates each column to convergence in its own thread would remove both; testing that, for example in DaCe directly, is the next measurement. Throughput on a data-centre GPU must be measured, not scaled from the RTX 3080's FP64 rate.

## 8. Proposed phases and task contracts

Estimates use coupling.md's measured 100–150 new context tokens per closely read/rewritten/validated C++ line plus roughly 10k fixed validation setup per independently delivered package. They are planning ranges, not measured tokens for this task. Ranges overlap at adapters; do not add overlapping raw line spans as independent production code counts.

| Task | C++ scope | Pure function inputs → outputs | Reference and dependencies | Estimate | Main risk |
| --- | --- | --- | --- | --- | --- |
| 1A ET extraction | :4969–5206; :9251–9451; :10747–10754 | soil parameters, theta, hgeo, root/coverage and remaining ET → theta, per-level ET, actual root extraction, ETR/EVAR | Verbatim helpers + bounded recorded-derived cases; exact soil lookup already available; independent after interface agreement | 55–85k | Ordered vegetation layers, saturated extraction redistribution, dry/wet thresholds |
| 1B Fluxes and restrictions | :4850–4965; :5461–5832; :6174–6558 | theta iterates, heads, geometry, conductivity, source/exchange arrays → ordered interface flows/interflow and restricted exchanges | Extraction with explicit context; trace deep-contact and both flux signs; can run alongside 1A/1C | 100–160k | In-place neighbour mutations, backward limiter iterations and a legacy assignment-in-condition |
| 1C State update / hydraulic heads | :5210–5264; :6563–6619 | theta iterates, fluxes, GW position, active bounds, history → theta, averages, heads, weight, error and initialized-range metadata | Phase-0 extraction/witness; resolve h_temp1 before claiming independent columns | 30–50k | Growing active range and retained history |
| 1D Surface, infiltration, ponds/macropores | :7320–7400; :9097–9250; :10586–10743; local :5269–5457 interface | precipitation/snow bypass, pond/lake/root state, soil/macropore parameters and explicit river availability → pond, infiltrated water, excess runoff, macro flux, exchange requests | Extraction for macropores; process snapshots and #8 contract for rivers; coordinate with 1A/1B | 65–100k | Do not duplicate #8 graph exchange or double-apply pond evaporation |
| 1E Lower boundary / groundwater | :6623–6703; :9538–9651; :10434–10576; :10890–10902 | previous lateral exchange, theta/fluxes, discrete+continuous GW state → updated GW state, recharge/up-rise/percolation, interflow/baseflow and S0 | Extraction + task-03 exact handovers; depends on 1B/1C field/sign contracts | 45–75k | Mixed mm/m/s units, layer crossing and groundwater depletion |
| 2 Assembled adaptive column | :7939–11325 remaining selected lines; glacier clear :11383–11409 | full ColumnState/forcing/parameters + explicit graph exchange → all states, fluxes, S0/ET handovers and split/iteration diagnostics | Whole-process dumps per true module substep, plus selected-cell call trace; depends on all 1A–1E and history/order decision | 200–350k | Frontend loop composition, exact stopping order, history and host/device dispatch |
| 3 Zonal runoff and couplings | :7732–7913 selected; :11898–12202 selected; :13105–13234; getters :4698–4835 | per-cell runoff + areas/FZT + persistent delay/reservoir state + glacier contribution → delayed QD/QI/QB, Qg and new zonal state | Ordered scalar reference and zone snapshots; depends on 2, #8/#10 handovers, #12/#14 interfaces | 70–120k | Non-associative weighted reductions and storage update timing |
| 4 Coverage completion / daily | Gaps in the retained line set | Variant controls/forcing and synthetic column states → exact per-helper and per-substep checks | #17 fixtures after 1 packages, full daily acceptance after 2/3 | 100–200k | Rare limiter/ET/GW states; correct daily forcing aggregation |

Total rough planning envelope is 0.7–1.2 M tokens, consistent with the issue's 0.8–1.5 M uncertainty after setup/trace work. The graph port (#8), heat/snow (#9), driver implementation (#14) and any separately approved legacy correction are not included. Extraction/call-capture infrastructure can add 15–35k; conditional-history resolution can add a further 20–50k. Parallel helper work is possible only after common array/unit/ownership contracts are fixed; assembled state must have one owner.

Acceptance should progress helper exactness → intermediate real call-state exactness → every active cell/output of every module substep → ordered zonal outputs. Always compare shapes, NaN, positive/negative infinity and nodata patterns explicitly, with a finite-value target of zero. Keep stage-specific conservation checks separate from the exact legacy oracle: a legacy quirk is not permission to change the arithmetic until a separate fix is approved.

## 9. Decisions before production porting

- **Maintainer (mmyrte):** Accept or reject the `h_temp1` fix (mmyrte/wasim-eth#3); reachability and order dependence are now established (§4). Agree the #7/#8 cut through calc_surface_interaction and limiter exchange restrictions.
- **Maintainer and #14 owner:** Agree physical K layout (0,1..L,L+1), table/horizon mapping, static vs evolving parameters, initialized-history metadata, ensemble shape, masks, immutable state/flux ownership, units at each handover and deterministic reduction order.
- **Port team:** Turn the phase-1 packages into individual tasks after the interface gate, select a trace/capture implementation for phase 2, and review the conservative line set. Keep method-3 macropores in scope despite the old estimate.
- **Project lead:** Decide the milestone budget and the go/no-go criteria at the #1 decision point: CPU bit identity first, target-domain performance next, and GPU tolerance only after device measurements. Confirm whether daily runs are required before that decision, and decide, after a per-column kernel test, whether the column is written in GT4Py or DaCe (§7).
- **Joint validation:** Keep the known Rhone full-model golden discrepancy separate from isolated helper/column checks. Do not overwrite goldens. Define how S0, baseflow phase units, river availability and zonal stores are compared and persisted.

No trace hook is proposed as an implemented diff; the trace specification above is a design for a later reviewed change.

## 10. Reproduction

This directory contains the document, the extraction and coverage scripts and the source-line documentation, without hydrological data. Build WaSiM single-threaded and without FMA (`make MPI=0 ARCH=-ffp-contract=off`) and run the Doubs and Rhone references as in `golden/README.md`.

```sh
# extract the unchanged helpers into a standalone C ABI and compile them
python3 gt4py/richards/extract_reference.py path/to/wasim-eth/src/unsatzon.cc /scratch/helpers.cc
clang++ -O3 -ffp-contract=off -shared -fPIC /scratch/helpers.cc -o /scratch/helpers.dylib

# line coverage of the in-scope set, from a no-OpenMP coverage build (--coverage, OMPFLAGS=) run on both references
python3 gt4py/richards/coverage.py path/to/wasim-eth/src <doubs coverage dir> <rhone coverage dir> <output dir>
```

Run coverage collections one after another: a coverage binary shares one counter directory. The extraction records function hashes of the source it read.
