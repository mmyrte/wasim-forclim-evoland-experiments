# Validating wasim-gt against legacy WaSiM

wasim-gt (the GT4Py reimplementation, separate repository) contains no legacy code and no legacy formats. Correctness against WaSiM is checked here. The CPU target is bit-identical results; GPU runs get a stated tolerance.

## Parts

| Part | State | Where |
| --- | --- | --- |
| Grid dump hook: WaSiM writes named grids before and after each process call | done | wasim-eth `src/modapuzr.cc` (`WASIM_DUMP_GRIDS`) |
| Per-process comparisons against the dumps | done for interception, snow evaporation, T-index snow, radiation correction, potential ET, table-based phenology; soil functions against WaSiM's C++ | `gt4py/*.py`, `gt4py/soil/` |
| Coverage of the WaSiM code a reference run executes | done | `wasim_coverage.py` |
| Harness running legacy WaSiM and wasim-gt on the same case from one TOML configuration, comparing process by process and at the outputs | planned | here; uses the legacy tools (`../legacy/`) to write the WaSiM control file |
| Reference cases | Doubs and Rhone (real); first synthetic variant: `../gt4py/doubs_et_june.wasim.in` (land-use parameters set to reach the surface-resistance branches of ET) | `../golden/`, below |

## Reference cases: real and synthetic

A port can only be checked on the code paths its reference run goes through. Doubs and Rhone cover the variants switched on in them, but not the rest of wasim-gt's scope. The simple permafrost model, for instance, is off in both because they lack its thaw-coefficient grid. They also miss rare paths within the processes they do use.

**Synthetic, small cases can cover most of it.** They test the translation, not hydrological realism, so they only need to be consistent enough for WaSiM to run them:

1. **Small grids cut from the real setups.** A headwater subbasin of 10×10 to 30×30 cells keeps the routing topology valid. For the soil model, a few columns suffice.
2. **Switch on the variant** in the control file and add the grids it needs with synthetic values in their documented ranges: an α grid for permafrost, a glacier fraction, drainage, clay layers, macropores.
3. **Synthetic forcing** that drives the regimes the code distinguishes: freeze and thaw, snowfall and melt, heavy rain up to saturation and ponding, drought down to the wilting point. Days to weeks are enough. Long-memory processes (groundwater, glacier dynamics with its annual update) need years, which is cheap on a small grid.
4. **Coverage as the acceptance criterion.** Build WaSiM with coverage (`wasim_coverage.py`, usage in its header), run the case, and check that every line of the ported branch ran. A case is good enough when nothing in the branch is left unexecuted.

On the Doubs golden run (2 days, hourly), coverage of the functions ported or next in line:

| Function | Lines run | Note |
| --- | --- | --- |
| interception, method-1/11 branch (`interzep.cc:1376–1465`) | 21 of 39 | all of method 1/11 ran; the rest is methods 2/12 and ArcVeg (not ported) |
| snow evaporation (`etpetr.cc`) | 38 of 103 | energy-balance branch only; the conventional branch needs the T-index case |
| Richards zone, `unsatzonmodellklasse::run` | 600 of 2,058 | includes heat transport, layered snow, tracers and irrigation, all off; the branch to port has to be measured separately |
| glacier dynamics (`GlacierRetreatAndAdvance`) | 0 of 88 | no glaciers in Doubs, and no year boundary in 2 days |

For the potential ET port, over four runs (Doubs golden, first 2 days; `../gt4py/doubs_et_june.wasim.in`; Rhone daily January and June, `../golden/rhone/daily_*.wasim.in`):

| Function | Lines run | Note |
| --- | --- | --- |
| `verdunstungsmodellklasse::getETP` | 109 of 110 | the missing line is net radiation as input (not ported) |
| `verdunstungsmodellklasse::potentiell`, up to the other ET methods (`etpetr.cc:1743–2111`) | 104 of 138 | missing: modelled radiation or sunshine duration (not ported), ArcVeg (not ported), the 1260 Wh/m² cap, the albedo without lake model, SSDTemp without a bright step (unreachable) |
| `expokorrklasse::run` | 178 of 193 | missing: interval counter > 1 (not ported), sun below 1° in a bright sub interval, temperature for `split` > 1, sunset after midnight (not at Swiss latitudes) |

**What still needs real setups:** performance and scaling, calibration, and any variant whose parameters are hard to invent consistently. A further real setup is most useful if it exercises something the Swiss setup will use and Doubs and Rhone do not, for example heat transport and permafrost in an alpine catchment.
