# Integrated land surface + land use change modelling

This repository is for experiments in coupling
[WASIM](https://www.wasim.ch/),
the Python version of [ForClim](https://ites-fe.ethz.ch/openaccess/products/forclim),
and [evoland-plus](https://github.com/ethzplus/evoland-plus/) 
at the Swiss national scale.
The goal is to describe landscape-level response and ecosystem services under climate change conditions, including some climate extremes.
The intention is that potential feedbacks between hydrology, forest structure, and land use could be described by a two-way coupling:
on the one hand, the biophysical models run at _monthly_ intervals
while on the other, the land use change model is evaluated at _decadal_ time steps.

- The biophysical, mechanistic models receive updated land cover data.
- The land use model receives explanatory variables to allow statistical estimation of land use transition potential.

The biophysical side is thus designed to inform land use change processes via metrics describing ecosystem state.
Ecosystem state can be useful as proxies for ecosystem services, i.e. ecosystem functions that are valued by society, and thus indicative of land use suitability.
Example state variables are soil moisture, streamflow, plant productivity, etc.
The monthly timestep is sufficient to reproduce drought stress events and plant productivity faithfully.
Monthly timesteps will not be sufficient to capture anything less than the most severe flooding events.

## To Dos

- evoland -> WASIM
  - [ ] How do land use classes map to land cover classes?
- WASIM -> evoland
  - [x] 500m resolution great
  - [x] monthly timestep great
  - [ ] Do we need to recalibrate WASIM for monthly timesteps? Compare monthly run to 2018 drought data?
  - [ ] Identify explanatory variates that could be useful
    - Plant available water
    - annual biomass yield
    - ...
  - [ ] Plant productivity: evaluate against yield data for crops?
  - [ ] Plant mortality from flooding - is there precedent for this in CH?
- WASIM + ForClim integration
  - Find interfacing options; avoid redundancies. ForClim has some hydro routines. WASIM has some dynamic plant phenologies.

## Layout

- `plan-switzerland.md`: working plan for the harmonized Swiss WaSiM setup.
- `golden/`: golden-test cases (Doubs, Rhone) and benchmarks for legacy WaSiM ([wasim-eth](https://github.com/mmyrte/wasim-eth)).
- `validation/`: checking wasim-gt (the reimplementation, separate repository; GT4Py so far, moving to DaCe and to be renamed wasim-dace) against legacy WaSiM. Includes reference cases, the coverage tool, performance comparisons and the CPU spike results.
- `gt4py/`: per-process checks of wasim-gt against WaSiM grid dumps, plus the analysis of WaSiM's module coupling; `gt4py/spike/`: the CPU spike comparing C++, DaCe, gt4py.cartesian and gt4py.next.
- `dace/known-issues/`: reproducers of DaCe defects met on the way, usable as regression tests.
- `legacy/`: tools that convert legacy WaSiM inputs, states and control files (planned).
