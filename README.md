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
