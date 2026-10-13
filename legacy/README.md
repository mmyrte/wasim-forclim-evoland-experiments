# Tools for legacy WaSiM inputs (planned)

wasim-gt reads only modern formats: rasters through GDAL, netCDF and zarr, TOML configuration. Everything that converts from or to legacy WaSiM lives here, so that the new codebase carries no legacy code.

| Tool | Direction | Needed for |
| --- | --- | --- |
| control file writer | TOML → WaSiM control file | running legacy WaSiM from the same configuration as wasim-gt: validation harness (`../validation/`), production runs until wasim-gt is complete; replaces the regex edits and the module adapter (`plan-switzerland.md`) |
| control file reader | WaSiM control file → TOML | importing parameters from the 17 inherited regional setups, if they are still needed |
| land-use and soil tables | WaSiM tables → TOML tables | parameters for wasim-gt |
| grids | WaSiM ASCII and binary grids → GeoTIFF or zarr | inputs; GDAL already reads the ASCII grids, the binary format needs a reader |
| station tables | WaSiM station files → netCDF or zarr | meteo input, observed discharge |
| states | WaSiM state grids and storage files (`StateIni_*`) → netCDF or zarr, and back | starting wasim-gt from WaSiM states; comparing states |

Each tool gets a round-trip test on the golden setups: converting and converting back must reproduce the original files, or the original WaSiM results.
