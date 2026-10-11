# Golden test cases for WaSiM

Control file templates for the golden tests in [wasim-eth](https://github.com/mmyrte/wasim-eth) (`tests/golden/`). Each case runs 2021-01-01 to 2021-01-11 at an hourly step.

| Case | Setup | Active cells | Run time, 1 thread | Peak memory |
| --- | --- | --- | --- | --- |
| `doubs/` | `sample-run-doubs.zip` (not uploaded yet, see below) | 6,022 | 20 s | 187 MB |
| `rhone/` | `sample-run-rhone.zip` | 21,584 | 49 s | 324 MB |

The templates use the placeholders `@CASE_DIR@`, `@GOLDEN_DIR@` and `@OUT_DIR@`, which `run_case.py` fills in. Their headers list what differs from the original setups. For Rhone, temperature and precipitation come from station regression, because the gridded inputs are not in the sample. Create the merged station files once:

```sh
python3 golden/rhone/merge_station_files.py path/to/sample-run-rhone
```

## Reference outputs

The reference outputs are not committed yet: this session could not upload LFS files. To recreate them, build wasim-eth at commit `a9b28ae` (11.01.01 sources plus a platform line in the banner), then run each case single-threaded:

```sh
python3 tests/golden/run_case.py --binary <a9b28ae build> --threads 1 \
    --template golden/rhone/control.wasim.in --case-dir path/to/sample-run-rhone --out golden/rhone/reference
```

The originals were built with `src/makefile_macos MPI=0 ARCHFLAGS=` on Linux, clang 18, `-O3`. Repeated single-threaded runs of that build were bit-identical, so a rebuild with the same compiler should reproduce them exactly. The wasim-eth branch `claude/amazing-allen-a17bpw` as of 2026-10-11 (bug fixes, unified Makefile) reproduces them exactly with both g++ 13 and clang 18.

## Benchmarks

`rhone/year_hourly.wasim.in` and `rhone/year_daily.wasim.in` run the Rhone sample for all of 2021, at hourly and daily time steps. The daily one needs daily inputs, aggregated from the hourly station files (means; sums for precipitation):

```sh
python3 golden/rhone/hourly_to_daily.py golden/rhone/daily \
    --sum golden/rhone/precipitation_merged.dat \
    --mean golden/rhone/temperature_merged.dat cases/sample-run-rhone/Input/timeseries/{windspeed_2025,radiation_2025,sunshine_2025,humidity_2025,inflow221}.dat
```

Measured on 4 cores (clang 18, `-O3`): hourly about 30 min on 1 thread (extrapolated from 10 days) and 17 min on 4 threads; daily 401 s on 1 thread and 119 s on 4 threads (`OMP_SCHEDULE=dynamic,4`). See `plan-switzerland.md`.
