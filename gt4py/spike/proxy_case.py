"""Synthetic case for the Richards proxy: columns with the Doubs soil table, rain and ET forcing.

The proxy follows the structure of WaSiM's Richards zone (unsatzon.cc, run, method 3):
per cell a Courant sub-step count, and per sub-step an iteration of heads, Darcy fluxes,
flux restrictions and new water contents with WaSiM's moving average, until the largest change
is below maxerror (5 times that after 20 iterations) or max_iterations is reached. It is not
WaSiM's Richards zone (no groundwater coupling, macropores, interflow, drainage, frost or heat),
but it does the same kind and amount of work per cell, so it can show how each formulation
handles the per-cell loops. See README.md.

Every implementation reads the case written here and writes its results, so they can be
compared: richards_ref.cpp (C++ reference), richards_dace.py, richards_next.py, richards_cartesian.py.

usage: proxy_case.py OUT_DIR [--cells N] [--steps T] [--seed S]
   writes OUT_DIR/case.npz and the same arrays as raw float64 files for the C++ reference
"""
import argparse
import os

import numpy as np

# Doubs [soil_table] (wasim-forclim-evoland-experiments golden/doubs): per horizon
# ksat (m/s), theta_sat, theta_res, alpha (1/m), n, layer thickness (m), number of layers
SOILS = [
    [(4.0e-7, 0.43, 0.078, 1.47, 1.25, 0.20, 10), (3.0e-7, 0.45, 0.067, 1.47, 1.25, 1.5, 10)],
    [(1.0e-5, 0.40, 0.080, 1.6, 1.42, 0.10, 10), (8.0e-6, 0.38, 0.110, 1.6, 1.42, 0.15, 10)],
    [(1.0e-5, 0.40, 0.080, 1.6, 1.42, 0.20, 10), (8.0e-6, 0.38, 0.110, 1.6, 1.42, 0.30, 10)],
    [(1.0e-5, 0.40, 0.080, 1.6, 1.42, 0.20, 10), (8.0e-6, 0.38, 0.110, 1.6, 1.42, 1.5, 10)],
    [(2.0e-5, 0.40, 0.080, 3.6, 2.1, 0.10, 10), (1.4e-5, 0.38, 0.110, 3.6, 2.1, 0.15, 10)],
    [(1.0e-4, 0.37, 0.080, 4.0, 3.0, 0.20, 10), (5.0e-5, 0.30, 0.10, 4.0, 3.0, 1.5, 10)],
    [(3.0e-7, 0.45, 0.067, 2.0, 1.41, 0.10, 9), (1.0e-7, 0.45, 0.067, 2.0, 1.41, 0.30, 2),
     (1.0e-7, 0.40, 0.095, 1.9, 1.31, 2.0, 9)],
]
LEVELS = 20
PARAMS = ("ks", "ths", "thr", "alpha", "inv_n", "m", "inv_m", "dz", "zc")


def soil_columns():
    """(n_soils, PARAMS, LEVELS) array."""
    out = np.zeros((len(SOILS), len(PARAMS), LEVELS))
    for s, horizons in enumerate(SOILS):
        rows = [h for h in horizons for _ in range(h[6])]
        assert len(rows) == LEVELS
        for k, (ks, ths, thr, alpha, n, dz, _) in enumerate(rows):
            m = 1.0 - 1.0 / n
            out[s, :8, k] = ks, ths, thr, alpha, 1.0 / n, m, 1.0 / m, dz
        dz = out[s, 7]
        out[s, 8] = -(np.cumsum(dz) - 0.5 * dz)  # centre of each level below the surface (m)
    return out


def field_capacity(cols, h=3.5):
    ks, ths, thr, alpha, inv_n, m = (cols[:, i] for i in range(6))
    return thr + (ths - thr) * (1.0 + (alpha * h) ** (1.0 / inv_n)) ** (-m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--cells", type=int, default=6022)
    ap.add_argument("--steps", type=int, default=48)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    n, t = args.cells, args.steps
    soils = soil_columns()
    soil = rng.choice(len(SOILS), size=n, p=[0.2, 0.15, 0.2, 0.15, 0.1, 0.1, 0.1])
    cols = soils[soil]                                     # (cells, PARAMS, LEVELS)
    theta = field_capacity(cols) * rng.uniform(0.8, 1.0, size=(n, LEVELS))
    theta = np.maximum(theta, cols[:, 2] + 0.01)
    # rain: two events with a spatial pattern (mm per hourly step), ET in daylight hours
    hours = np.arange(t)
    event = np.where((hours % 24 >= 3) & (hours % 24 < 9), 4.0, 0.0) + np.where((hours % 24 >= 15) & (hours % 24 < 17), 15.0, 0.0)
    event[hours >= 24] *= 0.5
    pattern = rng.uniform(0.0, 2.0, size=n) * (rng.uniform(size=n) < 0.8)
    rain = event[:, None] * pattern[None, :]
    et = np.where((hours % 24 >= 8) & (hours % 24 < 18), 0.3, 0.0)[:, None] * rng.uniform(0.5, 1.0, size=n)[None, :]
    case = dict(params=cols, theta=theta, rain=rain, et=et, dt_step=3600.0, max_err=1e-5, max_iter=100,
                min_step=1.0, root_levels=5)
    os.makedirs(args.out, exist_ok=True)
    np.savez(os.path.join(args.out, "case.npz"), **case)
    for k in ("params", "theta", "rain", "et"):
        np.ascontiguousarray(case[k], dtype=np.float64).tofile(os.path.join(args.out, k + ".f64"))
    with open(os.path.join(args.out, "meta.txt"), "w") as f:
        f.write(f"{n} {LEVELS} {t} {case['dt_step']!r} {case['max_err']!r} {case['max_iter']} {case['min_step']!r} "
                f"{case['root_levels']}\n")
    print(f"{n} cells, {LEVELS} levels, {t} steps; rain {rain.sum(0).mean():.1f} mm per cell in total")


if __name__ == "__main__":
    main()
