"""Numeric pattern mismatches must fail even if all finite differences are zero."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
spec = importlib.util.spec_from_file_location(
    "meteo_check", Path(__file__).resolve().parents[1] / "meteo_interpolation.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_nan_nodata_and_infinity_patterns():
    ref = np.array([[np.nan, -9999.0, np.inf, 1.0]])
    active = np.ones(ref.shape, dtype=bool)
    assert not module.compare(ref, ref.copy(), active)[1].any()
    for i, replacement in [(0, 0.0), (1, 0.0), (2, -np.inf), (3, np.nan)]:
        out = ref.copy()
        out[0, i] = replacement
        assert module.compare(ref, out, active)[1][0, i]


def test_finite_difference_and_inactive_cells():
    ref = np.array([[1.0, 2.0]])
    out = np.array([[np.nextafter(1.0, 2.0), 3.0]])
    maximum, bad = module.compare(ref, out, np.array([[True, False]]))
    assert maximum > 0 and bad[0, 0] and not bad[0, 1]
