"""Value and exceptional-pattern mismatches must fail groundwater validation."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location('groundwater_comparison', Path(__file__).parents[1] / 'groundwater.py')
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


@pytest.mark.parametrize('other', [1.1, np.nan, np.inf, -np.inf, -9999.0])
def test_mismatch_fails(other):
    expected = np.array([[1.0, np.nan, np.inf, -np.inf, -9999.0]])
    actual = expected.copy()
    actual[0, 0] = other
    result = comparison.compare_array(actual, expected, np.ones(expected.shape, dtype=bool))
    assert not result['ok']
    assert result['first']['column'] == 1


def test_matching_patterns_and_inactive_cells():
    expected = np.array([[np.nan, np.inf, -np.inf, -9999.0, 1.0]])
    actual = expected.copy()
    actual[0, -1] = np.nan
    active = np.array([[True, True, True, True, False]])
    assert comparison.compare_array(actual, expected, active)['ok']
    assert not comparison.compare_array(actual[:, :-1], expected, active)['ok']
