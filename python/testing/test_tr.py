"""
TEST_TR.PY: Unit tests for translate_result()
"""
import numpy as np
import pytest

from rcrbounds import translate_result


# Basic functionality
def test_tr_basic():
    """translate_result with default arguments"""
    mat = np.array([1., -1., -np.inf, np.inf, np.nan])
    expected_result = mat
    result = translate_result(mat)
    assert result == pytest.approx(expected_result, nan_ok=True)


# Optional arguments
def test_tr_withvalues():
    """translate_result with optional arguments"""
    mat = np.array([1., -1., -np.inf, np.inf, np.nan])
    expected_result = np.array([1., -1., -100., 100., 50.])
    result = translate_result(mat, inf=100.0, nan=50.0)
    assert result == pytest.approx(expected_result)
