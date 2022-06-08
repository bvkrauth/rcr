"""
TEST_LSF.PY: Unit tests for lambdastar()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import lambdastar, read_data


@pytest.fixture
def moment_vector():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testing/testin1.txt")
    return moment_vector


# Test with simple data
def test_true_ls_basic():
    mv1 = np.array([0, 0, 0, 1, 0.5, np.sqrt(0.2), 1, 0.5, 1.0])
    true_ls = 2.0
    ls = lambdastar(mv1)
    assert ls == pytest.approx(true_ls)


# Test with real data
def test_true_ls_realdata(moment_vector):
    true_ls = 12.310599093115798
    ls = lambdastar(moment_vector)
    assert ls == pytest.approx(true_ls)

# Special cases


# var(zhat) = 0 - should return infinity
def test_true_ls_zero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    true_ls = np.inf
    ls = lambdastar(mv1)
    assert ls == pytest.approx(true_ls)


# var(zhat) near 0 - should work normally
def test_true_ls_nearzero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    ls = lambdastar(mv1)
    assert np.isfinite(ls)


# var(zhat) > var(z) - should return zero
def test_true_ls_negvarz():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    ls = lambdastar(mv1)
    assert ls == 0.0
