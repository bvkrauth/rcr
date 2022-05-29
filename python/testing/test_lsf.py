"""
TEST_LSF.PY: Unit tests for lambdastar()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import lambdastar, read_data

tol = 1e-04


# Test with simple data
def test_lsf_basic():
    mv1 = np.array([0, 0, 0, 1, 0.5, np.sqrt(0.2), 1, 0.5, 1.0])
    lsf = 2.0
    assert abs(lambdastar(mv1) - lsf) < tol


# Test with real data
def test_lsf_realdata():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lsf = 12.310599093115798
    assert abs(lambdastar(moment_vector) - lsf) < tol

# Special cases


# var(zhat) = 0 - should return infinity
def test_lsf_zero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    assert lambdastar(mv1) == np.inf


# var(zhat) near 0 - should work normally
def test_lsf_nearzero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    assert np.isfinite(lambdastar(mv1))


# var(zhat) > var(z) - should return zero
def test_lsf_negvarz():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    assert lambdastar(mv1) == 0.0
