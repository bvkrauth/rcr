"""
TEST_TSF.PY: Unit tests for thetastar()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import thetastar, read_data


# Test with simple data
def test_tsf_basic():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    tsf = 1.0
    ts1 = thetastar(mv1)
    assert ts1 == pytest.approx(tsf)


# Test with real data
def test_tsf_realdata():
    (n_moments, n_theta, external_big_number, moment_vector,
        theta_range) = read_data("testin1.txt")
    tsf = 8.169709964904111
    ts1 = thetastar(moment_vector)
    assert ts1 == pytest.approx(tsf)


# Special cases


# var(zhat) = 0 - should return infinity
def test_tsf_zero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    assert np.isnan(thetastar(mv1))


# var(zhat) near 0 - should work normally
def test_tsf_nearzero():
    mv1 = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    assert np.isfinite(thetastar(mv1))
