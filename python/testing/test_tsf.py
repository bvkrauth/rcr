"""
TEST_TSF.PY: Unit tests for thetastar()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import thetastar, read_data


@pytest.fixture
def moment_vector():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testing/testin1.txt")
    return moment_vector


# Test with simple data
def test_true_ts_basic():
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    true_ts = 1.0
    ts = thetastar(moment_vector)
    assert ts == pytest.approx(true_ts)


# Test with real data
def test_true_ts_realdata(moment_vector):
    true_ts = 8.169709964904111
    ts = thetastar(moment_vector)
    assert ts == pytest.approx(true_ts)


# Special cases


# var(zhat) = 0 - should return infinity
def test_true_ts_zero():
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    ts = thetastar(moment_vector)
    assert np.isnan(ts)


# var(zhat) near 0 - should work normally
def test_true_ts_nearzero():
    moment_vector = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    ts = thetastar(moment_vector)
    assert np.isfinite(ts)
