"""
TEST_TSF.PY: Unit tests for thetastar()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import thetastar, read_data  # pylint: disable=wrong-import-position


@pytest.fixture
def moment_vector():
    """get moment vector from test data file"""
    mv = read_data("testing/testin1.txt")[3]
    return mv


# Basic functionality
def test_true_ts_basic():
    """thetastar from simple data"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    true_ts = 1.0
    ts = thetastar(moment_vector)
    assert ts == pytest.approx(true_ts)


def test_true_ts_realdata(moment_vector):
    """thetastar from real data"""
    true_ts = 8.169709964904111
    ts = thetastar(moment_vector)
    assert ts == pytest.approx(true_ts)


# Special cases
def test_true_ts_zero():
    """return inf if var(zhat) = 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    ts = thetastar(moment_vector)
    assert np.isnan(ts)


def test_true_ts_nearzero():
    """return a finite value if var(zhat) near 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    ts = thetastar(moment_vector)
    assert np.isfinite(ts)
