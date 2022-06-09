"""
TEST_TSF.PY: Unit tests for thetastar()
"""


import numpy as np
import pytest

from rcr import thetastar


# Basic functionality
def test_true_ts_basic():
    """thetastar from simple data"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    true_ts = 1.0
    theta_star = thetastar(moment_vector)
    assert theta_star == pytest.approx(true_ts)


def test_true_ts_realdata(moment_vector):
    """thetastar from real data"""
    true_ts = 8.169709964904111
    theta_star = thetastar(moment_vector)
    assert theta_star == pytest.approx(true_ts)


# Special cases
def test_true_ts_zero():
    """return inf if var(zhat) = 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    theta_star = thetastar(moment_vector)
    assert np.isnan(theta_star)


def test_true_ts_nearzero():
    """return a finite value if var(zhat) near 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    theta_star = thetastar(moment_vector)
    assert np.isfinite(theta_star)
