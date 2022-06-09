"""
TEST_EP.PY: Unit tests for estimate_parameter()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import estimate_parameter, thetastar, lambdastar, \
    simplify_moments  # pylint: disable=wrong-import-position


def vary(moment_vector):
    """calculate var(y) from moment vector"""
    sm = simplify_moments(moment_vector)
    return sm[0]


# Test with simple data
def test_ep_basic():
    """estimate parameter with simple data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    ep_true = np.array([1., 0, 0, 0, 0, 0, 0, 1, 0, 0])
    ep = estimate_parameter(vary, mv1)
    assert ep == pytest.approx(ep_true)


# Test with real data
def test_ep_realdata(moment_vector):
    """estimate parameter with real data"""
    ep_true = np.zeros(len(moment_vector)+1)
    ep_true[0] = 542.53831290783
    ep_true[7] = -102.89924396
    ep_true[42] = 1.0
    ep = estimate_parameter(vary, moment_vector)
    assert ep == pytest.approx(ep_true)


# Special cases


def test_ep_inf():
    """when the function is inf, derivative is zeros"""
    mv = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    ep = estimate_parameter(lambdastar, mv)
    assert ep[0] == np.inf and all(ep[1:] == 0.)


def test_ep_nan():
    """when the function is nan, derivative is zeros"""
    mv = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    ep = estimate_parameter(thetastar, mv)
    assert np.isnan(ep[0]) and all(ep[1:] == 0.)
