"""
TEST_EP.PY: Unit tests for estimate_parameter()
"""


import numpy as np
import pytest

from rcrbounds import estimate_parameter, effectinf, rcinf, \
    simplify_moments


def vary(moment_vector):
    """calculate var(y) from moment vector"""
    simplified_moments = simplify_moments(moment_vector)
    return simplified_moments[0]


# Test with simple data
def test_ep_basic():
    """estimate parameter with simple data"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    ep_true = np.array([1., 0, 0, 0, 0, 0, 0, 1, 0, 0])
    test_ep = estimate_parameter(vary, moment_vector)
    assert test_ep == pytest.approx(ep_true)


# Test with real data
def test_ep_realdata(moment_vector):
    """estimate parameter with real data"""
    ep_true = np.zeros(len(moment_vector)+1)
    ep_true[0] = 542.53831290783
    ep_true[7] = -102.89924396
    ep_true[42] = 1.0
    test_ep = estimate_parameter(vary, moment_vector)
    assert test_ep == pytest.approx(ep_true)


# Special cases


def test_ep_inf():
    """when the function is inf, derivative is zeros"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    test_ep = estimate_parameter(rcinf, moment_vector)
    assert test_ep[0] == np.inf and all(test_ep[1:] == 0.)


def test_ep_nan():
    """when the function is nan, derivative is zeros"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    test_ep = estimate_parameter(effectinf, moment_vector)
    assert np.isnan(test_ep[0]) and all(test_ep[1:] == 0.)
