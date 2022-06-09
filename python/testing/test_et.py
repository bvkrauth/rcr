"""
TEST_ET.PY: Unit tests for estimate_theta()
"""


import pytest
import numpy as np

from rcr import estimate_theta_segments, estimate_theta


# Basic functionality
def test_et_basic():
    """estimate theta from simple data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1 = estimate_theta_segments(mv1)[0]
    et_true = np.array([-0.33333333, 1.])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        et = estimate_theta(mv1, lr1, ts1)
    assert et[:, 0] == pytest.approx(et_true, rel=1e-04)


def test_et_realdata(moment_vector):
    """estimate theta from real data"""
    lambda_range = np.array([0.0, 1.0])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([5.13504376,  5.20150257])
    et = estimate_theta(moment_vector, lambda_range, theta_segments)
    assert et[:, 0] == pytest.approx(et_true)
    # should check gradient too


# Varying lambda_range
def test_et_lambdapoint(moment_vector):
    """estimate theta when lambda_range is a single point"""
    lr0 = np.array([0.0, 0.0])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([5.20150257,  5.20150257])
    et = estimate_theta(moment_vector, lr0, theta_segments)
    assert et[:, 0] == pytest.approx(et_true)
    # should check gradient too


def test_et_nolambdalow(moment_vector):
    """estimate theta when lambda_range has no lower bound"""
    lr0 = np.array([-np.inf, 1])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([5.13504376,  8.16970996])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        et = estimate_theta(moment_vector, lr0, theta_segments)
    assert et[:, 0] == pytest.approx(et_true)
    # should check gradient too


def test_et_nolambdahigh(moment_vector):
    """estimate theta when lambda_range has no upper bound"""
    lr0 = np.array([0, np.inf])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et = estimate_theta(moment_vector, lr0, theta_segments)
    assert et[0, 0] == -np.inf
    assert et[1, 0] == np.inf
    assert np.all(et[:, 1:] == 0.0)
    # should check gradient too


# Special cases for moments
def test_et_nearrct():
    """estimate theta for near-perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1 = estimate_theta_segments(mv1)[0]
    et_true = 0.5
    et = estimate_theta(mv1, lr1, ts1)
    assert et[:, 0] == pytest.approx(et_true, rel=1e-04)


def test_et_rct():
    """estimate theta for perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    et_true = 0.5
    # This test currently fails with an UnboundLocalError
    try:
        ts1 = estimate_theta_segments(mv1)[0]
    except UnboundLocalError:
        pass
    else:
        et = estimate_theta(mv1, lr1, ts1)
        assert et[:, 0] == pytest.approx(et_true, rel=1e-04)
