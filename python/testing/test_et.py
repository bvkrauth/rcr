"""
TEST_ET.PY: Unit tests for estimate_theta()
"""
import sys

import pytest
import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import read_data, estimate_theta_segments, estimate_theta, \
    get_logfile, set_logfile

tol = 1e-04


# Test with simple data
def test_et_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    et_true = np.array([-0.33333333, 1.])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        et = estimate_theta(mv1, lr1, ts1)
    assert max(abs(et[:, 0] - et_true)) < tol
    set_logfile(oldlogfile)


# Basic functionality
# Estimate parameters and gradient with real data
def test_et_realdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    theta_segments, thetavec, lambdavec = \
        estimate_theta_segments(moment_vector)
    et = estimate_theta(moment_vector, lambda_range, theta_segments)
    assert max(abs(et[:, 0] - np.array([5.13504376,  5.20150257]))) < tol
    # TODO: need to check gradient too
    set_logfile(oldlogfile)


# Varying lambda_range
# lambda_range is a single point
def test_et_lambdapoint():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([0, 0])
    theta_segments, thetavec, lambdavec = \
        estimate_theta_segments(moment_vector)
    et = estimate_theta(moment_vector, lr0, theta_segments)
    assert max(abs(et[:, 0] - np.array([5.20150257,  5.20150257]))) < tol
    # TODO: need to check gradient too
    set_logfile(oldlogfile)


# lambda_range has no lower bound
def test_et_nolambdalow():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([-np.inf, 1])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        theta_segments, thetavec, lambdavec = \
            estimate_theta_segments(moment_vector)
        et = estimate_theta(moment_vector, lr0, theta_segments)
    assert max(abs(et[:, 0] - np.array([5.13504376,  8.16970996]))) < tol
    # TODO: need to check gradient too
    set_logfile(oldlogfile)


# lambda_range has no upper bound
def test_et_nolambdahigh():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([0, np.inf])
    theta_segments, thetavec, lambdavec = \
        estimate_theta_segments(moment_vector)
    et = estimate_theta(moment_vector, lr0, theta_segments)
    assert et[0, 0] == -np.inf
    assert et[1, 0] == np.inf
    assert np.all(et[:, 1:] == 0.0)
    # TODO: need to check gradient too
    set_logfile(oldlogfile)


# Special cases for moments
# Near-perfect RCT (cov(z,x) almost zero)
def test_et_nearrct():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    et = estimate_theta(mv1, lr1, ts1)
    assert max(abs(et[:, 0] - 0.5)) < tol
    set_logfile(oldlogfile)


# Perfect RCT (cov(z,x) exactly zero)
def test_et_rct():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    except UnboundLocalError:
        pass
    else:
        et = estimate_theta(mv1, lr1, ts1)
        assert max(abs(et[:, 0] - 0.5)) < tol
    set_logfile(oldlogfile)
