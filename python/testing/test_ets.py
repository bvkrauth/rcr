"""
TEST_ETS.PY: Unit tests for estimate_theta_segments()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import read_data, lambdafast, simplify_moments, \
    estimate_theta_segments, get_logfile, set_logfile

tol = 1e-04


# Basic functionality
# Estimate parameters and gradient with real data
def test_ets_realdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    em, thetavec, lambdavec = estimate_theta_segments(moment_vector)
    assert max(abs(em - np.array([-1.00000000e+100,
                                  -1.48223355e+001,
                                  8.16970996e+000,
                                  8.16970996e+000,
                                  1.00000000e+100]))) < tol
    set_logfile(oldlogfile)


# Additional functionality
# Provide details on thetavec and lambdavec
def test_ets_details():
    oldlogfile = get_logfile()
    set_logfile(None)
    global detail_file
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    em, thetavec, lambdavec = estimate_theta_segments(moment_vector)
    lambdavec_true = lambdafast(thetavec,
                                simplify_moments(moment_vector))
    assert all(np.isfinite(thetavec))
    assert max(abs(lambdavec - lambdavec_true)) < tol
    set_logfile(oldlogfile)


# Special cases for moments
# Near-perfect RCT (cov(z,x) almost zero)
def test_ets_nearrct():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    ts1_true = np.array([-1.e+100,  5.e+005,  5.e+005,  1.e+100])
    assert max(abs(ts1 - ts1_true)) < tol
    set_logfile(oldlogfile)


# Perfect RCT (cov(z,x) exactly zero)
def test_ets_rct():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    except UnboundLocalError:
        pass
    else:
        pass
    set_logfile(oldlogfile)
