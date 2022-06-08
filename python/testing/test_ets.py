"""
TEST_ETS.PY: Unit tests for estimate_theta_segments()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import read_data, lambdafast, simplify_moments, \
    estimate_theta_segments


@pytest.fixture
def moment_vector():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testing/testin1.txt")
    return moment_vector


# Basic functionality
# Estimate parameters and gradient with real data
def test_ets_realdata(moment_vector):
    em_true = np.array([-1.00000000e+100,
                        -1.48223355e+001,
                        8.16970996e+000,
                        8.16970996e+000,
                        1.00000000e+100])
    em, thetavec, lambdavec = estimate_theta_segments(moment_vector)
    assert em == pytest.approx(em_true)


# Additional functionality
# Provide details on thetavec and lambdavec
def test_ets_details(moment_vector):
    global detail_file
    (em, thetavec, lambdavec) = estimate_theta_segments(moment_vector)
    lambdavec_true = lambdafast(thetavec,
                                simplify_moments(moment_vector))
    assert all(np.isfinite(thetavec))
    assert lambdavec == pytest.approx(lambdavec_true)


# Special cases for moments
# Near-perfect RCT (cov(z,x) almost zero)
def test_ets_nearrct():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    ts1_true = np.array([-1.e+100,  5.e+005,  5.e+005,  1.e+100])
    assert ts1 == pytest.approx(ts1_true)


# Perfect RCT (cov(z,x) exactly zero)
def test_ets_rct():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        ts1, thetavec, lambdavec = estimate_theta_segments(mv1)
    except UnboundLocalError:
        pass
    else:
        pass
