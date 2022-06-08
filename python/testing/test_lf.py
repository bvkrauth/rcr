"""
TEST_LF.PY: Unit tests for lambdafast()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import lambdafast, thetastar, lambdastar, simplify_moments, \
    read_data


@pytest.fixture
def moment_vector():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testing/testin1.txt")
    return moment_vector


# Basic functionality
# Takes a single value or array for theta, and a simpified moment vector
# Returns an array that is the lambda associated with each theta


# Test with simple data and a scalar theta
def test_lf_basic():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lf_true = 0.5773502691896257
    lf = lambdafast(0.0, simplify_moments(mv1))
    assert lf == pytest.approx(lf_true)


# Test with real data and an array theta
def test_lf_realdata(moment_vector):
    lf_true = np.array([28.93548917, 26.67790368])
    lf = lambdafast(np.array([0.0, 1.0, ]), simplify_moments(moment_vector))
    assert lf == pytest.approx(lf_true)


# Special cases

# The lambda(theta) function is undefined (NaN) for theta = thetastar
def test_lf_thetastar():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    theta_star = thetastar(mv1)
    lf = lambdafast(np.array([theta_star - 0.01,
                              theta_star,
                              theta_star + 0.01]), simplify_moments(mv1))
    assert np.isnan(lf[1]) and all(np.isfinite(lf[[0, 2]]))


# Large values of theta should produce something close to lambdastar
def test_lf_bigtheta(moment_vector):
    # This test fails for higher values of bignum
    bignum = 1.0e100
    lambda_star = lambdastar(moment_vector)
    lf = lambdafast(np.array([-bignum, bignum]),
                    simplify_moments(moment_vector))
    assert lf == pytest.approx(lambda_star)
