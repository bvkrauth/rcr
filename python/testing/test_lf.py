"""
TEST_LF.PY: Unit tests for lambdafast()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import lambdafast, thetastar, lambdastar, \
    simplify_moments  # pylint: disable=wrong-import-position


# Basic functionality
def test_lf_basic():
    """lambdafast for simple data and scalar theta"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lf_true = 0.5773502691896257
    lf = lambdafast(0.0, simplify_moments(mv1))
    assert lf == pytest.approx(lf_true)


def test_lf_realdata(moment_vector):
    """lambdafast for real data and array theta"""
    lf_true = np.array([28.93548917, 26.67790368])
    lf = lambdafast(np.array([0.0, 1.0, ]), simplify_moments(moment_vector))
    assert lf == pytest.approx(lf_true)


# Special cases
def test_lf_thetastar():
    """function is undefined (NaN) for theta = thetastar"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    theta_star = thetastar(mv1)
    lf = lambdafast(np.array([theta_star - 0.01,
                              theta_star,
                              theta_star + 0.01]), simplify_moments(mv1))
    assert np.isnan(lf[1]) and all(np.isfinite(lf[[0, 2]]))


def test_lf_bigtheta(moment_vector):
    """function is close to lambdastar for large theta"""
    # This test fails for higher values of bignum
    bignum = 1.0e100
    lambda_star = lambdastar(moment_vector)
    lf = lambdafast(np.array([-bignum, bignum]),
                    simplify_moments(moment_vector))
    assert lf == pytest.approx(lambda_star)
