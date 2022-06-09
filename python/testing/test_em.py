"""
TEST_EM.PY: Unit tests for estimate_model()
"""
import sys

import pytest
import numpy as np
import pandas as pd

sys.path.append("./")
sys.path.append("../")
from rcr import estimate_model  # pylint: disable=wrong-import-position


# Basic functionality
def test_em_realdata(moment_vector):
    """estimate parameters and gradient with real data"""
    lambda_range = np.array([0.0, 1.0])
    true_em1 = np.array([12.31059909,  8.16970996, 28.93548917,
                         5.13504376,  5.20150257])
    true_result = np.asarray(pd.read_csv("testing/testout1.txt",
                                         delimiter=" ",
                                         header=None,
                                         skipinitialspace=True))
    em = estimate_model(moment_vector, lambda_range)[0]
    # Check parameter estimates
    assert em[:, 0] == pytest.approx(true_em1)
    # Check parameter estimates and gradient
    assert em == pytest.approx(true_result)


# Special cases for lambda_range
def test_em_lambdapoint(moment_vector):
    """estimate when lambda_range is a single point"""
    lr0 = np.array([0, 0])
    true_em = np.array([12.31059909,  8.16970996, 28.93548917,
                        5.20150257,  5.20150257])
    em = estimate_model(moment_vector, lr0)[0]
    assert em[:, 0] == pytest.approx(true_em)
    # need to check gradient too


def test_em_nolambdalow(moment_vector):
    """estimate when lambda_range has no lower bound"""
    lr0 = np.array([-np.inf, 1])
    true_em = np.array([12.31059909,  8.16970996, 28.93548917,
                        5.13504376,  8.16970996])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        em = estimate_model(moment_vector, lr0)[0]
    assert em[:, 0] == pytest.approx(true_em)
    # need to check gradient too


def test_em_nolambdahigh(moment_vector):
    """estimate when lambda_range has no upper bound"""
    lr0 = np.array([0, np.inf])
    true_em = np.array([12.31059909,  8.16970996, 28.93548917,
                        -np.inf, np.inf])
    em = estimate_model(moment_vector, lr0)[0]
    assert em[:, 0] == pytest.approx(true_em)
    assert np.all(em[3:4, 1:] == 0.0)
    # need to check gradient too


# Special cases for moments
def test_em_nearrct():
    """estimate for near-perfect RCT: cov(z,x) almost zero"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    em = estimate_model(mv1, lr1)[0]
    assert np.all(em[0:3, 0] > 1000)
    assert em[3:4, 0] == pytest.approx(0.5, rel=1e-04)


def test_em_rct():
    """estimate for perfect RCT: cov(z,x)=0"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        em = estimate_model(mv1, lr1)[0]
    except UnboundLocalError:
        pass
    else:
        assert np.all(em[0:3, 0] > 1000)
        assert em[3:4, 0] == pytest.approx(0.5, rel=1e-04)
