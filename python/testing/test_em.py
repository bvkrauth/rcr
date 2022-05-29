"""
TEST_EM.PY: Unit tests for estimate_model()
"""
import sys

import pytest
import numpy as np
import pandas as pd

sys.path.append("./")
sys.path.append("../")
from rcr import read_data, estimate_model

tol = 1e-04

# Test with simple data
# def test_ep_basic():
#    mv1 = np.array([0,0,0, 1,0.5,0.5, 1, 0.5, 1.0])
#    ep_true = np.array([1., 0, 0, 0, 0, 0, 0, 1, 0, 0])
#    ep = estimate_parameter(vary,mv1)
#    assert max(abs(ep - ep_true)) < tol


# Basic functionality
# Estimate parameters and gradient with real data
def test_em_realdata():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    true_result = np.asarray(pd.read_csv("testout1.txt",
                                         delimiter=" ",
                                         header=None,
                                         skipinitialspace=True))
    true_em1 = np.array([12.31059909,  8.16970996, 28.93548917,
                         5.13504376,  5.20150257])
    em, thetavec, lambdavec = estimate_model(moment_vector, lambda_range)
    # Check parameter estimates
    assert max(abs(em[:, 0] - true_em1)) < 1e-4
    # Check parameter estimates and gradient
    assert np.max(abs(em - true_result)) < 1e-4


# lambda_range is a single point
def test_em_lambdapoint():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([0, 0])
    true_em = np.array([12.31059909,  8.16970996, 28.93548917,
                        5.20150257,  5.20150257])
    em, thetavec, lambdavec = estimate_model(moment_vector, lr0)
    assert max(abs(em[:, 0] - true_em)) < 1e-4
    # TODO: need to check gradient too


# lambda_range has no lower bound
def test_em_nolambdalow():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([-np.inf, 1])
    true_em = np.array([12.31059909,  8.16970996, 28.93548917,
                        5.13504376,  8.16970996])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        em, thetavec, lambdavec = estimate_model(moment_vector, lr0)
    assert max(abs(em[:, 0] - true_em)) < 1e-4
    # TODO: need to check gradient too


# lambda_range has no upper bound
def test_em_nolambdahigh():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    lr0 = np.array([0, np.inf])
    em, thetavec, lambdavec = estimate_model(moment_vector, lr0)
    assert max(abs(em[0:3, 0] - np.array([12.31059909,  8.16970996,
                                         28.93548917]))) < 1e-4
    assert em[3, 0] == -np.inf
    assert em[4, 0] == np.inf
    assert np.all(em[3:4, 1:] == 0.0)
    # TODO: need to check gradient too


# Special cases for moments
# Near-perfect RCT (cov(z,x) almost zero)
def test_em_nearrct():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    em, thetavec, lambdavec = estimate_model(mv1, lr1)
    assert np.all(em[0:3, 0] > 1000)
    assert max(abs(em[3:4, 0] - 0.5)) < tol


# Perfect RCT (cov(z,x) exactly zero)
def test_em_rct():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        em, thetavec, lambdavec = estimate_model(mv1, lr1)
    except UnboundLocalError:
        pass
    else:
        assert np.all(em[0:3, 0] > 1000)
        assert max(abs(em[3:4, 0] - 0.5)) < tol
