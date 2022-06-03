"""
TEST_RF.RC Unit tests for fit() method of RCR object
"""
import sys

import pytest
import numpy as np
import pandas as pd
import patsy

sys.path.append("./")
sys.path.append("../")
from rcr import RCR, RCR_results


@pytest.fixture
def dat():
    fname = "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
    return pd.read_stata(fname)


@pytest.fixture
def rcr_formula():
    rcr_left = "SAT + Small_Class ~ "
    rcr_right1 = "White_Asian + Girl + Free_Lunch + White_Teacher + "
    rcr_right2 = "Teacher_Experience + Masters_Degree"
    return rcr_left + rcr_right1 + rcr_right2


@pytest.fixture
def endog(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat)
    return endog


@pytest.fixture
def exog(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat)
    return exog


@pytest.fixture
def model(endog, exog):
    return RCR(endog, exog)


# Basic functionality
def test_rf_basic(model):
    results = model.fit()
    assert isinstance(results, RCR_results)
    assert isinstance(results.model, RCR)
    assert results.model == model
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             5.13504376,
                             5.20150257])
    truecov = np.asarray([[4.40273105e+00,  1.68091057e+00,  1.48603397e+01,
                           2.62163549e-02,  1.48105699e-02],
                          [1.68091057e+00,  9.36816074e+02, -3.30554494e+03,
                           -2.08604784e+01,  9.45995702e-02],
                          [1.48603397e+01, -3.30554494e+03,  1.17764763e+04,
                           7.63213528e+01,  2.09329548e+00],
                          [2.62163549e-02, -2.08604784e+01,  7.63213528e+01,
                           9.15729396e-01,  4.38565221e-01],
                          [1.48105699e-02,  9.45995702e-02,  2.09329548e+00,
                           4.38565221e-01,  4.30902711e-01]])
    assert results.params == pytest.approx(trueparams)
    # We will check values by checking other calculations
    assert results.cov_params.shape == (5, 5)
    assert results.cov_params == pytest.approx(truecov)
    # We will not check values here, though maybe we should
    assert results.details.shape == (2, 30000)
    assert results.param_names == ['lambdaInf',
                                   'betaxInf',
                                   'lambda0',
                                   'betaxL',
                                   'betaxH']
    assert results.cov_type == "nonrobust"
    assert results.vceadj == 1.0
    assert results.citype == "conservative"
    assert results.cilevel == 95


# Set lambda_range
def test_rf_lr(model, endog, exog):
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             5.20150257,
                             5.20150257])
    results = model.fit(lambda_range=np.asarray([0.0, 0.0]))
    assert results.params == pytest.approx(trueparams)
    model = RCR(endog, exog, lambda_range=np.asarray([0.0, 0.0]))
    results = model.fit()
    assert results.params == pytest.approx(trueparams)
    results = model.fit(lambda_range=np.asarray([0.0, 0.0]))
    assert results.params == pytest.approx(trueparams)


# lambda_range with no lower bound
def test_rf_lrnolb(model):
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             5.13504376,
                             8.16970997])
    # This will produce a warning
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        results = model.fit(lambda_range=np.asarray([-np.inf, 1]))
    assert results.params == pytest.approx(trueparams)


# lambda_range with no upper bound
def test_rf_lrnoub(model):
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             -np.inf,
                             np.inf])
    # for infinite values covariance should be NaN.  For
    # Stata compatibility it is zurrently zero.
    truecov = np.asarray([[4.40273105e+00,  1.68091057e+00,  1.48603397e+01,
                           0.00000000e+00,  0.00000000e+00],
                          [1.68091057e+00,  9.36816074e+02, -3.30554494e+03,
                           0.00000000e+00,  0.00000000e+00],
                          [1.48603397e+01, -3.30554494e+03,  1.17764763e+04,
                           0.00000000e+00,  0.00000000e+00],
                          [0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                           0.00000000e+00,  0.00000000e+00],
                          [0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                           0.00000000e+00,  0.00000000e+00]])
    results = model.fit(lambda_range=np.asarray([0, np.inf]))
    assert results.params == pytest.approx(trueparams)
    assert results.cov_params == pytest.approx(truecov)


# lambda_range is a 2-d array (should be 1-d)
def test_rf_lr2d(model):
    try:
        model.fit(lambda_range=np.zeros((2, 2)))
    except TypeError:
        pass
    else:
        raise AssertionError


# lambda_range has wrong number of elements
def test_rf_lr1e(model):
    try:
        model.fit(lambda_range=np.zeros(1))
    except TypeError:
        pass
    else:
        raise AssertionError


# lambda_range has NaNs
def test_rf_lrnan(model):
    try:
        model.fit(lambda_range=np.asarray([0, np.nan]))
    except ValueError:
        pass
    else:
        raise AssertionError


# lambda_range is not in ascending order
def test_rf_lrnotsorted(model):
    try:
        model.fit(lambda_range=np.asarray([1., 0.]))
    except ValueError:
        pass
    else:
        raise AssertionError


# Covariance matrix options
# cov_type optional argument
# Only default value ("nonrobust") is currently supported

# vceadj optional argument
def test_rf_vceadj(model):
    truecov = np.asarray([[2.20136553e+00,  8.40455283e-01,  7.43016985e+00,
                           1.31081775e-02,  7.40528497e-03],
                          [8.40455283e-01,  4.68408037e+02, -1.65277247e+03,
                           -1.04302392e+01,  4.72997851e-02],
                          [7.43016985e+00, -1.65277247e+03,  5.88823814e+03,
                           3.81606764e+01,  1.04664774e+00],
                          [1.31081775e-02, -1.04302392e+01,  3.81606764e+01,
                           4.57864698e-01,  2.19282610e-01],
                          [7.40528497e-03,  4.72997851e-02,  1.04664774e+00,
                           2.19282610e-01,  2.15451356e-01]])
    results = model.fit(cov_type="nonrobust", vceadj=0.5)
    assert results.cov_params == pytest.approx(truecov)


# cov_type is unsupported
def test_rf_ctunsup(model):
    try:
        model.fit(cov_type="unsupported")
    except ValueError:
        pass
    else:
        raise AssertionError


# vceadj is non-numeric
def test_rf_vcestr(model):
    try:
        model.fit(vceadj="a string")
    except TypeError:
        pass
    else:
        raise AssertionError


# vceadj is not a scalar
def test_rf_vcearr(model):
    try:
        model.fit(vceadj=np.asarray([1.0, 2.0]))
    except TypeError:
        pass
    else:
        raise AssertionError


# vceadj is negative
def test_rf_vceneg(model):
    try:
        model.fit(vceadj=-1.)
    except ValueError:
        pass
    else:
        raise AssertionError


# weights
def test_rf_weighted(endog, exog):
    wt = np.array(([1., 0.] * 2920)[0:5839])
    msk = (wt > 0.5)
    model0 = RCR(endog[msk], exog[msk])
    model1 = RCR(endog, exog, weights=wt)
    model2 = RCR(endog, exog)
    res0 = model0.fit()
    lf0 = res0._lambdafun(thetavals=np.zeros(1),
                          include_thetastar=False)
    res1 = model1.fit()
    res2 = model2.fit(weights=wt)
    assert res1.params == pytest.approx(res0.params)
    assert res1.cov_params == pytest.approx(res0.cov_params)
    lf1 = res1._lambdafun(thetavals=np.zeros(1),
                          include_thetastar=False)
    assert lf1 == pytest.approx(lf0)
    assert res2.params == pytest.approx(res0.params)
    assert res2.cov_params == pytest.approx(res0.cov_params)
    lf2 = res2._lambdafun(thetavals=np.zeros(1),
                          include_thetastar=False)
    assert lf2 == pytest.approx(lf0)
