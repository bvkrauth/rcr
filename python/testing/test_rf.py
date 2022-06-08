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
from rcr import RCR, RCRResults


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


@pytest.fixture
def weights(dat):
    wt = np.mod(dat["TCHID"], 2)
    return wt


@pytest.fixture
def cluster(dat):
    clust = dat["TCHID"]
    return clust


# Basic functionality
def test_rf_basic(model):
    results = model.fit()
    assert isinstance(results, RCRResults)
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
def test_rf_weighted(endog, exog, weights):
    msk = (weights > 0.5)
    model0 = RCR(endog[msk], exog[msk])
    model1 = RCR(endog, exog, weights=weights)
    model2 = RCR(endog, exog)
    res0 = model0.fit()
    res1 = model1.fit()
    res2 = model2.fit(weights=weights)
    assert res1.params == pytest.approx(res0.params)
    assert res1.cov_params == pytest.approx(res0.cov_params)
    assert res1.model.nobs == res0.model.nobs
    assert res1.nobs == res0.nobs
    assert res2.params == pytest.approx(res0.params)
    assert res2.cov_params == pytest.approx(res0.cov_params)
    assert res2.model.nobs == len(endog)
    assert res2.nobs == res0.nobs


# cluster-robust standard error
def test_rf_cluster(endog, exog, cluster):
    model = RCR(endog,
                exog,
                cov_type="cluster",
                groupvar=cluster)
    truecov = np.array([[6.91681472e+01,  1.26630806e+02, -1.21548954e+02,
                         -1.94406588e+00,  1.22940162e-01],
                        [1.26630806e+02,  1.90495752e+03, -6.12590889e+03,
                         -3.75141027e+01,  3.66381486e+00],
                        [-1.21548954e+02, -6.12590889e+03,  2.10764578e+04,
                         1.29096257e+02, -6.60769826e+00],
                        [-1.94406588e+00, -3.75141027e+01,  1.29096257e+02,
                         1.84604809e+00,  1.00506874e+00],
                        [1.22940162e-01,  3.66381486e+00, -6.60769826e+00,
                         1.00506874e+00,  1.06212946e+00]])
    result = model.fit()
    assert result.model.ngroups == 323
    assert result.cov_params == pytest.approx(truecov)


# clusters and weights
def test_rf_clust_and_wt(endog, exog, cluster, weights):
    model = RCR(endog,
                exog,
                cov_type="cluster",
                groupvar=cluster,
                weights=weights)
    truecov = np.array([[6.50239759e+02,  1.07844382e+02, -9.32223483e+02,
                         -1.31757653e-01, -2.38992069e+01],
                        [1.07844382e+02,  9.97990718e+03, -4.14671277e+03,
                         2.72173698e+01, -6.49984265e+01],
                        [-9.32223483e+02, -4.14671277e+03,  2.95964719e+03,
                         -1.76052641e+01,  5.18035063e+01],
                        [-1.31757653e-01,  2.72173698e+01, -1.76052641e+01,
                         2.18090187e+00,  1.96384166e+00],
                        [-2.38992069e+01, -6.49984265e+01,  5.18035063e+01,
                         1.96384166e+00,  3.39851221e+00]])
    result = model.fit()
    assert result.model.nobs == 3325
    assert result.model.ngroups == 184
    assert result.cov_params == pytest.approx(truecov)
