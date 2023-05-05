"""
TEST_RF.RC Unit tests for fit() method of RCR object
"""


import pytest
import numpy as np

from rcrbounds import RCR, RCRResults


# Basic functionality
def test_rf_basic(model):
    """RCR.fit with default options"""
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
    assert results.param_names == ['rcInf',
                                   'effectInf',
                                   'rc0',
                                   'effectL',
                                   'effectH']


# Set rc_range
def test_rf_lr(model, endog, exog):
    """RCR.fit with finite rc_range set"""
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             5.20150257,
                             5.20150257])
    results = model.fit(rc_range=np.asarray([0.0, 0.0]))
    assert results.params == pytest.approx(trueparams)
    model = RCR(endog, exog, rc_range=np.asarray([0.0, 0.0]))
    results = model.fit()
    assert results.params == pytest.approx(trueparams)
    results = model.fit(rc_range=np.asarray([0.0, 0.0]))
    assert results.params == pytest.approx(trueparams)


def test_rf_lrnolb(model):
    """RCR.fit with rc_range with no lower bound"""
    trueparams = np.asarray([12.31059909,
                             8.16970997,
                             28.93548917,
                             5.13504376,
                             8.16970997])
    # This will produce a warning
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        results = model.fit(rc_range=np.asarray([-np.inf, 1]))
    assert results.params == pytest.approx(trueparams)


def test_rf_lrnoub(model):
    """RCR.fit with rc_range with no upper bound"""
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
    results = model.fit(rc_range=np.asarray([0, np.inf]))
    assert results.params == pytest.approx(trueparams)
    assert results.cov_params == pytest.approx(truecov)


# Exceptions for rc_range
def test_rf_lr2d(model):
    """raise exception if rc_range is a 2-d array"""
    try:
        model.fit(rc_range=np.zeros((2, 2)))
    except TypeError:
        pass
    else:
        raise AssertionError


# rc_range has wrong number of elements
def test_rf_lr1e(model):
    """raise exception if rc_range has wrong # of elements"""
    try:
        model.fit(rc_range=np.zeros(1))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rf_lrnan(model):
    """raise exception if rc_range has NaN values"""
    try:
        model.fit(rc_range=np.asarray([0, np.nan]))
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rf_lrnotsorted(model):
    """raise exception if rc_range is out of order"""
    try:
        model.fit(rc_range=np.asarray([1., 0.]))
    except ValueError:
        pass
    else:
        raise AssertionError


# Covariance matrix options
def test_rf_vceadj(model):
    """use vceadj to modify covariance matrix"""
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


def test_rf_ctunsup(model):
    """raise exception if cov_type is unsupported"""
    try:
        model.fit(cov_type="unsupported")
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rf_vcestr(model):
    """raise exception if vceadj is non-numeric"""
    try:
        model.fit(vceadj="a string")
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rf_vcearr(model):
    """raise exception if vceadj is non-scalar"""
    try:
        model.fit(vceadj=np.asarray([1.0, 2.0]))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rf_vceneg(model):
    """raise exception if vceadj is negative"""
    try:
        model.fit(vceadj=-1.)
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rf_weighted(endog, exog, weights):
    """estimate with weights"""
    msk = weights > 0.5
    model0 = RCR(endog[msk], exog[msk])
    model1 = RCR(endog, exog, weights=weights)
    model2 = RCR(endog, exog)
    res0 = model0.fit()
    res1 = model1.fit()
    res2 = model2.fit(weights=weights)
    assert res1.params == pytest.approx(res0.params)
    assert res1.cov_params == pytest.approx(res0.cov_params)
    assert res1.model.nobs == res0.model.nobs
    assert res2.params == pytest.approx(res0.params)
    assert res2.cov_params == pytest.approx(res0.cov_params)
    assert res2.model.nobs == res0.model.nobs


def test_rf_cluster(endog, exog, clusters):
    """estimate with cluster-robust standard errors"""
    model = RCR(endog,
                exog,
                cov_type="cluster",
                groupvar=clusters)
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


def test_rf_nocluster(endog, exog):
    """estimate with cluster-robust standard errors but no groupvar"""
    model = RCR(endog,
                exog,
                cov_type="cluster")
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
    result = model.fit()
    assert result.cov_params == pytest.approx(truecov)


def test_rf_clust_and_wt(endog, exog, clusters, weights):
    """estimate with weights and clusters"""
    model = RCR(endog,
                exog,
                cov_type="cluster",
                groupvar=clusters,
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
