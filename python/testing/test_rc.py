"""
TEST_RC.RC Unit tests for RCR(), not including the fit method
"""


import numpy as np
import pandas as pd

from rcrbounds import RCR


# Basic functionality
def test_rc_patsy(endog, exog):
    """construct RCR model object using Patsy design matrices"""
    model = RCR(endog, exog)
    assert isinstance(model, RCR)
    assert isinstance(model.endog, np.ndarray)
    assert model.endog.shape == (5839, 2)
    assert isinstance(model.exog, np.ndarray)
    assert model.exog.shape == (5839, 7)
    assert model.endog_names == ["SAT", "Small_Class"]
    assert model.exog_names == ['Intercept',
                                'White_Asian',
                                'Girl',
                                'Free_Lunch',
                                'White_Teacher',
                                'Teacher_Experience',
                                'Masters_Degree']
    assert model.depvar == "SAT"
    assert model.treatvar == "Small_Class"
    control_vars = ("White_Asian Girl Free_Lunch White_Teacher " +
                    "Teacher_Experience Masters_Degree")
    assert model.controlvars == control_vars
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0
    assert model.cov_type == "nonrobust"
    assert model.vceadj == 1.0
    assert model.citype == "conservative"
    assert model.cilevel == 95
    assert model.weights is None


def test_rc_patsy_df(endog_df, exog_df):
    """construct RCR model object using Patsy data frames"""
    model = RCR(endog_df, exog_df)
    assert isinstance(model.endog, np.ndarray)
    assert model.endog.shape == (5839, 2)
    assert isinstance(model.exog, np.ndarray)
    assert model.exog.shape == (5839, 7)
    assert model.endog_names == ["SAT", "Small_Class"]
    assert model.exog_names == ['Intercept',
                                'White_Asian',
                                'Girl',
                                'Free_Lunch',
                                'White_Teacher',
                                'Teacher_Experience',
                                'Masters_Degree']
    assert model.depvar == "SAT"
    assert model.treatvar == "Small_Class"
    control_vars = ("White_Asian Girl Free_Lunch White_Teacher " +
                    "Teacher_Experience Masters_Degree")
    assert model.controlvars == control_vars
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0


def test_rc_dataframe(endog_df, exog_df):
    """construct RCR model object using data frames without design info"""
    model = RCR(pd.DataFrame(endog_df), pd.DataFrame(exog_df))
    assert isinstance(model.endog, np.ndarray)
    assert model.endog.shape == (5839, 2)
    assert isinstance(model.exog, np.ndarray)
    assert model.exog.shape == (5839, 7)
    assert model.endog_names == ["SAT", "Small_Class"]
    assert model.exog_names == ['Intercept',
                                'White_Asian',
                                'Girl',
                                'Free_Lunch',
                                'White_Teacher',
                                'Teacher_Experience',
                                'Masters_Degree']
    assert model.depvar == "SAT"
    assert model.treatvar == "Small_Class"
    cvars = ("White_Asian Girl Free_Lunch White_Teacher " +
             "Teacher_Experience Masters_Degree")
    assert model.controlvars == cvars
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0


def test_rc_array(endog, exog):
    """construct RCR model object using plain numpy arrays"""
    model = RCR(np.asarray(endog), np.asarray(exog))
    assert isinstance(model.endog, np.ndarray)
    assert model.endog.shape == (5839, 2)
    assert isinstance(model.exog, np.ndarray)
    assert model.exog.shape == (5839, 7)
    assert model.endog_names == ['y', 'treatment']
    assert model.exog_names == ['Intercept',
                                'x1',
                                'x2',
                                'x3',
                                'x4',
                                'x5',
                                'x6']
    assert model.depvar == "y"
    assert model.treatvar == "treatment"
    assert model.controlvars == 'x1 x2 x3 x4 x5 x6'
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0


# pylint: disable=duplicate-code
# lambda_range arguments
def test_rc_setlr(endog, exog):
    """set custom (finite) lambda range for RCR model object"""
    model = RCR(endog, exog, lambda_range=np.asarray([-1.0, 5.0]))
    assert model.lambda_range[0] == -1.0
    assert model.lambda_range[1] == 5.0


def test_rc_setlrinf(endog, exog):
    """set custom (non-finite) lambda range for RCR model object"""
    model = RCR(endog, exog, lambda_range=np.asarray([-np.inf, np.inf]))
    assert np.isneginf(model.lambda_range[0])
    assert np.isposinf(model.lambda_range[1])


# Exceptions
def test_rc_endog1d(endog_df, exog_df):
    """raise exception if endog is a 1-d array"""
    try:
        RCR(np.asarray(endog_df)[:, 1], np.asarray(exog_df))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_endog1c(endog_df, exog_df):
    """raise exception if endog has only 1 column"""
    try:
        RCR(np.asarray(endog_df)[:, 1:], np.asarray(exog_df))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_endog7c(exog_df):
    """raise exception if endog has more than 2 columns"""
    try:
        RCR(np.asarray(exog_df), np.asarray(exog_df))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_exog1d(endog_df, exog_df):
    """raise exception if exog is a 1-d array"""
    try:
        RCR(np.asarray(endog_df), np.asarray(exog_df)[:, 1])
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_exog1c(endog_df, exog_df):
    """raise exception if exog has only 1 column"""
    try:
        RCR(np.asarray(endog_df), np.asarray(exog_df)[:, :1])
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_nointercept(endog_df, exog_df):
    """raise exception if exog does not have an intercept"""
    try:
        RCR(np.asarray(endog_df), np.asarray(exog_df)[:, 1:])
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rc_nobsnoteq(endog_df, exog_df):
    """raise exception if exog and endog have different # of rows"""
    try:
        RCR(np.asarray(endog_df), np.asarray(exog_df)[0:100, ])
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_lr2d(endog_df, exog_df):
    """raise exception if lambda_range is a 2-d array"""
    try:
        RCR(endog_df, exog_df, lambda_range=np.asarray(endog_df))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_lr1e(endog_df, exog_df):
    """raise exception if lambda_range has wrong # of elements"""
    try:
        RCR(endog_df, exog_df, lambda_range=np.zeros(1))
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_lrnan(endog_df, exog_df):
    """raise exception if lambda_range includes NaN values (inf is OK)"""
    try:
        RCR(endog_df, exog_df, lambda_range=np.asarray([0, np.nan]))
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rc_lrnotsorted(endog_df, exog_df):
    """raise exception if lambda_range out of order"""
    try:
        RCR(endog_df, exog_df, lambda_range=np.asarray([1., 0.]))
    except ValueError:
        pass
    else:
        raise AssertionError


# Add clusters and weights
def test_rc_weights(endog, exog, weights):
    """construct model with weights"""
    model = RCR(endog, exog, weights=weights)
    assert all(model.weights == weights)
    assert model.weights_name == "TCHID_wt"
    assert model.nobs == 3325


def test_rc_weights_noname(endog, exog, weights):
    """construct model with weights (no variable name)"""
    model = RCR(endog, exog, weights=np.asarray(weights))
    assert all(model.weights == weights)
    assert model.weights_name == "(no name)"
    assert model.nobs == 3325


def test_rc_zeroweights(endog, exog):
    """estimate with weights"""
    weights = np.zeros(len(endog))
    try:
        RCR(endog, exog, weights=weights)
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rc_infweights(endog, exog):
    """estimate with weights"""
    weights = np.full((len(endog),), np.inf)
    try:
        RCR(endog, exog, weights=weights)
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rc_shortweights(endog, exog):
    """estimate with weights"""
    weights = np.ones((len(endog)-1,))
    try:
        RCR(endog, exog, weights=weights)
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_strweights(endog, exog):
    """estimate with weights"""
    weights = np.full(len(endog), "hello")
    try:
        RCR(endog, exog, weights=weights)
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_wrongdimweights(endog, exog):
    """estimate with weights"""
    weights = np.ones((len(endog), 2))
    try:
        RCR(endog, exog, weights=weights)
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rc_clusters(endog, exog, clusters):
    """construct model with clusters"""
    model = RCR(endog, exog, groupvar=clusters)
    assert all(model.groupvar == clusters)
    assert model.groupvar_name == "TCHID"
    assert model.ngroups == 323


def test_rc_clusters_noname(endog, exog, clusters):
    """construct model with clusters (no variable name)"""
    model = RCR(endog, exog, groupvar=np.asarray(clusters))
    assert all(model.groupvar == clusters)
    assert model.groupvar_name == "(no name)"
    assert model.ngroups == 323


def test_rc_clust_and_wt(endog, exog, clusters, weights):
    """construct model with weights and clusters"""
    model = RCR(endog, exog, weights=weights, groupvar=clusters)
    assert all(model.weights == weights)
    assert all(model.groupvar == clusters)
    assert model.weights_name == "TCHID_wt"
    assert model.groupvar_name == "TCHID"
    assert model.nobs == 3325
    assert model.ngroups == 184


def test_rc_badci(endog, exog):
    """estimate with weights"""
    try:
        RCR(endog, exog, citype="Some unsupported type")
    except ValueError:
        pass
    else:
        raise AssertionError
