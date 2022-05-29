"""
TEST_RC.RC Unit tests for RCR(), not including the fit method
"""
import sys

import numpy as np
import pandas as pd
import patsy
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import RCR


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


# Basic functionality
# Patsy design matrices
def test_rc_patsy(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
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
    cv = ("White_Asian Girl Free_Lunch White_Teacher " +
          "Teacher_Experience Masters_Degree")
    assert model.controlvars == cv
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0


# Data frames (with Patsy design_info)
def test_rc_patsy_df(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    model = RCR(endog, exog)
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
    cv = ("White_Asian Girl Free_Lunch White_Teacher " +
          "Teacher_Experience Masters_Degree")
    assert model.controlvars == cv
    assert model.nobs == 5839
    assert isinstance(model.lambda_range, np.ndarray)
    assert model.lambda_range.shape == (2, )
    assert model.lambda_range[0] == 0.0
    assert model.lambda_range[1] == 1.0


# Data frames (without Patsy design_info)
def test_rc_dataframe(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    model = RCR(pd.DataFrame(endog), pd.DataFrame(exog))
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


# Plain numpy arrays
def test_rc_array(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat)
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


# Set lambda_range
def test_rc_setlr(endog, exog):
    model = RCR(endog, exog, lambda_range=np.asarray([-1.0, 5.0]))
    assert model.lambda_range[0] == -1.0
    assert model.lambda_range[1] == 5.0


# Set lambda_range
def test_rc_setlrinf(endog, exog):
    model = RCR(endog, exog, lambda_range=np.asarray([-np.inf, np.inf]))
    assert np.isneginf(model.lambda_range[0])
    assert np.isposinf(model.lambda_range[1])


# Exceptions


# endog is a 1-d array, should be a 2-d array
def test_rc_endog1d(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(endog)[:, 1], np.asarray(exog))
    except TypeError:
        pass
    else:
        raise AssertionError


# endog has only 1 column, should have 2
def test_rc_endog1c(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(endog)[:, 1:], np.asarray(exog))
    except TypeError:
        pass
    else:
        raise AssertionError


# endog has more than 2 columns
def test_rc_endog7c(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(exog), np.asarray(exog))
    except TypeError:
        pass
    else:
        raise AssertionError


# exog is a 1-d array, should be a 2-d array
def test_rc_exog1d(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(endog), np.asarray(exog)[:, 1])
    except TypeError:
        pass
    else:
        raise AssertionError


# exog has only 1 column, should have at least 2
def test_rc_exog1c(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(endog), np.asarray(exog)[:, :1])
    except TypeError:
        pass
    else:
        raise AssertionError


# endog and exog do not have the same number of rows
def test_rc_nobsnoteq(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(np.asarray(endog), np.asarray(exog)[0:100, ])
    except TypeError:
        pass
    else:
        raise AssertionError


# lambda_range is a 2-d array, should be 1-d
def test_rc_lr2d(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(endog, exog, lambda_range=np.asarray(endog))
    except TypeError:
        pass
    else:
        raise AssertionError


# lambda_range has the wrong number of elements
def test_rc_lr1e(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(endog, exog, lambda_range=np.zeros(1))
    except TypeError:
        pass
    else:
        raise AssertionError


# lambda_range includes NaN values (inf values are OK)
def test_rc_lrnan(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(endog, exog, lambda_range=np.asarray([0, np.nan]))
    except ValueError:
        pass
    else:
        raise AssertionError


# lambda_range is not in weakly ascending order
def test_rc_lrnotsorted(dat, rcr_formula):
    endog, exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")
    try:
        RCR(endog, exog, lambda_range=np.asarray([1., 0.]))
    except ValueError:
        pass
    else:
        raise AssertionError
