"""
TEST_RP.PY: Unit tests for rcrplot() method
"""
import sys

import numpy as np
import pandas as pd
import pytest
import patsy

sys.path.append("./")
sys.path.append("../")
from rcr import RCR  # pylint: disable=wrong-import-position


@pytest.fixture
def dat():
    """get test data from web page"""
    fname = "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
    return pd.read_stata(fname)


@pytest.fixture
def rcr_formula():
    """construct formula for test example"""
    rcr_left = "SAT + Small_Class ~ "
    rcr_right1 = "White_Asian + Girl + Free_Lunch + White_Teacher + "
    rcr_right2 = "Teacher_Experience + Masters_Degree"
    return rcr_left + rcr_right1 + rcr_right2


@pytest.fixture
def endog(dat, rcr_formula):
    """get endogenous variables"""
    endog = patsy.dmatrices(rcr_formula, dat)[0]
    return endog


@pytest.fixture
def exog(dat, rcr_formula):
    """get endogenous variables"""
    exog = patsy.dmatrices(rcr_formula, dat)[1]
    return exog


@pytest.fixture
def model(endog, exog):
    """construct RCR mdodel"""
    return RCR(endog, exog)


@pytest.fixture
def results(model):
    """fit RCR mdodel"""
    return model.fit()


# Basic functionality
def test_rp_basic(results):
    """plot rcr function with default arguments"""
    lambdavals, thetavals = results.model.lambdavals(add_thetastar=True)
    ax = results.rcrplot()
    assert type(ax).__name__ == "AxesSubplot"
    assert ax.get_title() == ""
    assert ax.get_xlabel() == 'Effect ($\\beta_x$)'
    assert ax.get_ylabel() == 'Relative correlation ($\\lambda$)'
    flabel = '$\\lambda(\\beta_x)$ function'
    assert ax.get_legend_handles_labels()[1][0] == flabel
    assert ax.get_xlim() == (-55.0, 55.0)
    assert ax.get_ylim() == pytest.approx((-209.92807503735216,
                                           398.0573074081345))
    assert np.all(ax.get_lines()[0].get_xdata() == thetavals)
    assert np.all(np.logical_or(ax.get_lines()[0].get_ydata() == lambdavals,
                                np.isnan(lambdavals)))
    ax.clear()


def test_rp_xlim1(results):
    """plot rcr function with alternate xlim (pair)"""
    xlim = (0, 6)
    thetavals = np.linspace(xlim[0], xlim[1], 100)
    lambdavals, thetavals = results.model.lambdavals(thetavals,
                                                     add_thetastar=True)
    ax = results.rcrplot(xlim=xlim)
    assert np.all(ax.get_lines()[0].get_xdata() == thetavals)
    assert np.all(np.logical_or(ax.get_lines()[0].get_ydata() == lambdavals,
                                np.isnan(lambdavals)))
    ax.clear()


def test_rp_xlim2(results):
    """plot rcr function with alternate xlim (sequence)"""
    xlim = (0, 5, 6)
    thetavals = np.asarray(xlim)
    lambdavals, thetavals = results.model.lambdavals(thetavals,
                                                     add_thetastar=True)
    ax = results.rcrplot(xlim=xlim)
    assert np.all(ax.get_lines()[0].get_xdata() == thetavals)
    assert np.all(np.logical_or(ax.get_lines()[0].get_ydata() == lambdavals,
                                np.isnan(lambdavals)))
    ax.clear()


def test_rp_tsline(results):
    """plot rcr function with optional thetastar line"""
    ax = results.rcrplot(tsline=True)
    assert ax.get_legend_handles_labels()[1][1] == '$\\beta_x^{\\infty}$'
    assert ax.get_lines()[1].get_xdata() == [results.params[1]] * 2
    assert ax.get_lines()[1].get_ydata() == [0, 1]
    ax.clear()


def test_rp_lsline(results):
    """plot rcr function with optional lambdastar line"""
    ax = results.rcrplot(lsline=True)
    assert ax.get_legend_handles_labels()[1][1] == '$\\lambda^{\\infty}$'
    assert np.all(ax.get_lines()[1].get_xdata() == [0, 1])
    assert np.all(ax.get_lines()[1].get_ydata() == [results.params[0]] * 2)
    ax.clear()


def test_rp_idset(results):
    """plot rcr function with optional identified set polygon`"""
    ax = results.rcrplot(idset=True)
    # I don't know how to check polygons, so we will only make sure it runs
    ax.clear()


def test_rp_setlab(results):
    """plot rcr function with labels set by hand"""
    ax = results.rcrplot(tsline=True,
                         lsline=True,
                         idset=True,
                         title="TITLE",
                         xlabel="XLABEL",
                         ylabel="YLABEL",
                         flabel="FLABEL",
                         tslabel="TSLABEL",
                         lslabel="LSLABEL",
                         idlabels=("ID0", "ID1"))
    assert ax.get_title() == "TITLE"
    assert ax.get_xlabel() == "XLABEL"
    assert ax.get_ylabel() == "YLABEL"
    assert ax.get_legend_handles_labels()[1] == ['FLABEL',
                                                 'TSLABEL',
                                                 'LSLABEL',
                                                 'ID0',
                                                 'ID1']
    ax.clear()
