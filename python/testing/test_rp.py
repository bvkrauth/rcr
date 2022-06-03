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


@pytest.fixture
def model(endog, exog):
    return RCR(endog, exog)


@pytest.fixture
def results(model):
    return model.fit()


# Basic functionality
def test_rp_basic(results):
    """
    Test rcrplot with default arguments
    """
    thetavals, lambdavals = results._lambdafun()
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
    """
    Test rcrplot with alternate xlim pair
    """
    xlim = (0, 6)
    thetavals = np.linspace(xlim[0], xlim[1], 100)
    thetavals, lambdavals = results._lambdafun(thetavals)
    ax = results.rcrplot(xlim=xlim)
    assert np.all(ax.get_lines()[0].get_xdata() == thetavals)
    assert np.all(np.logical_or(ax.get_lines()[0].get_ydata() == lambdavals,
                                np.isnan(lambdavals)))
    ax.clear()


def test_rp_xlim2(results):
    """
    Test rcrplot with alternate xlim sequence
    """
    xlim = (0, 5, 6)
    thetavals = np.asarray(xlim)
    thetavals, lambdavals = results._lambdafun(thetavals)
    ax = results.rcrplot(xlim=xlim)
    assert np.all(ax.get_lines()[0].get_xdata() == thetavals)
    assert np.all(np.logical_or(ax.get_lines()[0].get_ydata() == lambdavals,
                                np.isnan(lambdavals)))
    ax.clear()


def test_rp_tsline(results):
    """
    Test rcrplot with added thetastar line
    """
    ax = results.rcrplot(tsline=True)
    assert ax.get_legend_handles_labels()[1][1] == '$\\beta_x^{\\infty}$'
    assert ax.get_lines()[1].get_xdata() == [results.params[1]] * 2
    assert ax.get_lines()[1].get_ydata() == [0, 1]
    ax.clear()


def test_rp_lsline(results):
    """
    Test rcrplot with added lambdastar line
    """
    ax = results.rcrplot(lsline=True)
    assert ax.get_legend_handles_labels()[1][1] == '$\\lambda^{\\infty}$'
    assert np.all(ax.get_lines()[1].get_xdata() == [0, 1])
    assert np.all(ax.get_lines()[1].get_ydata() == [results.params[0]] * 2)
    ax.clear()


def test_rp_idset(results):
    """
    Test rcrplot with added polygons showing the identified set
    """
    ax = results.rcrplot(idset=True)
    # I don't know how to check polygons, so we will only make sure it runs
    ax.clear()


def test_rp_setlab(results):
    """
    Test rcrplot with various labels set by hand
    """
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
