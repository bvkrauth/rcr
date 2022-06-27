"""
TEST_RP.PY: Unit tests for rcrplot() method
"""
import numpy as np
import pytest


# Basic functionality
def test_rp_basic(results):
    """plot rcr function with default arguments"""
    rcvals, effectvals = results.model.rcvals(add_effectinf=True)
    test_ax = results.rcrplot()
    assert type(test_ax).__name__ == "AxesSubplot"
    assert test_ax.get_title() == ""
    assert test_ax.get_xlabel() == 'Effect ($\\beta_x$)'
    assert test_ax.get_ylabel() == 'Relative correlation ($\\lambda$)'
    flabel = '$\\lambda(\\beta_x)$ function'
    assert test_ax.get_legend_handles_labels()[1][0] == flabel
    assert test_ax.get_xlim() == (-55.0, 55.0)
    assert test_ax.get_ylim() == pytest.approx((-209.92807503735216,
                                                398.0573074081345))
    assert np.all(test_ax.get_lines()[0].get_xdata() == effectvals)
    assert np.all(np.logical_or(test_ax.get_lines()[0].get_ydata() ==
                                rcvals,
                                np.isnan(rcvals)))
    assert test_ax.get_legend() is None
    test_ax.clear()


def test_rp_xlim1(results):
    """plot rcr function with alternate xlim (pair)"""
    xlim = (0, 6)
    effectvals = np.linspace(xlim[0], xlim[1], 100)
    rcvals, effectvals = results.model.rcvals(effectvals,
                                              add_effectinf=True)
    test_ax = results.rcrplot(xlim=xlim)
    assert np.all(test_ax.get_lines()[0].get_xdata() == effectvals)
    assert np.all(np.logical_or(test_ax.get_lines()[0].get_ydata() ==
                                rcvals,
                                np.isnan(rcvals)))
    test_ax.clear()


def test_rp_xlim2(results):
    """plot rcr function with alternate xlim (sequence)"""
    xlim = (0, 5, 6)
    effectvals = np.asarray(xlim)
    rcvals, effectvals = results.model.rcvals(effectvals,
                                              add_effectinf=True)
    test_ax = results.rcrplot(xlim=xlim)
    assert np.all(test_ax.get_lines()[0].get_xdata() == effectvals)
    assert np.all(np.logical_or(test_ax.get_lines()[0].get_ydata() ==
                                rcvals,
                                np.isnan(rcvals)))
    test_ax.clear()


def test_rp_ylim(results):
    """plot rcr function with alternate ylim"""
    ylim = np.asarray((0.0, 1.0))
    test_ax = results.rcrplot(ylim=ylim)
    assert np.all(test_ax.yaxis.get_view_interval() == ylim)
    test_ax.clear()


def test_rp_tsline(results):
    """plot rcr function with optional effectinf line"""
    test_ax = results.rcrplot(tsline=True)
    assert test_ax.get_legend_handles_labels()[1][1] == '$\\beta_x^{\\infty}$'
    assert test_ax.get_lines()[1].get_xdata() == [results.params[1]] * 2
    assert test_ax.get_lines()[1].get_ydata() == [0, 1]
    test_ax.clear()


def test_rp_lsline(results):
    """plot rcr function with optional rcinf line"""
    test_ax = results.rcrplot(lsline=True)
    assert test_ax.get_legend_handles_labels()[1][1] == '$\\lambda^{\\infty}$'
    assert np.all(test_ax.get_lines()[1].get_xdata() == [0, 1])
    assert np.all(test_ax.get_lines()[1].get_ydata() ==
                  [results.params[0]] * 2)
    test_ax.clear()


def test_rp_idset(results):
    """plot rcr function with optional identified set polygon`"""
    test_ax = results.rcrplot(idset=True)
    # I don't know how to check polygons, so we will only make sure it runs
    test_ax.clear()


def test_rp_setlab(results):
    """plot rcr function with labels set by hand"""
    test_ax = results.rcrplot(tsline=True,
                              lsline=True,
                              idset=True,
                              title="TITLE",
                              xlabel="XLABEL",
                              ylabel="YLABEL",
                              flabel="FLABEL",
                              tslabel="TSLABEL",
                              lslabel="LSLABEL",
                              idlabels=("ID0", "ID1"))
    assert test_ax.get_title() == "TITLE"
    assert test_ax.get_xlabel() == "XLABEL"
    assert test_ax.get_ylabel() == "YLABEL"
    assert test_ax.get_legend_handles_labels()[1] == ['FLABEL',
                                                      'TSLABEL',
                                                      'LSLABEL',
                                                      'ID0',
                                                      'ID1']
    test_ax.clear()


def test_rp_legend(results):
    """plot rcr function with labels set by hand"""
    test_ax = results.rcrplot(legend=True)
    assert test_ax.get_legend() is not None
    test_ax.clear()
