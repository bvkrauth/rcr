"""
TEST_RF.RR Unit tests for RCRResults object and its methods
"""
import sys

import numpy as np
import pandas as pd
import patsy
import pytest

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


@pytest.fixture
def results(model):
    return model.fit()


# Basic functionality
def test_rr_basic(model):
    results = model.fit()
    assert isinstance(results, RCR_results)


# se() method
def test_rr_se(results):
    truese = np.asarray([2.09826858,  30.60745128, 108.51947421,   0.95693751,
                         0.6564318])
    se = results.se()
    assert se == pytest.approx(truese)


# z() method
def test_rr_z(results):
    truez = np.asarray([5.86702731, 0.26691899, 0.26663868,
                        5.36612236, 7.92390398])
    z = results.z()
    assert z == pytest.approx(truez)


# pz() method
def test_rr_pz(results):
    truepz = np.asarray([4.43677606e-09, 7.89531535e-01, 7.89747372e-01,
                         8.04473756e-08, 2.22044605e-15])
    pz = results.pz()
    assert pz == pytest.approx(truepz)


# ci() method, default options
def test_rr_ci(results):
    trueci = np.asarray([[8.19806824, -51.8197922, -183.75877191, 3.25948071,
                          3.91491988],
                        [16.42312995, 68.15921213, 241.62975025, 7.01060682,
                         6.48808526]])
    ci = results.ci()
    assert ci == pytest.approx(trueci)


# ci() method, optional cilevel = 90
def test_rr_ci90(results):
    trueci = np.asarray([[8.8592544, -42.17506728, -149.56316158, 3.56102163,
                          4.12176834],
                         [15.76194378, 58.51448721, 207.43413992, 6.7090659,
                          6.2812368]])
    ci = results.ci(cilevel=90)
    assert ci == pytest.approx(trueci)


# ci() method, invalid cilevel (string)
def test_rr_cistr(results):
    try:
        results.ci(cilevel="this should be a number")
    except TypeError:
        pass
    else:
        raise AssertionError


# ci() method, invalid cilevel (bad number)
def test_rr_cineg(results):
    try:
        results.ci(cilevel=-50)
    except ValueError:
        pass
    else:
        raise AssertionError


# betaxCI_conservative() method, default options
def test_rr_bciconservative(results):
    trueci = np.asarray([3.25948071, 6.48808526])
    ci = results.betaxCI_conservative()
    assert ci == pytest.approx(trueci)


# betaxCI_upper() method, default options
def test_rr_bciupper(results):
    trueci = np.asarray([3.56102163, np.inf])
    ci = results.betaxCI_upper()
    assert ci == pytest.approx(trueci)


# betaxCI_upper() method, default options
def test_rr_bcilower(results):
    trueci = np.asarray([-np.inf, 6.281236804882139])
    ci = results.betaxCI_lower()
    assert ci == pytest.approx(trueci)


# betaxCI_imbensmanski() method, default options
def test_rr_bciimbensmanski(results):
    trueci = np.asarray([3.29158006, 6.46606603])
    ci = results.betaxCI_imbensmanski()
    assert ci == pytest.approx(trueci)


# handling when identified set is (-inf, inf)
def test_rr_noid(model):
    results = model.fit(lambda_range=np.asarray([0.0, np.inf]))
    assert np.isneginf(results.params[3])
    assert np.isposinf(results.params[4])
    msk = np.full((5, 5), True)
    msk[0:3, 0:3] = False
    assert all(results.cov_params[msk] == 0.)
    assert all(results.se()[3:] == 0)
    assert np.isneginf(results.z()[3])
    assert np.isposinf(results.z()[4])
    assert all(results.pz()[3:] == 0.0)
    assert all(np.isneginf(results.ci()[:, 3]))
    assert all(np.isposinf(results.ci()[:, 4]))
    assert np.isneginf(results.betaxCI_conservative()[0])
    assert np.isposinf(results.betaxCI_conservative()[1])
    assert np.isneginf(results.betaxCI_upper()[0])
    assert np.isposinf(results.betaxCI_upper()[1])
    assert np.isneginf(results.betaxCI_lower()[0])
    assert np.isposinf(results.betaxCI_lower()[1])
    assert np.isneginf(results.betaxCI_imbensmanski()[0])
    assert np.isposinf(results.betaxCI_imbensmanski()[1])
