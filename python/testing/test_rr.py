"""
TEST_RF.RR Unit tests for RCRResults object and its methods
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import RCRResults  # pylint: disable=wrong-import-position


# Basic functionality
def test_rr_basic(model):
    """check that fit produces RCRResults object"""
    results = model.fit()
    assert isinstance(results, RCRResults)


# Methods
def test_rr_se(results):
    """calculate standard errors with the se() method"""
    truese = np.asarray([2.09826858,  30.60745128, 108.51947421,   0.95693751,
                         0.6564318])
    se = results.se()
    assert se == pytest.approx(truese)


def test_rr_z(results):
    """calculate z-statistics with the z() method"""
    truez = np.asarray([5.86702731, 0.26691899, 0.26663868,
                        5.36612236, 7.92390398])
    z = results.z()
    assert z == pytest.approx(truez)


def test_rr_pz(results):
    """calculate p-values with the pz() method"""
    truepz = np.asarray([4.43677606e-09, 7.89531535e-01, 7.89747372e-01,
                         8.04473756e-08, 2.22044605e-15])
    pz = results.pz()
    assert pz == pytest.approx(truepz)


def test_rr_ci(results):
    """calculate confidence intervals with default options"""
    trueci = np.asarray([[8.19806824, -51.8197922, -183.75877191, 3.25948071,
                          3.91491988],
                        [16.42312995, 68.15921213, 241.62975025, 7.01060682,
                         6.48808526]])
    ci = results.ci()
    assert ci == pytest.approx(trueci)


def test_rr_ci90(results):
    """calculate confidence intervals with optional cilevel = 90"""
    trueci = np.asarray([[8.8592544, -42.17506728, -149.56316158, 3.56102163,
                          4.12176834],
                         [15.76194378, 58.51448721, 207.43413992, 6.7090659,
                          6.2812368]])
    ci = results.ci(cilevel=90)
    assert ci == pytest.approx(trueci)


def test_rr_cistr(results):
    """raise exception if cilevel is non-numeric"""
    try:
        results.ci(cilevel="this should be a number")
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rr_cineg(results):
    """raise exception if cilevel is out of range"""
    try:
        results.ci(cilevel=-50)
    except ValueError:
        pass
    else:
        raise AssertionError


def test_rr_bciconservative(results):
    """conservative confidence interval, default options"""
    trueci = np.asarray([3.25948071, 6.48808526])
    ci1 = results.betax_ci_conservative()
    ci2 = results.betax_ci(citype="conservative")
    assert ci1 == pytest.approx(trueci)
    assert ci2 == pytest.approx(trueci)


def test_rr_bciupper(results):
    """upper confidence interval, default options"""
    trueci = np.asarray([3.56102163, np.inf])
    ci1 = results.betax_ci_upper()
    ci2 = results.betax_ci(citype="upper")
    assert ci1 == pytest.approx(trueci)
    assert ci2 == pytest.approx(trueci)


def test_rr_bcilower(results):
    """lower confidence interval, default options"""
    trueci = np.asarray([-np.inf, 6.281236804882139])
    ci1 = results.betax_ci_lower()
    ci2 = results.betax_ci(citype="lower")
    assert ci1 == pytest.approx(trueci)
    assert ci2 == pytest.approx(trueci)


def test_rr_bciimbensmanski(results):
    """Imbens-Manski confidence interval, default options"""
    trueci = np.asarray([3.29158006, 6.46606603])
    ci1 = results.betax_ci_imbensmanski()
    ci2 = results.betax_ci(citype="Imbens-Manski")
    assert ci1 == pytest.approx(trueci)
    assert ci2 == pytest.approx(trueci)


def test_rr_testbetax(results):
    """test_betax() method with default options"""
    t0 = results.test_betax()
    t1 = results.test_betax(0.)
    t2 = results.test_betax(5.2)
    assert t0 == pytest.approx(1.1920928955078125e-07)
    assert t1 == pytest.approx(1.1920928955078125e-07)
    assert t2 == 1.0


def test_rr_summary(results):
    """summary() method with default options"""
    summary = results.summary()
    assert type(summary).__name__ == "Summary"
    assert type(summary.tables) == list
    assert len(summary.tables) == 3
    assert len(summary.extra_txt) > 0


def test_rr_noid(model):
    """handle when identified set is (-inf, inf)"""
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
    assert np.isneginf(results.betax_ci_conservative()[0])
    assert np.isposinf(results.betax_ci_conservative()[1])
    assert np.isneginf(results.betax_ci_upper()[0])
    assert np.isposinf(results.betax_ci_upper()[1])
    assert np.isneginf(results.betax_ci_lower()[0])
    assert np.isposinf(results.betax_ci_lower()[1])
    assert np.isneginf(results.betax_ci_imbensmanski()[0])
    assert np.isposinf(results.betax_ci_imbensmanski()[1])
