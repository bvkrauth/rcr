"""
TEST_RF.RR Unit tests for RCRResults object and its methods
"""
import numpy as np
import pytest

from rcrbounds import RCRResults


# Basic functionality
def test_rr_basic(model):
    """check that fit produces RCRResults object"""
    results = model.fit()
    assert isinstance(results, RCRResults)


# Methods
# pylint: disable=duplicate-code
def test_rr_se(results):
    """calculate standard errors with the se() method"""
    truese = np.asarray([2.09826858,  30.60745128, 108.51947421,   0.95693751,
                         0.6564318])
    test_se = results.params_se()
    assert test_se == pytest.approx(truese)


def test_rr_z(results):
    """calculate z-statistics with the z() method"""
    truez = np.asarray([5.86702731, 0.26691899, 0.26663868,
                        5.36612236, 7.92390398])
    test_z = results.params_z()
    assert test_z == pytest.approx(truez)


def test_rr_pz(results):
    """calculate p-values with the pz() method"""
    truepz = np.asarray([4.43677606e-09, 7.89531535e-01, 7.89747372e-01,
                         8.04473756e-08, 2.22044605e-15])
    test_pz = results.params_pvalue()
    assert test_pz == pytest.approx(truepz)


def test_rr_ci(results):
    """calculate confidence intervals with default options"""
    trueci = np.asarray([[8.19806824, -51.8197922, -183.75877191, 3.25948071,
                          3.91491988],
                        [16.42312995, 68.15921213, 241.62975025, 7.01060682,
                         6.48808526]])
    test_ci = results.params_ci()
    assert test_ci == pytest.approx(trueci)


def test_rr_ci90(results):
    """calculate confidence intervals with optional cilevel = 90"""
    trueci = np.asarray([[8.8592544, -42.17506728, -149.56316158, 3.56102163,
                          4.12176834],
                         [15.76194378, 58.51448721, 207.43413992, 6.7090659,
                          6.2812368]])
    test_ci = results.params_ci(cilevel=90)
    assert test_ci == pytest.approx(trueci)


def test_rr_cistr(results):
    """raise exception if cilevel is non-numeric"""
    try:
        results.params_ci(cilevel="this should be a number")
    except TypeError:
        pass
    else:
        raise AssertionError


def test_rr_cineg(results):
    """raise exception if cilevel is out of range"""
    try:
        results.params_ci(cilevel=-50)
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


def test_rr_bcibad(results):
    """incorrectly-named confidence interval, should return NaN"""
    ci1 = results.betax_ci(citype="Some Unsupported Type")
    assert np.all(np.isnan(ci1))


def test_rr_testbetax(results):
    """test_betax() method with default options"""
    test_t0 = results.test_betax()
    test_t1 = results.test_betax(0.)
    test_t2 = results.test_betax(5.2)
    assert test_t0 == pytest.approx(1.1920928955078125e-07)
    assert test_t1 == pytest.approx(1.1920928955078125e-07)
    assert test_t2 == 1.0


def test_rr_summary(results):
    """summary() method with default options"""
    # pylint: disable=unidiomatic-typecheck
    test_result = results.summary()
    assert type(test_result).__name__ == "Summary"
    assert type(test_result.tables) == list
    assert len(test_result.tables) == 3
    assert len(test_result.extra_txt) > 0


def test_rr_summary_weights(model, weights):
    """summary() method with weights"""
    test_result = model.fit(weights=weights).summary()
    assert isinstance(test_result.tables, list)
    assert len(test_result.tables) == 3
    assert len(test_result.extra_txt) > 0


@pytest.mark.skip(reason="not yet working")
def test_rr_summary_cluster(model, clusters):
    """summary() method with weights"""
    test_result = model.fit(groupvar=clusters, cov_type="cluster").summary()
    assert isinstance(test_result.tables, list)
    assert len(test_result.tables) == 3
    assert len(test_result.extra_txt) > 0


def test_rr_noid(model):
    """handle when identified set is (-inf, inf)"""
    results = model.fit(lambda_range=np.asarray([0.0, np.inf]))
    assert np.isneginf(results.params[3])
    assert np.isposinf(results.params[4])
    msk = np.full((5, 5), True)
    msk[0:3, 0:3] = False
    assert all(results.cov_params[msk] == 0.)
    assert all(results.params_se()[3:] == 0)
    assert np.isneginf(results.params_z()[3])
    assert np.isposinf(results.params_z()[4])
    assert all(results.params_pvalue()[3:] == 0.0)
    assert all(np.isneginf(results.params_ci()[:, 3]))
    assert all(np.isposinf(results.params_ci()[:, 4]))
    assert np.isneginf(results.betax_ci_conservative()[0])
    assert np.isposinf(results.betax_ci_conservative()[1])
    assert np.isneginf(results.betax_ci_upper()[0])
    assert np.isposinf(results.betax_ci_upper()[1])
    assert np.isneginf(results.betax_ci_lower()[0])
    assert np.isposinf(results.betax_ci_lower()[1])
    assert np.isneginf(results.betax_ci_imbensmanski()[0])
    assert np.isposinf(results.betax_ci_imbensmanski()[1])
