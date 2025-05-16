"""
TEST_LF.PY: Unit tests for rcfast()
"""
import numpy as np
import pytest

from rcrbounds import rcfast, effectinf, rcinf, simplify_moments
from rcrbounds import scalar_rcfast, negative_rcfast

# Basic functionality
def test_lf_basic():
    """rcfast for simple data and scalar effect"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lf_true = 0.5773502691896257
    test_lf = rcfast(0.0, simplify_moments(mv1))
    assert test_lf == pytest.approx(lf_true)
    assert type(test_lf) is np.ndarray


def test_lf_realdata(moment_vector):
    """rcfast for real data and array effect"""
    lf_true = np.array([28.93548917, 26.67790368])
    test_lf = rcfast(np.array([0.0, 1.0, ]),
                     simplify_moments(moment_vector))
    assert test_lf == pytest.approx(lf_true)
    assert type(test_lf) is np.ndarray


def test_lf_scalar():
    """rcfast for simple data and scalar effect"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lf_true = 0.5773502691896257
    test_lf = scalar_rcfast(0.0, simplify_moments(mv1))
    assert test_lf == pytest.approx(lf_true)
    assert type(test_lf) is float


def test_lf_negative():
    """rcfast for simple data and scalar effect"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lf_true = -0.5773502691896257
    test_lf = negative_rcfast(0.0, simplify_moments(mv1))
    assert test_lf == pytest.approx(lf_true)
    assert type(test_lf) is float


# Special cases
def test_lf_effectinf():
    """function is undefined (NaN) for effect = effectinf"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    effect_inf = effectinf(mv1)
    test_lf = rcfast(np.array([effect_inf - 0.01,
                               effect_inf,
                               effect_inf + 0.01]),
                     simplify_moments(mv1))
    assert np.isnan(test_lf[1]) and all(np.isfinite(test_lf[[0, 2]]))


def test_lf_bigeffect(moment_vector):
    """function is close to rcinf for large effect"""
    # This test fails for higher values of bignum
    bignum = 1.0e100
    rc_inf = rcinf(moment_vector)
    test_lf = rcfast(np.array([-bignum, bignum]),
                     simplify_moments(moment_vector))
    assert test_lf == pytest.approx(rc_inf)
