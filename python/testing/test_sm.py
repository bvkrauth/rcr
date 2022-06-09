"""
TEST_SM.PY: Unit tests for simplify_moments()
"""


import numpy as np
import pytest

from rcr import simplify_moments


# Basic functionality
def test_sm_basic():
    """get simple moments from simple test data"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    true_sm = np.array([1, 1, .5, .25, .25, .25])
    test_sm = simplify_moments(moment_vector)
    assert test_sm == pytest.approx(true_sm)
    assert test_sm[5] == (np.sign(test_sm[5]) *
                          np.sqrt(test_sm[3] * test_sm[4]))


def test_sm_realdata(moment_vector):
    """get simple moments from real test data"""
    true_sm = np.array([5.42538313e+02, 2.05839484e-01,
                        1.07467966e+00, 4.47643916e+01,
                        1.34931719e-03, 1.10235301e-02])
    test_sm = simplify_moments(moment_vector)
    assert test_sm == pytest.approx(true_sm)
    assert test_sm[5] <= (np.sign(test_sm[5]) *
                          np.sqrt(test_sm[3] * test_sm[4]))


# Special cases
def test_sm_wronglen():
    """raise exception if wrong length"""
    moment_vector = np.zeros(8)
    try:
        simplify_moments(moment_vector)
    except AssertionError:
        pass
    else:
        raise AssertionError("simplify_moments has accepted invalid input")


def test_sm_allzeros():
    """return mix of zeros and NaN if input is all zeros"""
    moment_vector = np.zeros(9)
    true_sm = np.array([0, 0, 0, np.nan, np.nan, np.nan])
    test_sm = simplify_moments(moment_vector)
    assert test_sm == pytest.approx(true_sm, nan_ok=True)


def test_sm_integer():
    """handle integer data"""
    moment_vector = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    true_sm = np.array([2, 2, 1, .5, .5, .5])
    test_sm = simplify_moments(moment_vector)
    assert test_sm == pytest.approx(true_sm)
    assert test_sm[5] == (np.sign(test_sm[5]) *
                          np.sqrt(test_sm[3] * test_sm[4]))


def test_sm_varx0():
    """return NaN if var(x)=0"""
    moment_vector = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    true_sm = np.array([1, 1, 0.5, np.nan, np.nan, np.nan])
    test_sm = simplify_moments(moment_vector)
    assert test_sm == pytest.approx(true_sm, nan_ok=True)


def test_sm_varxneg():
    """return numeric values if var(x)<0"""
    moment_vector = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    test_sm = simplify_moments(moment_vector)
    assert all(np.isfinite(test_sm))


def test_sm_varyneg():
    """return numeric values if var(y) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    test_sm = simplify_moments(moment_vector)
    assert all(np.isfinite(test_sm))
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    test_sm = simplify_moments(moment_vector)
    assert all(np.isfinite(test_sm))


def test_sm_varzneg():
    """return numeric values if var(z) < 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    test_sm = simplify_moments(moment_vector)
    assert all(np.isfinite(test_sm))
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    test_sm = simplify_moments(moment_vector)
    assert all(np.isfinite(test_sm))
