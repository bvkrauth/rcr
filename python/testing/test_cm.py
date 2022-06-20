"""
TEST_CM.PY: Unit tests for check_moments()
"""
import numpy as np
import pytest

from rcrbounds import check_moments


# Basic functionality
def test_cm_basic():
    """check moments with simple valid data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    assert check_moments(mv1) == (True, True)


def test_cm_realdata(moment_vector):
    """check moments from real data"""
    assert check_moments(moment_vector) == (True, True)


# Special cases
def test_cm_wronglen():
    """raise an exception if invalid length"""
    try:
        check_moments(np.zeros(8))
    except AssertionError:
        pass
    else:
        raise AssertionError("check_moments has accepted invalid input")


def test_cm_allzeros():
    """invalid if all zeros"""
    moment_vector = np.zeros(9)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_integer():
    """should accept integer data"""
    moment_vector = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    result = check_moments(moment_vector)
    assert result == (True, True)


def test_cm_varx0():
    """invalid if var(x) = 0"""
    moment_vector = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varxneg():
    """invalid if var(x) < 0"""
    moment_vector = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varyneg():
    """invalid if var(y) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varzneg():
    """invalid if var(z) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varyhatzero():
    """valid but not identified if var(yhat) = 0"""
    moment_vector = np.array([0, 0, 0, 1, 0., 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Model not identified:"):
        result = check_moments(moment_vector)
    assert result == (True, False)


def test_cm_noresid():
    """valid but not identified if corr(x,y) = 1"""
    moment_vector = np.array([0, 0, 0, 1, 1, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Model not identified:"):
        result = check_moments(moment_vector)
    assert result == (True, False)
