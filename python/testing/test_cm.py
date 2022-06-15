"""
TEST_CM.PY: Unit tests for check_moments()
"""
import numpy as np
import pytest

from rcrbounds import check_moments, read_data


# Basic functionality
def test_cm_basic():
    """check moments with simple valid data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    assert check_moments(mv1) == (True, True)


def test_cm_realdata():
    """check moments from real data"""
    moment_vector = read_data("testing/testin1.txt")[3]
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
    """return mics of zero and NaN if all zeros"""
    moment_vector = np.zeros(9)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_integer():
    """check moments with integer data"""
    moment_vector = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    result = check_moments(moment_vector)
    assert result == (True, True)


def test_cm_varx0():
    """return NaN if var(x) = 0"""
    moment_vector = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varxneg():
    """return numeric values if var(x) < 0"""
    moment_vector = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varyneg():
    """return numeric values if var(y) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varzneg():
    """return numeric values if var(z) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
