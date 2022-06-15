"""
TEST_RD.PY Unit tests for read_data()
"""
import os


import pytest

from rcrbounds import read_data


# Basic functionality
def test_rd_basic():
    """read data from the specified text file"""
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    assert n_moments == 44
    assert n_lambda == 1
    assert external_big_number.dtype == "float"
    assert moment_vector.shape == (44, )
    assert lambda_range.shape == (2, )

# Exceptions


def test_rd_nonexistent():
    """raise an exception if file does not exist"""
    assert ~os.path.exists("nonexistent-file")
    try:
        read_data("nonexistent-file")
    except RuntimeError:
        pass
    else:
        raise AssertionError

# File cannot be opened for reading
# NOT YET TESTED


def test_rd_badfile():
    """raise an exception if data in wrong format"""
    try:
        read_data("testing/out.txt")
    except RuntimeError:
        pass
    else:
        raise AssertionError


def test_rd_badnmoments():
    """reset n_moments and warn if it does not match moment_vector"""
    with pytest.warns(UserWarning, match="n_moments reset"):
        n_moments = read_data("testing/badin1.txt")[0]
    assert n_moments == 44


def test_rd_badnlambda():
    """reset n_lmambda and warn if it does not match lambda_range"""
    with pytest.warns(UserWarning, match="n_lambda reset"):
        n_lambda = read_data("testing/badin2.txt")[1]
    assert n_lambda == 1


def test_rd_badmoments():
    """raise exception if n_moments is invalid"""
    try:
        read_data("testing/badin3.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError


@pytest.mark.skip(reason="not yet implemented")
def test_rd_badlambda():
    """raise exception if n_lambda is invalid"""
    try:
        read_data("testing/badin4.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError


def test_rd_badbignum():
    """raise exception if extarnal_big_number is invalid"""
    try:
        read_data("testing/badin5.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError
