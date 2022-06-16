"""
TEST_RD.PY Unit tests for read_data()
"""
import os


import pytest

from rcrbounds import read_data


# Basic functionality
def test_rd_basic(infile):
    """read data from the specified text file"""
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data(infile)
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


def test_rd_badfile(read_only_file):
    """raise an exception if data in wrong format"""
    try:
        read_data(read_only_file)
    except RuntimeError:
        pass
    else:
        raise AssertionError


def test_rd_badnmoments(badin1):
    """reset n_moments and warn if it does not match moment_vector"""
    with pytest.warns(UserWarning, match="n_moments reset"):
        n_moments = read_data(badin1)[0]
    assert n_moments == 44


def test_rd_badnlambda(badin2):
    """reset n_lmambda and warn if it does not match lambda_range"""
    with pytest.warns(UserWarning, match="n_lambda reset"):
        n_lambda = read_data(badin2)[1]
    assert n_lambda == 1


def test_rd_badmoments(badin3):
    """raise exception if n_moments is invalid"""
    try:
        read_data(badin3)
    except AssertionError:
        pass
    else:
        raise AssertionError


@pytest.mark.skip(reason="not yet implemented")
def test_rd_badlambda(badin4):
    """raise exception if n_lambda is invalid"""
    try:
        read_data(badin4)
    except AssertionError:
        pass
    else:
        raise AssertionError


def test_rd_badbignum(badin5):
    """raise exception if extarnal_big_number is invalid"""
    try:
        read_data(badin5)
    except AssertionError:
        pass
    else:
        raise AssertionError
