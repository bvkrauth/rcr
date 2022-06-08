"""
TEST_RD.PY Unit tests for read_data()
"""
import os
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import read_data


# Basic functionality
# Should read data from the specified text file with the correct format
def test_rd_basic():
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    assert n_moments == 44
    assert n_lambda == 1
    assert external_big_number.dtype == "float"
    assert moment_vector.shape == (44, )
    assert lambda_range.shape == (2, )

# Exceptions


# Nonexistent file
# Should raise an exception
def test_rd_nonexistent():
    assert ~os.path.exists("nonexistent-file")
    try:
        read_data("nonexistent-file")
    except RuntimeError:
        pass
    else:
        raise AssertionError

# File cannot be opened for reading
# NOT YET TESTED


# File with data in the wrong format
# Should raise an exception
def test_rd_badfile():
    try:
        read_data("testing/out.txt")
    except RuntimeError:
        pass
    else:
        raise AssertionError


# n_moments does not match the actual number of moments
# Should reset n_moments and issue a warning
def test_rd_badnmoments():
    with pytest.warns(UserWarning, match="n_moments reset"):
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("testing/badin1.txt")
    assert n_moments == 44


# n_lambda does not match the actual number of lambdas
# reset n_lambda
def test_rd_badnlambda():
    with pytest.warns(UserWarning, match="n_lambda reset"):
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("testing/badin2.txt")
    assert n_lambda == 1


# (calculated) n_moments is not a valid value
# Should raise an exception
def test_rd_badmoments():
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("testing/badin3.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError


# (calculated) n_lambda is not a valid value (currently, anything but 1)
# Should raise an exception
@pytest.mark.skip(reason="not yet implemented")
def test_rd_badlambda():
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
           lambda_range = read_data("testing/badin4.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError


# external_big_number is not a valid value
# Should raise an exception
def test_rd_badbignum():
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("testing/badin5.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError
