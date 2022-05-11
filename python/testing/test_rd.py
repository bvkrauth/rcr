"""
TEST_RD.PY Unit tests for read_data()
"""
import os
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import get_logfile, set_logfile, read_data


# Basic functionality
# Should read data from the specified text file with the correct format
def test_rd_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testin1.txt")
    assert n_moments == 44
    assert n_lambda == 1
    assert external_big_number.dtype == "float"
    assert moment_vector.shape == (44, )
    assert lambda_range.shape == (2, )
    set_logfile(oldlogfile)

# Exceptions


# Nonexistent file
# Should raise an exception
def test_rd_nonexistent():
    oldlogfile = get_logfile()
    set_logfile(None)
    assert ~os.path.exists("nonexistent-file")
    try:
       read_data("nonexistent-file")
    except RuntimeError:
        pass
    else:
        raise AssertionError
    set_logfile(oldlogfile)

# File cannot be opened for reading
# NOT YET TESTED


# File with data in the wrong format
# Should raise an exception
def test_rd_badfile():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        read_data("out.txt")
    except RuntimeError:
        pass
    else:
        raise AssertionError
    set_logfile(oldlogfile)


# n_moments does not match the actual number of moments
# Should reset n_moments and issue a warning
def test_rd_badnmoments():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="n_moments reset"):
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("badin1.txt")
    assert n_moments == 44
    set_logfile(oldlogfile)


# n_lambda does not match the actual number of lambdas
# reset n_lambda
def test_rd_badnlambda():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="n_lambda reset"):
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("badin2.txt")
    assert n_lambda == 1
    set_logfile(oldlogfile)


# (calculated) n_moments is not a valid value
# Should raise an exception
def test_rd_badmoments():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("badin3.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError
    set_logfile(oldlogfile)


# (calculated) n_lambda is not a valid value (currently, anything but 1)
# Should raise an exception
@pytest.mark.skip(reason="not yet implemented")
def test_rd_badlambda():
    set_logfile(None)
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
           lambda_range = read_data("badin4.txt")
    except AssertionError:
        pass
    else:
       raise AssertionError


# external_big_number is not a valid value
# Should raise an exception
def test_rd_badbignum():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        n_moments, n_lambda, external_big_number, moment_vector, \
            lambda_range = read_data("badin5.txt")
    except AssertionError:
        pass
    else:
        raise AssertionError
    set_logfile(oldlogfile)
