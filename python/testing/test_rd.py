"""
TEST_RD.PY Unit tests for read_data()
"""
import os


import pytest

from rcrbounds import read_data


# Basic functionality
def test_rd_basic(infile):
    """read data from the specified text file"""
    n_moments, n_rc, external_big_number, moment_vector, \
        rc_range = read_data(infile)
    assert n_moments == 44
    assert n_rc == 1
    assert external_big_number == pytest.approx(8.98846567e+306)
    assert moment_vector.shape == (44, )
    assert rc_range.shape == (2, )

# Exceptions


def test_rd_nonexistent():
    """raise an exception if file does not exist"""
    assert not os.path.exists("nonexistent-file")
    try:
        read_data("nonexistent-file")
    except RuntimeError:
        pass
    else:
        raise AssertionError


def test_rd_notastring():
    """raise exception if filename is not a string"""
    try:
        read_data([1, 10])
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


def test_rd_badnrc(badin2):
    """reset n_lmambda and warn if it does not match rc_range"""
    with pytest.warns(UserWarning, match="n_rc reset"):
        n_rc = read_data(badin2)[1]
    assert n_rc == 1


def test_rd_badmoments(badin3):
    """raise exception if n_moments is invalid"""
    try:
        read_data(badin3)
    except AssertionError:
        pass
    else:
        raise AssertionError


@pytest.mark.skip(reason="not yet implemented")
def test_rd_badrc(badin4):
    """raise exception if n_rc is invalid"""
    try:
        read_data(badin4)
    except AssertionError:
        pass
    else:
        raise AssertionError


def test_rd_badbignum1(badin5):
    """raise exception if extarnal_big_number is invalid"""
    try:
        read_data(badin5)
    except AssertionError:
        pass
    else:
        raise AssertionError


def test_rd_badbignum2(badin6):
    """raise warning if extarnal_big_number is too big"""
    with pytest.warns(UserWarning, match="Largest Python real"):
        read_data(badin6)
