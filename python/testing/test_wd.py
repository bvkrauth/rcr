"""
TEST_WD.PY Unit tests for write_details()
"""
import os
import tempfile


import pytest
import numpy as np

from rcr import write_details


# Basic functionality
def test_wd_basic():
    """write the specified arrays to the specified test file"""
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with tempfile.TemporaryDirectory() as tmp:
        detfile = os.path.join(tmp, 'pdet.txt')
        assert ~os.path.exists(detfile)
        write_details(thetavec, lambdavec, detfile)
        assert os.path.exists(detfile)


# Exceptions to handle
def test_wd_nofile():
    """do nothing if file name is blank"""
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    write_details(thetavec, lambdavec, "")


def test_wd_readonly():
    """warn and continue if read-only file"""
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "testing/read-only-file.txt")


def test_wd_badfolder():
    """warn and continue if non-existend folder"""
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "nonexistent-path-name/pout.txt")


def test_wd_illegalname():
    """warn and continue if illegal file name"""
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "?")
