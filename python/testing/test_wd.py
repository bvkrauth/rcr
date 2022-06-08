"""
TEST_WD.PY Unit tests for write_details()
"""
import os
import tempfile
import sys

import pytest
import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import write_details


# Basic fnctionality
# Write the specified arrays to the specified text file
def test_wd_basic():
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with tempfile.TemporaryDirectory() as tmp:
        detfile = os.path.join(tmp, 'pdet.txt')
        assert ~os.path.exists(detfile)
        write_details(thetavec, lambdavec, detfile)
        assert os.path.exists(detfile)


# Exceptions to handle


# File name is blank ("")
# Should do nothing
def test_wd_nofile():
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    write_details(thetavec, lambdavec, "")


# Read-only file
# Should issue a warning and continue
def test_wd_readonly():
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "testing/read-only-file.txt")


# Nonexistent folder name
# Should issue a warning and continue
def test_wd_badfolder():
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "nonexistent-path-name/pout.txt")


# Illegal file name
# Should issue a warning and continue
def test_wd_illegalname():
    thetavec = np.zeros(3)
    lambdavec = np.zeros(3)
    with pytest.warns(UserWarning, match="Cannot write"):
        write_details(thetavec, lambdavec, "?")
