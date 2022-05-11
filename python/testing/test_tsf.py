"""
TEST_TSF.PY: Unit tests for thetastar()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import thetastar, read_data, set_logfile, get_logfile

tol = 1e-04


# Test with simple data
def test_tsf_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    tsf = 1.0
    assert abs(thetastar(mv1) - tsf) < tol
    set_logfile(oldlogfile)


# Test with real data
def test_tsf_realdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_theta, external_big_number, moment_vector,
        theta_range) = read_data("testin1.txt")
    tsf = 8.169709964904111
    assert abs(thetastar(moment_vector) - tsf) < tol
    set_logfile(oldlogfile)


# Special cases


# var(zhat) = 0 - should return infinity
def test_tsf_zero():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    assert np.isnan(thetastar(mv1))
    set_logfile(oldlogfile)


# var(zhat) near 0 - should work normally
def test_tsf_nearzero():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    assert np.isfinite(thetastar(mv1))
    set_logfile(oldlogfile)
