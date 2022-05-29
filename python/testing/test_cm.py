"""
TEST_CM.PY: Unit tests for check_moments()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import check_moments, read_data, set_logfile, get_logfile

tol = 1e-04


# Test with simple data
def test_cm_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    assert check_moments(mv1) == (True, True)
    set_logfile(oldlogfile)


# Test with real data
def test_cm_realdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    assert check_moments(moment_vector) == (True, True)
    set_logfile(oldlogfile)


def random_mv():
    oldlogfile = get_logfile()
    set_logfile(None)
    exyz = np.random.randn(3)
    vxyz = abs(np.random.randn(3))
    corxyz = np.random.uniform(-1, 1, 3)
    ex = exyz[0]
    ey = exyz[1]
    ez = exyz[2]
    vx = vxyz[0]
    vy = vxyz[1]
    vz = vxyz[2]
    cxy = corxyz[0] * np.sqrt(vx * vy)
    cxz = corxyz[1] * np.sqrt(vx * vz)
    cyz = corxyz[2] * np.sqrt(vy * vz)
    exx = vx + ex ** 2
    eyy = vy + ey ** 2
    ezz = vz + ez ** 2
    exy = cxy + ex * ey
    exz = cxz + ex * ez
    eyz = cyz + ey * ez
    vyhat = (cxy ** 2) / vx
    vzhat = (cxz ** 2) / vx
    covhat = cxy * cxz / vx
    mvi = np.array([ex, ey, ez, exx, exy, exz, eyy, eyz, ezz])
    smi = np.array([vy, vz, cyz, vyhat, vzhat, covhat])
    set_logfile(oldlogfile)
    return mvi, smi


def test_cm_randomdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    # Test with randomly-generated data
    for i in range(10):
        mvi, true_smi = random_mv()
        assert check_moments(mvi) == (True, True)
    set_logfile(oldlogfile)


# Special cases
# Invalid length
# This should throw an exception
def test_cm_wronglen():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        check_moments(np.zeros(8))
    except AssertionError:
        pass
    else:
        raise AssertionError("check_moments has accepted invalid input")
    set_logfile(oldlogfile)


# All zeros
# This should return a mix of zero and NaN
def test_cm_allzeros():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        assert check_moments(np.zeros(9)) == (False, False)
    set_logfile(oldlogfile)


# Integer data
def test_cm_integer():
    oldlogfile = get_logfile()
    set_logfile(None)
    assert check_moments(np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])) == (True, True)
    set_logfile(oldlogfile)


# var(x) = 0
# Should return NaN
def test_cm_varx0():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        assert check_moments(mv) == (False, False)
    set_logfile(oldlogfile)


# var(x) < 0
# Should return numeric values
def test_cm_varxneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        assert check_moments(mv) == (False, False)
    set_logfile(oldlogfile)


# var(y) <= 0
# Should return numeric values
def test_cm_varyneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="Invalid data:"):
        mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
        assert(check_moments(mv0)) == (False, False)
    with pytest.warns(UserWarning, match="Invalid data:"):
        mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
        assert(check_moments(mv0)) == (False, False)
    set_logfile(oldlogfile)


# var(z) <= 0
# Should return numeric values
def test_cm_varzneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="Invalid data:"):
        mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
        assert(check_moments(mv0)) == (False, False)
    with pytest.warns(UserWarning, match="Invalid data:"):
        mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
        assert(check_moments(mv0)) == (False, False)
    set_logfile(oldlogfile)
