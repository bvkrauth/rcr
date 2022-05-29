"""
TEST_SM.PY: Unit tests for simplify_moments()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import simplify_moments, read_data, set_logfile, get_logfile

tol = 1e-04

# SIMPLIFY_MOMENTS


# Test with simple data
def test_sm_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    sm1 = simplify_moments(mv1)
    assert sum(abs(sm1 - np.array([1, 1, .5, .25, .25, .25]))) < tol
    assert sm1[5] == (np.sign(sm1[5]) * np.sqrt(sm1[3] * sm1[4]))
    set_logfile(oldlogfile)


# Test with real data
def test_sm_realdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testin1.txt")
    sm = np.array([5.42538313e+02, 2.05839484e-01,
                   1.07467966e+00, 4.47643916e+01,
                   1.34931719e-03, 1.10235301e-02])
    assert max(abs(simplify_moments(moment_vector) - sm)) < tol
    assert sm[5] <= (np.sign(sm[5]) * np.sqrt(sm[3] * sm[4]))
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


def test_sm_randomdata():
    oldlogfile = get_logfile()
    set_logfile(None)
    # Test with randomly-generated data
    for i in range(10):
        mvi, true_smi = random_mv()
        smi = simplify_moments(mvi)
        assert max(abs(smi - true_smi)) < tol
        assert smi[5] == (np.sign(smi[5]) * np.sqrt(smi[3] * smi[4]))
    set_logfile(oldlogfile)

# Special cases


# Invalid length
# This should throw an exception
def test_sm_wronglen():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        mv0 = np.zeros(8)
        simplify_moments(mv0)
    except AssertionError:
        pass
    else:
        raise AssertionError("simplify_moments has accepted invalid input")
    set_logfile(oldlogfile)


# All zeros
# This should return a mix of zero and NaN
def test_sm_allzeros():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv0 = np.zeros(9)
    sm0 = simplify_moments(mv0)
    assert all(sm0[0:3] == 0.0) and all(np.isnan(sm0[3:]))
    set_logfile(oldlogfile)


# Integer data
def test_sm_integer():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv1 = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    sm1 = simplify_moments(mv1)
    assert sum(abs(sm1 - np.array([2, 2, 1, .5, .5, .5]))) < tol
    assert sm1[5] == (np.sign(sm1[5]) * np.sqrt(sm1[3] * sm1[4]))
    set_logfile(oldlogfile)


# var(x) = 0
# Should return NaN
def test_sm_varx0():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv0 = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    sm0 = simplify_moments(mv0)
    assert all(np.isnan(sm0[3:]))
    set_logfile(oldlogfile)


# var(x) < 0
# Should return numeric values
def test_sm_varxneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv0 = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    sm0 = simplify_moments(mv0)
    assert all(np.isfinite(sm0))
    set_logfile(oldlogfile)


# var(y) <= 0
# Should return numeric values
def test_sm_varyneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    sm0 = simplify_moments(mv0)
    assert all(np.isfinite(sm0))
    mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    sm0 = simplify_moments(mv0)
    assert all(np.isfinite(sm0))
    set_logfile(oldlogfile)


# var(z) <= 0
# Should return numeric values
def test_sm_varzneg():
    oldlogfile = get_logfile()
    set_logfile(None)
    mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    sm0 = simplify_moments(mv0)
    assert all(np.isfinite(sm0))
    mv0 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    sm0 = simplify_moments(mv0)
    assert all(np.isfinite(sm0))
    set_logfile(oldlogfile)
