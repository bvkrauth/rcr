"""
TEST_CM.PY: Unit tests for check_moments()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import check_moments, read_data


# Test with simple data
def test_cm_basic():
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    assert check_moments(mv1) == (True, True)


# Test with real data
def test_cm_realdata():
    (n_moments, n_lambda, external_big_number, moment_vector,
        lambda_range) = read_data("testing/testin1.txt")
    assert check_moments(moment_vector) == (True, True)


def random_mv():
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
    return mvi, smi


def test_cm_randomdata():
    # Test with randomly-generated data
    for i in range(10):
        mvi, true_smi = random_mv()
        result = check_moments(mvi)
        assert result == (True, True)


# Special cases
# Invalid length
# This should throw an exception
def test_cm_wronglen():
    try:
        check_moments(np.zeros(8))
    except AssertionError:
        pass
    else:
        raise AssertionError("check_moments has accepted invalid input")


# All zeros
# This should return a mix of zero and NaN
def test_cm_allzeros():
    mv = np.zeros(9)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(mv)
    assert result == (False, False)


# Integer data
def test_cm_integer():
    mv = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    result = check_moments(mv)
    assert result == (True, True)


# var(x) = 0
# Should return NaN
def test_cm_varx0():
    mv = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(mv)
    assert result == (False, False)


# var(x) < 0
# Should return numeric values
def test_cm_varxneg():
    mv = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)


# var(y) <= 0
# Should return numeric values
def test_cm_varyneg():
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)


# var(z) <= 0
# Should return numeric values
def test_cm_varzneg():
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
