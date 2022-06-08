"""
TEST_CM.PY: Unit tests for check_moments()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import check_moments, \
    read_data  # pylint: disable=wrong-import-position


def test_cm_basic():
    """check moments with simple valid data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    assert check_moments(mv1) == (True, True)


# Test with real data
def test_cm_realdata():
    """check moments from real data"""
    moment_vector = read_data("testing/testin1.txt")[3]
    assert check_moments(moment_vector) == (True, True)


def random_mv():
    """generate random moment_vector"""
    # pylint: disable=too-many-locals
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
    """ check moments with randomly-generated data"""
    for _ in range(10):
        mvi = random_mv()[0]
        result = check_moments(mvi)
        assert result == (True, True)


# Special cases
def test_cm_wronglen():
    """raise an exception if invalid length"""
    try:
        check_moments(np.zeros(8))
    except AssertionError:
        pass
    else:
        raise AssertionError("check_moments has accepted invalid input")


def test_cm_allzeros():
    """return mics of zero and NaN if all zeros"""
    mv = np.zeros(9)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(mv)
    assert result == (False, False)


def test_cm_integer():
    """check moments with integer data"""
    mv = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    result = check_moments(mv)
    assert result == (True, True)


def test_cm_varx0():
    """return NaN if var(x) = 0"""
    mv = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(mv)
    assert result == (False, False)


def test_cm_varxneg():
    """return numeric values if var(x) < 0"""
    mv = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)


def test_cm_varyneg():
    """return numeric values if var(y) <= 0"""
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)


def test_cm_varzneg():
    """return numeric values if var(z) <= 0"""
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
    mv = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(mv)
    assert result == (False, False)
