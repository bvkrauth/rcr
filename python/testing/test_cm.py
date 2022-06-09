"""
TEST_CM.PY: Unit tests for check_moments()
"""
import numpy as np
import pytest

from rcr import check_moments, read_data


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
    e_xyz = np.random.randn(3)
    v_xyz = abs(np.random.randn(3))
    cor_xyz = np.random.uniform(-1, 1, 3)
    e_x = e_xyz[0]
    e_y = e_xyz[1]
    e_z = e_xyz[2]
    v_x = v_xyz[0]
    v_y = v_xyz[1]
    v_z = v_xyz[2]
    cov_xy = cor_xyz[0] * np.sqrt(v_x * v_y)
    cov_xz = cor_xyz[1] * np.sqrt(v_x * v_z)
    cov_yz = cor_xyz[2] * np.sqrt(v_y * v_z)
    e_xx = v_x + e_x ** 2
    e_yy = v_y + e_y ** 2
    e_zz = v_z + e_z ** 2
    e_xy = cov_xy + e_x * e_y
    e_xz = cov_xz + e_x * e_z
    e_yz = cov_yz + e_y * e_z
    v_yhat = (cov_xy ** 2) / v_x
    v_zhat = (cov_xz ** 2) / v_x
    cov_hat = cov_xy * cov_xz / v_x
    moment_vector = np.array([e_x, e_y, e_z,
                              e_xx, e_xy, e_xz,
                              e_yy, e_yz, e_zz])
    simplified_moments = np.array([v_y, v_z, cov_yz,
                                   v_yhat, v_zhat, cov_hat])
    return moment_vector, simplified_moments


def test_cm_randomdata():
    """ check moments with randomly-generated data"""
    for _ in range(10):
        moment_vector = random_mv()[0]
        result = check_moments(moment_vector)
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
    moment_vector = np.zeros(9)
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_integer():
    """check moments with integer data"""
    moment_vector = np.array([0, 0, 0, 2, 1, 1, 2, 1, 2])
    result = check_moments(moment_vector)
    assert result == (True, True)


def test_cm_varx0():
    """return NaN if var(x) = 0"""
    moment_vector = np.array([0, 0, 0, 0, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data: nonsingular"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varxneg():
    """return numeric values if var(x) < 0"""
    moment_vector = np.array([0, 0, 0, -1, 0.5, 0.5, 1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varyneg():
    """return numeric values if var(y) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, -1, 0.5, 1.0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)


def test_cm_varzneg():
    """return numeric values if var(z) <= 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 0])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    with pytest.warns(UserWarning, match="Invalid data:"):
        result = check_moments(moment_vector)
    assert result == (False, False)
