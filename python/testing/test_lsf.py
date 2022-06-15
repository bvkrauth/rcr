"""
TEST_LSF.PY: Unit tests for lambdastar()
"""


import numpy as np
import pytest

from rcrbounds import lambdastar


# Basic functionality
def test_true_ls_basic():
    """lambdastar from simple data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, np.sqrt(0.2), 1, 0.5, 1.0])
    true_ls = 2.0
    test_ls = lambdastar(mv1)
    assert test_ls == pytest.approx(true_ls)


def test_true_ls_realdata(moment_vector):
    """lambdastar from real data"""
    true_ls = 12.310599093115798
    test_ls = lambdastar(moment_vector)
    assert test_ls == pytest.approx(true_ls)


# Special cases
def test_true_ls_zero():
    """return inf if var(zhat) = 0"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    true_ls = np.inf
    test_ls = lambdastar(mv1)
    assert test_ls == pytest.approx(true_ls)


def test_true_ls_nearzero():
    """work normally if var(zhat) near 0"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    test_ls = lambdastar(mv1)
    assert np.isfinite(test_ls)


def test_true_ls_negvarz():
    """return zero if var(zhat) > var(z)"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, -1])
    test_ls = lambdastar(mv1)
    assert test_ls == 0.0
