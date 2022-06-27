"""
TEST_TSF.PY: Unit tests for effectinf()
"""


import numpy as np
import pytest

from rcrbounds import effectinf


# Basic functionality
def test_true_ts_basic():
    """effectinf from simple data"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    true_ts = 1.0
    effect_inf = effectinf(moment_vector)
    assert effect_inf == pytest.approx(true_ts)


def test_true_ts_realdata(moment_vector):
    """effectinf from real data"""
    true_ts = 8.169709964904111
    effect_inf = effectinf(moment_vector)
    assert effect_inf == pytest.approx(true_ts)


# Special cases
def test_true_ts_zero():
    """return inf if var(zhat) = 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 0, 1, 0.5, 1])
    effect_inf = effectinf(moment_vector)
    assert np.isnan(effect_inf)


def test_true_ts_nearzero():
    """return a finite value if var(zhat) near 0"""
    moment_vector = np.array([0, 0, 0, 1, 0.5, 1e-100, 1, 0.5, 1])
    effect_inf = effectinf(moment_vector)
    assert np.isfinite(effect_inf)
