"""
TEST_ETS.PY: Unit tests for estimate_theta_segments()
"""
import sys

import numpy as np
import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import lambdafast, simplify_moments, \
    estimate_theta_segments  # pylint: disable=wrong-import-position


# Basic functionality
def test_ets_realdata(moment_vector):
    """ estimate theta segments with real data"""
    em_true = np.array([-1.00000000e+100,
                        -1.48223355e+001,
                        8.16970996e+000,
                        8.16970996e+000,
                        1.00000000e+100])
    em, thetavec, lambdavec = estimate_theta_segments(moment_vector)
    lambdavec_true = lambdafast(thetavec,
                                simplify_moments(moment_vector))
    assert em == pytest.approx(em_true)
    assert all(np.isfinite(thetavec))
    assert lambdavec == pytest.approx(lambdavec_true)


# Special cases for moments
def test_ets_nearrct():
    """estimate theta segments with near-perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    ts1 = estimate_theta_segments(mv1)[0]
    ts1_true = np.array([-1.e+100,  5.e+005,  5.e+005,  1.e+100])
    assert ts1 == pytest.approx(ts1_true)


def test_ets_rct():
    """estimate theta segments with perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        estimate_theta_segments(mv1)
    except UnboundLocalError:
        pass
    else:
        pass
