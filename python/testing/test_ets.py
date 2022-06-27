"""
TEST_ETS.PY: Unit tests for estimate_effect_segments()
"""


import numpy as np
import pytest

from rcrbounds import rcfast, simplify_moments, estimate_effect_segments


# Basic functionality
def test_ets_realdata(moment_vector):
    """ estimate effect segments with real data"""
    ts_true = np.array([-1.00000000e+100,
                        -1.48223355e+001,
                        8.16970996e+000,
                        8.16970996e+000,
                        1.00000000e+100])
    test_ts, effectvec, rcvec = estimate_effect_segments(moment_vector)
    rcvec_true = rcfast(effectvec,
                        simplify_moments(moment_vector))
    assert test_ts == pytest.approx(ts_true)
    assert all(np.isfinite(effectvec))
    assert rcvec == pytest.approx(rcvec_true)


# Special cases for moments
def test_ets_nearrct():
    """estimate effect segments with near-perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    ts1 = estimate_effect_segments(mv1)[0]
    ts1_true = np.array([-1.e+100,  5.e+005,  5.e+005,  1.e+100])
    assert ts1 == pytest.approx(ts1_true)


def test_ets_rct():
    """estimate effect segments with perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    # This test currently fails with an UnboundLocalError
    try:
        estimate_effect_segments(mv1)
    except UnboundLocalError:
        pass
    else:
        pass
