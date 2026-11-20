"""
TEST_PB.PY: Unit tests for pole_bounds()
"""
import numpy as np

from rcrbounds import pole_bounds, simplify_moments, effectinf


def test_pb_basic(moment_vector):
    """identify bounds that are next to effect_inf"""
    simplified_moments = simplify_moments(moment_vector)
    effect_inf = effectinf(moment_vector)
    sign_below = np.sign(simplified_moments[2] -
                         simplified_moments[5] *
                         simplified_moments[1] /
                         simplified_moments[4])
    below = np.array([effect_inf - 10.0, effect_inf - 1e-8])
    above = np.array([effect_inf + 1e-8, effect_inf + 10.0])
    unbounded_below = np.array([-np.inf, 0.0])
    unbounded_above = np.array([0.0, np.inf])
    if sign_below < 0:
        range_below, range_above = unbounded_below, unbounded_above
    else:
        range_below, range_above = unbounded_above, unbounded_below
    assert pole_bounds(below, range_below, effect_inf,
                       simplified_moments) == [False, True]
    assert pole_bounds(above, range_above, effect_inf,
                       simplified_moments) == [True, False]
    # Wrong side of effect_inf
    assert pole_bounds(above, range_below, effect_inf,
                       simplified_moments) == [False, False]
    # Bounded rc_range
    assert pole_bounds(below, np.array([-1e6, 0.0]), effect_inf,
                       simplified_moments) == [False, False]


def test_pb_nonfinite(moment_vector):
    """never identify infinite bounds or non-finite effect_inf"""
    simplified_moments = simplify_moments(moment_vector)
    effect_inf = effectinf(moment_vector)
    rc_range = np.array([-np.inf, np.inf])
    assert pole_bounds(np.array([-np.inf, np.inf]), rc_range, effect_inf,
                       simplified_moments) == [False, False]
    assert pole_bounds(np.array([0.0, 1.0]), rc_range, np.nan,
                       simplified_moments) == [False, False]
