"""
TEST_PB.PY: Unit tests for pole_bounds()
"""
import numpy as np

from rcrbounds import (pole_bounds, simplify_moments, effectinf,
                       estimate_effect_segments, rcfast)


def get_bracket(moment_vector):
    """get effect_segments and rc at the values bracketing effect_inf"""
    simplified_moments = simplify_moments(moment_vector)
    effect_inf = effectinf(moment_vector)
    effect_segments = estimate_effect_segments(moment_vector)[0]
    bracket = np.array([np.max(effect_segments[effect_segments < effect_inf]),
                        np.min(effect_segments[effect_segments > effect_inf])])
    rc_bracket = rcfast(bracket, simplified_moments)
    return simplified_moments, effect_inf, effect_segments, rc_bracket


def test_pb_basic(moment_vector):
    """identify bounds that are next to effect_inf"""
    simplified_moments, effect_inf, effect_segments, rc_bracket = \
        get_bracket(moment_vector)
    below = np.array([effect_inf - 10.0, effect_inf - 1e-8])
    above = np.array([effect_inf + 1e-8, effect_inf + 10.0])

    def pb(bounds, rc_range):
        return pole_bounds(bounds, np.array(rc_range), effect_inf,
                           effect_segments, simplified_moments)
    # rc_range includes rc at the lower bracketing value
    range_below = (min(rc_bracket[0], 0.0), max(rc_bracket[0], 0.0))
    assert pb(below, range_below) == [False, True]
    assert pb(above, range_below) == [False, False]
    # rc_range includes rc at the upper bracketing value
    range_above = (min(rc_bracket[1], 0.0), max(rc_bracket[1], 0.0))
    assert pb(above, range_above) == [True, False]
    assert pb(below, range_above) == [False, False]
    # Unbounded rc_range on the same side as rc at the bracketing value
    unbounded = (-np.inf, 0.0) if rc_bracket[0] < 0 else (0.0, np.inf)
    assert pb(below, unbounded) == [False, True]
    # Finite rc_range that stops short of rc at the bracketing value
    short = tuple(np.sort([0.0, rc_bracket[0] / 1e3]))
    assert pb(below, short) == [False, False]


def test_pb_nonfinite(moment_vector):
    """never identify infinite bounds or non-finite effect_inf"""
    simplified_moments, effect_inf, effect_segments, _ = \
        get_bracket(moment_vector)
    rc_range = np.array([-np.inf, np.inf])
    assert pole_bounds(np.array([-np.inf, np.inf]), rc_range, effect_inf,
                       effect_segments, simplified_moments) == [False, False]
    assert pole_bounds(np.array([0.0, 1.0]), rc_range, np.nan,
                       effect_segments, simplified_moments) == [False, False]
