"""
TEST_ADT.PY: Unit tests for add_detail_tails()
"""
import numpy as np

from rcrbounds import estimate_model, add_detail_tails, rcfun


# Basic functionality
def test_adt_basic(moment_vector):
    """add tail points beyond +/-50 with correct rc values"""
    rc_range = np.array([0.0, 1.0])
    effectvec, rcvec = estimate_model(moment_vector, rc_range)[1:]
    new_effectvec, new_rcvec = add_detail_tails(effectvec,
                                                rcvec,
                                                moment_vector)
    # Original points are all kept, and the result is sorted
    assert np.all(np.isin(effectvec, new_effectvec))
    assert np.all(np.diff(new_effectvec) > 0)
    assert len(new_effectvec) == len(new_rcvec)
    # Same end points as before
    assert new_effectvec[0] == effectvec[0]
    assert new_effectvec[-1] == effectvec[-1]
    # Gaps between 50 and 1e6 (in absolute value) are at most about 1%
    for sign in (1.0, -1.0):
        mid = np.sort(sign * new_effectvec[(sign * new_effectvec >= 50.0) &
                                           (sign * new_effectvec <= 1.0e6)])
        assert mid[0] == 50.0
        assert mid[-1] == 1.0e6
        assert np.max(mid[1:] / mid[:-1]) < 1.0101
    # rc values at the new points match rcfun
    tails = np.abs(new_effectvec) > 50.0
    assert np.allclose(new_rcvec[tails],
                       rcfun(moment_vector, new_effectvec[tails]),
                       equal_nan=True)
    # rc values at the original points are unchanged
    keep = np.isin(new_effectvec, effectvec)
    assert np.array_equal(new_rcvec[keep], rcvec, equal_nan=True)


def test_adt_small_effectmax(moment_vector):
    """add nothing if the grid does not extend beyond +/-50"""
    effectvec = np.linspace(-50.0, 50.0, 11)
    rcvec = rcfun(moment_vector, effectvec)
    new_effectvec, new_rcvec = add_detail_tails(effectvec,
                                                rcvec,
                                                moment_vector)
    assert np.array_equal(new_effectvec, effectvec)
    assert np.array_equal(new_rcvec, rcvec, equal_nan=True)
