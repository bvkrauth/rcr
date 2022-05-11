"""
TEST_DW.PY Unit tests for die() and warn()
"""
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import warn, die, set_logfile, get_logfile


# Basic functionality of warn()
# Write to the logfile and issue a warning
def test_warn():
    oldlogfile = get_logfile()
    set_logfile(None)
    with pytest.warns(UserWarning, match="Test warning"):
        warn("Test warning")
    set_logfile(oldlogfile)


# Basic functionality of die()
# Issue an exception
def test_die():
    oldlogfile = get_logfile()
    set_logfile(None)
    try:
        die("Test dying")
    except RuntimeError:
        pass
    else:
        raise AssertionError
    set_logfile(oldlogfile)
