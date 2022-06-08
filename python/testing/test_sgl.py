"""
TEST_SGL.PY Unit tests for set_logfile() and get_logfile()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import get_logfile, \
    set_logfile  # pylint: disable=wrong-import-position


# set_logfile(str) sets the global variable logfile to str
# get_logfile() retrieves the value of the global variable logfile
def test_sgl():
    """get and set the global variable logfile"""
    tmp = get_logfile()
    set_logfile("any string")
    assert get_logfile() == "any string"
    set_logfile("another string")
    assert get_logfile() == "another string"
    set_logfile(tmp)


def test_sgl_none():
    """set_logfile accepts strings or None"""
    tmp = get_logfile()
    set_logfile(None)
    assert get_logfile() is None
    set_logfile(tmp)


# set_logfile should ignore all other inputs
def test_sgl_bool():
    """ignore logical arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(True)
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_int():
    """ignore integer arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(0)
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_real():
    """ignore real arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(0.)
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_list():
    """ignore list arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(["another string"])
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_tuple():
    """ignore tuple arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(("another string", "a third string"))
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_array():
    """ignore array arguments"""
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(np.zeros(2))
    assert get_logfile() == "a string"
    set_logfile(tmp)


def test_sgl_undef():
    """initial value is None"""
    assert get_logfile() is None
