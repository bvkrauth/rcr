"""
TEST_SGL.PY Unit tests for set_logfile() and get_logfile()
"""
import sys

import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import get_logfile, set_logfile


# set_logfile(str) sets the global variable logfile to str
# get_logfile() retrieves the value of the global variable logfile
def test_sgl():
    tmp = get_logfile()
    set_logfile("any string")
    assert get_logfile() == "any string"
    set_logfile("another string")
    assert get_logfile() == "another string"
    set_logfile(tmp)


# set_logfile should accept a string or None
def test_sgl_none():
    tmp = get_logfile()
    set_logfile(None)
    assert get_logfile() is None
    set_logfile(tmp)


# set_logfile should ignore all other inputs
# Boolean
# Should leave logfile unchanged
def test_sgl_bool():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(True)
    assert get_logfile() == "a string"
    set_logfile(tmp)


# Integer
# Should leave logfile unchanged
def test_sgl_int():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(0)
    assert get_logfile() == "a string"
    set_logfile(tmp)


# Real
# Should leave logfile unchanged
def test_sgl_real():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(0.)
    assert get_logfile() == "a string"
    set_logfile(tmp)


# List
# Should leave logfile unchanged
def test_sgl_list():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(["another string"])
    assert get_logfile() == "a string"
    set_logfile(tmp)


# Tuple
# Should leave logfile unchanged
def test_sgl_tuple():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(("another string", "a third string"))
    assert get_logfile() == "a string"
    set_logfile(tmp)


# Array
# Should leave logfile unchanged
def test_sgl_array():
    tmp = get_logfile()
    set_logfile("a string")
    set_logfile(np.zeros(2))
    assert get_logfile() == "a string"
    set_logfile(tmp)


# When set_logfile() has not been called yet, logfile is undefined
# get_logfile() should return None
def test_sgl_undef():
    assert get_logfile() is None
