"""
TEST_WTL.PY Unit tests for write_to_logfile()
"""
import os
import tempfile

import pytest

from rcrbounds import write_to_logfile, set_logfile, get_logfile


# Basic functionality
def test_wtl_basic():
    """create the log file and write to it"""
    oldlogfile = get_logfile()
    with tempfile.TemporaryDirectory() as tmp:
        logfile = os.path.join(tmp, 'log.txt')
        set_logfile(logfile)
        assert not os.path.exists(logfile)
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
        assert os.path.exists(logfile)
        # would be nice if we could check the contents of the file too
    set_logfile(oldlogfile)


def test_wtl_none():
    """do nothing if logfile = None"""
    oldlogfile = get_logfile()
    set_logfile(None)
    write_to_logfile("Line 1\n", mode="w")
    write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


# Exceptions to handle
def test_wtl_readonly(read_only_file):
    """warn and continue if file is read-only"""
    oldlogfile = get_logfile()
    set_logfile(read_only_file)
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


def test_wtl_badfolder():
    """warn and continue if folder does not exist"""
    oldlogfile = get_logfile()
    set_logfile("nonexistent-folder/log.txt")
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


def test_wtl_badfilename():
    """warn and continue if file name is illegal"""
    oldlogfile = get_logfile()
    set_logfile("?/:")
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)
