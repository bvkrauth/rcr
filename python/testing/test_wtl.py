"""
TEST_WTL.PY Unit tests for write_to_logfile()
"""
import os
import tempfile
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import write_to_logfile, set_logfile, get_logfile


# Basic functionality - create the log file (if mode = "w") and write to it
def test_wtl_basic():
    oldlogfile = get_logfile()
    with tempfile.TemporaryDirectory() as tmp:
        logfile = os.path.join(tmp, 'log.txt')
        set_logfile(logfile)
        assert ~os.path.exists(logfile)
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
        assert os.path.exists(logfile)
        # would be nice if we could check the contents of the file too
    set_logfile(oldlogfile)


# If the logfile is None, it should do nothing
def test_wtl_none():
    oldlogfile = get_logfile()
    set_logfile(None)
    write_to_logfile("Line 1\n", mode="w")
    write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


# Exceptions to handle
# Read-only file
# Issue a warning but continue
def test_wtl_readonly():
    oldlogfile = get_logfile()
    set_logfile("testing/read-only-file.txt")
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


# Nonexistent folder name
# Issue a warning but continue
def test_wtl_badfolder():
    oldlogfile = get_logfile()
    set_logfile("nonexistent-folder/log.txt")
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)


# Illegal file name
# Issue a warning but continue
def test_wtl_badfilename():
    oldlogfile = get_logfile()
    set_logfile("?/:")
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        write_to_logfile("Line 1\n", mode="w")
        write_to_logfile("Line 2\n")
    set_logfile(oldlogfile)
