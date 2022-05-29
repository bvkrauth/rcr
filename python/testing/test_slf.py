"""
TEST_SLF.PY Unit tests for start_logfile()
"""
import os
import tempfile
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import get_logfile, set_logfile, start_logfile


# Basic functionality: set logfile name to its argument,
# open it for writing, and write the first line
def test_slf_basic():
    oldlogfile = get_logfile()
    with tempfile.TemporaryDirectory() as tmp:
        logfile = os.path.join(tmp, 'log.txt')
        assert ~os.path.exists(logfile)
        start_logfile(logfile)
        assert os.path.exists(logfile)
    set_logfile(oldlogfile)


# logfile set to None
# Do nothing
def test_slf_none():
    oldlogfile = get_logfile()
    start_logfile(None)
    set_logfile(oldlogfile)


# Exceptions

# Read-only file
# Issue a warning but continue
def test_slf_readonly():
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("read-only-file.txt")
    set_logfile(oldlogfile)


# Nonexistent folder name
# Issue a warning but continue
def test_slf_badfolder():
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("nonexistent-folder/log.txt")
    set_logfile(oldlogfile)


# Illegal file name
# Issue a warning but continue
def test_slf_badfilename():
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("?/:")
    set_logfile(oldlogfile)


# Non-string argument
# Do nothing
def test_slf_notastring():
    oldlogfile = get_logfile()
    start_logfile(1.0)
    start_logfile(True)
    set_logfile(oldlogfile)
