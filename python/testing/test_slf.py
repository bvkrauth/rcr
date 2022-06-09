"""
TEST_SLF.PY Unit tests for start_logfile()
"""
import os
import tempfile


import pytest

from rcr import get_logfile, set_logfile, start_logfile


# Basic functionality
def test_slf_basic():
    """set logfile name, open for writing, write the first line"""
    oldlogfile = get_logfile()
    with tempfile.TemporaryDirectory() as tmp:
        logfile = os.path.join(tmp, 'log.txt')
        assert ~os.path.exists(logfile)
        start_logfile(logfile)
        assert os.path.exists(logfile)
    set_logfile(oldlogfile)


def test_slf_none():
    """do nothing if logfile is set to None"""
    oldlogfile = get_logfile()
    start_logfile(None)
    set_logfile(oldlogfile)


# Exceptions
def test_slf_readonly():
    """warn and continue if read-only file"""
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("testing/read-only-file.txt")
    set_logfile(oldlogfile)


def test_slf_badfolder():
    """warn and continue if folder does not exist"""
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("nonexistent-folder/log.txt")
    set_logfile(oldlogfile)


def test_slf_badfilename():
    """warn and continue if illegal file name"""
    oldlogfile = get_logfile()
    with pytest.warns(UserWarning, match="Cannot write to logfile"):
        start_logfile("?/:")
    set_logfile(oldlogfile)


def test_slf_notastring():
    """do nothing if non-string file name"""
    oldlogfile = get_logfile()
    start_logfile(1.0)
    start_logfile(True)
    set_logfile(oldlogfile)
