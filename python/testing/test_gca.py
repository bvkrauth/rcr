"""
TEST_GCA.PY Unit tests for get_command_arguments()
"""
import sys

import pytest
import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import get_command_arguments  # pylint: disable=wrong-import-position

# get_command_arguments() takes a list of 0+ strings
# and returns a list of 4 strings.  Any strings not
# provided are replaced with default values, extra
# strings are discarded.


# Base case
def test_gca_noargs():
    """use defaults if no command arguments"""
    args = ["program name"]
    assert get_command_arguments(args) == ("in.txt",
                                           "pout.txt",
                                           "plog.txt",
                                           "")


def test_gca_infile():
    """change infile if one argument"""
    args = ["program name", "alt_infile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "pout.txt",
                                           "plog.txt",
                                           "")


def test_gca_outfile():
    """change infile, outfile if two arguments"""
    args = ["program name", "alt_infile", "alt_outfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "plog.txt",
                                           "")


def test_gca_logfile():
    """change infile, outfile, logfile if three arguments"""
    args = ["program name", "alt_infile", "alt_outfile", "alt_logfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "alt_logfile",
                                           "")


def test_gca_detfile():
    """change ...detail_file if four arguments"""
    args = ["program name", "alt_infile", "alt_outfile",
            "alt_logfile", "alt_detailfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "alt_logfile",
                                           "alt_detailfile")


def test_gca_extra():
    """issue a warning if five+ arguments"""
    args = ["program name", "alt_infile", "alt_outfile",
            "alt_logfile", "alt_detailfile", "EXTRA JUNK"]
    with pytest.warns(UserWarning, match="Unused program arguments"):
        tmp = get_command_arguments(args)
    assert tmp == ("alt_infile", "alt_outfile",
                   "alt_logfile", "alt_detailfile")


# Blank/whitespace arguments
# This will be an invalid file name so maybe we should issue a warning
def test_gca_blank():
    """return a blank string if blank/whitespace arguments"""
    args = ["program name", "", "    "]
    assert get_command_arguments(args) == ("",
                                           "",
                                           "plog.txt",
                                           "")


def test_gca_nonarray():
    """warn and return defaults if non-array arguments"""
    msg = "Invalid command arguments, using defaults"
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(None)
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments("string")
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(True)
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(1.0)
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")


def test_gca_nonstring():
    """warn and return defaults if non-string arguments"""
    msg = "Invalid command arguments, using defaults"
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(np.array([True, False]))
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(np.zeros(2))
    assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
