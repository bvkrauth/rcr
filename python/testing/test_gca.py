"""
TEST_GCA.PY Unit tests for get_command_arguments()
"""
import sys

import pytest
import numpy as np

sys.path.append("./")
sys.path.append("../")
from rcr import get_command_arguments

# get_command_arguments() takes a list of 0+ strings
# and returns a list of 4 strings.  Any strings not
# provided are replaced with default values, extra
# strings are discarded.


# Base case: No arguments
def test_gca_noargs():
    args = ["program name"]
    assert get_command_arguments(args) == ("in.txt",
                                           "pout.txt",
                                           "plog.txt",
                                           "")


# One argument (infile)
def test_gca_infile():
    args = ["program name", "alt_infile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "pout.txt",
                                           "plog.txt",
                                           "")


# Two arguments (infile outfile)
def test_gca_outfile():
    args = ["program name", "alt_infile", "alt_outfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "plog.txt",
                                           "")


# Three arguments (infile outfile logfile)
def test_gca_logfile():
    args = ["program name", "alt_infile", "alt_outfile", "alt_logfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "alt_logfile",
                                           "")


# Four arguments (infile outfile logfile detail_file)
def test_gca_detfile():
    args = ["program name", "alt_infile", "alt_outfile",
            "alt_logfile", "alt_detailfile"]
    assert get_command_arguments(args) == ("alt_infile",
                                           "alt_outfile",
                                           "alt_logfile",
                                           "alt_detailfile")


# Five+ arguments (infile outfile logfile detail_file)
# Should issue a warning
def test_gca_extra():
    args = ["program name", "alt_infile", "alt_outfile",
            "alt_logfile", "alt_detailfile", "EXTRA JUNK"]
    with pytest.warns(UserWarning, match="Unused program arguments"):
        tmp = get_command_arguments(args)
        assert tmp == ("alt_infile", "alt_outfile",
                       "alt_logfile", "alt_detailfile")


# Blank/whitespace arguments
# Returns a blank string
# This will be an invalid file name so maybe we should issue a warning
def test_gca_blank():
    args = ["program name", "", "    "]
    assert get_command_arguments(args) == ("",
                                           "",
                                           "plog.txt",
                                           "")


# Non-array argument
# Should issue a warning and return the defaults
def test_gca_nonarray():
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


# Non-string array argument
# Should issue a warning and return the defaults
def test_gca_nonstring():
    msg = "Invalid command arguments, using defaults"
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(np.array([True, False]))
        assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
    with pytest.warns(UserWarning, match=msg):
        tmp = get_command_arguments(np.zeros(2))
        assert tmp == ("in.txt", "pout.txt", "plog.txt", "")
