"""
TEST_RD.PY Unit tests for write_results()
"""
import os
import tempfile
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import write_results, read_data, get_logfile, set_logfile


# Basic fnctionality
# Write the specified array to the specified text file
def test_wr_basic():
    oldlogfile = get_logfile()
    set_logfile(None)
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testin1.txt")
    with tempfile.TemporaryDirectory() as tmp:
        outfile = os.path.join(tmp, 'pout.txt')
        write_results(moment_vector, outfile)
    set_logfile(oldlogfile)


# Exceptions to handle
# Read-only file
# Issue a warning
def test_wr_readonly():
    oldlogfile = get_logfile()
    set_logfile(None)
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "read-only-file.txt")
    set_logfile(oldlogfile)


# Nonexistent folder name
# Issue a warning
def test_wr_badfolder():
    oldlogfile = get_logfile()
    set_logfile(None)
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "nonexistent-path-name/pout.txt")
    set_logfile(oldlogfile)


# Illegal file name
# Issue a warning
def test_wr_illegalname():
    oldlogfile = get_logfile()
    set_logfile(None)
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "?")
    set_logfile(oldlogfile)
