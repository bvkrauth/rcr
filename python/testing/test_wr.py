"""
TEST_WR.PY Unit tests for write_results()
"""
import os
import tempfile
import sys

import pytest

sys.path.append("./")
sys.path.append("../")
from rcr import write_results, read_data


# Basic fnctionality
# Write the specified array to the specified text file
def test_wr_basic():
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    with tempfile.TemporaryDirectory() as tmp:
        outfile = os.path.join(tmp, 'pout.txt')
        write_results(moment_vector, outfile)


# Exceptions to handle
# Read-only file
# Issue a warning
def test_wr_readonly():
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "testing/read-only-file.txt")


# Nonexistent folder name
# Issue a warning
def test_wr_badfolder():
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "nonexistent-path-name/pout.txt")


# Illegal file name
# Issue a warning
def test_wr_illegalname():
    n_moments, n_lambda, external_big_number, moment_vector, \
        lambda_range = read_data("testing/testin1.txt")
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "?")
