"""
TEST_WR.PY Unit tests for write_results()
"""
import os
import tempfile

import pytest

from rcrbounds import write_results, read_data


# Basic functionality
def test_wr_basic():
    """write the specified array to the specified text file"""
    moment_vector = read_data("testing/testin1.txt")[3]
    with tempfile.TemporaryDirectory() as tmp:
        outfile = os.path.join(tmp, 'pout.txt')
        write_results(moment_vector, outfile)


# Exceptions to handle
def test_wr_readonly():
    """warn and continue if file is read-only"""
    moment_vector = read_data("testing/testin1.txt")[3]
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "testing/read-only-file.txt")


def test_wr_badfolder():
    """warn and continue if folder does not exist"""
    moment_vector = read_data("testing/testin1.txt")[3]
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "nonexistent-path-name/pout.txt")


def test_wr_illegalname():
    """warn and continue if file name is illegal"""
    moment_vector = read_data("testing/testin1.txt")[3]
    with pytest.warns(UserWarning, match="Cannot write"):
        write_results(moment_vector, "?")
