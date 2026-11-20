"""
TEST_SE.PY: Unit tests for stata_exe()
"""
import os

import numpy as np

from rcrbounds import stata_exe


# Basic functionality
def test_se_basic(tmp_path, infile):
    """translate_result with default arguments"""
    outfile = os.path.join(tmp_path, 'out.txt')
    logfile = os.path.join(tmp_path, 'log.txt')
    argv = ["rcrbounds", infile, outfile, logfile]
    assert os.path.exists(infile)
    assert not os.path.exists(outfile)
    assert not os.path.exists(logfile)
    stata_exe(argv)
    assert os.path.exists(outfile)
    assert os.path.exists(logfile)


def test_se_details(tmp_path, infile):
    """translate_result with default arguments"""
    outfile = os.path.join(tmp_path, 'out.txt')
    logfile = os.path.join(tmp_path, 'log.txt')
    detfile = os.path.join(tmp_path, 'det.txt')
    argv = ["rcrbounds", infile, outfile, logfile, detfile]
    assert os.path.exists(infile)
    assert not os.path.exists(outfile)
    assert not os.path.exists(logfile)
    assert not os.path.exists(detfile)
    stata_exe(argv)
    assert os.path.exists(outfile)
    assert os.path.exists(logfile)
    assert os.path.exists(detfile)


def test_se_details_range(tmp_path, infile):
    """details file covers effect values well beyond +/-50 (see #13)"""
    outfile = os.path.join(tmp_path, 'out.txt')
    logfile = os.path.join(tmp_path, 'log.txt')
    detfile = os.path.join(tmp_path, 'det.txt')
    argv = ["rcrbounds", infile, outfile, logfile, detfile]
    stata_exe(argv)
    effectvec = np.loadtxt(detfile, delimiter=",", skiprows=1)[:, 0]
    assert np.sum((effectvec > 50.0) & (effectvec < 1000.0)) > 100
    assert np.sum((effectvec < -50.0) & (effectvec > -1000.0)) > 100


def test_se_unbounded_rc(tmp_path, infile):
    """Stata's -maxdouble() for missing lambda_L is treated as -infinity,
    so betaxH equals betaxInf (issue #10)"""
    with open(infile, encoding="utf-8") as in_file:
        lines = in_file.read().splitlines()
    lines[2] = " -1.797693134862e+308\t 1"
    infile2 = os.path.join(tmp_path, 'in.txt')
    with open(infile2, mode="w", encoding="utf-8") as in_file:
        in_file.write("\n".join(lines) + "\n")
    outfile = os.path.join(tmp_path, 'out.txt')
    logfile = os.path.join(tmp_path, 'log.txt')
    stata_exe(["rcrbounds", infile2, outfile, logfile])
    results = np.loadtxt(outfile)
    assert np.array_equal(results[4, :], results[1, :])
