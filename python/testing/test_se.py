"""
TEST_SE.PY: Unit tests for stata_exe()
"""
import os

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
