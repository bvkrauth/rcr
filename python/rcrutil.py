"""
RCRUTIL.PY Function library for RCR program

Details TBA.
"""
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from warnings import warn


def get_command_arguments(args):
    """Retrieve command arguments"""
    global logfile
    # bug here - if you give "  " as an argument it will crash
    infile = args[1].strip() if len(args) > 1 else "in.txt"
    outfile = args[2].strip() if len(args) > 2 else "pout.txt"
    logfile = args[3].strip() if len(args) > 3 else "plog.txt"
    detail_file = args[4].strip() if len(args) > 4 else ""
    print("infile = {0}, outfile = {1}".format(infile, outfile))
    print("logfile = {0}, detail_file = {1}".format(logfile, detail_file))
    write_to_logfile("Log file {0} for RCR version 1.0\n".format(logfile),
                     mode="w")
    if (len(args) > 5):
        msg = "WARNING: Unused program arguments {0}".format(args[5:])
        write_to_logfile(msg)
        warn(msg)
    return infile, outfile, detail_file


def read_data(infile):
    """Read RCR data from infile"""
    start_time = datetime.now().strftime("%H:%M on %m/%d/%y")
    write_to_logfile("Run at {0}.\n".format(start_time))
    write_to_logfile("Reading data from input file {0}.\n".format(infile))
    try:
        line1 = pd.read_csv(infile,
                            sep="\s+",
                            skiprows=[1, 2],
                            header=None).values[0, ]
        n_moments, n_lambda, external_big_number = tuple(line1)
    except FileNotFoundError:
        msg = "infile {0} not found.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    except ValueError:
        msg = "Incorrect format in line 1 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    except:
        msg = "Unknown problem with line 1 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    else:
        msg = "Line 1: n_moments = {0}, n_lambda = {1}, external_big_number = {2}.\n".format(n_moments, n_lambda, external_big_number)
        write_to_logfile(msg)
    try:
        moment_vector = pd.read_csv(infile,
                                    sep="\s+",
                                    skiprows=[0, 2],
                                    header=None).values[0, ].astype(np.float64)
    except ValueError:
        msg = "Incorrect format in line 2 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    except:
        msg = "Unknown problem with line 2 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    else:
        msg = "Line 2: moment_vector = a vector of length {0}.\n".format(len(moment_vector))
        write_to_logfile(msg)
    try:
        lambda_range = pd.read_csv(infile,
                                   sep="\s+",
                                   skiprows=[0, 1],
                                   header=None).values[0, ].astype(np.float64)
    except ValueError:
        msg = "Incorrect format in line 3 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    except:
        msg = "Unknown problem with line 3 of infile {0}.\n".format(infile)
        write_to_logfile("FATAL ERROR: " + msg)
        return sys.exit()
    else:
        write_to_logfile("Line 3: lambda_range = {0}.\n".format(lambda_range))
        write_to_logfile("For calculations, lambda_range,...\n")
    write_to_logfile("Data successfully loaded from file {0}\n".format(infile))
    # Check to make sure n_lambda is a valid (i.e., positive) value
    n_lambda = int(n_lambda)
    #   1. It should be twice as long as lambda_range. if not, just reset it
    if len(lambda_range) != 2*n_lambda:
        msg = "n_lambda reset from {0} to len(lambda_range)/2 = {1}.".format(n_lambda, int(len(lambda_range)/2))
        warn(msg)
        n_lambda = int(len(lambda_range) / 2)
    #   2. It should be positive.
    assert n_lambda > 0
    #   3. For now, it should be one.
    assert n_lambda == 1
    # Check to make sure n_moments is a valid value
    n_moments = int(n_moments)
    #   1. It should be the same as the length of moment_vector.  if not,
    #      just reset it.
    if n_moments != len(moment_vector):
        msg = "n_moments reset from {0} to len(moment_vector) = {1}.".format(n_moments, len(moment_vector))
        warn(msg)
        n_moments = len(moment_vector)
    #   2. It must be at least 9 (i.e., there must be at least one explanatory
    #      variable)
    assert n_moments >= 9
    #   3. The number of implied explanatory variables must be an integer
    k = int((np.sqrt(9 + 8 * n_moments) - 1) / 2)
    assert (2 * (n_moments + 1)) == int(k ** 2 + k)
    # Check to make sure external_big_number is a valid value
    assert external_big_number > 0.0
    # If external_big_number is bigger than sys.float_info.max, then issue a
    # warning but don't stop program.
    # TODO: I'm not satisfied with this.
    if (external_big_number > sys.float_info.max):
        write_to_logfile("Warning: largest real number in executable ({0}) is less than largest in Stata {1}".format(sys.float_info.max, external_big_number))
        warn("Largest real number in executable ({0}) is less than largest in Stata {1}".format(sys.float_info.max, external_big_number))
    return n_moments, n_lambda, external_big_number, moment_vector, lambda_range


def write_results(result_matrix, outfile):
    """Write results to outfile."""
    write_to_logfile("Writing results to output file {0}.\n".format(outfile))
    write_to_logfile("Actual results = ...\n")
    try:
        np.savetxt(outfile, result_matrix, delimiter=" ")
    except:
        write_to_logfile("Cannot write to output file {0}.".format(outfile))
    else:
        write_to_logfile("RCR successfully concluded.\n")


def write_to_logfile(str, mode="a"):
    """Write a note to the log file."""
    global logfile
    # print(str)
    try:
        lf = open(logfile, mode)
        lf.write(str)
        lf.close()
    except:
        warn("Cannot write to logfile {0}.".format(logfile))
        pass
