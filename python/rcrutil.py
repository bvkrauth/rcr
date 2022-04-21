"""
RCRUTIL.PY Function library for RCR program

Details TBA.
"""
# Standard library imports
import sys
import warnings
from datetime import datetime

# Third party imports
import numpy as np
import pandas as pd

# Local application imports (none)


def get_command_arguments(args):
    """Retrieve command arguments, usually from sys.argv."""
    # ARGS should be a list of 1 to 5 strings like sys.argv
    if isinstance(args, list) and all([type(item) == str for item in args]):
        if (len(args) > 5):
            msg = "Unused program arguments {0}".format(args[5:])
            warnings.warn(msg)
            pass
    else:
        msg = "Invalid command arguments, using defaults: {0}".format(args)
        warnings.warn(msg)
        args = []
    infile = args[1].strip() if len(args) > 1 else "in.txt"
    outfile = args[2].strip() if len(args) > 2 else "pout.txt"
    logfile = args[3].strip() if len(args) > 3 else "plog.txt"
    detail_file = args[4].strip() if len(args) > 4 else ""
    return infile, outfile, logfile, detail_file


def start_logfile(logfile):
    """Start the log file."""
    set_logfile(logfile)
    # TODO: Add any warnings that have already been thrown
    write_to_logfile("Log file {0} for RCR version 1.0\n".format(logfile),
                     mode="w")
    start_time = datetime.now().strftime("%H:%M on %m/%d/%y")
    write_to_logfile("Run at {0}.\n".format(start_time))


def read_data(infile):
    """Read RCR data from input file"""
    write_to_logfile("Reading data from input file {0}.\n".format(infile))
    # infile argument should be a single string
    if not isinstance(infile, str):
        msg = "Infile should be a single string"
        warn(msg)
    # Line 1 should be three whitespace delimited numbers
    try:
        line1 = pd.read_csv(infile,
                            delim_whitespace=True,
                            skiprows=[1, 2],
                            header=None).values[0, ]
        n_moments, n_lambda, external_big_number = tuple(line1)
    except FileNotFoundError:
        msg = "infile {0} not found.\n".format(infile)
        die(msg)
    except ValueError:
        msg = "Incorrect format in line 1 of infile {0}.\n".format(infile)
        die(msg)
    except:
        msg = "Unknown problem with line 1 of infile {0}.\n".format(infile)
        die(msg)
    else:
        msg1 = "Line 1: n_moments = {0}, n_lambda = {1}".format(n_moments,
                                                                n_lambda)
        msg2 = "external_big_number = {0}.\n".format(external_big_number)
        write_to_logfile(msg1 + ", " + msg2)
    # Line 2 should be n_moments whitespace delimited numbers
    try:
        moment_vector = pd.read_csv(infile,
                                    delim_whitespace=True,
                                    skiprows=[0, 2],
                                    header=None).values[0, ].astype(np.float64)
    except ValueError:
        msg = "Incorrect format in line 2 of infile {0}.\n".format(infile)
        die(msg)
    except:
        msg = "Unknown problem with line 2 of infile {0}.\n".format(infile)
        die(msg)
    else:
        msg = "Line 2: moment_vector = a vector of length {0}.\n"
        write_to_logfile(msg.format(len(moment_vector)))
    # Lines 3+ should be two whitespace delimited numbers each
    try:
        lambda_range = pd.read_csv(infile,
                                   delim_whitespace=True,
                                   skiprows=[0, 1],
                                   header=None).values[0, ].astype(np.float64)
    except ValueError:
        msg = "Incorrect format in line 3 of infile {0}.\n".format(infile)
        die(msg)
    except:
        msg = "Unknown problem with line 3 of infile {0}.\n".format(infile)
        die(msg)
    else:
        write_to_logfile("Line 3: lambda_range = {0}.\n".format(lambda_range))
        write_to_logfile("For calculations, lambda_range,...\n")
    write_to_logfile("Data successfully loaded from file {0}\n".format(infile))
    # Check to make sure n_lambda is a valid (i.e., positive) value
    n_lambda = int(n_lambda)
    #   1. It should be twice as long as lambda_range. if not, just reset it
    if len(lambda_range) != 2*n_lambda:
        msg = "n_lambda reset from {0} to len(lambda_range)/2 = {1}."
        warn(msg.format(n_lambda, int(len(lambda_range)/2)))
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
        msg = "n_moments reset from {0} to len(moment_vector) = {1}."
        warn(msg.format(n_moments, len(moment_vector)))
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
        msg = "Largest Python real ({0}) is less than largest in Stata {1}"
        warn(msg.format(sys.float_info.max, external_big_number))
    return n_moments, n_lambda, external_big_number, \
        moment_vector, lambda_range


def write_results(result_matrix, outfile):
    """Write results to outfile."""
    write_to_logfile("Writing results to output file {0}.\n".format(outfile))
    write_to_logfile("Actual results = ...\n")
    try:
        with np.printoptions(threshold=np.inf, linewidth=np.inf):
            np.savetxt(outfile, result_matrix, delimiter=" ")
    except:
        msg = "Cannot write to output file {0}.".format(outfile)
        warn(msg)
    else:
        write_to_logfile("RCR successfully concluded.\n")


def write_to_logfile(str, mode="a"):
    """Write a note to the log file."""
    logfile = get_logfile()
    if logfile is None:
        return
    try:
        with open(logfile, mode) as lf:
            lf.write(str)
    except:
        msg = "Cannot write to logfile {0}.".format(logfile)
        warnings.warn(msg)
        pass


def warn(msg):
    """Issue warning (to logfile and python warning system) but continue."""
    write_to_logfile("WARNING: " + msg + "\n")
    warnings.warn(msg)


def die(msg):
    """Fatal error - write message to log file and then shut down."""
    write_to_logfile("FATAL ERROR: " + msg)
    return sys.exit(msg)


def set_logfile(fname):
    """Set name of log file"""
    global logfile
    if isinstance(fname, str) or fname is None:
        logfile = fname
    else:
        pass


def get_logfile():
    """Retrieve name of log file.  If undefined, return None"""
    global logfile
    if "logfile" not in globals():
        logfile = None
    return logfile
