"""
TEST_RP.PY: Unit tests for rcrplot() method
"""
import os

import numpy as np
import pandas as pd
import pytest
import patsy

# sys.path.append("./")
# sys.path.append("../")
from rcrbounds import RCR, read_data  # pylint: disable=wrong-import-position

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture(name="badin1")
def fixture_badin1():
    """get bad data"""
    infile = os.path.join(FIXTURE_DIR, "badin1.txt")
    return str(infile)


@pytest.fixture(name="badin2")
def fixture_badin2():
    """get bad data"""
    infile = os.path.join(FIXTURE_DIR, "badin2.txt")
    return str(infile)


@pytest.fixture(name="badin3")
def fixture_badin3():
    """get bad data"""
    infile = os.path.join(FIXTURE_DIR, "badin3.txt")
    return str(infile)


@pytest.fixture(name="badin4")
def fixture_badin4():
    """get bad data"""
    infile = os.path.join(FIXTURE_DIR, "badin4.txt")
    return str(infile)


@pytest.fixture(name="badin5")
def fixture_badin5():
    """get bad data"""
    infile = os.path.join(FIXTURE_DIR, "badin5.txt")
    return str(infile)


@pytest.fixture(name="true_result")
def fixture_true_result():
    """get true results for test data"""
    outfile = os.path.join(FIXTURE_DIR, "testout1.txt")
    result = np.asarray(pd.read_csv(outfile,
                                    delimiter=" ",
                                    header=None,
                                    skipinitialspace=True))
    return result


@pytest.fixture(name="read_only_file")
def fixture_read_only_file(tmp_path):
    """get read only file"""
    infile = tmp_path / "read-only-file.txt"
    infile.write_text("This file is read only")
    # mode = infile.stat().st_mode
    infile.chmod(0o444)
    return str(infile)


@pytest.fixture(name="infile")
def fixture_infile():
    """get test data"""
    infile = os.path.join(FIXTURE_DIR, "testin1.txt")
    return str(infile)


@pytest.fixture(name="moment_vector")
def fixture_moment_vector(infile):
    """get moment vector from test data"""
    return read_data(infile)[3]


@pytest.fixture(name="dat")
def fixture_dat():
    """get test data from web page"""
    fname = "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
    return pd.read_stata(fname)


@pytest.fixture(name="rcr_formula")
def fixture_rcr_formula():
    """construct formula for test example"""
    rcr_left = "SAT + Small_Class ~ "
    rcr_right1 = "White_Asian + Girl + Free_Lunch + White_Teacher + "
    rcr_right2 = "Teacher_Experience + Masters_Degree"
    return rcr_left + rcr_right1 + rcr_right2


@pytest.fixture(name="endog")
def fixture_endog(dat, rcr_formula):
    """get endogenous variables"""
    # pylint: disable=no-member
    endog = patsy.dmatrices(rcr_formula, dat)[0]
    return endog


@pytest.fixture(name="exog")
def fixture_exog(dat, rcr_formula):
    """get endogenous variables"""
    # pylint: disable=no-member
    exog = patsy.dmatrices(rcr_formula, dat)[1]
    return exog


@pytest.fixture(name="endog_df")
def fixture_endog_df(dat, rcr_formula):
    """get endogenous variables"""
    # pylint: disable=no-member
    endog = patsy.dmatrices(rcr_formula,
                            dat,
                            return_type="dataframe")[0]
    return endog


@pytest.fixture(name="exog_df")
def fixture_exog_df(dat, rcr_formula):
    """get endogenous variables"""
    # pylint: disable=no-member
    exog = patsy.dmatrices(rcr_formula, dat, return_type="dataframe")[1]
    return exog


@pytest.fixture(name="model")
def fixture_model(endog, exog):
    """construct RCR mdodel"""
    return RCR(endog, exog)


@pytest.fixture(name="results")
def fixture_results(model):
    """fit RCR mdodel"""
    return model.fit()


@pytest.fixture(name="weights")
def fixture_weights(dat):
    """get weights"""
    weight = np.mod(dat["TCHID"], 2)
    weight.name = "TCHID_wt"
    return weight


@pytest.fixture(name="clusters")
def fixture_clusters(dat):
    """get cluster IDs"""
    clust = dat["TCHID"]
    return clust
