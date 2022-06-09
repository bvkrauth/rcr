"""
TEST_RP.PY: Unit tests for rcrplot() method
"""
import sys

import numpy as np
import pandas as pd
import pytest
import patsy

sys.path.append("./")
sys.path.append("../")
from rcr import RCR, read_data  # pylint: disable=wrong-import-position


@pytest.fixture
def moment_vector():
    """get moment vector from test data"""
    mv = read_data("testing/testin1.txt")[3]
    return mv


@pytest.fixture
def dat():
    """get test data from web page"""
    fname = "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
    return pd.read_stata(fname)


@pytest.fixture
def rcr_formula():
    """construct formula for test example"""
    rcr_left = "SAT + Small_Class ~ "
    rcr_right1 = "White_Asian + Girl + Free_Lunch + White_Teacher + "
    rcr_right2 = "Teacher_Experience + Masters_Degree"
    return rcr_left + rcr_right1 + rcr_right2


@pytest.fixture
def endog(dat, rcr_formula):
    """get endogenous variables"""
    endog = patsy.dmatrices(rcr_formula, dat)[0]
    return endog


@pytest.fixture
def exog(dat, rcr_formula):
    """get endogenous variables"""
    exog = patsy.dmatrices(rcr_formula, dat)[1]
    return exog


@pytest.fixture
def model(endog, exog):
    """construct RCR mdodel"""
    return RCR(endog, exog)


@pytest.fixture
def results(model):
    """fit RCR mdodel"""
    return model.fit()


@pytest.fixture
def weights(dat):
    """get weights"""
    wt = np.mod(dat["TCHID"], 2)
    wt.name = "TCHID_wt"
    return wt


@pytest.fixture
def clusters(dat):
    """get cluster IDs"""
    clust = dat["TCHID"]
    return clust
