"""
A set of examples you can run as a Python script

By default, only the basic examples run.  If you want
to run all examples, call this script with the optional
argument "all".  For example, in iPython this would be:

%run rcr_example all

"""
# pylint: disable=no-member
# Built in modules
import sys

# Third party modules
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
try:
    import matplotlib.pyplot as plt
except ImportError:
    PLT = False
else:
    PLT = True

# The rcrbounds module
from rcrbounds import RCR

# Setup

# check for "all" argument
SHOW_ALL = any(x.lower() == "all" for x in sys.argv)
# Load example data from the paper - Project STAR
dat = pd.read_stata("http://www.sfu.ca/~bkrauth/code/rcr_example.dta")
clust = dat["TCHID"]
wt = np.mod(clust, 2)

# statsmodels OLS example
#
#   The rcrbounds package is patterned after the OLS estimation procedures
#   in statsmodels.  The example below shows how those procedures work
#
# STEP 1: Construct the matrices of outcome and explanatory variables.
# I am using patsy to build them from a formula, or you can just
# build them with numpy or pandas
BASIC_FORMULA = ('SAT ~ Small_Class + White_Asian+  Girl ' +
                 ' + Free_Lunch + White_Teacher + Teacher_Experience' +
                 ' + Masters_Degree')
ols_endog, ols_exog = patsy.dmatrices(BASIC_FORMULA, dat)
# STEP 2: Construct the OLS model. The statsmodels.OLS class is a class
# defining an OLS regression model. It takes two main arguments:
# the outcome variable (n x 1 array-like object) and the
# explanatory variables (n x K array-like object) which include
# an intercept.
ols_model = sm.OLS(ols_endog, ols_exog)
# STEP 3: Fit the model. This can be done using the fit() method, which
# produces an object of class RegressionResults
ols_results = ols_model.fit()
# STEP 4: Display results. The RegressionResults object has a summary()
# method which can be used to display a table summarizing the results.
ols_summary = ols_results.summary()
print("sm.OLS(ols_endog, ols_exog).fit().summary(): \n \n",
      ols_summary,
      "\n\n")


# RCR example - basics

#   Similarly, you follow the same steps to estimate RCR bounds.
#   Note that the rcrbounds model has *two* endogenous varibles,
#   the outcome and the causal variable of interest
#
# STEP 1: construct the matrices.  endog is an n x 2 array of
# endogenous variables (outcome in the first column, causal variable
# in the second column).  As in the OLS example, you can use patsy
# to construct them from a formula, or use numpy/pandas to
# construct the arrays.
#
RCR_FORMULA = ('SAT + Small_Class ~ White_Asian+  Girl + Free_Lunch' +
               '+ White_Teacher + Teacher_Experience + Masters_Degree')
endog, exog = patsy.dmatrices(RCR_FORMULA, dat)
# STEP 2: Construct the RCR model using the RCR command.  It has
# two required arguments (endog and exog).  Optional arguments are
# covered in the examples below.
rcr_model = RCR(endog, exog)
# STEP 3: Fit the model using the RCR object's fit() method, which
# creates an RCRResult object.
rcr_results = rcr_model.fit()
# STEP 4: Print a summary table using the RCRResult object's summary()
# method.
rcr_summary = rcr_results.summary()
print("Basic RCR bounds example\n",
      "  RCR(endog, exog).fit().summary():\n\n",
      rcr_summary,
      "\n\n")

# Model specification options
if SHOW_ALL:
    # lambda_range:
    #
    # By default, bounds will be estimated for a relative correlation
    # parameter lambda betwen 0 and 1.  You can specify the lambda range
    # with the lambda_range optional argument:
    #
    lambda_example = RCR(endog, exog, lambda_range=(-1.0, 1.0)).fit().summary()
    print("lambda_range example\n",
          "  RCR(endog, exog, lambda_range=(-1.0, 1.0)).fit().summary()\n\n",
          lambda_example,
          "\n\n")
    #
    # weights:
    #
    # By default, all observations are weighted equally.  You
    # can provide observation weights using the weights argument.
    #
    weights_example = RCR(endog, exog, weights=wt).fit().summary()
    print("weights example\n\n",
          "  RCR(endog, exog, weights=wt).fit().summary()\n\n",
          weights_example,
          "\n\n")
    #
    # clustering:
    #
    # By default, standard errors are calculated under the independence
    # assumption. You can request cluster-robust standard errors
    # providing the clustering variable with groupvar:
    #
    cluster_example = RCR(endog, exog,
                          cov_type="cluster",
                          groupvar=clust).fit().summary()
    print("clustering example\n\n",
          "  RCR(endog, exog,\n",
          '      cov_type="cluster",\n',
          "      groupvar=clust).fit().summary()\n\n",
          cluster_example,
          "\n\n")

# Model fit/estimation options
if SHOW_ALL:
    #
    # cov_type:
    #
    # By default, cluster-robust standard errors are calculated if a
    # groupvar has been provided.  You can override this setting by
    # by setting cov_type to "nonrobust"
    #
    cov_type_example = RCR(endog, exog,
                           cov_type="nonrobust",
                           groupvar=clust).fit().summary()
    print("cov_type example\n\n",
          "  RCR(endog, exog,\n",
          '      cov_type="nonrobust",\n',
          "      groupvar=clust).fit().summary()\n\n",
          cluster_example,
          "\n\n")
    #
    # vceadj
    #
    # By default, the program calculates asymptotic covariance matrices
    # are calculated. You can add a finite sample degrees of freedom correction
    # using the vceadj option. vceadj takes a single numeric value
    # (default 1.0) and the estimated covariance matrix is simply rescaled
    # by the value of vceadj.  For example, suppose you want to apply the
    # standard OLS DOF correction of n/(n-K):
    #
    (n, k) = exog.shape
    vceadj_example = RCR(endog, exog,
                         vceadj=n/(n - k)).fit().summary()
    print("vceadj example\n\n",
          "  RCR(endog, exog,\n",
          "      vceadj=n/(n - k)).fit().summary()\n\n",
          vceadj_example,
          "\n\n")
    #
    # cilevel
    #
    # The default confidence intervals are 95%, but you can change that
    # with the cilevel option
    #
    cilevel_example = RCR(endog, exog, cilevel=90).fit().summary()
    print("cilevel example\n\n",
          "  RCR(endog, exog, cilevel=90).fit().summary()\n\n",
          cilevel_example,
          "\n\n")
    #
    # citype
    #
    # The default confidence interval for the parameter of
    # interest (betax_ci) is conservative in the sense that it does not
    # account for the width of the identified set. The Imbens-Manski
    # confidence interval is also available, and accounts for the width
    # of the identified set.  You can also produce one-tailed confidence
    # intervals with citype="lower" or citype="upper"
    citype_example = RCR(endog, exog).fit().summary(citype="Imbens-Manski")
    print("citype example\n\n",
          '  RCR(endog, exog).fit().summary(citype="Imbens-Manski")\n\n',
          citype_example,
          "\n\n")

# Postestimation options
if SHOW_ALL:
    #
    # Object properties
    #
    # The RCRResult object includes various properties you can view directly
    #
    print("nobs is the number of observations:\n",
          rcr_results.model.nobs)
    print("params is the vector of parameter estimates:\n",
          rcr_results.params)
    print("param_names is the parameter names:\n",
          rcr_results.param_names)
    print("cov_params is the estimated covariance matrix:\n",
          rcr_results.cov_params)
    print("The results object also includes the original model object:",
          rcr_results.model,
          "\nwhich you can use to re-fit the model with other settings:\n\n",
          rcr_results.model.fit(lambda_range=(0, 2)).summary(),
          "\n\n")
    #
    # Simple object methods
    #
    # The RCRResult object also includes various methods you can use to
    # calculate additional statistics or produceS tables and plots.
    #
    print("the params_se() method calculates standard errors:\n",
          rcr_results.params_se())
    print("the params_z() method calculates z-statistics",
          "(for the null parameter value of zero):\n",
          rcr_results.params_z())
    print("the params_pvalue() method calculates p-values",
          "(for the null parameter value of zero):\n",
          rcr_results.params_pvalue())
    print("the betax_ci() method calculates a confidence interval",
          "for the parameter of interest:\n",
          rcr_results.betax_ci(citype="Imbens-Manski"))
    print("the params_ci() method calculates confidence intervals",
          "for the other (point-identified) parameters:\n",
          rcr_results.params_ci())
    print("the test_betax() method performs hypothesis tests",
          "for the parameter of interest:\n",
          "P-value for H0: betax = 0:", rcr_results.test_betax(), "\n",
          "P-value for H0: betax = 5:", rcr_results.test_betax(5), "\n")
    #
    # summary() method
    #
    # The summary() method produces a statsmodels Summary object, which
    # includes as_csv, as_html, and as_latex methods in adddition
    # to the default as_text display.
    #
    print("Printing a summary of results as a CSV table:\n\n",
          rcr_summary.as_csv(),
          "\n\n")
    #
    # rcrplot() method
    #
    # The rcrplot() method creates a matplotlib plot. Unlike the other
    # libraries used in this package, matplotlib is not pre-installed
    # in the Anaconda distribution of Python.  So you may need to install
    # it if you want to create plots.
    #
    if PLT:
        # The base plot shows the lambda(betax) function over the range
        # from -50 to 50
        fig, axes = plt.subplots()
        ax = rcr_results.rcrplot(ax=axes)
        plt.show(block=False)
        # You can add various other plot elements, including a legend
        fig, axes = plt.subplots()
        ax = rcr_results.rcrplot(ax=axes,
                                 title="Example: added plot elements",
                                 tsline=True,
                                 lsline=True,
                                 idset=True,
                                 legend=True)
        plt.show(block=False)
        # The default range is very wide, so you may want to show a smaller
        # range using the xlim option
        fig, axes = plt.subplots()
        ax = rcr_results.rcrplot(ax=axes,
                                 title="Example: setting xlim",
                                 xlim=(0.0, 6.0),
                                 tsline=True,
                                 lsline=True,
                                 idset=True,
                                 legend=True)
        plt.show(block=False)
        # You can also add or change titles, change the y-axis range,
        # and adjust line styles and colors.  See the method's
        # docstring for details.
