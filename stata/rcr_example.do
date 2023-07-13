capture log close
log using "rcr_example", replace
/*******************************************************************
* RCR_EXAMPLE.DO
*
* This is a do-file designed to provide a demonstration of the RCR
* program and its main features.
*
*
********************************************************************/
clear
set mem 10m
set more off

/* This example uses the kindergarten group from the Project STAR data set used in the paper. */
/* The data set has already been subject to the fixed effects transformation */
/* Use the local version of this file if available */
local fname "rcr_example.dta"
capture confirm file `fname'
/* Otherwise use the web version */
if (_rc != 0) {
    local fname "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
}
di "Using data file `fname'"
use "`fname'", clear

/* This is just an ordinary regression, similar to that in Table 2 of the paper. */
reg SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree
/* This is RCR analysis of the same regression, with default options */
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree
/* The results of the RCR command are stored in e() and so can be saved and used in tables */
estimates store default
estimates table default , stats(betaxCI_L betaxCI_H N)

/* The LAMBDA option allows for the user to control the range of lambda values */
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, lambda(0 2)
/* Use missing (.) to denote absence of an upper or lower bound. */
/* Note that in this particular example, the bounds on betax are (-infinity,infinity).  Since Stata
 * has no facility for manipulating infinity, infinity is reported as just a very large number. */
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, lambda(0 .)

/* By default confidence intervals are estimated conservatively.  The CITYPE option allows the choice of other methods. */
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, citype("Imbens-Manski")

/* The VCEADJ allows for the user to adjust the covariance matrix estimates by a multiplicative factor.  This can be 
 * used to implement degrees-of-freedom corrections when the data has been transformed (e.g., normalization or
 * fixed-effects transformations). */
quietly duplicates report SCHID
local dofcorrection = r(N) / (r(N) - r(unique_value))
di "There are " r(N) " observations, but the data have received the fixed-effects transformation." _newline /*
*/ "That is, each variable is expressed in terms of its deviation from the school-level average." _newline  /*
*/ "There are " r(unique_value) " schools in the data, so the effective number of observations is " r(N)-r(unique_value) "." _newline /*
*/ "So we should multiply the covariance matrix by " r(N) "/(" r(N) "-" r(unique_value) ") = `dofcorrection'."
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, vceadj(`dofcorrection')

/* The LEVEL and CLUSTER options work the same as they do for most Stata commands (e.g., REG) */
/* IF, IN, and weights (all 4 types) are supported as well. */
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, level(90)
rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, cluster(TCHID)

/* The ESTAT postestimation command is supported */
estat vce
estat summarize

/* For hypothesis testing of point-identified parameters (including betaxH and betaxL), the TESTNL and NLCOM postestimation commands can be used */
testnl _b[betaxL] = 0
nlcom _b[betaxH] - _b[betaxL]

/* There is also a postestimation command test_betax for testing hypotheses about the interval-identified parameter betax */
/* You can test any point null hypothesis.  This line of code tests H0: betax = 1.25 */
test_betax = 1.25
/* If you call test_betax without any arguments, the default is to test H0: betax = 0. */
test_betax

log close
