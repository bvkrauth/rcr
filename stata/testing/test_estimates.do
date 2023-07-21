run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test ESTIMATES postestimation command
*******************************************************************
rcr SAT Small_Class ${controls}
* ESTIMATES STORE
estimates store basic
* ESTIMATES STATS
estimates stats basic
* ESTIMATES TABLE
estimates table basic, b(%9.3f) se p title("Results from RCR") stats(N)
