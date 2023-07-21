run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test DETAILS option
*******************************************************************
* Call RCR with DETAILS specified. This also produces a plot
rcr SAT Small_Class ${controls}, details
* Check to make sure the DETAILS data has not changed
assert _N == 30000
summarize lambda
assert reldif( r(mean)   , -83505016503.90866) < ${tol}
summarize betax
assert reldif( r(mean)   , -4.66921103756e-08) < ${tol}
