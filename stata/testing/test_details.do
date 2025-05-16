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
summarize lambda , detail
assert reldif(r(p75), 44.78824776422209 ) <  ${tol}
assert reldif(r(p50), 36.85231368869719 ) <  ${tol}
assert reldif(r(p25), 33.43887260299969 ) <  ${tol}
summarize betax , detail
assert reldif(r(p75), 25.00250025002501 ) <  ${tol}
assert r(p50) == 0
assert reldif(r(p25), -25.002500250025  ) <  ${tol}
