run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test TESTNL postestimation command
*******************************************************************
rcr SAT Small_Class ${controls}
* Testing TESTNL
testnl _b[betaxH] = 0
assert reldif( r(chi2)  , 62.78825430596226 ) <  ${tol}
* What happens when the parameter is not identified?
rcr SAT Small_Class ${controls}, lambda(. .)
testnl _b[betaxH] = 0
assert         r(p)    == .
assert         r(chi2) == 0
assert         r(df)   == 0
