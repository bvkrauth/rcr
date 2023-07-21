run config_test `0'
which rcr
which test_betax
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test TEST_BETAX postestimation command
*******************************************************************
rcr SAT Small_Class ${controls}
* Default is to test betax = 0
test_betax
assert reldif( r(p)  , 6.75341735534e - 08 ) <  ${tol}
* Can also use =exp
test_betax =  3.29158
assert reldif( r(p)  , .0500000163998828 ) <  ${tol}
