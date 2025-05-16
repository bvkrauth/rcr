run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test LAMBDA option
*******************************************************************
* Default lambda, for comparison
rcr SAT Small_Class ${controls}
savedresults save basic e()
* Specified at default (should be identical)
rcr SAT Small_Class ${controls}, lambda(0 1)
savedresults compare basic e()
* A single point (should give the same value for betaxH and betaxL
rcr SAT Small_Class ${controls}, lambda(0 0)
assert reldif( e(betaxCI_H)  , 6.488085264868137 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.914919882302702 ) <  ${tol}
* Going from - infinity
rcr SAT Small_Class ${controls}, lambda(. 0)
assert reldif( e(betaxCI_L)  , 3.914919882302702 ) <  ${tol}
if inlist("${exe}","python") & inlist(c(os),"Windows") {
	assert reldif(_b[betaxH], 8.169709965020822) <  ${tol}
}
else if inlist("${exe}","python") & inlist(c(os),"Unix") {
    assert reldif( e(betaxCI_H)  , 13.94704306786475 ) <  ${tol}
}
else if inlist("${exe}","windows-fortran") {
    * assert reldif( e(betaxCI_H)  , 12.30252920699166 ) <  ${tol}
}
else if inlist("${exe}","unix-fortran") {
    * assert reldif( e(betaxCI_H)  , 12.30252920699166 ) <  ${tol}
}
* Going to infinity
rcr SAT Small_Class ${controls}, lambda(0 .)
assert         e(betaxCI_H) > 8.99e305
assert         e(betaxCI_L) < -8.99e305
* Going from - infinity to infinity
rcr SAT Small_Class ${controls}, lambda(. .)
assert         e(betaxCI_H) > 8.99e305
assert         e(betaxCI_L) < -8.99e305
* Just below lambdaInf
rcr SAT Small_Class ${controls}, lambda(0 12.3)
assert reldif( e(betaxCI_H)  , 6.488085264868137 ) <  ${tol}
assert reldif( e(betaxCI_L)  , -18.08208794130959) <  ${tol}
* Just above lambdaInf
rcr SAT Small_Class ${controls}, lambda(12 12.32)
assert         e(betaxCI_H) > 8.99e305
assert         e(betaxCI_L) < -8.99e305
