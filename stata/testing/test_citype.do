run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test CITYPE option
*******************************************************************
*Default options
rcr SAT Small_Class ${controls}
assert reldif( e(betaxCI_H)  , 6.488085264868137 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.259480713112168 ) <  ${tol}
savedresults save basic e()
* Specify the same citype as default (result should be the same)
rcr SAT Small_Class ${controls}, citype("conservative")
savedresults compare basic e()
* Specify Imbens - Manski (should be narrower than default)
rcr SAT Small_Class ${controls}, citype("imbens-manski")
assert reldif( e(betaxCI_H)  , 6.466066180253211 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.291579840373562 ) <  ${tol}
* Specify Lower.
rcr SAT Small_Class ${controls}, citype("lower")
assert reldif( r(betaxCI_H)  , 6.281236804837625 ) <  ${tol}
assert         r(betaxCI_L) < -8.0e306
* Specify Upper.
rcr SAT Small_Class ${controls}, citype("upper")
assert         r(betaxCI_H) > 8.0e306
assert reldif( r(betaxCI_L)  , 3.561021633566327 ) <  ${tol}
* Specify NONSENSE (or any other unsupported type).  Should generate a
* warning message and a missing CI
rcr SAT Small_Class ${controls}, citype("NONSENSE")
assert         e(betaxCI_L) == .
assert         e(betaxCI_H) == .
