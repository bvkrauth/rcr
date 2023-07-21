run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test IF option
*******************************************************************
* This specification should work
rcr SAT Small_Class ${controls} if  Free_Lunch > 0
assert reldif( e(betaxCI_H)  , 8.284396403084376 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 4.126516220297534 ) <  ${tol}
* This specification excludes all observations, so it should (and does)
* produce an error message
rcof "noisily rcr SAT Small_Class ${controls} if Free_Lunch > 2" == 2000
