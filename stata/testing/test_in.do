run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test IN option
*******************************************************************
rcr SAT Small_Class ${controls} in 1/100
assert reldif( e(betaxCI_H)  , 27.26163444874422 ) <  ${tol}
assert reldif( e(betaxCI_L)  , .9203357905110554 ) <  ${tol}
* These two specifications exclude too many observations.
* They should produce an error message.
rcof "noisily rcr SAT Small_Class ${controls} in 1/1" == 2001
rcof "noisily rcr SAT Small_Class ${controls} in 1/3" == 2001
