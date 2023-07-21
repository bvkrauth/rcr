run config_test `0'
which rcr
which rcr_predict
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test PREDICT postestimation command
*******************************************************************
rcr SAT Small_Class ${controls}
* The predict command is not supported, so this should produce an error
* message
rcof "noisily predict" == 321
