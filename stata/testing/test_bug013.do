run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*****************************************************************************
* Test BUG013
*****************************************************************************
* Description of bug: Standard error results are inaccurate when
*                     lambda_L = -infty.
rcr SAT Small_Class White_Asian , lambda(. 0)
