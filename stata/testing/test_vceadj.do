run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test VCEADJ option
*******************************************************************
* No adjustment, for comparison
rcr SAT Small_Class ${controls}
* Quadruple the covariance matrix (should double the standard errors)
rcr SAT Small_Class ${controls}, vceadj(4.0)
assert reldif( e(betaxCI_H)  , 7.774667956150855 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 1.383917661240516 ) <  ${tol}
* Zero.  Should make the standard errors zero.
rcr SAT Small_Class ${controls}, vceadj(0.0)
tempname T_V
mat `T_V' = J(5,5,0)
tempname C_V
matrix `C_V' = e(V)
assert mreldif( `C_V' , `T_V' ) < ${tol}
* Negative.  Should produce an error message
rcof "noisily rcr SAT Small_Class ${controls}, vceadj(-6)" == 111
