run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test NLCOM postestimation command
*******************************************************************
rcr SAT Small_Class ${controls}
* Testing NLCOM for the width of the identified set
nlcom _b[betaxH] - _b[betaxL]
tempname b v
matrix `b' = r(b)
assert reldif( `b'[1,1] , .0664588086015998) < ${tol}
matrix `v' = r(V)
assert reldif( `v'[1,1] , .4695016651521872) < ${tol}
mat drop `b' `v'
* When the parameter is not identified, the width is infinite and so is its
* standard error. Stata does not handle infinity well, so the width is a very
* large number and the standard error is zero (which Stata displays as a
* missing value).
rcr SAT Small_Class ${controls}, lambda(. .)
nlcom _b[betaxH] - _b[betaxL]
matrix `b' = r(b)
assert `b'[1,1] > 8e305
matrix `v' = r(V)
assert `v'[1,1] == 0
mat drop `b' `v'
