run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test basic regression with different explanatory variables
*******************************************************************
quietly gen zero = 0
if (c(version) >= 14) {
    local rng = c(rng)
    set rng kiss32
}
set seed 53
forvalues vnum = 1 / 20 {
    quietly gen x`vnum' = runiform()
}
if (c(version) >= 14) {
    set rng `rng'
}
***** Things that should produce an error message
* No treatment or control variable. Correctly produces an error message.
rcof "noisily rcr SAT" == 102
* No control variable. Correctly produces an error message.
rcof "noisily rcr SAT Small_Class" == 102
* Outcome variable doesn't vary. Correctly produces an error message.
rcof "noisily rcr zero Small_Class White_Asian" == 1
* Treatment variable doesn't vary. Correctly produces an error message.
rcof "noisily rcr SAT zero White_Asian" == 1
* Control variable doesn't vary. Correctly produces an error message.
rcof "noisily rcr SAT Small_Class zero" == 1
* More than 25 control variables. Correctly produces an error message.
scalar max_controls = floor(sqrt(c(max_matsize))) - 3
if max_controls <= 25 {
    rcof "noisily rcr SAT Small_Class ${controls} x1-x20" == 103
}
else {
    rcof "noisily rcr SAT Small_Class ${controls} x1-x20" == 0
}

***** Things that should produce a warning message but give it a try
* Outcome and treatment are identical (collinear)
rcof "noisily rcr SAT SAT White_Asian" == 0
* Treatment and control are identical (collinear)
rcof "noisily rcr SAT Small_Class Small_Class" == 0
* Outcome and control are identical (collinear)
* This leads to an error in Unix but not in Windows
if inlist(c(os),"Windows") {
    rcof "noisily rcr SAT Small_Class SAT" == 0
}
else if inlist(c(os),"Unix") {
    rcof "noisily rcr SAT Small_Class SAT" > 0
}
* Outcome and control are exactly unrelated
* This leads to an error in the external program RCR.EXE. It would be
* cleaner to catch it in the ado file
preserve
quietly expand 2, generate(unrelated_variable)
if inlist("${exe}","python") {
    rcof "noisily rcr SAT Small_Class unrelated_variable" == 0
}
else if inlist("${exe}","windows-fortran","unix-fortran") {
    rcof "noisily rcr SAT Small_Class unrelated_variable" == 1
}
restore

***** Things that should work
* Just one control variable
rcr SAT Small_Class White_Asian
assert reldif( e(betaxCI_H)  , 6.494755100453821 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 2.893150619179314 ) <  ${tol}
* Up to 25 control variables
rcr SAT Small_Class ${controls} x1 - x19
assert reldif( e(betaxCI_H)  , 7.341531127530879 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.976272660323884 ) <  ${tol}
* Collinearity among control variables
* Note the behavior here has changed.  Previously, the program issued an
* error message when control variables were collinear. Now it just omits
* variables as needed, just like REG does.
rcr SAT Small_Class ${controls} zero
assert reldif( e(betaxCI_H)  , 6.488085265277207 ) <  ${tol}
rcr SAT Small_Class ${controls} Girl
assert reldif( e(betaxCI_H)  , 6.488085265277207 ) <  ${tol}
