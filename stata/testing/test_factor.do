run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test factor variables
*******************************************************************
label define studcat ///
    1 "White/Asian boy" ///
    2 "White/Asian girl" ///
    3 "NWA boy" ///
    4 "NWA girl", replace
gen studcat:studcat = 1 if (White_Asian > 0) & (Girl < 0)
replace studcat = 2 if (White_Asian > 0) & (Girl > 0)
replace studcat = 3 if (White_Asian < 0) & (Girl < 0)
replace studcat = 4 if (White_Asian < 0) & (Girl > 0)
label variable studcat "Student category"
gen wt = (studcat == 1)
local fvcontrols "Free_Lunch White_Teacher Teacher_Experience Masters_Degree "
* Issue an error if either of the first two variables are factor variables
rcof "noisily rcr i.studcat Small_Class ${controls}" == 198
rcof "noisily rcr SAT i.studcat ${controls}" == 198
* Issue an error if there is no variation in the control variables
rcof "noisily rcr SAT Small_Class i.studcat if studcat == 1" == 1
* Otherwise, expand out factor variables just like in REG
reg SAT Small_Class i.studcat Free_Lunch White_Teacher Teacher_Experience
scalar reg_coef = _b[Small_Class]
* The RCR version should produce similar results
rcr SAT Small_Class i.studcat Free_Lunch White_Teacher Teacher_Experience
assert reldif(_b[betaxL] , reg_coef ) < ${tol}
* Omit factor variables that do not vary given IF/IN/WEIGHTS
reg SAT Small_Class i.studcat `fvcontrols' if studcat == 1
scalar reg_coef = _b[Small_Class]
rcr SAT Small_Class i.studcat `fvcontrols' if studcat == 1
assert reldif(_b[betaxL] , reg_coef ) < ${tol}
rcr SAT Small_Class i.studcat `fvcontrols' [pw = wt]
assert reldif(_b[betaxL] , reg_coef ) < ${tol}
