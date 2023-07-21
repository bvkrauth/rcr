run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test time series operators
*******************************************************************
gen timevar1 = _n
sort TCHID timevar1
by TCHID : gen timevar2 = _n
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
local tscontrols "White_Teacher Teacher_Experience Masters_Degree"

* Should issue an error message if tsset has not been called
tsset, clear
rcof "noisily rcr SAT Small_Class i.studcat L.Free_Lunch `tscontrols'" == 111

* Should allow time series operators for dependent, treatment and control,
* variables, possibly in combination with factor variables
tsset timevar1
reg L.SAT L.Small_Class i.studcat L.Free_Lunch `tscontrols' 
scalar reg_coef = _b[L.Small_Class]
rcr L.SAT L.Small_Class i.studcat L.Free_Lunch `tscontrols'
assert min(reldif(_b[betaxL], reg_coef), reldif(_b[betaxH], reg_coef)) < ${tol}

* Should also work with panel data
tsset TCHID timevar2
reg L.SAT L.Small_Class i.studcat L.Free_Lunch `tscontrols'
scalar reg_coef = _b[L.Small_Class]
rcr L.SAT L.Small_Class i.studcat L.Free_Lunch `tscontrols'
assert min(reldif(_b[betaxL], reg_coef), reldif(_b[betaxH], reg_coef)) < ${tol}
