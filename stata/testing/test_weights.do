run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test WEIGHTS option
*******************************************************************
quietly gen one = 1
quietly gen two = 2
* Unweighted, for comparison
rcr SAT Small_Class ${controls}
savedresults save unweighted e()
* FW = frequency weights, should lead to double the number of observations
* and smaller SE's
rcr SAT Small_Class ${controls} [fw = two]
assert reldif( e(betaxCI_H)  , 6.111214963474101 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.808877201524171 ) <  ${tol}
* PW = probability weights, should lead to same number of observations
* and same SE's
rcr SAT Small_Class ${controls} [pw = two]
assert reldif( e(betaxCI_H)  , 6.488085264868134 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.259480713112181 ) <  ${tol}
* AW = analytic weights, should lead to same number of observations and same
* SE's
rcr SAT Small_Class ${controls} [aw = two]
assert reldif( e(betaxCI_H)  , 6.488085264868134 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.259480713112181 ) <  ${tol}
* IW = importance weights, should lead to double the number of observations
* and smaller SE's
rcr SAT Small_Class ${controls} [iw = two]
assert reldif( e(betaxCI_H)  , 6.111214963474101 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.808877201524171 ) <  ${tol}
* With weights, but all weights set to one.  This should generate the same
* results as unweighted
quietly rcr SAT Small_Class ${controls} [fw = one]
savedresults compare unweighted e()
quietly rcr SAT Small_Class ${controls} [aw = one]
savedresults compare unweighted e()
