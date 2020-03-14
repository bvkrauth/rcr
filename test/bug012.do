/* BUG012.DO
*
* Description of bug: Standard error results are not invariant to adding constants to data.
*                     This s because the algorithm has trouble getting the derivatives
*                     right when the moments are on substantially different scales.  This
*		          issue is also addressed by bug010.do.  There the problem was that
*			    the program didn't issue a warning message when the standard errors
*		          were bad.  Now it issues a warning message, but it would be nice if 
*			    the derivative algorithm were more robust.
*/
clear
discard
set mem 400m
set more off
capture log close
log using "bug012.log", replace

use rcr_example, clear

/* I've constructed an example to be super-simple */

/* Normalize everything */
foreach vname in SAT Small_Class White_Asian {
	egen tmp1 = mean(`vname')
	egen tmp2 = sd(`vname')
	replace `vname' = (`vname' - tmp1)/tmp2
	drop tmp1 tmp2
}

replace SAT = SAT + 4095.999998
gen SATbad = SAT + 0.000001
rcr SAT Small_Class White_Asian 
est store good
rcr SATbad Small_Class White_Asian 
est store bad
est table _all, b(%9.2f) se keep(lambda0)

log close
