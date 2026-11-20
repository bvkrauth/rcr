run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*****************************************************************************
* Test BUG012
*****************************************************************************
* Description of bug: Standard error results are not invariant to adding
*                     constants to data. This is because the algorithm has
*                     trouble getting the derivatives right when the moments
*                     are on substantially different scales.  This issue is
*                     also addressed by bug010.do.  There the problem was
*                     that the program didn't issue a warning message when
*                     the standard errors were bad.
*
*                     Fix (issue #11): rcr now always centers every variable
*                     at its mean before calculating the moment vector. The
*                     estimates and the delta-method covariance matrix are
*                     invariant to this, but the numerical derivatives are
*                     much more accurate. This do-file now checks that
*                     results (including V) are the same with and without
*                     large constants added to the data.
*
* I've constructed an example to be super simple
foreach vname in SAT Small_Class White_Asian {
    egen m = mean(`vname')
    egen sd = sd(`vname')
    replace `vname' = (`vname' - m)/sd
    drop m sd
}
replace SAT = SAT + 4095.999998
* Use double so SATbad is exactly SAT plus a constant
gen double SATbad = SAT + 0.000001
gen double SATbig = SAT + 1e6
gen double Small_Class_big = Small_Class + 1e4

* Baseline results
rcr SAT Small_Class White_Asian
estimates store good
savedresults save good e()
* Results should be the same after adding constants
rcr SATbad Small_Class White_Asian, rescale(no)
estimates store bad
savedresults compare good e(), tol(1e-6) ///
    exclude(macros: depvar _estimates_name)
rcr SATbig Small_Class_big White_Asian, rescale(no)
estimates store big
savedresults compare good e(), tol(1e-6) ///
    exclude(macros: depvar treatvar _estimates_name)
rcr SATbad Small_Class White_Asian, rescale(yes)
estimates store rescaled
savedresults compare good e(), tol(1e-6) ///
    exclude(macros: depvar _estimates_name)
* These should all produce the same result
estimates table _all, b(%9.2f) se keep(lambda0)
