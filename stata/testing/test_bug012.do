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
*                     the standard errors were bad.  Now it issues a warning
*                     message, but it would be nice if the derivative
*                     algorithm were more robust.
*
* I've constructed an example to be super simple
foreach vname in SAT Small_Class White_Asian {
    egen m = mean(`vname')
    egen sd = sd(`vname')
    replace `vname' = (`vname' - m)/sd
    drop m sd
}
replace SAT = SAT + 4095.999998
gen SATbad = SAT + 0.000001

* Good results
rcr SAT Small_Class White_Asian
estimates store good
savedresults save good e()
* Bad results with large standard errors
rcr SATbad Small_Class White_Asian, rescale(no)
estimates store bad
* The point estimates are similar but the covariance matrix is much larger
savedresults compare good e(), tol(1e-3) ///
    exclude(macros: depvar _estimates_name matrix: V)
* This bug has been partially solved (Jul 2023) by adding a new RESCALE
* option. This option allows the user to automatically rescale the data to
* zero mean and unit variance. Currently, rescaling only occurs if the user
* specifically requests it, but it may be made the default at a later date.
rcr SATbad Small_Class White_Asian, rescale(yes)
estimates store rescaled
savedresults compare good e(), tol(1e-3) exclude(macros: depvar _estimates_name)
* These should all produce the same result, but note that the "bad"
* regression has much larger standard errors.
estimates table _all, b(%9.2f) se keep(lambda0)
