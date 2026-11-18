run config_test `0'
which rcr
which rcrplot
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*****************************************************************************
* Test BUG014
*****************************************************************************
* Description of bug: rcrplot doesn't work with xrange outside of (-50,50).
*                     The problem was that rcr, details only calculated
*                     lambda(betax) for a grid of betax values between
*                     -50 and 50 (plus +/-1e100).
*
*                     Fix (issue #13): the Python version of rcr now adds
*                     log-spaced tail values of betax (about 1% apart out to
*                     |betax| = 1e6, and sparser beyond) to the details data.
*                     The deprecated Fortran version still has this bug.
rcr SAT Small_Class ${controls}, details
* The lines on this graph should now run all the way from -100 to 100
rcrplot , xrange(-100 100)
if "${exe}" == "python" {
    quietly count if betax > 50 & betax < 100
    assert r(N) > 50
    quietly count if betax < -50 & betax > -100
    assert r(N) > 50
}
