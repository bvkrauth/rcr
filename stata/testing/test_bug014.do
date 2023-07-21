run config_test `0'
which rcr
which rcrplot
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*****************************************************************************
* Test BUG014
*****************************************************************************
* Description of bug: rcrplot doesn't work with xrange outside of (-50,50).
*                     The problem is that rcr, details only calculates
*                     lambda(betax) for a grid of betax values between
*                     -50 and 50.
*
*                     Possible solutions:
*                       (1) Don't accept values of xrange outside that range
*                       (2) Change rcr to calculate lambda(betax) for some
*                           tail values and interpolate between them.
rcr SAT Small_Class ${controls}, details
rcrplot , xrange(-100 100)
