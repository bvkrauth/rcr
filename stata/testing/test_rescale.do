run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test RESCALE option
*******************************************************************
* RESCALE is optional
rcr SAT Small_Class ${controls}
savedresults save basic e()
* RESCALE(NO) says do not rescale.  It is the default
rcr SAT Small_Class ${controls}, rescale("no")
savedresults compare basic e()
* RESCALE(YES) means that variables should be standardized before
* analysis begins. Coefficients and standard errors will be rescaled
* to match the scale of the original variables.
rcr SAT Small_Class ${controls}, rescale("yes")
* Compare against basic results.  They do not need to be identical, but
* should be roughly similar for this case
savedresults compare basic e(), tol(1e-4)
* Rescale options are not case sensitive
rcr SAT Small_Class ${controls}, rescale("NO")
savedresults compare basic e()
* Only YES and NO are allowed, other options produce an error message
rcof "noisily rcr SAT Small_Class ${controls}, rescale(BAD_OPTION)" == 111
