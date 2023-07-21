run config_test `0'
which rcr
which rcrplot
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test RCRPLOT postestimation command
*******************************************************************
* Use a consistent and simple scheme
set scheme s2mono
* Call RCRPLOT without having specified DETAILS.  Should generate an error
* message
rcr SAT Small_Class ${controls}
rcof "noisily rcrplot" == 321

* Call RCR with DETAILS specified. This also produces a plot
rcr SAT Small_Class ${controls}, details
* Check to make sure that the plot has not changed. This is done by saving
* the plot as a postscript file, and using checksums to compare it to a
* previously-saved version of this plot. Minor changes in Stata will
* sometimes change the exact contents of the postscript file, so a failed
* comparison here does not necessarily mean anything substantive has
* changed. In the event of a failed checksum comparison, visually compare
* the current and previously-saved versions of the plot:
*     - If the two versions look different, find and fix the bug.
*     - If the two versions look the same, delete the previously-saved
*       version and re-run this test to update the file.
capture confirm file "rcrplot1.ps"
* Create comparison file if it does not exist
if (_rc != 0) {
    graph export rcrplot1.ps, as(ps)
}
* Create current plot file
tempfile rcrplot
quietly graph export `rcrplot', as(ps) replace
* Compare the two files
checksum "rcrplot1.ps"
local checksum = r(checksum)
checksum `rcrplot'
if (r(checksum) != `checksum') {
    di as error "Current graph does not exactly match rcrplot1.ps. " ///
        "Please visually compare the two graphs. " ///
        "If they match, delete the current rcrplot1.ps and rerun tests. " ///
        "If they do not match, the test has failed."
    error 1
}

* Call RCRPLOT with arguments
rcrplot, xrange(-20 20) yrange(-40 40)
* Create comparison file if it does not exist
capture confirm file "rcrplot2.ps"
if (_rc != 0) {
    quietly graph export rcrplot2.ps, as(ps)
}
* Create current plot file
quietly graph export `rcrplot', as(ps) replace
* Compare the two files
checksum "rcrplot2.ps"
local checksum = r(checksum)
checksum `rcrplot'
if (r(checksum) != `checksum') {
    di as error "Current graph does not exactly match rcrplot2.ps. " ///
        "Please visually compare the two graphs. " ///
        "If they match, delete the current rcrplot2.ps and rerun tests. " ///
        "If they do not match, the test has failed."
    error 1
}
