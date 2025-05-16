******************************************************************************
* Configuration script for Stata tests
*
* This do-file should run before each unit test in the Stata/testing folder
******************************************************************************
* Set environment
clear all
adopath ++ ../

* Get test data
* Use local version if it exists
global fname "rcr_example.dta"
capture confirm file ${fname}
* Otherwise use the web version
if (_rc != 0) {
    global fname "http://www.sfu.ca/~bkrauth/code/rcr_example.dta"
}
di "Using data file ${fname}"
use "${fname}", clear
global controls = "White_Asian Girl Free_Lunch White_Teacher " + ///
    "Teacher_Experience Masters_Degree"

* Determine which version is supported on this system
global os = c(os)
capture rcr_config
if (_rc == 0) {
    global exe = r(default_version)
}
else {
    if "${os}" == "Windows" {
        global exe "windows-fortran"
    }
    else if "${os}" == "Unix" {
        global exe "unix-fortran"
    }
}

* Set tolerances for numeric comparisons.  The Windows Python version
* is the "base" version for regular testing.  So we need to catch
* any changes to that even if fairly small.  For the other platforms,
* we only need to make sure the results are close to those for the
* Windows Python version.  So we use a higher relative tolerance.
if c(os) == "Windows" & "${exe}" == "python" {
    global tol "1E-7"
}
else {
    global tol "1E-4"
}

di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
