run config_test `0'
which rcr_config
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test RCR_CONFIG command
*******************************************************************
* The result of rcr_config depends on the system setup. So i have
* added some undocumented optional arguments that force a particular
* result for testing purposes.
rcr_config
if r(default_version) == "python" {
    * missing python modules, use fortran if possible
    quietly rcr_config , forceversion(17) forceos("Windows") forcereq("junk")
    assert r(default_version) == "windows-fortran"
    quietly rcr_config , forceversion(17) forceos("Unix") forcereq("junk")
    assert r(default_version) == "unix-fortran"
    quietly rcr_config , forceversion(17) forceos("MacOS") forcereq("junk")
    assert r(default_version) == "none"
}
* version < 16; use fortran if possible
quietly rcr_config , forceversion(12) forceos("Windows")
assert r(default_version) == "windows-fortran"
quietly rcr_config , forceversion(12) forceos("Unix")
assert r(default_version) == "unix-fortran"
quietly rcr_config , forceversion(12) forceos("MacOS")
assert r(default_version) == "none"
* no python installation, use fortran if possible
quietly rcr_config , forceversion(17) forceos("Windows") nopython
assert r(default_version) == "windows-fortran"
quietly rcr_config , forceversion(17) forceos("Unix") nopython
assert r(default_version) == "unix-fortran"
quietly rcr_config , forceversion(17) forceos("MacOS") nopython
assert r(default_version) == "none"
