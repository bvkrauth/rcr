run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test LEVEL option
*******************************************************************
* Default level (i.e. 95)
rcr SAT Small_Class ${controls}
assert reldif( e(betaxCI_H)  , 6.488085264868137 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.259480713112168 ) <  ${tol}
savedresults save basic e()
* Specify same level as default (should get same result)
rcr SAT Small_Class ${controls}, level(95)
savedresults compare basic e()
* Now replay with a different level
* Doing this will leave the current results in e() unchanged, but report new
* levels and put them in r()
rcr, level(90)
savedresults compare basic e()
assert reldif( r(betaxCI_H)  , 6.281236804837625 ) <  ${tol}
assert reldif( r(betaxCI_L)  , 3.561021633566327 ) <  ${tol}
* Estimate with that level. Should get same result, but now stored in e().
rcr SAT Small_Class ${controls}, level(90)
assert reldif( e(betaxCI_H)  , 6.281236804837625 ) <  ${tol}
assert reldif( e(betaxCI_L)  , 3.561021633566327 ) <  ${tol}
