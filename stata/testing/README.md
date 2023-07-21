# rcr/test: Files for testing Stata RCR package

This folder contains files that can be used to test the
Stata RCR package.  

Tests for the Python module are located in 
[../python/tests](../python/tests).

## Tests available

### Testing the full Stata package 

*rcr_certification.do* is a Stata do-file that thoroughly tests
the entire Stata package.  To execute the tests, open
that file in Stata and then execute the Stata command:
```stata
do rcr_certification
```
The do-file will run through various test cases, and will issue an
error if any of these cases produce an unexpected result.

### Bugs/issues

The following Stata do-files reproduce a current bug/issue:

 - *test_bug012.do*: Standard error results are not invariant to
   adding constants to data.
 - *test_bug013.do*: Standard error results are inaccurate when 
    lambda_L = -infty.
 - *test_bug014.do*: rcrplot doesn't work with xrange outside of 
   (-50,50).



