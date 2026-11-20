run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*****************************************************************************
* Test BUG013
*****************************************************************************
* Description of bug: Standard error results are inaccurate when
*                     lambda_L = -infty.
*
*                     Fix (issue #10): when the identified set approaches
*                     betaxInf (where lambda(betax) goes to -infinity), the
*                     bound is now exactly betaxInf, with the same gradient.
*                     Previously the bound was a grid point next to betaxInf
*                     and its standard error was badly understated.
*                     The deprecated Fortran version still has this bug.
rcr SAT Small_Class White_Asian , lambda(. 0)
if "${exe}" == "python" {
    assert _b[betaxH] == _b[betaxInf]
    assert reldif(_se[betaxH], _se[betaxInf]) < ${tol}
}
rcr SAT Small_Class ${controls} , lambda(. 1)
if "${exe}" == "python" {
    assert _b[betaxH] == _b[betaxInf]
    assert reldif(_se[betaxH], _se[betaxInf]) < ${tol}
}
