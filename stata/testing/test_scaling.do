run config_test `0'
which rcr
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test effect of variable scale
*******************************************************************
* Generate rescaled variables
quietly gen SAT1000 = SAT * 1000
quietly gen Small_Class1000 = Small_Class * 1000
quietly gen SAT0001 = SAT * 0.0001
quietly gen Small_Class0001 = Small_Class * 0.0001

* Rescaling both, should leave results roughly unchanged
rcr SAT1000 Small_Class1000 ${controls}
assert reldif( e(betaxCI_H)  , 6.488093355120499 ) <  ${tol}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , 3.259476611208612 ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , 3.259476500967557 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , 3.259476791718824 ) <  ${tol}
}
rcr SAT0001 Small_Class0001 ${controls}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , 3.259111536946617 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6.488083259649406 ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , 3.259481216660567 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6.488085194193078 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , 3.25948034557183  ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6.488085407673287 ) <  ${tol}
}

* Scaling outcome up or treatment down, either should multiply coefficient
* by 1000.
rcr SAT1000 Small_Class ${controls}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , 3259.480923495696 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6488.085335644481 ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , 3259.481232785231 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6488.085284518956 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , 3259.481445798481 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 6488.085263438398 ) <  ${tol}
}

* This example produces bad standard errors
rcr SAT Small_Class0001 ${controls}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , 32599.0544917943  ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 64878.0395553112  ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , 32607.22424757308 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 64871.83283612684 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , 32599.04261071111 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , 64871.8357401337  ) <  ${tol}
}

* Scaling outcome down or treatment up, either should multiply coefficient
* by 0.0001
rcr SAT0001 Small_Class ${controls}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , .0003259480936204 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .00064880852586   ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , .000325948074516  ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .0006488085303737 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , .00032594807454   ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .0006488085303706 ) <  ${tol}
}
rcr SAT Small_Class1000 ${controls}
if inlist("`exe'","python") {
    assert reldif( e(betaxCI_L)  , .0032594806861852 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .0064880852480392 ) <  ${tol}
}
else if inlist("`exe'","windows-fortran") {
    assert reldif( e(betaxCI_L)  , .0032594807011477 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .0064880852473763 ) <  ${tol}
}
else if inlist("`exe'","unix-fortran") {
    assert reldif( e(betaxCI_L)  , .0032594806689177 ) <  ${tol}
    assert reldif( e(betaxCI_H)  , .00648808524812   ) <  ${tol}
}
