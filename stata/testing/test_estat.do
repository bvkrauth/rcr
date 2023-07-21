run config_test `0'
which rcr
which rcr_estat
di "Data from file ${fname}" _newline ///
   "Parameter values for tests: os = ${os}, exe = ${exe}, tol = ${tol}"
*******************************************************************
* Test ESTAT postestimation commands
*******************************************************************
rcr SAT Small_Class ${controls}

* ESTAT VCE
estat vce
assert mreldif( e(V) , r(V)) < ${tol}
* ESTAT VCE with correlation option
estat vce, correlation
qui {
mat T_V = J(5,5,0)
mat T_V[1,1] =                  1
mat T_V[1,2] =  .0261731723888593
mat T_V[1,3] =  .0652619421662122
mat T_V[1,4] =  .0130565284494662
mat T_V[1,5] =  .0107527886530741
mat T_V[2,1] =  .0261731723888593
mat T_V[2,2] =                  1
mat T_V[2,3] = -.9951950739017771
mat T_V[2,4] = -.7122189402928949
mat T_V[2,5] =  .0047083896415582
mat T_V[3,1] =  .0652619421662122
mat T_V[3,2] = -.9951950739017771
mat T_V[3,3] =                  1
mat T_V[3,4] =  .7349449348545146
mat T_V[3,5] =  .0293855113145831
mat T_V[4,1] =  .0130565284494662
mat T_V[4,2] = -.7122189402928949
mat T_V[4,3] =  .7349449348545146
mat T_V[4,4] =                  1
mat T_V[4,5] =  .6981696999117305
mat T_V[5,1] =  .0107527886530741
mat T_V[5,2] =  .0047083896415582
mat T_V[5,3] =  .0293855113145831
mat T_V[5,4] =  .6981696999117305
mat T_V[5,5] =                  1
}
matrix C_V = r(V)
assert mreldif( C_V , T_V ) < ${tol}
_assert_streq `"`: rowfullnames C_V'"' `"lambdaInf betaxInf lambda0 betaxL betaxH"'
_assert_streq `"`: colfullnames C_V'"' `"lambdaInf betaxInf lambda0 betaxL betaxH"'
mat drop C_V T_V

* ESTAT SUMMARIZE
estat summarize
tempname T_stats
mat `T_stats' = J(8,4,0)
mat `T_stats'[1,1] =  51.44962197971362
mat `T_stats'[1,2] =  23.29444665859316
mat `T_stats'[1,3] = -15.76239002946902
mat `T_stats'[1,4] =    127.28907462949
mat `T_stats'[2,1] =  .3024490494947765
mat `T_stats'[2,2] =  .4537342203017396
mat `T_stats'[2,3] = - .1323335592008756
mat `T_stats'[2,4] =   1.13789208746946
mat `T_stats'[3,1] =  .6723754067477308
mat `T_stats'[3,2] =  .2394896046268605
mat `T_stats'[3,3] = - .3176245932522692
mat `T_stats'[3,4] =  1.656502390874715
mat `T_stats'[4,1] =  .4874122281212536
mat `T_stats'[4,2] =  .4965590599853063
mat `T_stats'[4,3] = - .1860571596338484
mat `T_stats'[4,4] =  1.187412228121254
mat `T_stats'[5,1] =  .4836444596677513
mat `T_stats'[5,2] =  .4164020368805063
mat `T_stats'[5,3] = - .4973079212846296
mat `T_stats'[5,4] =  1.468492944516236
mat `T_stats'[6,1] =  .8405548895358794
mat `T_stats'[6,2] =  .2504450623319503
mat `T_stats'[6,3] = - .0344451104641206
mat `T_stats'[6,4] =  1.662777111758102
mat `T_stats'[7,1] =  9.266997773591369
mat `T_stats'[7,2] =  5.152943131556667
mat `T_stats'[7,3] = - 2.399668893075297
mat `T_stats'[7,4] =  25.05647145780189
mat `T_stats'[8,1] =  .3519438259976023
mat `T_stats'[8,2] =  .3908864114709845
mat `T_stats'[8,3] = - .4980561740023977
mat `T_stats'[8,4] =  1.242700128518611
tempname C_stats
matrix `C_stats' = r(stats)
assert mreldif( `C_stats' , `T_stats' ) < ${tol}
_assert_streq `"`: rowfullnames `C_stats''"' `"SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree"'
_assert_streq `"`: colfullnames `C_stats''"' `"mean sd min max"'
mat drop `C_stats' `T_stats'

* ESTAT IC
* This should produce an error message
rcof "noisily estat ic" == 321
