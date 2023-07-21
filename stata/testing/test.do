capture log close
set more off
local tmp = c(current_date)
log using "rcr_certification `tmp'", text replace
*******************************************************************
* RCR_CERTIFICATION.DO
*
* This is a do-file designed to provide a detailed test of the RCR
* program, including its optional arguments and postestimation
* commands. The individual tests below can also be run separately.
********************************************************************
cscript "Testing script for RCR and its postestimation commands"

* Test basic functionality with default options 
do test_basic `0'
* Test functionality with different variables
do test_variables `0'
* Test with factor variables
do test_factor `0'
* Test with time series operators
do test_timeseries `0'
* Test sensitivity to variable scale
do test_scaling `0'
* Test IF option
do test_if `0'
* Test IN option
do test_in `0'
* Test WEIGHTS option
do test_weights `0'
* Test CLUSTER and VCE options
do test_cluster `0'
* Test VCEADJ option
do test_vceadj `0'
* Test LAMBDA option
do test_lambda `0'
* Test LEVEL option
do test_level `0'
* Test CITYPE option
do test_citype `0'
* Test SAVE option
do test_save `0'
* Test DETAILS option
do test_details `0'
* Test RESCALE option
do test_rescale `0'
* Test NLCOM postestimation command
do test_nlcom `0'
* Test TESTNL postestimation command
do test_testnl `0'
* Test ESTAT postestimation commands
do test_estat `0'
* Test ESTIMATES postestimation commands
do test_estimates `0'
* Test PREDICT postestimation command
do test_predict `0'
* Test RCRPLOT postestimation command
do test_rcrplot `0'
* Test TEST_BETAX postestimation command
do test_betax `0'
* Test RCR_CONFIG command
do test_rcr_config `0'

* Test/demonstrate outstanding bugs
do test_bug012 `0'
do test_bug013 `0'
do test_bug014 `0'

log close
