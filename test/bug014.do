/* BUG014.DO
*
* Description of bug: rcrplot doesn't work with xrange outside of (-50,50).
*                     The problem is that rcr, details only calculates lambda(betax)
*                     for a grid of betax values between -50 and 50. 
*
*                     Possible solutions:
*                       (1) Don't accept values of xrange outside that range
*                       (2) Change rcr to calculate lambda(betax) for some 
*                           tail values and interpolate between them.
*/
clear
discard
set mem 400m
set more off
capture log close
adopath ++ ../stata
log using "bug014.log", replace

use "http://www.sfu.ca/~bkrauth/code/rcr_example.dta", clear

rcr SAT Small_Class White_Asian Girl Free_Lunch White_Teacher Teacher_Experience Masters_Degree, details
rcrplot , xrange(-100 100)


log close
