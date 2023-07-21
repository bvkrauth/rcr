* BUG013.DO
*
* Description of bug: Standard error results are inaccurate when
* lambda_L = -infty.
clear
discard
set mem 400m
set more off
capture log close
adopath ++ ../stata
log using "bug013.log", replace

use "http://www.sfu.ca/~bkrauth/code/rcr_example.dta", clear

* I've constructed an example to be super simple

rcr SAT Small_Class White_Asian , lambda(. 0)

log close
