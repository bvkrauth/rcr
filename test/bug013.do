/* BUG013.DO
*
* Description of bug: Standard error results are inaccurate when lambda_L = -infty.
*/
clear
discard
set mem 400m
set more off
capture log close
log using "bug013.log", replace

use rcr_example, clear


/* I've constructed an example to be super-simple */

rcr SAT Small_Class White_Asian , lambda(. 0)

log close
