# rcr/test: Files for testing Stata RCR package

This folder contains files that can be used for testing the Stata RCR package

## Testing the RCR executable

To perform a basic test of the executable without using Stata, 
just open a command window in this folder and execute the command
   ../bin/rcr
This will read in the file in.txt, perform the rcr calculations, and output the results to out.txt

If this concludes without error, this means the rcr executable exists and runs.

## Testing the full Stata package 

The file rcr_certification.do is a Stata do-file that thoroughly tests the entire Stata package.

## Bugs/issues

The following Stata do-files reproduce a current bug/issue:

 - bug012.do
 - bug013.do



