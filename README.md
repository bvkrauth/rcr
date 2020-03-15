# rcr - Relative correlation restrictions (Krauth 2016)

This repository contains the RCR Stata package to implement the estimator described in my paper
"Bounding a linear causal effect using relative correlation restrictions" (Journal of Econometric Methods 2016).

It also contains the Fortran source code, and support files for compiling and testing.

## Installing the package in Windows

The easiest way to install the current release of the package is to execute the
Stata command:

  net install rcr, from("http://www.sfu.ca/~bkrauth/code")

I have not yet set the package up to be distributed via github.

## Installing the package in Linux

Tyler Ransom has ported the code to Linux. It is available on his site at https://github.com/tyleransom/rcr.

## Folder structure

The folder structure is

  - stata: Stata package
  - src: Fortran source 
  - bin: binary executables for Windows and Linux
  - test: files for testing
  
Please see the README.md file in each folder for additional details.
