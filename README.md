# rcr - Relative correlation restrictions (Krauth 2016)

This repository contains the RCR Stata package to implement the estimator described in my paper
"Bounding a linear causal effect using relative correlation restrictions" (Journal of Econometric Methods 2016).

It also contains the Fortran source code, and support files for compiling and testing.

## Installing the package in Windows

### Current general release

The current general release of the package can be obtained by executing the
Stata command:

```stata
net install rcr, from("http://www.sfu.ca/~bkrauth/code")
```

To see how the command works, you can call `help rcr`.  You can
also download and execute a sample file by executing

```stata
net install rcr, from("http://www.sfu.ca/~bkrauth/code")
do rcr_example
```

### Current developmental version

The current developmental version can be obtained from this site by executing the Stata command:

```stata
net instalrcr, from("https://raw.githubusercontent.com/bvkrauth/rcr/master/stata/")
```

To see how the command works, you can call `help rcr`.  You can
also download and execute a sample file by executing

```stata
net install rcr, from("https://raw.githubusercontent.com/bvkrauth/rcr/master/stata/")
do rcr_example
```

## Installing the package in Linux

Tyler Ransom has ported the code to Linux. It is available on his site at https://github.com/tyleransom/rcr.

## Folder structure

The folder structure is

  - stata: Stata package
  - src: Fortran source 
  - bin: binary executables for Windows and Linux
  - test: files for testing
  
Please see the README.md file in each folder for additional details.
