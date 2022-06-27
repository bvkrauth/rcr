# rcr: Bounding a linear causal effect using relative correlation restrictions

This repository contains the Stata and Python packages to implement the estimator described in my paper
"Bounding a linear causal effect using relative correlation restrictions" (Journal of Econometric Methods 2016).

- The Stata package is called **rcr** and is available in the 
  [stata](stata) folder.
- The Python module is called **rcrbounds**  and is available
  in the [python](python) folder.

The repository also contains Fortran source code for an older version of the
Stata package, and various support files for testing.

## Folder structure

The folder structure is:

  - [python](python): Python module
  - [stata](stata): Stata package
  - [test](test): support files for testing Stata package
  - [bin](bin): binary executables for Fortran version (deprecated)
  - [src](src): Fortran source code (deprecated)
  
Please see the *README.md* file in each of these folders for additional
details.

## Support

Please feel free to email me at <bkrauth@sfu.ca> with questions,
bugs, or feature requests.  You can also add bugs or feature
requests as [Github Issues](https://github.com/bvkrauth/rcr/issues).
