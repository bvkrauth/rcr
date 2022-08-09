# rcr/r: RCRBOUNDS package for R

This folder contains the R **rcrbounds** package.  The package
is currently under development, and can change before general
release.  I do not anticipate major changes.

## Installation

### Current (development) release

You can install this package by executing:
```
remotes::install_github("bvkrauth/rcr/r/rcrbounds")
```
from the R command line.

## Usage

Once installed, the **rcrbounds** package can be loaded in R
with the statement:
```
library(rcrbounds)
```

The first time you use this package, you should execute the following 
command to ensure that you have all required Python modules:
```
rcrbounds::install_rcrpy()
```

The RCR model can then be estimated using the `rcr()` function. Use:
```
? rcr
```
to get help on this function.


## Support

Please feel free to email me at <bkrauth@sfu.ca> with questions,
bugs, or feature requests.  You can also add bugs or feature
requests as [Github Issues](https://github.com/bvkrauth/rcr/issues).
