# rcr/python: RCRBOUNDS package for Python

This folder contains the Python **rcrbounds** module.

## Installation

### Current general release

The current general release is available in the [Python Package Index](https://pypi.org/)
and can be installed by executing:
```
pip install rcrbounds
```
from the command line.

## Usage

Once installed, the **rcrbounds** module can be imported in Python
with the statement:
```
import rcrbounds
```

RCR estimation can be done in three steps:

1. Use the `RCR()` function to create an RCR object.
2. Use the `fit()` method of the RCR object to
   fit the model and create an RCRResults object.
3. View parameter estimates with the `params` property
   of the RCRResults object, use the `summary()` method
   to view a full table of results, or use the 
   `rcrplot()` method to view a plot of the results.

See `rcr_example.py` for an example, or execute the Python
statement:
```
rcrbounds.RCR?
```
to see the docstring for the RCR class.

## Support

Please feel free to email me at <bkrauth@sfu.ca> with questions,
bugs, or feature requests.  You can also add bugs or feature
requests as [Github Issues](https://github.com/bvkrauth/rcr/issues).
