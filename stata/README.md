# rcr/stata: RCR package for Stata

This folder contains the Stata **rcr** package.

## Installation

### Current general release

The current general release of the package can be obtained by executing the
Stata command:

```stata
net install rcr, from("http://www.sfu.ca/~bkrauth/code")
```

Future general releases will also be posted to that address, so
you will be able to update using Stata's `adoupdate` command.

### Development releases

The current developmental version can be obtained from this site by executing the Stata command:

```stata
net install rcr, from("https://raw.githubusercontent.com/bvkrauth/rcr/master/stata/")
```

Alternatively, you can copy the contents of this folder to 
your Stata working directory.

### Earlier releases

Earlier releases can be found at https://github.com/bvkrauth/rcr/releases.

You can install any of them by unzipping the associated zip file
and copying the contents of the *Stata* folder to your working
directory or any folder in Stata's ado-path (use `adopath` to
find a list of such folders).

## Usage

To see how the command works, you can call `help rcr`.  You can
also download and execute a sample file by executing the
Stata commands:

```stata
net get rcr, from("http://www.sfu.ca/~bkrauth/code")
do rcr_example
```

## Technical notes

The `rcr` command is a wrapper to an external (non-Stata) program
that does the key calculations. There are two versions of this
external program: a Python module, and an older Fortran executable.
When possible, the command will use the Python module.

### Windows

The Python module *rcrbounds.py* will be used if you have:

- Stata version 16 or higher.
- Python installed (and Stata can find it).
- All required Python modules installed.

The Fortran executable *rcr.exe* will be used otherwise.

Use the Stata command `rcr_config` to check your
configuration and get advice on how to fix it.

### Linux

The Python module *rcrbounds.py* will be used if you have:

- Stata version 16 or higher.
- Python installed (and Stata can find it).
- All required Python modules installed.

The Fortran executable *rcr* will be used otherwise. Thanks
to Tyler Ransom for porting the Fortran code to Linux.

If the supplied executable does not work, you may need to
compile from source or configure libraries. Please
see the https://github.com/bvkrauth/rcr/tree/master/src
folder for details and advice.

Use the Stata command `rcr_config` to check your
configuration and get advice on how to fix it.

### MacOS

The Python module *rcrbounds.py* will be used if you have:

- Stata version 16 or higher.
- Python installed (and Stata can find it).
- All required Python modules installed.

Unfortunately, there is no MacOS executable of the Fortran version,
so the requirements above must be satisfied to use the `rcr` command.

Use the Stata command `rcr_config` to check your
configuration and get advice on how to fix it.


## Support

Please feel free to email me at <bkrauth@sfu.ca> with questions,
bugs, or feature requests.  You can also add bugs or feature
requests as [Github Issues](https://github.com/bvkrauth/rcr/issues).
