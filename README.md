# rcr
This directory contains the source code for the RCR program, along with supporting files for compiling and testing.

Description of files:

## The current program, as downloaded by (Windows) users
RCR.EXE			This is the most recently compiled copy of the executable
RCR.ADO			These are the Stata files implementing the RCR command and its relatives.			
RCR_ESTAT.ADO
RCR_PREDICT.ADO
RCRPLOT.ADO
TEST_THETA.ADO
RCR.HLP
RCR_POSTESTIMATION.HLP
RCRPLOT.HLP
TEST_THETA.HLP
STATA.TOC
RCR.PKG			This Stata file contains the version number and date.
			These must be updated by hand when the version is updated.
RCR_EXAMPLE.DO		This is an example do-file and data set that are provided to users when the RCR package is downloaded.
RCR_EXAMPLE.DTA

## Files needed to compile the Fortran executable in Windows
SETUP.BAT		Sets environment variables needed to run compile.bat
COMPILE.BAT		Compiles the program
RCR.F90			The source code (main program)  
RCRUTIL.F90		The source code (utility functions)
RCRLIB_IFC.F90		The source code (compiler-specific functions for the Intel Fortran Compiler)

## Testing and version control
ARCHIVE.BAT		Once a new version is complete, run (for example, if it's version 0.9.45) 
                        >  ARCHIVE 0.9.45
			and this batch file will create a source code archive and update the local 
			copy of the website (updates to the local copy still need to be uploaded by Dreamweaver).

RCR_CERTIFICATION.DO	This certification script should always be run after any changes to the program.  Read the file to see how to use.

CHANGES_LOG.TXT		Log of all changes made to the program since version 0.9.


BUGnnn.DO		Any do-files with this kind of name demonstrate currently outstanding bugs in the program.  


## Miscellaneous files
RCRLIB_AAA.F90		The source code (compiler-specific functions for user-defined compiler).  Not currently functional.
RCRLIB_IBM.F90		The source code (compiler-specific functions for the IBM compiler).  Not currently functional.
COMPILE.LINUX		Script for compiling program in Linux.  Not currently functional.
CPS.DTA			Not sure if this is used.
INSTALLATION.DO		This do-file shows how to download and update the program from the web.
ZEROTEST.TXT		This can be used to test the executable's handling of zeros.  
			Should be worked into a full bug report and/or certification item, but haven't done that yet.
README.TXT		This file.

There are also various log files, temporary files, and SMCL files that have been generated in testing.  These can usually be deleted.




# A note on linking, and new installations of the compiler

This program links in all of its libraries statically so that everything is in a single executable.  

This means that we need some library files mkl_blas95.lib and mkl_lapack95.lib that do not exist when the Intel Fortran Compiler is installed.  They need to be generated (this only needs to be done once for a given installation of the compiler).  Instructions for doing so are copied below from the MKL User Guide.


## Fortran 95 Interfaces and Wrappers to LAPACK and BLAS

Fortran 95 interfaces are provided for pure procedures and along with wrappers are
delivered as sources. (For more information, see Compiler-dependent Functions and
Fortran 90 Modules). The simplest way to use them is building corresponding libraries and
linking them as user's libraries. 

To do this, you must have administrator rights. [More specifically you must be running the command window as administrator].

Provided the product directory is open for writing, the procedure is simple:

1. Go to the respective directory <mkl_directory>\interfaces\blas95 or
<mkl_directory>\interfaces\lapack95

2. Type one of the following commands:

  nmake PLAT=win32 lib - for IA-32 architecture
  nmake PLAT=win32e lib - for Intel® 64 architecture
  nmake PLAT=win64 lib - for IA-64 architecture.

As a result, the required library and a respective .mod file will be built and installed in the
standard catalog of the release.


# Compiling in Linux using open-source versions of LAPACK and BLAS

Additional files have been made available for use in compiling on Linux systems using the open-source compiler `gfortran` and open-source version of the LAPACK and BLAS libraries.

## BLAS & LAPACK libraries (OpenBLAS)

To install the open-source BLAS and LAPACK libraries, do the following:

1. Clone the OpenBLAS GitHub repository: `git clone https://github.com/xianyi/OpenBLAS`

2. Type `make FC=gfortran CC=gcc USE_OPENMP=1` in the directory where you cloned the repository

3. Then change to the OpenBLAS directory and type `make FC=gfortran CC=gcc USE_OPENMP=1 PREFIX=/path/where/you/installed/it/ install`

## BLAS & LAPACK libraries (Netlib)

To install the open-source BLAS and LAPACK libraries, do the following:

1. Type (from your home directory): `wget http://www.netlib.org/lapack/lapack.tgz`

2. Decompress the file: `gunzip -c lapack.tgz | tar xvf -`

3. Remove the compressed file: `rm -f lapack.tgz`

4. Copy and edit the file `lapack-3.7.1/make.inc.example` to `lapack-3.7.1/make.inc`

5. Edit the file `lapack-3.7.1/Makefile` and type `make`. (And wait for 5-10 minutes)

6. (Optional) move the compiled libraries to a different folder, e.g. `~/lib/LAPACK/`

## Test files

To test the validity of the build, execute the script `test.compile` and then execute the file `randomsys1`. Verify that the numerical errors are small (on the order of 1e-10 or smaller)

