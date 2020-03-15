# rcr/src: Fortran source code for RCR executable

This directory contains the Fortran source code for the RCR package executable.

## Compiling the Windows executable (RCR.EXE) using the Intel Fortran Compiler

To compile the Windows executable using the Intel Fortran Compiler (IFC) and Math Kernel Library (MKL):

 1. Generate any needed MKL library files (see below)
 2. Open a command window with the appropriate build environment (start menu item "Fortran Build Environment for applications running on IA-32")
 3. Run COMPILE.BAT to compile and link.  The source code files to be compiled are:
    - RCR.F90		The source code (main program)  
    - RCRUTIL.F90		The source code (utility functions)
    - RCRLIB_IFC.F90	The source code (compiler-specific functions for the Intel Fortran Compiler)
 4. You can run a *very* basic test after compilation by simply typing "rcr".
    If no error message, then the program exists and runs.

### Generating the needed MKL library files

This program links in all of its libraries statically so that everything can be distributed as a single executable.  

This means that we need some library files mkl_blas95.lib and mkl_lapack95.lib that do not exist when IFC/MKL is installed.  
They need to be generated. This only needs to be done once for a given installation of the compiler.
Instructions for doing so are copied below from the MKL User Guide.

#### Fortran 95 Interfaces and Wrappers to LAPACK and BLAS (copied from Intel MKL User Guide)

Fortran 95 interfaces are provided for pure procedures and along with wrappers are
delivered as sources. (For more information, see Compiler-dependent Functions and
Fortran 90 Modules). The simplest way to use them is building corresponding libraries and
linking them as user's libraries. 

To do this, you must have administrator rights. [More specifically you must be running the command window as administrator].

Provided the product directory is open for writing, the procedure is simple:

1. Go to the respective directory <mkl_directory>\interfaces\blas95 or
<mkl_directory>\interfaces\lapack95

2. Type one of the following commands:

  - nmake PLAT=win32 lib - for IA-32 architecture
  - nmake PLAT=win32e lib - for Intel® 64 architecture
  - nmake PLAT=win64 lib - for IA-64 architecture.

As a result, the required library and a respective .mod file will be built and installed in the
standard catalog of the release.

## Compiling in Linux using open-source versions of LAPACK and BLAS

Tyler Ransom has ported the RCR code to Linux using open-source tools. 

I have incorporated his changes into my own code and reproduced his documentation below.
You can also access his code directly at https://github.com/tyleransom/rcr. 

### Files

Additional files have been made available for use in compiling on Linux systems using the open-source compiler `gfortran` and 
open-source version of the LAPACK and BLAS libraries.

 - COMPILE.LINUX	Script for compiling program in Linux.
 - RCR_GNU.F90		The source code (main program)  
 - RCRUTIL_GNU.F90	The source code (utility functions)
 - RCRLIB_GNU.F90	The source code (compiler-specific functions for GCC/GFortran compiler and open-source versions of BLAS and LAPACK)

### BLAS & LAPACK libraries (OpenBLAS)

To install the open-source BLAS and LAPACK libraries, do the following:

1. Clone the OpenBLAS GitHub repository: `git clone https://github.com/xianyi/OpenBLAS`

2. Type `make FC=gfortran CC=gcc USE_OPENMP=1` in the directory where you cloned the repository

3. Then change to the OpenBLAS directory and type `make FC=gfortran CC=gcc USE_OPENMP=1 PREFIX=/path/where/you/installed/it/ install`

### BLAS & LAPACK libraries (Netlib)

To install the open-source BLAS and LAPACK libraries, do the following:

1. Type (from your home directory): `wget http://www.netlib.org/lapack/lapack.tgz`

2. Decompress the file: `gunzip -c lapack.tgz | tar xvf -`

3. Remove the compressed file: `rm -f lapack.tgz`

4. Copy and edit the file `lapack-3.7.1/make.inc.example` to `lapack-3.7.1/make.inc`

5. Edit the file `lapack-3.7.1/Makefile` and type `make`. (And wait for 5-10 minutes)

6. (Optional) move the compiled libraries to a different folder, e.g. `~/lib/LAPACK/`

### Test files

To test the validity of the build, execute the script `test.compile` and then execute the file `randomsys1`. Verify that the numerical errors are small (on the order of 1e-10 or smaller)

### Other notes for Linux users

- You will need a version of the GCC compiler 5.0 or above. On some machines, multiple version of the GNU C and FORTRAN compilers are stored in "modules" (e.g. `module load GCC/5.4.0`). You will need to load the appropriate module before executing the `rcr` command.
- You will need to edit the `rcr.ado` file to contain the appropriate paths to the GCC and OpenBLAS libraries. This is handled in the macro `path_to_libs` in the code. An example of the contents of this macro is: "LD_LIBRARY_PATH=/opt/apps/rhel7/gcc-5.4.0/lib64:/lib64/:/hpchome/econ/tmr17/lib/OpenBLAS/lib/"
- You will also need to put the *.ado and *.hlp files in the appropriate places on Stata's search path. N.B. `test_betax.ado` should be put in the `/t/` ado subfolder, while `rcr.ado` should be placed in the `/r/` subfolder.

## Compiling in other settings

If you have access to different tools, or are in a different operating system, you will need to:
  1. Rewrite the RCRLIB module to fit your environment.
     - RCRLIB_AAA.F90		A template for your new RCRLIB module.
     - RCRLIB_IBM.F90		An RCRLIB module I wrote for the IBM compiler.  It is not maintained.
  2. Adapt the compilation script to fit your environment.

 