@rem This is the script for compiling the RCR program in Windows

@rem \progfile\g95\bin\g95 rcr.f90 -o rcr.exe
@rem ifort -o rcr.exe rcr.f90
@rem copy rcr.exe "C:\Work\Research\CMU\Stata program\rcr.exe"


@rem ifort /MT rcr.f90 mkl_blas95.lib mkl_intel_c.lib mkl_core.lib libguide.lib
ifort /MT /fltconsistency -o rcr.exe rcrlib_ifc.f90 rcrutil.f90 rcr.f90 mkl_blas95.lib mkl_lapack95.lib mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libguide.lib 

del *.mod *.obj
