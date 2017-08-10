! System-specific code for RCR program.
! This version is for use with the Intel fortran compiler (IFC) and the Intel Math Kernel library

include "mkl_lapack.f90"

module rcrlib
	use mkl95_precision, only : SP,DP
	use mkl95_lapack, only : sysv
	use ifport, only : sortqq, srt$real8
	implicit none
	integer, parameter :: stderr=0
	public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr, initialize


contains

subroutine initialize(infty,nan)
	real(kind=DP), intent(out) :: infty, nan
	infty = 1.0_dp/0.0_dp
	nan = 0.0_dp/0.0_dp
end subroutine initialize

elemental function is_finite(x)
	real(kind=DP), intent(in) :: x
	logical :: is_finite
	is_finite = (fp_class(x) > 3) ! TOO LITERAL
end function is_finite 

elemental function is_positive_normal(x)
	real(kind=DP), intent(in) :: x
	logical :: is_positive_normal
	is_positive_normal = (fp_class(x) == fp_class(1.0_dp))
end function is_positive_normal

function sort(array) 
	real(kind=DP), dimension(:), intent(in) :: array
	real(kind=DP), dimension(size(array)) :: sort
	integer :: n
	sort = array
	n=size(array)
	if (DP == 8) then
		call sortqq(loc(sort),n,SRT$REAL8)
		if (n /= size(array)) then 
			write (*,*) "Error trying to sort"
		end if
	else
		write (*,*) "Wrong type"
	end if
end function sort

! SYSV is in MKL_LAPACK

end module rcrlib
