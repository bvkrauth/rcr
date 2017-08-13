! System-specific code for RCR program.
! This version is for use with the GNU fortran compiler (gfortran)

! include "lapack.f90" ... 
! Do I need something like this, or is gfortran smart enough to find the dgesv() function on its own???

! module rcrlib_gnu
	! implicit none
	! integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99) 
	! integer, parameter :: stderr=0
	! public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr, initialize


module rcrlib_gnu
	! use, intrinsic :: ieee_arithmetic ! IBM-specific call
	implicit none
	integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99)
	integer, parameter :: stderr=0
	public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr

contains

! subroutine initialize(infty,nan)
	! real(kind=DP), intent(out) :: infty, nan
	! infty = 1.0_dp/0.0_dp
	! nan = 0.0_dp/0.0_dp
! end subroutine initialize

subroutine sysv(a,b) 
  real(kind=DP), dimension(:,:), intent(inout) :: a,b
  integer :: n, nrhs, lda, ldb, ipvt(size(b,1)), info
  n = size(b,1)
  nrhs = size(b,2)
  lda = n
  ipvt = 0
  ldb = n
  info = -1
  call dgesv(n,nrhs,a,lda,ipvt,b,ldb,info) ! Should call LAPACK dgesv which is included in gfortran compiler
end subroutine sysv

elemental function is_finite(x)
	real(kind=DP), intent(in) :: x
	logical :: is_finite
	double precision :: infinity
	! This function returns .TRUE. if x is finite and .FALSE. if it is not.
	infinity = huge(1.0_dp)
	if (x > infinity) then
		is_finite = .false.
	else
		is_finite = .true.
	end if
end function is_finite 

elemental function is_positive_normal(x)
	real(kind=DP), intent(in) :: x
	logical :: is_positive_normal
	is_positive_normal = .true.
end function is_positive_normal

function sort(array) 
	real(kind=DP), dimension(:), intent(in) :: array
	real(kind=DP), dimension(size(array)) :: sort
	integer :: n
	sort = array
end function sort

end module rcrlib_gnu
