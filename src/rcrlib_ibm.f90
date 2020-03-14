! System-specific code for RCR program.
! This version is for use with the IBM fortran compiler and IBM ESSL library

module rcrlib
	use, intrinsic :: ieee_arithmetic ! IBM-specific call
	implicit none
	integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99)
	integer, parameter :: stderr=0
	public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr

contains

subroutine sysv(a,b) 
  real(kind=DP), dimension(:,:), intent(inout) :: a,b
  integer :: n, nrhs, lda, ldb, ipvt(size(b,1)), info
  n = size(b,1)
  nrhs = size(b,2)
  lda = n
  ipvt = 0
  ldb = n
  info = -1
  call dgesv(n,nrhs,a,lda,ipvt,b,ldb,info) ! IBM-specific call
end subroutine sysv

elemental function is_finite(x)
	real(kind=DP), intent(in) :: x
	logical :: is_finite
    is_finite = ieee_is_finite(x) ! IBM-specific call
end function is_finite 

elemental function is_positive_normal(x)
	real(kind=DP), intent(in) :: x
	logical :: is_positive_normal
	is_positive_normal = (ieee_class(x)== IEEE_POSITIVE_NORMAL)	! IBM-specific call
end function is_positive_normal

function sort(array) 
	real(kind=DP), dimension(:), intent(in) :: array
	real(kind=DP), dimension(size(array)) :: sort
	integer :: n
	sort = array
	n=size(array)
	if (DP == 8) then
        call qsort_up(sort,n,DP) ! IBM-specific call
		if (n /= size(array)) then 
			write (*,*) "Error trying to sort"
		end if
	else
		write (*,*) "Wrong type"
	end if
end function sort


end module rcrlib
