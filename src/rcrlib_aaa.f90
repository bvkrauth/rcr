! RCRLIB module
! 
! System-specific code for RCR program.
!
! This file contains all system-specific procedures and data elements in the 
! RCR program.  This particular version (RCRLIB_AAA.F90) is a template
! for users who might need to compile their own versions of the program.
! 

module rcrlib
	implicit none
	integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99) 
	integer, parameter :: stderr=0
	public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr

contains

subroutine sysv(a,b) 
  real(kind=DP), dimension(:,:), intent(inout) :: a,b
end subroutine sysv

elemental function is_finite(x)
	real(kind=DP), intent(in) :: x
	logical :: is_finite
	! This function returns .TRUE. if x is finite and .FALSE. if it is not.
    is_finite = .true. 
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

end module rcrlib
