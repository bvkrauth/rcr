module rcrlib_gnu
    use, intrinsic :: ieee_arithmetic ! requires gfortran version 5.0 or higher
    implicit none
    integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99)
    integer, parameter :: stderr=0
    public SP, DP, is_finite, stderr, initialize

contains

subroutine initialize(infty,nan)
    real(kind=DP), intent(out) :: infty, nan
    infty = 1.0_dp/0.0_dp ! huge(1.0_dp)
    nan = 0.0_dp/0.0_dp
end subroutine initialize

elemental function is_finite(x)
    real(kind=DP), intent(in) :: x
    logical :: is_finite
    is_finite = ieee_is_finite(x) ! This call requires "ieee_arithmetic"
end function is_finite 

end module rcrlib_gnu