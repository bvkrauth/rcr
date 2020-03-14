! System-specific code for RCR program.
! This version is for use with the GNU fortran compiler (gfortran)
! Requires gfortran version 5.0 or higher


! Import qsort function from C
module gqsort
  use, intrinsic :: iso_c_binding
  implicit none
  private

  interface
    subroutine qsort(base, nel, width, compar) bind(c, name='qsort')
      import c_size_t, c_int
      implicit none
      type(*), intent(inout) :: base(*)
      integer(c_size_t), value :: nel
      integer(c_size_t), value :: width
      abstract interface
        function compar_iface(a, b) bind(c)
          import c_int, c_ptr
          implicit none
          integer(c_int) compar_iface
          type(c_ptr), value :: a, b
        end function
      end interface
      procedure(compar_iface) compar
    end subroutine
  end interface

  interface gnu_qsort
    module procedure gnu_qsort_int4
    module procedure gnu_qsort_int8
    module procedure gnu_qsort_real4
    module procedure gnu_qsort_real8
  end interface
  public gnu_qsort
contains
  subroutine gnu_qsort_int4(a, nel)
    integer(c_int), intent(inout) :: a(*)
    integer(4), value :: nel
    call qsort(a, int(nel, c_size_t), c_sizeof(a(1)), less_int4)
  end subroutine

  subroutine gnu_qsort_int8(a, nel)
    integer(c_long_long), intent(inout) :: a(*)
    integer(4), value :: nel
    call qsort(a, int(nel, c_size_t), c_sizeof(a(1)), less_int8)
  end subroutine

  subroutine gnu_qsort_real4(a, nel)
    real(c_float), intent(inout) :: a(*)
    integer(4), value :: nel
    call qsort(a, int(nel, c_size_t), c_sizeof(a(1)), less_real4)
  end subroutine

  subroutine gnu_qsort_real8(a, nel)
    real(c_double), intent(inout) :: a(*)
    integer(4), value :: nel
    call qsort(a, int(nel, c_size_t), c_sizeof(a(1)), less_real8)
  end subroutine

  function less_int4(a, b) result(result)
    integer(c_int) result
    type(c_ptr), value :: a, b
    integer(c_int), pointer :: ap, bp
    call c_f_pointer(a, ap)
    call c_f_pointer(b, bp)
    result = int(ap - bp, c_int)
  end function

  function less_int8(a, b) result(result)
    integer(c_int) result
    type(c_ptr), value :: a, b
    integer(c_long_long), pointer :: ap, bp
    call c_f_pointer(a, ap)
    call c_f_pointer(b, bp)
    result = int(ap - bp, c_int)
  end function

  function less_real4(a, b) result(result)
    integer(c_int) result
    type(c_ptr), value :: a, b
    real(c_float), pointer :: ap, bp
    call c_f_pointer(a, ap)
    call c_f_pointer(b, bp)
    result = int(ap - bp, c_int)
  end function

  function less_real8(a, b) result(result)
    integer(c_int) result
    type(c_ptr), value :: a, b
    real(c_double), pointer :: ap, bp
    call c_f_pointer(a, ap)
    call c_f_pointer(b, bp)
    result = int(ap - bp, c_int)
  end function
end module



! functions required for RCR procedure
module rcrlib_gnu
    use, intrinsic :: ieee_arithmetic ! requires gfortran version 5.0 or higher
    use gqsort
    implicit none
    integer, parameter :: SP=kind(1.0), DP=selected_real_kind(9,99)
    integer, parameter :: stderr=0
    public SP, DP, sysv, is_finite, sort, is_positive_normal, stderr, initialize

contains

subroutine initialize(infty,nan)
    real(kind=DP), intent(out) :: infty, nan
    infty = huge(1.0_dp)+100
    nan = infty-infty
end subroutine initialize

subroutine sysv(a,b)
  real(kind=DP), dimension(:,:), intent(inout) :: a,b
  integer :: n, nrhs, lda, ldb, ipvt(size(b,1)), info
  n = size(b,1)
  nrhs = size(b,2)
  lda = n
  ipvt = 0
  ldb = n
  info = -1
  call dgesv(n,nrhs,a,lda,ipvt,b,ldb,info) ! Should call LAPACK dgesv which needs to be built from open source
end subroutine sysv

elemental function is_finite(x)
    real(kind=DP), intent(in) :: x
    logical :: is_finite
    is_finite = ieee_is_finite(x) ! This call requires "ieee_arithmetic"
end function is_finite 

elemental function is_positive_normal(x)
    real(kind=DP), intent(in) :: x
    logical :: is_positive_normal
    is_positive_normal = (ieee_class(x)== IEEE_POSITIVE_NORMAL) ! This call requires "ieee_arithmetic"
end function is_positive_normal

function sort(array)
    real(kind=DP), dimension(:), intent(in) :: array
    real(kind=DP), dimension(size(array)) :: sort
    integer :: n
    sort = array
    n=size(array)
    if (DP == 8) then
        call gnu_qsort(sort,n) ! call from C-imported module above
        if (n /= size(array)) then 
            write (*,*) "Error trying to sort"
        end if
    else
        write (*,*) "Wrong type"
    end if
end function sort

end module rcrlib_gnu
