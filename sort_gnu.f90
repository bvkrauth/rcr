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
          import c_int
          implicit none
          integer(c_int) compar_iface
          type(*), intent(in) :: a, b
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
    integer(c_int), intent(in) :: a, b
    result = a - b
  end function

  function less_int8(a, b) result(result)
    integer(c_int) result
    integer(c_long_long), intent(in) :: a, b
    result = a - b
  end function

  function less_real4(a, b) result(result)
    integer(c_int) result
    real(c_float), intent(in) :: a, b
    result = a - b
  end function

  function less_real8(a, b) result(result)
    integer(c_int) result
    real(c_double), intent(in) :: a, b
    result = a - b
  end function
end module
