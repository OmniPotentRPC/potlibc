! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause
program test_fortran_api
  use iso_c_binding, only: c_char, c_double, c_int, c_null_char, c_ptr, c_size_t
  implicit none
  interface
    integer(c_int) function potlibc_create(id, handle) bind(c)
      import c_char, c_int, c_ptr
      character(c_char), intent(in) :: id(*)
      type(c_ptr), intent(out) :: handle
    end function
    subroutine potlibc_destroy(handle) bind(c)
      import c_ptr
      type(c_ptr), value :: handle
    end subroutine
    integer(c_int) function potlibc_evaluate(handle,count,x,e,f) bind(c)
      import c_double, c_int, c_ptr, c_size_t
      type(c_ptr), value :: handle
      integer(c_size_t), value :: count
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: e(*), f(*)
    end function
  end interface
  character(c_char), parameter :: id(8) = [ &
    'c','h','5','2','0','0','8',c_null_char ]
  real(c_double), parameter :: x(18) = [ &
    2.5448986937601746d0,3.0319878740728727d0,2.603957609783651d0, &
    1.1100980488026544d0,1.3613545004038685d0,1.1486133528262976d0, &
    0.013880929393421192d0,0.2086334261632169d0,2.457516169583444d0, &
    2.4611712350631576d0,0.21486419033998002d0,0.07550700261319757d0, &
    -0.09665148160819977d0,2.4937031206965017d0,-0.0782634154328256d0, &
    3.3853544918540113d0,4.100244654765713d0,3.5175358525837166d0 ]
  type(c_ptr) :: handle
  real(c_double) :: energy(1), force(18)
  if (potlibc_create(id,handle) /= 0) error stop "create"
  if (potlibc_evaluate(handle,1_c_size_t,x,energy,force) /= 0) &
    error stop "evaluate"
  if (energy(1) /= 0.023947511165499313_c_double) &
    error stop "CBE parity"
  call potlibc_destroy(handle)
end program
