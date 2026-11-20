! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause

! Integer powers of a jet at zero are exact through third order, including the
! exponents whose derivative factors vanish at the origin.
program test_jet3_powers
  use jet3_mod, only: jet3, operator(**)
  use iso_c_binding, only: c_double
  implicit none
  type(jet3) :: a, c
  real(c_double) :: actual(12)
  real(c_double), parameter :: expected(12) = [ &
    1d0, 0d0, 0d0, 0d0, &
    0d0, 1d0, 5d0, 7d0, &
    0d0, 0d0, 2d0, 30d0]

  a = jet3(0d0, 1d0, 5d0, 7d0)
  c = a ** 0
  actual(1:4) = [c%re, c%d1, c%d2, c%d3]
  c = a ** 1
  actual(5:8) = [c%re, c%d1, c%d2, c%d3]
  c = a ** 2
  actual(9:12) = [c%re, c%d1, c%d2, c%d3]
  if (any(actual /= expected)) then
    print *, actual
    error stop 1
  end if
end program
