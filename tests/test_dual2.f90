! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause
program test_dual2
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dual2_mod
  implicit none

  type(dual2) :: x, y, z

  x = dual2(0d0, 1d0, -2d0, 3d0)
  z = sqrt(x)
  call expect(z, 0d0, 0d0, 0d0, 0d0, "sqrt zero subgradient")

  x = dual2(-0d0, 1d0, -2d0, 3d0)
  z = abs(x)
  call expect(z, 0d0, 0d0, 0d0, 0d0, "abs signed-zero subgradient")

  x = dual2(-2d0, 3d0, 5d0, 7d0)
  y = dual2(1d0, 0d0, 0d0, 0d0)
  z = sign(x, y)
  call expect(z, 2d0, -3d0, -5d0, -7d0, "sign positive selector")
  y%re = -1d0
  z = sign(x, y)
  call expect(z, -2d0, 3d0, 5d0, 7d0, "sign negative selector")

  x = dual2(-2d0, 1d0, 2d0, 3d0)
  z = x**2d0
  call expect(z, 4d0, -4d0, -8d0, -8d0, "integral real power")
  x = dual2(0d0, 1d0, 2d0, 3d0)
  z = x**0
  call expect(z, 1d0, 0d0, 0d0, 0d0, "zero power")

  x = dual2(1d0, 1d0, 2d0, 3d0)
  z = min(1d0, x)
  call expect(z, 1d0, 0d0, 0d0, 0d0, "min constant-bound tie")
  x%re = -1d0
  z = max(-1d0, x)
  call expect(z, -1d0, 0d0, 0d0, 0d0, "max constant-bound tie")

  x = dual2(1000d0, 0d0, 0d0, 0d0)
  z = exp(x)
  if (ieee_is_finite(z%re)) error stop "exp overflow primal hidden"
  if (z%a /= 0d0 .or. z%b /= 0d0 .or. z%ab /= 0d0) &
    error stop "exp constant tangent is nonzero"

contains

  subroutine expect(value, re, a, b, ab, label)
    type(dual2), intent(in) :: value
    double precision, intent(in) :: re, a, b, ab
    character(len=*), intent(in) :: label

    if (value%re /= re .or. value%a /= a .or. value%b /= b .or. &
        value%ab /= ab) then
      write(*, '(a,4(1x,es24.16))') trim(label), value%re, value%a, &
        value%b, value%ab
      error stop "dual edge-case mismatch"
    end if
  end subroutine expect
end program test_dual2
