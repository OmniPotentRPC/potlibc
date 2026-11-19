! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause
! Forward-mode jet through third derivative.
!
! A value carries the function and its first three derivatives with respect
! to one parameter. Arithmetic is the product rule, so the derivatives are
! exact: there is no step. d1 of a gradient is H v, d2 is D³V(v, v) and d3
! is D⁴V(v, v, v).
module jet3_mod
  use iso_c_binding, only: c_int64_t, c_double
  implicit none
  integer(c_int64_t), save :: jet3_ops = 0_c_int64_t
  type :: jet3
    sequence
    real(c_double) :: re
    real(c_double) :: d1
    real(c_double) :: d2
    real(c_double) :: d3
  end type
  interface operator(*)
    module procedure mul_dd, mul_dr, mul_rd, mul_d4, mul_4d
  end interface
  interface operator(/)
    module procedure div_dd, div_dr, div_rd, div_d4, div_4d
  end interface
  interface operator(+)
    module procedure add_dd, add_dr, add_rd, add_d4, add_4d
  end interface
  interface operator(-)
    module procedure sub_dd, sub_dr, sub_rd, sub_d4, sub_4d, neg_d
  end interface
  interface operator(<)
    module procedure lt_dd, lt_dr, lt_rd, lt_d4, lt_4d, lt_di, lt_id
  end interface
  interface operator(<=)
    module procedure le_dd, le_dr, le_rd, le_d4, le_4d, le_di, le_id
  end interface
  interface operator(>)
    module procedure gt_dd, gt_dr, gt_rd, gt_d4, gt_4d, gt_di, gt_id
  end interface
  interface operator(>=)
    module procedure ge_dd, ge_dr, ge_rd, ge_d4, ge_4d, ge_di, ge_id
  end interface
  interface operator(==)
    module procedure eq_dd, eq_dr, eq_rd, eq_d4, eq_4d, eq_di, eq_id
  end interface
  interface operator(/=)
    module procedure ne_dd, ne_dr, ne_rd, ne_d4, ne_4d, ne_di, ne_id
  end interface
  interface operator(**)
    module procedure pow_di
  end interface
  interface assignment(=)
    module procedure assign_real, assign_int
  end interface
  interface real
    module procedure jet3_real
  end interface
  interface sqrt
    module procedure jet3_sqrt
  end interface
  interface exp
    module procedure jet3_exp
  end interface
  interface abs
    module procedure jet3_abs
  end interface
  interface cos
    module procedure jet3_cos
  end interface
  interface sin
    module procedure jet3_sin
  end interface
  interface tan
    module procedure jet3_tan
  end interface
  interface tanh
    module procedure jet3_tanh
  end interface
  interface cosh
    module procedure jet3_cosh
  end interface
  interface acos
    module procedure jet3_acos
  end interface
  interface atan
    module procedure jet3_atan
  end interface
  interface max
    module procedure jet3_max
  end interface
contains
  subroutine jet3_reset()
    jet3_ops = 0_c_int64_t
  end subroutine

  integer(c_int64_t) function jet3_operations()
    jet3_operations = jet3_ops
  end function

  impure elemental subroutine tally()
    jet3_ops = jet3_ops + 1_c_int64_t
  end subroutine

  impure elemental subroutine zero_deriv(c)
    type(jet3), intent(inout) :: c
    c%d1 = 0d0
    c%d2 = 0d0
    c%d3 = 0d0
  end subroutine

  impure elemental function add_dd(a, b) result(c)
    type(jet3), intent(in) :: a, b
    type(jet3) :: c
    call tally()
    c%re = a%re + b%re
    c%d1 = a%d1 + b%d1
    c%d2 = a%d2 + b%d2
    c%d3 = a%d3 + b%d3
  end function

  impure elemental function add_dr(a, b) result(c)
    type(jet3), intent(in) :: a
    real(c_double), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re + b
    c%d1 = a%d1
    c%d2 = a%d2
    c%d3 = a%d3
  end function

  impure elemental function add_rd(a, b) result(c)
    real(c_double), intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a + b%re
    c%d1 = b%d1
    c%d2 = b%d2
    c%d3 = b%d3
  end function

  impure elemental function add_d4(a, b) result(c)
    type(jet3), intent(in) :: a
    real, intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re + b
    c%d1 = a%d1
    c%d2 = a%d2
    c%d3 = a%d3
  end function

  impure elemental function add_4d(a, b) result(c)
    real, intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a + b%re
    c%d1 = b%d1
    c%d2 = b%d2
    c%d3 = b%d3
  end function

  impure elemental function sub_dd(a, b) result(c)
    type(jet3), intent(in) :: a, b
    type(jet3) :: c
    call tally()
    c%re = a%re - b%re
    c%d1 = a%d1 - b%d1
    c%d2 = a%d2 - b%d2
    c%d3 = a%d3 - b%d3
  end function

  impure elemental function sub_dr(a, b) result(c)
    type(jet3), intent(in) :: a
    real(c_double), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re - b
    c%d1 = a%d1
    c%d2 = a%d2
    c%d3 = a%d3
  end function

  impure elemental function sub_rd(a, b) result(c)
    real(c_double), intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a - b%re
    c%d1 = -b%d1
    c%d2 = -b%d2
    c%d3 = -b%d3
  end function

  impure elemental function sub_d4(a, b) result(c)
    type(jet3), intent(in) :: a
    real, intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re - b
    c%d1 = a%d1
    c%d2 = a%d2
    c%d3 = a%d3
  end function

  impure elemental function sub_4d(a, b) result(c)
    real, intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a - b%re
    c%d1 = -b%d1
    c%d2 = -b%d2
    c%d3 = -b%d3
  end function

  impure elemental function neg_d(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    call tally()
    c%re = -a%re
    c%d1 = -a%d1
    c%d2 = -a%d2
    c%d3 = -a%d3
  end function

  impure elemental function mul_dd(a, b) result(c)
    type(jet3), intent(in) :: a, b
    type(jet3) :: c
    call tally()
    c%re = a%re * b%re
    c%d1 = a%re * b%d1 + a%d1 * b%re
    c%d2 = a%d2 * b%re + 2d0 * a%d1 * b%d1 + a%re * b%d2
    c%d3 = a%d3 * b%re + 3d0 * a%d2 * b%d1 + 3d0 * a%d1 * b%d2 + a%re * b%d3
  end function

  impure elemental function mul_dr(a, b) result(c)
    type(jet3), intent(in) :: a
    real(c_double), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re * b
    c%d1 = a%d1 * b
    c%d2 = a%d2 * b
    c%d3 = a%d3 * b
  end function

  impure elemental function mul_rd(a, b) result(c)
    real(c_double), intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a * b%re
    c%d1 = a * b%d1
    c%d2 = a * b%d2
    c%d3 = a * b%d3
  end function

  impure elemental function mul_d4(a, b) result(c)
    type(jet3), intent(in) :: a
    real, intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re * b
    c%d1 = a%d1 * b
    c%d2 = a%d2 * b
    c%d3 = a%d3 * b
  end function

  impure elemental function mul_4d(a, b) result(c)
    real, intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a * b%re
    c%d1 = a * b%d1
    c%d2 = a * b%d2
    c%d3 = a * b%d3
  end function

  impure elemental function div_dd(a, b) result(c)
    type(jet3), intent(in) :: a, b
    type(jet3) :: c
    type(jet3) :: inv
    call tally()
    inv = reciprocal(b)
    c%re = a%re / b%re
    c%d1 = (a%d1 * b%re - a%re * b%d1) / (b%re * b%re)
    c%d2 = a%d2 * inv%re + 2d0 * a%d1 * inv%d1 + a%re * inv%d2
    c%d3 = a%d3 * inv%re + 3d0 * a%d2 * inv%d1 + 3d0 * a%d1 * inv%d2 + a%re * inv%d3
  end function

  impure elemental function div_dr(a, b) result(c)
    type(jet3), intent(in) :: a
    real(c_double), intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re / b
    c%d1 = a%d1 / b
    c%d2 = a%d2 / b
    c%d3 = a%d3 / b
  end function

  impure elemental function div_rd(a, b) result(c)
    real(c_double), intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c, inv
    call tally()
    inv = reciprocal(b)
    c%re = a / b%re
    c%d1 = -a * b%d1 / (b%re * b%re)
    c%d2 = a * inv%d2
    c%d3 = a * inv%d3
  end function

  impure elemental function div_d4(a, b) result(c)
    type(jet3), intent(in) :: a
    real, intent(in) :: b
    type(jet3) :: c
    call tally()
    c%re = a%re / b
    c%d1 = a%d1 / b
    c%d2 = a%d2 / b
    c%d3 = a%d3 / b
  end function

  impure elemental function div_4d(a, b) result(c)
    real, intent(in) :: a
    type(jet3), intent(in) :: b
    type(jet3) :: c, inv
    call tally()
    inv = reciprocal(b)
    c%re = a / b%re
    c%d1 = -a * b%d1 / (b%re * b%re)
    c%d2 = a * inv%d2
    c%d3 = a * inv%d3
  end function

  impure elemental function reciprocal(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: inv, inv2, inv3, inv4
    inv = 1d0 / a%re
    inv2 = inv * inv
    inv3 = inv2 * inv
    inv4 = inv2 * inv2
    c%re = inv
    c%d1 = -a%d1 * inv2
    c%d2 = -a%d2 * inv2 + 2d0 * a%d1 * a%d1 * inv3
    c%d3 = -a%d3 * inv2 + 6d0 * a%d1 * a%d2 * inv3 - 6d0 * a%d1 * a%d1 * a%d1 * inv4
  end function

  impure elemental function pow_di(a, n) result(c)
    type(jet3), intent(in) :: a
    integer, intent(in) :: n
    type(jet3) :: c
    real(c_double) :: base, p2, p3
    if (n == 0) then
      c%re = 1d0
      c%d1 = 0d0
      c%d2 = 0d0
      c%d3 = 0d0
      return
    end if
    if (n == 1) then
      c = a
      return
    end if
    if (n == 2) then
      c = a * a
      return
    end if
    call tally()
    base = a%re ** (n - 1)
    c%re = a%re * base
    c%d1 = real(n, c_double) * base * a%d1
    p2 = real(n * (n - 1), c_double) * a%re ** (n - 2)
    p3 = real(n * (n - 1) * (n - 2), c_double) * a%re ** (n - 3)
    c%d2 = p2 * a%d1 * a%d1 + real(n, c_double) * base * a%d2
    c%d3 = p3 * a%d1 * a%d1 * a%d1 + 3d0 * p2 * a%d1 * a%d2 &
      + real(n, c_double) * base * a%d3
  end function

  impure elemental subroutine assign_real(a, b)
    type(jet3), intent(out) :: a
    real(c_double), intent(in) :: b
    call tally()
    a%re = b
    a%d1 = 0d0
    a%d2 = 0d0
    a%d3 = 0d0
  end subroutine

  impure elemental subroutine assign_int(a, b)
    type(jet3), intent(out) :: a
    integer, intent(in) :: b
    call tally()
    a%re = real(b, c_double)
    a%d1 = 0d0
    a%d2 = 0d0
    a%d3 = 0d0
  end subroutine

  impure elemental function jet3_real(a) result(c)
    type(jet3), intent(in) :: a
    real(c_double) :: c
    c = a%re
  end function

  impure elemental function jet3_sqrt(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: s, inv
    call tally()
    s = sqrt(a%re)
    inv = 1d0 / s
    c%re = s
    c%d1 = a%d1 / (2d0 * c%re)
    c%d2 = 0.5d0 * a%d2 * inv - c%d1 * c%d1 * inv
    c%d3 = 0.5d0 * a%d3 * inv - 0.5d0 * a%d2 * c%d1 * inv * inv &
      - 2d0 * c%d1 * c%d2 * inv + c%d1 * c%d1 * c%d1 * inv * inv
  end function

  impure elemental function jet3_exp(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    call tally()
    c%re = exp(a%re)
    c%d1 = c%re * a%d1
    c%d2 = c%re * (a%d1 * a%d1 + a%d2)
    c%d3 = c%re * (a%d1 * (a%d1 * a%d1 + a%d2) + 2d0 * a%d1 * a%d2 + a%d3)
  end function

  impure elemental function jet3_abs(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: s
    call tally()
    s = sign(1d0, a%re)
    c%re = abs(a%re)
    c%d1 = s * a%d1
    c%d2 = s * a%d2
    c%d3 = s * a%d3
  end function

  impure elemental function jet3_sin(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: s, co
    call tally()
    s = sin(a%re)
    co = cos(a%re)
    c%re = s
    c%d1 = co * a%d1
    c%d2 = co * a%d2 - s * a%d1 * a%d1
    c%d3 = co * a%d3 - co * a%d1 * a%d1 * a%d1 - 3d0 * s * a%d1 * a%d2
  end function

  impure elemental function jet3_cos(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: s, co
    call tally()
    s = sin(a%re)
    co = cos(a%re)
    c%re = co
    c%d1 = -s * a%d1
    c%d2 = -co * a%d1 * a%d1 - s * a%d2
    c%d3 = -s * a%d3 - 3d0 * co * a%d1 * a%d2 + s * a%d1 * a%d1 * a%d1
  end function

  impure elemental function jet3_tan(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: co, sec2, t
    call tally()
    t = tan(a%re)
    co = cos(a%re)
    sec2 = 1d0 / (co * co)
    c%re = t
    c%d1 = a%d1 / (co * co)
    c%d2 = sec2 * a%d2 + 2d0 * t * sec2 * a%d1 * a%d1
    c%d3 = sec2 * a%d3 + 6d0 * t * sec2 * a%d1 * a%d2 + (2d0 * sec2 * sec2 + 4d0 * t * t * sec2) * a%d1 * a%d1 * a%d1
  end function

  impure elemental function jet3_tanh(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: t, u
    call tally()
    t = tanh(a%re)
    u = 1d0 - t * t
    c%re = t
    c%d1 = (1d0 - c%re * c%re) * a%d1
    c%d2 = u * a%d2 - 2d0 * t * u * a%d1 * a%d1
    c%d3 = u * a%d3 - 6d0 * t * u * a%d1 * a%d2 + (4d0 * t * t * u - 2d0 * u * u) * a%d1 * a%d1 * a%d1
  end function

  impure elemental function jet3_cosh(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: ch, sh
    call tally()
    ch = cosh(a%re)
    sh = sinh(a%re)
    c%re = ch
    c%d1 = sh * a%d1
    c%d2 = ch * a%d1 * a%d1 + sh * a%d2
    c%d3 = sh * a%d1 * a%d1 * a%d1 + 3d0 * ch * a%d1 * a%d2 + sh * a%d3
  end function

  impure elemental function jet3_acos(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: s, inv
    call tally()
    s = sqrt(1d0 - a%re * a%re)
    inv = -1d0 / s
    c%re = acos(a%re)
    c%d1 = -a%d1 / sqrt(1d0 - a%re * a%re)
    c%d2 = inv * a%d2 + inv * inv * inv * a%re * a%d1 * a%d1
    c%d3 = inv * a%d3 + 3d0 * inv * inv * inv * a%re * a%d1 * a%d2 &
      + (inv * inv * inv + 3d0 * inv * inv * inv * inv * inv * a%re * a%re) * a%d1 * a%d1 * a%d1
  end function

  impure elemental function jet3_atan(a) result(c)
    type(jet3), intent(in) :: a
    type(jet3) :: c
    real(c_double) :: den, den2
    call tally()
    den = 1d0 + a%re * a%re
    den2 = den * den
    c%re = atan(a%re)
    c%d1 = a%d1 / (1d0 + a%re * a%re)
    c%d2 = a%d2 / den - 2d0 * a%re * a%d1 * a%d1 / den2
    c%d3 = a%d3 / den - 6d0 * a%re * a%d1 * a%d2 / den2 &
      + (8d0 * a%re * a%re / (den2 * den) - 2d0 / den2) * a%d1 * a%d1 * a%d1
  end function

  impure elemental function jet3_max(a, b) result(c)
    type(jet3), intent(in) :: a, b
    type(jet3) :: c
    call tally()
    if (a%re >= b%re) then
      c = a
    else
      c = b
    end if
  end function
  elemental logical function lt_dd(x, y)
    type(jet3), intent(in) :: x, y
    lt_dd = x%re < y%re
  end function

  elemental logical function lt_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    lt_dr = x%re < y
  end function

  elemental logical function lt_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    lt_rd = x < y%re
  end function

  elemental logical function lt_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    lt_d4 = x%re < y
  end function

  elemental logical function lt_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    lt_4d = x < y%re
  end function

  elemental logical function lt_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    lt_di = x%re < y
  end function

  elemental logical function lt_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    lt_id = x < y%re
  end function

  elemental logical function le_dd(x, y)
    type(jet3), intent(in) :: x, y
    le_dd = x%re <= y%re
  end function

  elemental logical function le_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    le_dr = x%re <= y
  end function

  elemental logical function le_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    le_rd = x <= y%re
  end function

  elemental logical function le_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    le_d4 = x%re <= y
  end function

  elemental logical function le_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    le_4d = x <= y%re
  end function

  elemental logical function le_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    le_di = x%re <= y
  end function

  elemental logical function le_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    le_id = x <= y%re
  end function

  elemental logical function gt_dd(x, y)
    type(jet3), intent(in) :: x, y
    gt_dd = x%re > y%re
  end function

  elemental logical function gt_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    gt_dr = x%re > y
  end function

  elemental logical function gt_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    gt_rd = x > y%re
  end function

  elemental logical function gt_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    gt_d4 = x%re > y
  end function

  elemental logical function gt_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    gt_4d = x > y%re
  end function

  elemental logical function gt_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    gt_di = x%re > y
  end function

  elemental logical function gt_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    gt_id = x > y%re
  end function

  elemental logical function ge_dd(x, y)
    type(jet3), intent(in) :: x, y
    ge_dd = x%re >= y%re
  end function

  elemental logical function ge_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    ge_dr = x%re >= y
  end function

  elemental logical function ge_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    ge_rd = x >= y%re
  end function

  elemental logical function ge_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    ge_d4 = x%re >= y
  end function

  elemental logical function ge_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    ge_4d = x >= y%re
  end function

  elemental logical function ge_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    ge_di = x%re >= y
  end function

  elemental logical function ge_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    ge_id = x >= y%re
  end function

  elemental logical function eq_dd(x, y)
    type(jet3), intent(in) :: x, y
    eq_dd = x%re == y%re
  end function

  elemental logical function eq_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    eq_dr = x%re == y
  end function

  elemental logical function eq_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    eq_rd = x == y%re
  end function

  elemental logical function eq_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    eq_d4 = x%re == y
  end function

  elemental logical function eq_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    eq_4d = x == y%re
  end function

  elemental logical function eq_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    eq_di = x%re == y
  end function

  elemental logical function eq_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    eq_id = x == y%re
  end function

  elemental logical function ne_dd(x, y)
    type(jet3), intent(in) :: x, y
    ne_dd = x%re /= y%re
  end function

  elemental logical function ne_dr(x, y)
    type(jet3), intent(in) :: x
    real(c_double), intent(in) :: y
    ne_dr = x%re /= y
  end function

  elemental logical function ne_rd(x, y)
    real(c_double), intent(in) :: x
    type(jet3), intent(in) :: y
    ne_rd = x /= y%re
  end function

  elemental logical function ne_d4(x, y)
    type(jet3), intent(in) :: x
    real, intent(in) :: y
    ne_d4 = x%re /= y
  end function

  elemental logical function ne_4d(x, y)
    real, intent(in) :: x
    type(jet3), intent(in) :: y
    ne_4d = x /= y%re
  end function

  elemental logical function ne_di(x, y)
    type(jet3), intent(in) :: x
    integer, intent(in) :: y
    ne_di = x%re /= y
  end function

  elemental logical function ne_id(x, y)
    integer, intent(in) :: x
    type(jet3), intent(in) :: y
    ne_id = x /= y%re
  end function
end module
