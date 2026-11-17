! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause
! Two-direction forward-mode arithmetic with the mixed epsilon_a*epsilon_b term.
module dual2_mod
  use iso_c_binding, only: c_double
  implicit none

  type :: dual2
    sequence
    real(c_double) :: re, a, b, ab
  end type

  interface assignment(=)
    module procedure assign_r, assign_s, assign_i
  end interface
  interface operator(+)
    module procedure add_dd, add_dr, add_rd, add_ds, add_sd, add_di, add_id, pos_d
  end interface
  interface operator(-)
    module procedure sub_dd, sub_dr, sub_rd, sub_ds, sub_sd, sub_di, sub_id, neg_d
  end interface
  interface operator(*)
    module procedure mul_dd, mul_dr, mul_rd, mul_ds, mul_sd, mul_di, mul_id
  end interface
  interface operator(/)
    module procedure div_dd, div_dr, div_rd, div_ds, div_sd, div_di, div_id
  end interface
  interface operator(**)
    module procedure pow_di, pow_dr, pow_ds, pow_dd, pow_rd, pow_sd
  end interface
  interface operator(<)
    module procedure lt_dd, lt_dr, lt_rd, lt_ds, lt_sd, lt_di, lt_id
  end interface
  interface operator(<=)
    module procedure le_dd, le_dr, le_rd, le_ds, le_sd, le_di, le_id
  end interface
  interface operator(>)
    module procedure gt_dd, gt_dr, gt_rd, gt_ds, gt_sd, gt_di, gt_id
  end interface
  interface operator(>=)
    module procedure ge_dd, ge_dr, ge_rd, ge_ds, ge_sd, ge_di, ge_id
  end interface
  interface operator(==)
    module procedure eq_dd, eq_dr, eq_rd, eq_ds, eq_sd, eq_di, eq_id
  end interface
  interface operator(/=)
    module procedure ne_dd, ne_dr, ne_rd, ne_ds, ne_sd, ne_di, ne_id
  end interface
  interface sqrt
    module procedure sqrt_d
  end interface
  interface exp
    module procedure exp_d
  end interface
  interface log
    module procedure log_d
  end interface
  interface log10
    module procedure log10_d
  end interface
  interface abs
    module procedure abs_d
  end interface
  interface sin
    module procedure sin_d
  end interface
  interface cos
    module procedure cos_d
  end interface
  interface tan
    module procedure tan_d
  end interface
  interface tanh
    module procedure tanh_d
  end interface
  interface sinh
    module procedure sinh_d
  end interface
  interface cosh
    module procedure cosh_d
  end interface
  interface acos
    module procedure acos_d
  end interface
  interface asin
    module procedure asin_d
  end interface
  interface atan
    module procedure atan_d
  end interface
  interface atan2
    module procedure atan2_dd, atan2_dr, atan2_rd, atan2_ds, atan2_sd
  end interface
  interface erf
    module procedure erf_d
  end interface
  interface erfc
    module procedure erfc_d
  end interface
  interface min
    module procedure min_dd, min_dr, min_rd, min_ds, min_sd, min_di, min_id
  end interface
  interface max
    module procedure max_dd, max_dr, max_rd, max_ds, max_sd, max_di, max_id
  end interface
  interface minval
    module procedure minval_d
  end interface
  interface maxval
    module procedure maxval_d
  end interface
  interface sign
    module procedure sign_dd, sign_dr, sign_rd, sign_ds, sign_sd, sign_di, sign_id
  end interface
  interface mod
    module procedure mod_dd, mod_dr, mod_rd, mod_ds, mod_sd, mod_di, mod_id
  end interface
  interface nint
    module procedure nint_d
  end interface
  interface int
    module procedure int_d
  end interface
  interface dble
    module procedure dble_d
  end interface
  interface dot_product
    module procedure dot_product_dd
  end interface
  interface matmul
    module procedure matmul_dd_22, matmul_dd_21, matmul_dd_12
  end interface
  interface real
    module procedure dble_d
  end interface
  interface conjg
    module procedure conjg_d
  end interface

contains
  function minval_d(values) result(value)
    type(dual2), intent(in) :: values(:)
    type(dual2) :: value
    integer :: index
    value = values(1)
    do index = 2, size(values)
      value = min_dd(value, values(index))
    end do
  end function

  function maxval_d(values) result(value)
    type(dual2), intent(in) :: values(:)
    type(dual2) :: value
    integer :: index
    value = values(1)
    do index = 2, size(values)
      value = max_dd(value, values(index))
    end do
  end function

  function dot_product_dd(left, right) result(value)
    type(dual2), intent(in) :: left(:), right(:)
    type(dual2) :: value
    integer :: index
    value = constant_d(0d0)
    do index = 1, size(left)
      value = value + left(index) * right(index)
    end do
  end function

  function matmul_dd_22(left, right) result(value)
    type(dual2), intent(in) :: left(:,:), right(:,:)
    type(dual2) :: value(size(left,1),size(right,2))
    integer :: row, column, inner
    do column = 1, size(right,2)
      do row = 1, size(left,1)
        value(row,column) = constant_d(0d0)
        do inner = 1, size(left,2)
          value(row,column) = value(row,column) + left(row,inner) * right(inner,column)
        end do
      end do
    end do
  end function

  function matmul_dd_21(left, right) result(value)
    type(dual2), intent(in) :: left(:,:), right(:)
    type(dual2) :: value(size(left,1))
    integer :: row, inner
    do row = 1, size(left,1)
      value(row) = constant_d(0d0)
      do inner = 1, size(left,2)
        value(row) = value(row) + left(row,inner) * right(inner)
      end do
    end do
  end function

  function matmul_dd_12(left, right) result(value)
    type(dual2), intent(in) :: left(:), right(:,:)
    type(dual2) :: value(size(right,2))
    integer :: column, inner
    do column = 1, size(right,2)
      value(column) = constant_d(0d0)
      do inner = 1, size(left)
        value(column) = value(column) + left(inner) * right(inner,column)
      end do
    end do
  end function

  elemental function pos_d(value) result(result_value)
    type(dual2), intent(in) :: value
    type(dual2) :: result_value
    result_value = value
  end function


  elemental subroutine assign_r(x, value)
    type(dual2), intent(out) :: x
    real(c_double), intent(in) :: value
    x = dual2(value, 0.0_c_double, 0.0_c_double, 0.0_c_double)
  end subroutine

  elemental subroutine assign_s(x, value)
    type(dual2), intent(out) :: x
    real, intent(in) :: value
    x = dual2(real(value, c_double), 0.0_c_double, 0.0_c_double, 0.0_c_double)
  end subroutine

  elemental subroutine assign_i(x, value)
    type(dual2), intent(out) :: x
    integer, intent(in) :: value
    x = dual2(real(value, c_double), 0.0_c_double, 0.0_c_double, 0.0_c_double)
  end subroutine

  elemental function constant_d(value) result(x)
    real(c_double), intent(in) :: value
    type(dual2) :: x
    x = dual2(value, 0.0_c_double, 0.0_c_double, 0.0_c_double)
  end function

  elemental function unary_d(x, value, first, second) result(y)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: value, first, second
    type(dual2) :: y
    y%re = value
    y%a = first*x%a
    y%b = first*x%b
    y%ab = first*x%ab + second*x%a*x%b
  end function

  elemental function add_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    z = dual2(x%re+y%re, x%a+y%a, x%b+y%b, x%ab+y%ab)
  end function

  elemental function add_dr(x, y) result(z)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: y
    type(dual2) :: z
    z = x; z%re = z%re + y
  end function

  elemental function add_rd(x, y) result(z)
    real(c_double), intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = y; z%re = x + z%re
  end function

  elemental function add_ds(x, y) result(z)
    type(dual2), intent(in) :: x
    real, intent(in) :: y
    type(dual2) :: z
    z = add_dr(x, real(y, c_double))
  end function

  elemental function add_sd(x, y) result(z)
    real, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = add_rd(real(x, c_double), y)
  end function

  elemental function add_di(x, y) result(z)
    type(dual2), intent(in) :: x
    integer, intent(in) :: y
    type(dual2) :: z
    z = add_dr(x, real(y, c_double))
  end function

  elemental function add_id(x, y) result(z)
    integer, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = add_rd(real(x, c_double), y)
  end function

  elemental function sub_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    z = dual2(x%re-y%re, x%a-y%a, x%b-y%b, x%ab-y%ab)
  end function

  elemental function sub_dr(x, y) result(z)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: y
    type(dual2) :: z
    z = x; z%re = z%re - y
  end function

  elemental function sub_rd(x, y) result(z)
    real(c_double), intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = dual2(x-y%re, -y%a, -y%b, -y%ab)
  end function

  elemental function sub_ds(x, y) result(z)
    type(dual2), intent(in) :: x
    real, intent(in) :: y
    type(dual2) :: z
    z = sub_dr(x, real(y, c_double))
  end function

  elemental function sub_sd(x, y) result(z)
    real, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = sub_rd(real(x, c_double), y)
  end function

  elemental function sub_di(x, y) result(z)
    type(dual2), intent(in) :: x
    integer, intent(in) :: y
    type(dual2) :: z
    z = sub_dr(x, real(y, c_double))
  end function

  elemental function sub_id(x, y) result(z)
    integer, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = sub_rd(real(x, c_double), y)
  end function

  elemental function neg_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = dual2(-x%re, -x%a, -x%b, -x%ab)
  end function

  elemental function mul_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    z%re = x%re*y%re
    z%a = x%a*y%re + x%re*y%a
    z%b = x%b*y%re + x%re*y%b
    z%ab = x%ab*y%re + x%a*y%b + x%b*y%a + x%re*y%ab
  end function

  elemental function mul_dr(x, y) result(z)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: y
    type(dual2) :: z
    z = dual2(x%re*y, x%a*y, x%b*y, x%ab*y)
  end function

  elemental function mul_rd(x, y) result(z)
    real(c_double), intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = mul_dr(y, x)
  end function

  elemental function mul_ds(x, y) result(z)
    type(dual2), intent(in) :: x
    real, intent(in) :: y
    type(dual2) :: z
    z = mul_dr(x, real(y, c_double))
  end function

  elemental function mul_sd(x, y) result(z)
    real, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = mul_dr(y, real(x, c_double))
  end function

  elemental function mul_di(x, y) result(z)
    type(dual2), intent(in) :: x
    integer, intent(in) :: y
    type(dual2) :: z
    z = mul_dr(x, real(y, c_double))
  end function

  elemental function mul_id(x, y) result(z)
    integer, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = mul_dr(y, real(x, c_double))
  end function

  elemental function inverse_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, 1.0_c_double/x%re, -1.0_c_double/x%re**2, 2.0_c_double/x%re**3)
  end function

  elemental function div_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    z = mul_dd(x, inverse_d(y))
    z%re = x%re/y%re
    z%a = (x%a*y%re-x%re*y%a)/(y%re*y%re)
    z%b = (x%b*y%re-x%re*y%b)/(y%re*y%re)
  end function

  elemental function div_dr(x, y) result(z)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: y
    type(dual2) :: z
    z = mul_dr(x, 1.0_c_double/y)
    z%re = x%re/y
    z%a = x%a/y
    z%b = x%b/y
    z%ab = x%ab/y
  end function

  elemental function div_rd(x, y) result(z)
    real(c_double), intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = mul_dr(inverse_d(y), x)
    z%re = x/y%re
    z%a = -x*y%a/(y%re*y%re)
    z%b = -x*y%b/(y%re*y%re)
  end function

  elemental function div_ds(x, y) result(z)
    type(dual2), intent(in) :: x
    real, intent(in) :: y
    type(dual2) :: z
    z = div_dr(x, real(y, c_double))
  end function

  elemental function div_sd(x, y) result(z)
    real, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = div_rd(real(x, c_double), y)
  end function

  elemental function div_di(x, y) result(z)
    type(dual2), intent(in) :: x
    integer, intent(in) :: y
    type(dual2) :: z
    z = div_dr(x, real(y, c_double))
  end function

  elemental function div_id(x, y) result(z)
    integer, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = div_rd(real(x, c_double), y)
  end function

  elemental function pow_di(x, n) result(z)
    type(dual2), intent(in) :: x
    integer, intent(in) :: n
    type(dual2) :: z
    if (n == 0) then
      z = constant_d(1.0_c_double)
    else if (n == 1) then
      z = x
    else
      z = unary_d(x, x%re**n, real(n,c_double)*x%re**(n-1), &
                  real(n,c_double)*real(n-1,c_double)*x%re**(n-2))
    end if
  end function

  elemental function pow_dr(x, n) result(z)
    type(dual2), intent(in) :: x
    real(c_double), intent(in) :: n
    type(dual2) :: z
    if (n == 0.0_c_double) then
      z = constant_d(1.0_c_double)
    else if (x%re == 0.0_c_double .and. n > 0.0_c_double) then
      z = constant_d(0.0_c_double)
    else if (n == 1.0_c_double) then
      z = x
    else if (n == anint(n)) then
      z = pow_di(x, int(n))
    else
      z = unary_d(x, x%re**n, n*x%re**(n-1.0_c_double), &
                  n*(n-1.0_c_double)*x%re**(n-2.0_c_double))
    end if
  end function

  elemental function pow_ds(x, n) result(z)
    type(dual2), intent(in) :: x
    real, intent(in) :: n
    type(dual2) :: z
    z = pow_dr(x, real(n, c_double))
  end function

  elemental function pow_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    if (y%a == 0d0 .and. y%b == 0d0 .and. y%ab == 0d0) then
      z = pow_dr(x, y%re)
    else
      z = exp_d(mul_dd(y, log_d(x)))
    end if
  end function

  elemental function pow_rd(x, y) result(z)
    real(c_double), intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    real(c_double) :: value, lx
    value = x**y%re; lx = log(x)
    z = unary_d(y, value, value*lx, value*lx*lx)
  end function

  elemental function pow_sd(x, y) result(z)
    real, intent(in) :: x
    type(dual2), intent(in) :: y
    type(dual2) :: z
    z = pow_rd(real(x, c_double), y)
  end function

  elemental function sqrt_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: value
    if (x%re == 0d0) then
      z = constant_d(0d0)
      return
    end if
    value = sqrt(x%re)
    z = unary_d(x, value, 0.5_c_double/value, -0.25_c_double/value**3)
    z%a = x%a/(2.0_c_double*value)
    z%b = x%b/(2.0_c_double*value)
  end function

  elemental function exp_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: value
    value = exp(x%re)
    if (x%a == 0d0 .and. x%b == 0d0 .and. x%ab == 0d0) then
      z = constant_d(value)
    else
      z = unary_d(x, value, value, value)
    end if
  end function

  elemental function log_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, log(x%re), 1.0_c_double/x%re, -1.0_c_double/x%re**2)
    z%a = x%a/x%re
    z%b = x%b/x%re
  end function

  elemental function log10_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double), parameter :: ln10 = log(10.0_c_double)
    z = unary_d(x, log10(x%re), 1.0_c_double/(x%re*ln10), -1.0_c_double/(x%re**2*ln10))
  end function

  elemental function abs_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    if (x%re > 0.0_c_double) then
      z = x
    else if (x%re < 0.0_c_double) then
      z = neg_d(x)
    else
      z = constant_d(0d0)
    end if
  end function

  elemental function sin_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, sin(x%re), cos(x%re), -sin(x%re))
  end function

  elemental function cos_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, cos(x%re), -sin(x%re), -cos(x%re))
  end function

  elemental function tan_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: value, sec2
    value = tan(x%re); sec2 = 1.0_c_double/cos(x%re)**2
    z = unary_d(x, value, sec2, 2.0_c_double*value*sec2)
    z%a = x%a/(cos(x%re)**2)
    z%b = x%b/(cos(x%re)**2)
  end function

  elemental function tanh_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: value, first
    value = tanh(x%re); first = 1.0_c_double-value**2
    z = unary_d(x, value, first, -2.0_c_double*value*first)
  end function

  elemental function sinh_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, sinh(x%re), cosh(x%re), sinh(x%re))
  end function

  elemental function cosh_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = unary_d(x, cosh(x%re), sinh(x%re), cosh(x%re))
  end function

  elemental function acos_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: den
    den = sqrt(1.0_c_double-x%re**2)
    z = unary_d(x, acos(x%re), -1.0_c_double/den, -x%re/den**3)
    z%a = -x%a/den
    z%b = -x%b/den
  end function

  elemental function asin_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: den
    den = sqrt(1.0_c_double-x%re**2)
    z = unary_d(x, asin(x%re), 1.0_c_double/den, x%re/den**3)
    z%a = x%a/den
    z%b = x%b/den
  end function

  elemental function atan_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: den
    den = 1.0_c_double+x%re**2
    z = unary_d(x, atan(x%re), 1.0_c_double/den, -2.0_c_double*x%re/den**2)
    z%a = x%a/den
    z%b = x%b/den
  end function

  elemental function atan2_dd(y, x) result(z)
    type(dual2), intent(in) :: y, x
    type(dual2) :: z
    real(c_double) :: r2, r4, fy, fx, fyy, fyx, fxx
    r2 = x%re**2 + y%re**2; r4 = r2**2
    fy = x%re/r2; fx = -y%re/r2
    fyy = -2.0_c_double*x%re*y%re/r4
    fyx = (y%re**2-x%re**2)/r4
    fxx = 2.0_c_double*x%re*y%re/r4
    z%re = atan2(y%re, x%re)
    z%a = fy*y%a + fx*x%a
    z%b = fy*y%b + fx*x%b
    z%ab = fy*y%ab + fx*x%ab + fyy*y%a*y%b + &
           fyx*(y%a*x%b+y%b*x%a) + fxx*x%a*x%b
  end function

  elemental function atan2_dr(y, x) result(z)
    type(dual2), intent(in) :: y
    real(c_double), intent(in) :: x
    type(dual2) :: z
    z = atan2_dd(y, constant_d(x))
  end function

  elemental function atan2_rd(y, x) result(z)
    real(c_double), intent(in) :: y
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = atan2_dd(constant_d(y), x)
  end function

  elemental function atan2_ds(y, x) result(z)
    type(dual2), intent(in) :: y
    real, intent(in) :: x
    type(dual2) :: z
    z = atan2_dr(y, real(x, c_double))
  end function

  elemental function atan2_sd(y, x) result(z)
    real, intent(in) :: y
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = atan2_rd(real(y, c_double), x)
  end function

  elemental function erf_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: first
    first = 2.0_c_double/sqrt(acos(-1.0_c_double))*exp(-x%re**2)
    z = unary_d(x, erf(x%re), first, -2.0_c_double*x%re*first)
  end function

  elemental function erfc_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    real(c_double) :: first
    first = -2.0_c_double/sqrt(acos(-1.0_c_double))*exp(-x%re**2)
    z = unary_d(x, erfc(x%re), first, -2.0_c_double*x%re*first)
  end function

  elemental function min_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    if (x%re < y%re) then
      z=x
    else if (y%re < x%re) then
      z=y
    else
      z=x
    end if
  end function

  elemental function min_dr(x, y) result(z)
    type(dual2), intent(in) :: x; real(c_double), intent(in) :: y; type(dual2) :: z
    z=min_dd(x,constant_d(y))
  end function
  elemental function min_rd(x, y) result(z)
    real(c_double), intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=min_dd(constant_d(x),y)
  end function
  elemental function min_ds(x, y) result(z)
    type(dual2), intent(in) :: x; real, intent(in) :: y; type(dual2) :: z
    z=min_dr(x,real(y,c_double))
  end function
  elemental function min_sd(x, y) result(z)
    real, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=min_rd(real(x,c_double),y)
  end function
  elemental function min_di(x, y) result(z)
    type(dual2), intent(in) :: x; integer, intent(in) :: y; type(dual2) :: z
    z=min_dr(x,real(y,c_double))
  end function
  elemental function min_id(x, y) result(z)
    integer, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=min_rd(real(x,c_double),y)
  end function

  elemental function max_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    if (x%re > y%re) then
      z=x
    else if (y%re > x%re) then
      z=y
    else
      z=x
    end if
  end function

  elemental function max_dr(x, y) result(z)
    type(dual2), intent(in) :: x; real(c_double), intent(in) :: y; type(dual2) :: z
    z=max_dd(x,constant_d(y))
  end function
  elemental function max_rd(x, y) result(z)
    real(c_double), intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=max_dd(constant_d(x),y)
  end function
  elemental function max_ds(x, y) result(z)
    type(dual2), intent(in) :: x; real, intent(in) :: y; type(dual2) :: z
    z=max_dr(x,real(y,c_double))
  end function
  elemental function max_sd(x, y) result(z)
    real, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=max_rd(real(x,c_double),y)
  end function
  elemental function max_di(x, y) result(z)
    type(dual2), intent(in) :: x; integer, intent(in) :: y; type(dual2) :: z
    z=max_dr(x,real(y,c_double))
  end function
  elemental function max_id(x, y) result(z)
    integer, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=max_rd(real(x,c_double),y)
  end function

  elemental function sign_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    z = abs_d(x)
    if (y%re < 0.0_c_double) z = neg_d(z)
  end function
  elemental function sign_dr(x, y) result(z)
    type(dual2), intent(in) :: x; real(c_double), intent(in) :: y; type(dual2) :: z
    z=sign_dd(x,constant_d(y))
  end function
  elemental function sign_rd(x, y) result(z)
    real(c_double), intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=sign_dd(constant_d(x),y)
  end function
  elemental function sign_ds(x, y) result(z)
    type(dual2), intent(in) :: x; real, intent(in) :: y; type(dual2) :: z
    z=sign_dr(x,real(y,c_double))
  end function
  elemental function sign_sd(x, y) result(z)
    real, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=sign_rd(real(x,c_double),y)
  end function
  elemental function sign_di(x, y) result(z)
    type(dual2), intent(in) :: x; integer, intent(in) :: y; type(dual2) :: z
    z=sign_dr(x,real(y,c_double))
  end function
  elemental function sign_id(x, y) result(z)
    integer, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=sign_rd(real(x,c_double),y)
  end function

  elemental function mod_dd(x, y) result(z)
    type(dual2), intent(in) :: x, y
    type(dual2) :: z
    real(c_double) :: q
    q = aint(x%re/y%re)
    z = sub_dd(x, mul_dr(y, q))
  end function
  elemental function mod_dr(x, y) result(z)
    type(dual2), intent(in) :: x; real(c_double), intent(in) :: y; type(dual2) :: z
    z=mod_dd(x,constant_d(y))
  end function
  elemental function mod_rd(x, y) result(z)
    real(c_double), intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=mod_dd(constant_d(x),y)
  end function
  elemental function mod_ds(x, y) result(z)
    type(dual2), intent(in) :: x; real, intent(in) :: y; type(dual2) :: z
    z=mod_dr(x,real(y,c_double))
  end function
  elemental function mod_sd(x, y) result(z)
    real, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=mod_rd(real(x,c_double),y)
  end function
  elemental function mod_di(x, y) result(z)
    type(dual2), intent(in) :: x; integer, intent(in) :: y; type(dual2) :: z
    z=mod_dr(x,real(y,c_double))
  end function
  elemental function mod_id(x, y) result(z)
    integer, intent(in) :: x; type(dual2), intent(in) :: y; type(dual2) :: z
    z=mod_rd(real(x,c_double),y)
  end function

  elemental integer function nint_d(x)
    type(dual2), intent(in) :: x
    nint_d = nint(x%re)
  end function
  elemental integer function int_d(x)
    type(dual2), intent(in) :: x
    int_d = int(x%re)
  end function
  elemental real(c_double) function dble_d(x)
    type(dual2), intent(in) :: x
    dble_d = x%re
  end function
  elemental function conjg_d(x) result(z)
    type(dual2), intent(in) :: x
    type(dual2) :: z
    z = x
  end function

  elemental logical function lt_dd(x,y)
    type(dual2),intent(in)::x,y; lt_dd=x%re<y%re
  end function
  elemental logical function lt_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; lt_dr=x%re<y
  end function
  elemental logical function lt_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; lt_rd=x<y%re
  end function
  elemental logical function lt_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; lt_ds=x%re<y
  end function
  elemental logical function lt_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; lt_sd=x<y%re
  end function
  elemental logical function lt_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; lt_di=x%re<y
  end function
  elemental logical function lt_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; lt_id=x<y%re
  end function

  elemental logical function le_dd(x,y)
    type(dual2),intent(in)::x,y; le_dd=x%re<=y%re
  end function
  elemental logical function le_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; le_dr=x%re<=y
  end function
  elemental logical function le_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; le_rd=x<=y%re
  end function
  elemental logical function le_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; le_ds=x%re<=y
  end function
  elemental logical function le_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; le_sd=x<=y%re
  end function
  elemental logical function le_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; le_di=x%re<=y
  end function
  elemental logical function le_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; le_id=x<=y%re
  end function

  elemental logical function gt_dd(x,y)
    type(dual2),intent(in)::x,y; gt_dd=x%re>y%re
  end function
  elemental logical function gt_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; gt_dr=x%re>y
  end function
  elemental logical function gt_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; gt_rd=x>y%re
  end function
  elemental logical function gt_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; gt_ds=x%re>y
  end function
  elemental logical function gt_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; gt_sd=x>y%re
  end function
  elemental logical function gt_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; gt_di=x%re>y
  end function
  elemental logical function gt_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; gt_id=x>y%re
  end function

  elemental logical function ge_dd(x,y)
    type(dual2),intent(in)::x,y; ge_dd=x%re>=y%re
  end function
  elemental logical function ge_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; ge_dr=x%re>=y
  end function
  elemental logical function ge_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; ge_rd=x>=y%re
  end function
  elemental logical function ge_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; ge_ds=x%re>=y
  end function
  elemental logical function ge_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; ge_sd=x>=y%re
  end function
  elemental logical function ge_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; ge_di=x%re>=y
  end function
  elemental logical function ge_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; ge_id=x>=y%re
  end function

  elemental logical function eq_dd(x,y)
    type(dual2),intent(in)::x,y; eq_dd=x%re==y%re
  end function
  elemental logical function eq_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; eq_dr=x%re==y
  end function
  elemental logical function eq_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; eq_rd=x==y%re
  end function
  elemental logical function eq_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; eq_ds=x%re==y
  end function
  elemental logical function eq_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; eq_sd=x==y%re
  end function
  elemental logical function eq_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; eq_di=x%re==y
  end function
  elemental logical function eq_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; eq_id=x==y%re
  end function

  elemental logical function ne_dd(x,y)
    type(dual2),intent(in)::x,y; ne_dd=x%re/=y%re
  end function
  elemental logical function ne_dr(x,y)
    type(dual2),intent(in)::x; real(c_double),intent(in)::y; ne_dr=x%re/=y
  end function
  elemental logical function ne_rd(x,y)
    real(c_double),intent(in)::x; type(dual2),intent(in)::y; ne_rd=x/=y%re
  end function
  elemental logical function ne_ds(x,y)
    type(dual2),intent(in)::x; real,intent(in)::y; ne_ds=x%re/=y
  end function
  elemental logical function ne_sd(x,y)
    real,intent(in)::x; type(dual2),intent(in)::y; ne_sd=x/=y%re
  end function
  elemental logical function ne_di(x,y)
    type(dual2),intent(in)::x; integer,intent(in)::y; ne_di=x%re/=y
  end function
  elemental logical function ne_id(x,y)
    integer,intent(in)::x; type(dual2),intent(in)::y; ne_id=x/=y%re
  end function

  subroutine potlibc_dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
    character, intent(in) :: jobz,uplo
    integer, intent(in) :: n,lda,lwork
    type(dual2), intent(inout) :: a(lda,*),work(*)
    type(dual2), intent(out) :: w(*)
    integer, intent(out) :: info
    type(dual2), allocatable :: original(:,:)
    real(c_double), allocatable :: matrix(:,:),values(:),workspace(:)
    real(c_double) :: da,db,dab,left_a,left_b,denominator
    integer :: i,j,p,q,real_lwork

    if (lwork == -1) then
      work(1)=real(max(1,3*n-1),c_double)
      info=0
      return
    end if
    allocate(original(n,n),matrix(n,n),values(n))
    do q=1,n
      do p=1,n
        if ((uplo == 'U' .or. uplo == 'u') .eqv. (p <= q)) then
          original(p,q)=a(p,q)
        else
          original(p,q)=a(q,p)
        end if
        matrix(p,q)=original(p,q)%re
      end do
    end do
    real_lwork=max(1,lwork)
    allocate(workspace(real_lwork))
    call dsyev(jobz,uplo,n,matrix,n,values,workspace,real_lwork,info)
    if (info /= 0) return
    do i=1,n
      da=0.0_c_double
      db=0.0_c_double
      dab=0.0_c_double
      do q=1,n
        do p=1,n
          da=da+matrix(p,i)*original(p,q)%a*matrix(q,i)
          db=db+matrix(p,i)*original(p,q)%b*matrix(q,i)
          dab=dab+matrix(p,i)*original(p,q)%ab*matrix(q,i)
        end do
      end do
      do j=1,n
        if (j == i) cycle
        denominator=values(i)-values(j)
        if (denominator == 0.0_c_double) cycle
        left_a=0.0_c_double
        left_b=0.0_c_double
        do q=1,n
          do p=1,n
            left_a=left_a+matrix(p,j)*original(p,q)%a*matrix(q,i)
            left_b=left_b+matrix(p,j)*original(p,q)%b*matrix(q,i)
          end do
        end do
        dab=dab+2.0_c_double*left_a*left_b/denominator
      end do
      w(i)=dual2(values(i),da,db,dab)
    end do
    do q=1,n
      do p=1,n
        a(p,q)=dual2(matrix(p,q),0.0_c_double,0.0_c_double,0.0_c_double)
      end do
    end do
  end subroutine

end module dual2_mod
