! SPDX-FileCopyrightText: 2026 Rohit Goswami
! SPDX-License-Identifier: BSD-3-Clause
! One-direction forward-mode dual arithmetic used by generated CBE sources.
module dual_mod
  use iso_c_binding, only: c_double
  implicit none
  type :: dual
    sequence
    real(c_double) :: re, du
  end type
  interface assignment(=)
    module procedure assign_d, assign_s, assign_i
  end interface
  interface operator(+)
    module procedure add_dd, add_dr, add_rd, add_ds, add_sd, add_di, add_id
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
    module procedure pow_di, pow_dd, pow_ds, pow_dr, pow_sd, pow_rd
  end interface
  interface operator(<)
    module procedure lt_dd, lt_di, lt_id, lt_ds, lt_sd, lt_dr, lt_rd
  end interface
  interface operator(<=)
    module procedure le_dd, le_di, le_id, le_ds, le_sd, le_dr, le_rd
  end interface
  interface operator(>)
    module procedure gt_dd, gt_di, gt_id, gt_ds, gt_sd, gt_dr, gt_rd
  end interface
  interface operator(>=)
    module procedure ge_dd, ge_di, ge_id, ge_ds, ge_sd, ge_dr, ge_rd
  end interface
  interface operator(==)
    module procedure eq_dd, eq_di, eq_id, eq_ds, eq_sd, eq_dr, eq_rd
  end interface
  interface operator(/=)
    module procedure ne_dd, ne_di, ne_id, ne_ds, ne_sd, ne_dr, ne_rd
  end interface
  interface real
    module procedure real_d
  end interface
  interface int
    module procedure int_d
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
  interface cosh
    module procedure cosh_d
  end interface
  interface acos
    module procedure acos_d
  end interface
  interface atan
    module procedure atan_d
  end interface
  interface max
    module procedure max_dd
  end interface
  interface min
    module procedure min_dd
  end interface
  interface sign
    module procedure sign_dd
  end interface
contains
  elemental subroutine assign_d(a,b)
    type(dual),intent(out)::a
    real(c_double),intent(in)::b
    a%re=b; a%du=0d0
  end subroutine
  elemental subroutine assign_i(a,b)
    type(dual),intent(out)::a
    integer,intent(in)::b
    a%re=real(b,c_double); a%du=0d0
  end subroutine
  elemental subroutine assign_s(a,b)
    type(dual),intent(out)::a
    real,intent(in)::b
    a%re=real(b,c_double); a%du=0d0
  end subroutine
  elemental function add_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    c%re=a%re+b%re; c%du=a%du+b%du
  end function
  elemental function add_dr(a,b) result(c)
    type(dual),intent(in)::a; real(c_double),intent(in)::b; type(dual)::c
    c%re=a%re+b; c%du=a%du
  end function
  elemental function add_rd(a,b) result(c)
    real(c_double),intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a+b%re; c%du=b%du
  end function
  elemental function add_ds(a,b) result(c)
    type(dual),intent(in)::a; real,intent(in)::b; type(dual)::c
    c%re=a%re+b; c%du=a%du
  end function
  elemental function add_sd(a,b) result(c)
    real,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a+b%re; c%du=b%du
  end function
  elemental function add_di(a,b) result(c)
    type(dual),intent(in)::a; integer,intent(in)::b; type(dual)::c
    c%re=a%re+b; c%du=a%du
  end function
  elemental function add_id(a,b) result(c)
    integer,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a+b%re; c%du=b%du
  end function
  elemental function sub_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    c%re=a%re-b%re; c%du=a%du-b%du
  end function
  elemental function sub_dr(a,b) result(c)
    type(dual),intent(in)::a; real(c_double),intent(in)::b; type(dual)::c
    c%re=a%re-b; c%du=a%du
  end function
  elemental function sub_rd(a,b) result(c)
    real(c_double),intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a-b%re; c%du=-b%du
  end function
  elemental function sub_ds(a,b) result(c)
    type(dual),intent(in)::a; real,intent(in)::b; type(dual)::c
    c%re=a%re-b; c%du=a%du
  end function
  elemental function sub_sd(a,b) result(c)
    real,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a-b%re; c%du=-b%du
  end function
  elemental function neg_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=-a%re; c%du=-a%du
  end function
  elemental function sub_di(a,b) result(c)
    type(dual),intent(in)::a; integer,intent(in)::b; type(dual)::c
    c%re=a%re-b; c%du=a%du
  end function
  elemental function sub_id(a,b) result(c)
    integer,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a-b%re; c%du=-b%du
  end function
  elemental function mul_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    c%re=a%re*b%re; c%du=a%re*b%du+a%du*b%re
  end function
  elemental function mul_dr(a,b) result(c)
    type(dual),intent(in)::a; real(c_double),intent(in)::b; type(dual)::c
    c%re=a%re*b; c%du=a%du*b
  end function
  elemental function mul_rd(a,b) result(c)
    real(c_double),intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a*b%re; c%du=a*b%du
  end function
  elemental function mul_ds(a,b) result(c)
    type(dual),intent(in)::a; real,intent(in)::b; type(dual)::c
    c%re=a%re*b; c%du=a%du*b
  end function
  elemental function mul_sd(a,b) result(c)
    real,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a*b%re; c%du=a*b%du
  end function
  elemental function mul_di(a,b) result(c)
    type(dual),intent(in)::a; integer,intent(in)::b; type(dual)::c
    c%re=a%re*b; c%du=a%du*b
  end function
  elemental function mul_id(a,b) result(c)
    integer,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a*b%re; c%du=a*b%du
  end function
  elemental function div_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    c%re=a%re/b%re; c%du=(a%du*b%re-a%re*b%du)/(b%re*b%re)
  end function
  elemental function div_dr(a,b) result(c)
    type(dual),intent(in)::a; real(c_double),intent(in)::b; type(dual)::c
    c%re=a%re/b; c%du=a%du/b
  end function
  elemental function div_rd(a,b) result(c)
    real(c_double),intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a/b%re; c%du=-a*b%du/(b%re*b%re)
  end function
  elemental function div_ds(a,b) result(c)
    type(dual),intent(in)::a; real,intent(in)::b; type(dual)::c
    c%re=a%re/b; c%du=a%du/b
  end function
  elemental function div_sd(a,b) result(c)
    real,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a/b%re; c%du=-a*b%du/(b%re*b%re)
  end function
  elemental function div_di(a,b) result(c)
    type(dual),intent(in)::a; integer,intent(in)::b; type(dual)::c
    c%re=a%re/b; c%du=a%du/b
  end function
  elemental function div_id(a,b) result(c)
    integer,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a/b%re; c%du=-a*b%du/(b%re*b%re)
  end function
  elemental function pow_di(a,n) result(c)
    type(dual),intent(in)::a; integer,intent(in)::n; type(dual)::c
    c%re=a%re**n; c%du=real(n,c_double)*a%re**(n-1)*a%du
  end function
  elemental function pow_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    c%re=a%re**b%re
    c%du=c%re*(b%du*log(a%re)+b%re*a%du/a%re)
  end function
  elemental function pow_ds(a,b) result(c)
    type(dual),intent(in)::a; real,intent(in)::b; type(dual)::c
    c%re=a%re**b
    c%du=b*a%re**(b-1)*a%du
  end function
  elemental function pow_dr(a,b) result(c)
    type(dual),intent(in)::a; real(c_double),intent(in)::b; type(dual)::c
    c%re=a%re**b
    c%du=b*a%re**(b-1)*a%du
  end function
  elemental function pow_sd(a,b) result(c)
    real,intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a**b%re
    c%du=c%re*log(real(a,c_double))*b%du
  end function
  elemental function pow_rd(a,b) result(c)
    real(c_double),intent(in)::a; type(dual),intent(in)::b; type(dual)::c
    c%re=a**b%re
    c%du=c%re*log(a)*b%du
  end function
  elemental real(c_double) function real_d(a)
    type(dual),intent(in)::a
    real_d=a%re
  end function
  elemental integer function int_d(a)
    type(dual),intent(in)::a
    int_d=int(a%re)
  end function
  elemental function sqrt_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=sqrt(a%re); c%du=a%du/(2d0*c%re)
  end function
  elemental function exp_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=exp(a%re); c%du=c%re*a%du
  end function
  elemental function log_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=log(a%re); c%du=a%du/a%re
  end function
  elemental function abs_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=abs(a%re); c%du=sign(1d0,a%re)*a%du
  end function
  elemental function sin_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=sin(a%re); c%du=cos(a%re)*a%du
  end function
  elemental function cos_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=cos(a%re); c%du=-sin(a%re)*a%du
  end function
  elemental function tan_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=tan(a%re); c%du=a%du/(cos(a%re)**2)
  end function
  elemental function tanh_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=tanh(a%re); c%du=(1d0-c%re*c%re)*a%du
  end function
  elemental function cosh_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=cosh(a%re); c%du=sinh(a%re)*a%du
  end function
  elemental function acos_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=acos(a%re); c%du=-a%du/sqrt(1d0-a%re*a%re)
  end function
  elemental function atan_d(a) result(c)
    type(dual),intent(in)::a; type(dual)::c
    c%re=atan(a%re); c%du=a%du/(1d0+a%re*a%re)
  end function
  elemental function max_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    if (a%re >= b%re) then; c=a; else; c=b; end if
  end function
  elemental function min_dd(a,b) result(c)
    type(dual),intent(in)::a,b; type(dual)::c
    if (a%re <= b%re) then; c=a; else; c=b; end if
  end function
  elemental logical function lt_dd(a,b)
    type(dual),intent(in)::a,b
    lt_dd=a%re < b%re
  end function
  elemental logical function lt_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    lt_di=a%re < b
  end function
  elemental logical function lt_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    lt_id=a < b%re
  end function
  elemental logical function lt_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    lt_ds=a%re < b
  end function
  elemental logical function lt_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    lt_sd=a < b%re
  end function
  elemental logical function lt_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    lt_dr=a%re < b
  end function
  elemental logical function lt_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    lt_rd=a < b%re
  end function
  elemental logical function le_dd(a,b)
    type(dual),intent(in)::a,b
    le_dd=a%re <= b%re
  end function
  elemental logical function le_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    le_di=a%re <= b
  end function
  elemental logical function le_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    le_id=a <= b%re
  end function
  elemental logical function le_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    le_ds=a%re <= b
  end function
  elemental logical function le_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    le_sd=a <= b%re
  end function
  elemental logical function le_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    le_dr=a%re <= b
  end function
  elemental logical function le_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    le_rd=a <= b%re
  end function
  elemental logical function gt_dd(a,b)
    type(dual),intent(in)::a,b
    gt_dd=a%re > b%re
  end function
  elemental logical function gt_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    gt_di=a%re > b
  end function
  elemental logical function gt_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    gt_id=a > b%re
  end function
  elemental logical function gt_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    gt_ds=a%re > b
  end function
  elemental logical function gt_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    gt_sd=a > b%re
  end function
  elemental logical function gt_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    gt_dr=a%re > b
  end function
  elemental logical function gt_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    gt_rd=a > b%re
  end function
  elemental logical function ge_dd(a,b)
    type(dual),intent(in)::a,b
    ge_dd=a%re >= b%re
  end function
  elemental logical function ge_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    ge_di=a%re >= b
  end function
  elemental logical function ge_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    ge_id=a >= b%re
  end function
  elemental logical function ge_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    ge_ds=a%re >= b
  end function
  elemental logical function ge_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    ge_sd=a >= b%re
  end function
  elemental logical function ge_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    ge_dr=a%re >= b
  end function
  elemental logical function ge_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    ge_rd=a >= b%re
  end function
  elemental logical function eq_dd(a,b)
    type(dual),intent(in)::a,b
    eq_dd=a%re == b%re
  end function
  elemental logical function eq_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    eq_di=a%re == b
  end function
  elemental logical function eq_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    eq_id=a == b%re
  end function
  elemental logical function eq_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    eq_ds=a%re == b
  end function
  elemental logical function eq_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    eq_sd=a == b%re
  end function
  elemental logical function eq_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    eq_dr=a%re == b
  end function
  elemental logical function eq_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    eq_rd=a == b%re
  end function
  elemental logical function ne_dd(a,b)
    type(dual),intent(in)::a,b
    ne_dd=a%re /= b%re
  end function
  elemental logical function ne_di(a,b)
    type(dual),intent(in)::a; integer,intent(in)::b
    ne_di=a%re /= b
  end function
  elemental logical function ne_id(a,b)
    integer,intent(in)::a; type(dual),intent(in)::b
    ne_id=a /= b%re
  end function
  elemental logical function ne_ds(a,b)
    type(dual),intent(in)::a; real,intent(in)::b
    ne_ds=a%re /= b
  end function
  elemental logical function ne_sd(a,b)
    real,intent(in)::a; type(dual),intent(in)::b
    ne_sd=a /= b%re
  end function
  elemental logical function ne_dr(a,b)
    type(dual),intent(in)::a; real(c_double),intent(in)::b
    ne_dr=a%re /= b
  end function
  elemental logical function ne_rd(a,b)
    real(c_double),intent(in)::a; type(dual),intent(in)::b
    ne_rd=a /= b%re
  end function
  elemental function sign_dd(a,b) result(c)
    type(dual),intent(in)::a,b
    type(dual)::c
    c%re=sign(a%re,b%re)
    c%du=sign(a%du,b%re)
  end function
end module
