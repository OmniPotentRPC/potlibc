
program main
use pes_shell
implicit none

  real :: v, v0
  logical :: havegrad
  character (len=1) :: symb
  integer :: i, j, natm, numero
  real, dimension (:,:), allocatable :: xyz
  real, dimension (:), allocatable :: grad, g0
  character(len=32)::fname, bname
  real,parameter::aucm = 219474.63
  real,parameter::auang= 0.5291772083d0

  havegrad = .true.

  call getarg(1,fname)
  i=index(fname,'.',.true.)
  bname=fname(1:i-1)

  open(12,file=trim(fname),status='old')
  open(13,file=trim(bname)//'.out',status='unknown')

  call pes_init()
  read (12,*) natm
  read (12,*) v0

  allocate(xyz(3,natm))
  if (havegrad) allocate(grad(3*natm),g0(3*natm))

  do i=1,natm
     if (havegrad) then
        read (12,*)symb, xyz(:,i), g0(3*i-2:3*i)
     else
        read (12,*)symb, xyz(:,i)
     end if
  end do

  xyz=xyz/auang
  v=f(xyz)
  write(13,'(A)') "# E_ai (hartree)   E_fit (hartree)"
  write(13,'(2F16.8)') v0, v

  if (havegrad) then
     grad=g(xyz)
     write(13,'(A)') "# Grad_ai (hartree/bohr)  Grad_fit:"
     do i=1,size(grad)
        write(13,'(2F16.8)') g0(i), grad(i)
     end do
     deallocate(g0,grad)
  end if

  deallocate(xyz)
  close (12)
  close (13)

end program
