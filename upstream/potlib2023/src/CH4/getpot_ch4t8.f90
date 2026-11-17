subroutine getpot(v,xx)
  implicit none
  integer ::  ifirst=0
  double precision :: v0
  double precision :: v, energy_hartree
  double precision :: xx(5,3), z(10)
!
! This routine takes a set of Cartesian coordinates (xx), converts it to Radau coordinates, 
! and calls the T8 PES of Schwenke and Partridge returning the energy as V.  
! The potential is re-zeroed so that the minimum configuration has an energy of 0 Hartree 
!
! On input the C is the first atom.
!
  if(ifirst.eq.0)then
    ifirst=1
    write(6,*)'Calling T8 potential for CH4 of Schwenke and Partridge'
    write(6,*)'D. W. Schwenke and H. Partridge, Spectrochim. Acta A 57, (2001) 887'
  endif

  call methaneradau(xx,z)
  call vibpot(z, energy_hartree, 1)

  v0=-28.11193019d0*4.556335d-6
  v = energy_hartree-v0

end subroutine getpot


subroutine methaneradau(x,z)
  implicit none
  integer :: i, j, ifirst=0
  double precision :: x(5,3), r(4,3), s(4,3), shift(3), z(10), dist(3), dnorm
  double precision :: mass_c=12.d0, mass_h=1.007825d0, mass_tot, alpha0
  double precision, parameter :: unitconvm = 5.48579903d-04
  save :: ifirst, alpha0

! get the Radau coordinates for the special case of methane needed to call the T8 methane PES
! x are the cartesian coordinates, z are the radau coordinates
! Assumes that the C coordinates are first
!
  if(ifirst.eq.0)then
! The potential used NUCLEAR masses so substract the electron mass....
! DO NOT change masses here when using isotopic sustitution of the CH4 PES
    mass_c=mass_c/unitconvm  -6d0
    mass_h=mass_h/unitconvm  -1d0
    ifirst=1
    alpha0=sqrt(mass_c/(mass_c+4*mass_h))
  endif

  do i=1,4
    r(i,:)=x(i+1,:)-x(1,:)
  enddo

  shift(:)=0.25d0*(1d0-alpha0)*(r(1,:)+r(2,:)+r(3,:)+r(4,:))

  do i=1,4
    s(i,:)=r(i,:)-shift(:)
  enddo

  z(1)=sqrt(s(1,1)**2+s(1,2)**2+s(1,3)**2)
  z(2)=sqrt(s(2,1)**2+s(2,2)**2+s(2,3)**2)
  z(3)=sqrt(s(3,1)**2+s(3,2)**2+s(3,3)**2)
  z(4)=sqrt(s(4,1)**2+s(4,2)**2+s(4,3)**2)

  dist(:)=r(1,:)-r(2,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(5)=acos( (z(1)*z(1)+z(2)*z(2)-dnorm*dnorm)/(2d0*z(1)*z(2)))    !H1CH2 angle

  dist(:)=r(1,:)-r(3,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(6)=acos( (z(1)*z(1)+z(3)*z(3)-dnorm*dnorm)/(2d0*z(1)*z(3)))    !H1CH3 angle

  dist(:)=r(1,:)-r(4,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(7)=acos( (z(1)*z(1)+z(4)*z(4)-dnorm*dnorm)/(2d0*z(1)*z(4)))    !H1CH4 angle

  dist(:)=r(2,:)-r(3,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(8)=acos( (z(2)*z(2)+z(3)*z(3)-dnorm*dnorm)/(2d0*z(2)*z(3)))    !H2CH3 angle

  dist(:)=r(2,:)-r(4,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(9)=acos( (z(2)*z(2)+z(4)*z(4)-dnorm*dnorm)/(2d0*z(2)*z(4)))    !H2CH4 angle

  dist(:)=r(3,:)-r(4,:)
  dnorm=sqrt( dist(1)*dist(1)+dist(2)*dist(2)+dist(3)*dist(3) )
  z(10)=acos( (z(3)*z(3)+z(4)*z(4)-dnorm*dnorm)/(2d0*z(3)*z(4)))    !H3CH4 angle

end subroutine methaneradau
