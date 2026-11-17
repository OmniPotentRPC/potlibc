!**********************************************************************
!   System:                     N3
!   Functional form:
!   Common name:		N3
!   Number of derivatives:      0
!   Number of bodies:           3
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References: doi 10.1063/1.4983813
!               Author:  Tapan K. Mankodi, Upendra V. Bhandarkar,
!               Bhalchandra P. Puranik
!
!   Notes:
!
!   Data file:  coeff_n3_o3.txt
!
***********************************************************************

!***********************************************************************!
! Global Potential energy surface for N3 system using CASSCF+CASPT2  	!
! Permutational Invariant Global Least Square fit                       !
! Input Variable: Interatomic distance R1-R2-R3 in Angstrom             !
! Output Variable: Energy in kcal/mol                                   !	
! Author of the Code: Tapan K. Mankodi                                  !
! 		      Indian Institute of Technology Bombay             !
! 		      Mumbai, India                                     !
! Ref: http://aip.scitation.org/doi/abs/10.1063/1.4983813               !
! Author:  Tapan K. Mankodi, Upendra V. Bhandarkar,                     !
!    	   Bhalchandra P. Puranik                                       ! 
! Corresponding Address: tapan.mankodi@iitb.ac.in                       !
!						 bhandarkar@iitb.ac.in  !					
!***********************************************************************!



program pes

	implicit none

	integer							::	mpoints
	double precision				::	coeff(52)
	double precision				::	indices(52,3)
	double precision				::	ae,re

	integer							::	i,j,k
	double precision				::	a

	re = 1.098d0
	ae = 0.9d0

		! Read the Coefficient file for the 3d PES
	open(unit=701,file='coeff_n3_o3.txt')
	j=0
	mpoints = 0
	do j=0,51
		mpoints = mpoints + 1
		read(701,*) i,indices(mpoints,1:3),coeff(mpoints),a
	end do
	close(701)
	
	print*, n3pes(1.098d0, 2.0d0, 2.0d0)
	
	
	
	
contains

double precision function n3pes(rx,ry,rz)

double precision :: rx,ry,rz,r4
double precision :: xrx,xry,xrz
integer			::	ii

xrx = exp(-(rx-re)/ae)
xry = exp(-(ry-re)/ae)
xrz = exp(-(rz-re)/ae)

n3pes = 0.0
do ii=1,mpoints
	n3pes = n3pes + coeff(ii)*(xrx**(indices(ii,1))*xry**(indices(ii,2))*xrz**(indices(ii,3)) + &
	xrx**(indices(ii,1))*xry**(indices(ii,3))*xrz**(indices(ii,2)) + &
	xrx**(indices(ii,2))*xry**(indices(ii,3))*xrz**(indices(ii,1)) + &
	xrx**(indices(ii,2))*xry**(indices(ii,1))*xrz**(indices(ii,3)) + &
	xrx**(indices(ii,3))*xry**(indices(ii,1))*xrz**(indices(ii,2)) + &
	xrx**(indices(ii,3))*xry**(indices(ii,2))*xrz**(indices(ii,1)))
end do	
end function

	

end program




 
