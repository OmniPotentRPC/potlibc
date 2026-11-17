!***********************************************************************
!   System:                     O4
!   Functional form:            Permutational Invariant Global Least Square
!   Common name:                O4 singlet
!   Number of derivatives:      0
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:  doi 10.1063/1.5008973
!	 Tapan K. Mankodi, Upendra V. Bhandarkar,
!        Bhalchandra P. Puranik
!
!   Notes:    
!
!   Data file: coeff_o4.txt
!
!***********************************************************************

!***********************************************************************!
! Global Potential energy surface for 									!
! singlet O4 system using CASSCF+CASPT2  								!	
! Permutational Invariant Global Least Square fit						!
! Input Variable: Interatomic distance R1-R2-R3-R4-R5-R6 in Angstrom	!
! Output Variable: Energy in kcal/mol									!	
! Author of the Code: Tapan K. Mankodi									!
! 					  Indian Institute of Technology Bombay				!
! 					  Mumbai, India										!
! Ref: http://aip.scitation.org/doi/abs/10.1063/1.5008973			!
! Author:  Tapan K. Mankodi, Upendra V. Bhandarkar, 					!
!    	   Bhalchandra P. Puranik										! 
! Corresponding Address: tapan.mankodi@iitb.ac.in						!
!						 bhandarkar@iitb.ac.in							!					
!***********************************************************************!



program pes

	implicit none

	integer							::	mpoints
	double precision				::	coeff(306)
	double precision				::	indices(306,6)
	integer							::	permutation(24,6)
	double precision				::	ae,re
	double precision				::	r1(6)
	integer							::	i,j,k
	double precision				::	a

	re = 1.208d0
	ae = 0.7d0
	mpoints = 306

		! Reading data from the coefficient files
	open(unit=701,file='coeff_o4.txt')
	do j=1,306
		read(701,*) i,indices(j,1:6),coeff(j)
	end do
	close(701)
	
		! Reading data for permutational invariance
	open(unit=702,file='permutation.txt')
	do j=1,24
		read(702,*) permutation(j,:)
	end do
	close(702)

	r1 = (/1.208, 8.87100487063, 9.96282078241, 9.82694933693, 10.9808487407, 10.6246822707/)
	print*, o4pes(r1)
	

	
contains

double precision function o4pes(r)

	double precision		::		r(6),xr(6)
	integer					::		ci,cj,ck,cl
	double precision		::		mult
	
	do ci=1,6
		xr(ci) = exp((1.208d0-r(ci))/0.7d0)
	end do
	
	o4pes = 0d0	
	do ci=1,306
		do cj = 1,24
			mult = 1.0
			do ck =1,6
				mult = mult*(xr(ck)**indices(ci,permutation(cj,ck)+1))
			end do	
			o4pes = o4pes + coeff(ci)*mult
		end do
	end do

end function	

end program




 
