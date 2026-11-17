MODULE pot_mod

CONTAINS
	SUBROUTINE pot_initialize
		write(*,*)
		write(*,*) 'POT : CZLT SURFACE FOR H + C2H6 ----> C2H5 + H2 '
	END SUBROUTINE pot_initialize


END MODULE pot_mod


!!******************************************************************************************!!
SUBROUTINE getpot(V,NATOM,XX)
!!******************************************************************************************!!
! Geometry in Bohr
! Energy in kcal/mol
!
! Geometry format:
! 
!     H2 H6
!     |  |
!  H4-C1-C5-H8.....H9
!     |  |
!     H3 H7
!     
implicit none
integer,intent(in) :: NATOM
double precision,intent(in)  :: XX(NATOM,3)
double precision,intent(out) :: V

double precision, parameter :: HARTREE_IN_KCALMOL = 627.510d0
double precision, parameter :: BOHR_IN_ANG = 0.52917720830d0

double precision,parameter :: V0=-335.3578466555570690d0
double precision,parameter :: VR=5.072024009677989E-007
double precision,parameter :: VP=0.4276623940976380d0
double precision :: eng_kcal
double precision :: geo_ang(NATOM,3)

    geo_ang = XX * BOHR_IN_ANG
    call c2h7pot(geo_ang,eng_kcal)

! reactant as zero of eneergy
   V = ((eng_kcal-V0)-VR)/HARTREE_IN_KCALMOL

! product as zero of energy
!  V = (((eng_kcal-V0)-VR)+VP)/HARTREE_IN_KCALMOL

END SUBROUTINE getpot
