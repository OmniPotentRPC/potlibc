      program test 
***********************************************************************
*  Program to calculate the six interatomic distance 
*  by reading XYZ coordinate
***********************************************************************
      implicit none

      integer i
      double precision E
      double precision X(4),Y(4),Z(4),dEdX(4),dEdY(4),dEdZ(4)
      

C Read coordinates from "n4pes.dat"
      open(unit=20,file="n4pes.dat")

***********************************************************************
* Read the Cartesian coordinate of 4 atoms from input
* Cartesian coordinates:
* N1: X1, Y1, Z1
* N2: X2, Y2, Z2
* N3: X3, Y3, Z3
* N4: X4, Y4, Z4
***********************************************************************

      do i=1,4
         read(20,*) X(i),Y(i),Z(i)
      enddo
      close(unit=20)
C
C Now call pot to calculate V and dVdX
C
      call pot(X,Y,Z,E,dEdX,dEdY,dEdZ)

C
C Print Energy and Gradient to n4pes.out
C
      open(unit=30,file="test.out")
      write(30,9999) E
 9999 Format(' The energy is ',F18.10,' hartree.')
      write(30,*) 'The Cartesian gradients (in hartree/bohr):'
      do i=1,4
        write(30,9998) dEdX(i),dEdY(i),dEdZ(i)
      enddo
 9998 Format(2X,3F12.6)
      close(unit=30)

      end
