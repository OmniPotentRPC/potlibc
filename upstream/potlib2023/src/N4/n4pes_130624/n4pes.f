C   System:                     N4
C   Functional form:            permutation-invariant polynomials
C   Common name:                N4(adiabatic ground state)
C   Number of derivatives:      1
C   Number of bodies:           4
C   Number of electronic surfaces: 1
C   Interface: Section-2
C
C   References:: Yuliya, Paukku, Ke R. Yang, Zoltan Varga,
C   and Donald G. Truhlar
C      J. Chem. Phys., submitted
C
C   Notes:    PES of N4 with special emphasize for
C             N2 + N2 --> N2 + N + N
C
C
C     N1--N2
C
C     N3--N4
C
C
C
C   Input: X(4),Y(4),Z(4)               in unit of bohr
C   Output: E                           in unit of hartree
C   Output: dEdX(4),dEdY(4),dEdZ(4)     hartree/bohr
C
      subroutine pot(X,Y,Z,E,dEdX,dEdY,dEdZ)

      implicit double precision (a-h,o-z)

C
C Convert factors
C Cconv: bohr to Angstrom
C Econv: kcal/mol to hartree
C Gconv: kcal/(mol*Angstrom) to hartree/bohr
C
      double precision Cconv
      double precision Econv
      double precision Gconv
      parameter(Cconv=0.5291772108d0)
      parameter(Econv=0.159360145d-2)
      parameter(Gconv=0.843297571d-3)
C
C Reference energy of infinite seperated N2 + N2 in hartree (taken from
C Yuliya in)
C
      double precision Eref
      parameter(Eref=-218.40801323d0)

      integer i
      double precision E,V
      double precision X(4),Y(4),Z(4),dEdX(4),dEdY(4),dEdZ(4)
      double precision Xcart(12),dVdX(12)

C Convert to local variables
      do i=1,4
        Xcart(3*i-2)=X(i)*Cconv
        Xcart(3*i-1)=Y(i)*Cconv
        Xcart(3*i)=Z(i)*Cconv
      enddo

      call n4pes(Xcart,v,dVdX,1)

C Convet local output to the ones ANT wants

      E = v*Econv + Eref
      do i=1,4
        dEdX(i)=dVdX(3*i-2)*Gconv
        dEdY(i)=dVdX(3*i-1)*Gconv
        dEdZ(i)=dVdX(3*i)*Gconv
      enddo

      end

C
C Local variables used in the N4 PES subroutine
C input coordinate matrix: X(12)in Ang
C                          N1: X(1),X(2),X(3)
C                          N2: X(4),X(5),X(6)
C                          N3: X(7),X(8),X(9)
C                          N4: X(10),X(11),X(12)
C input flag: igrad     igrad=0 energy-only calculation
C                       igrad=1 energy + gradient
C output potential energy:      v    in kcal/mol
C output gradient:              dVdX in kcal/(mol*Ang)
C
      subroutine n4pes(X,v,dVdX,igrad)
***********************************************************************
* Subroutine to calculate the potential energy V and gradient dVdX
* for given Cartesian coordinate X(12)  
* R:		Interatomic bond distance (6)
* V:		Calculated potential energy
* dVdX:		The derivative of V w.r.t. X, dim(12)
* dVdR:		The derivative of V w.r.t. R, dim(6) 
* dPdR:		The derivative of basis functions w.r.t. R
*		dim(6*306)
* dMdR:		The derivative of monomials w.r.t. R
*		dim(6*112)
* dRdX:		The derivative of R w.r.t. X, dim(6*12)
***********************************************************************
      
      implicit double precision (a-h,o-z) 

      integer i,igrad,j,nob,k
      double precision V
      double precision dVdX(12),X(12) 

C      common /coord/    R(6)
C      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
C      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
C     $                  dRdX(6,12),dBdR(6,276)
C      common /msprmt/   a,re

C Read cartesian coordinate from input file
      call coord_convt(X)

      if (igrad .le. 1) then
C Call subroutine EvV to evaluate potential energy V
        call evv(V)

        if (igrad .eq. 1) then
C Call EvdVdX to evaluate the derivatives of V w.r.t. X
          call evdvdx(X,dVdX)
        endif
      else
        write (*,*) 'Only igrad = 0, 1 is allowed!'
      endif

      end 
      subroutine coord_convt(X)
***********************************************************************
*  Program to calculate the six interatomic distance 
*  by reading XYZ coordinate
***********************************************************************
      implicit double precision (a-h,o-z) 

      integer i
      double precision X(12)

      common /coord/    R(6)
      
***********************************************************************
*  Now, calculate the inter-atomic distance
*  r1 = r(N1N2)    r2 = r(N1N3)
*  r3 = r(N1N4)    r4 = r(N2N3)
*  r5 = r(N2N4)    r6 = r(N3N4)       
***********************************************************************

      R(1)=Sqrt((X(4)-X(1))**2 + (X(5)-X(2))**2 + (X(6)-X(3))**2)
      R(2)=Sqrt((X(7)-X(1))**2 + (X(8)-X(2))**2 + (X(9)-X(3))**2)
      R(3)=Sqrt((X(10)-X(1))**2 + (X(11)-X(2))**2 + (X(12)-X(3))**2)
      R(4)=Sqrt((X(4)-X(7))**2 + (X(5)-X(8))**2 + (X(6)-X(9))**2)
      R(5)=Sqrt((X(4)-X(10))**2 + (X(5)-X(11))**2 + (X(6)-X(12))**2)
      R(6)=Sqrt((X(7)-X(10))**2 + (X(8)-X(11))**2 + (X(9)-X(12))**2)
 
      return

      end
      subroutine EvV(V)
***********************************************************************
* Subroutine to evaluate V for giving R 
* V(R) = C*P
* C:		Coefficients, stored in 'dim.inc' 
* P:		Basis functions evaluated for giving R
* rMs:		rMs(6), six Morse terms
* a:		Nonlinear parameters in Morse terrms(Angstrom)
* re:		Equilibrium bond length(Angstrom)
* nop:		number of points
* nom:  	number of monomials
* nob:  	number of basis functions(polynomials)
* rM(0:111):	Array to store monomials
* P(0:305):	Array to store polynomials
* B(1:276):     Array to store basis functions
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j,k
      double precision dist,dv2dr,V,V2

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)

C Calculate the six Morse terms for each point
      call evmorse

C Calculate the monomials for each point by using six Morse terms
      call evmono

C Calculate the polynomials (basis functions) by using monomials
      call evpoly 

C Calculate the basis functions by removing unconnected and 2-body terms
      call evbas

C Initialized v to be 2De 2*228.7 kcal/mol
      v=457.4d0
C Evaluate 2-body interactions
      do i=1,6
        dist=r(i)
        call ev2gm2(dist,v2,dv2dr,1,0)
        v=v+v2
      enddo

C Evaluate V by taken the product of C and Basis function array
      do i=1,276
        v=v + c(i)*b(i)
      enddo

C      Write(*,9999) V 
C 9999 Format('The potential energy is ',F20.14,' kcal/mol')

      return

      end 
      subroutine EvdVdX(X,dVdX)
***********************************************************************
* Subroutine to evaluate dRdX for giving R and X 
* R:		R(6), 6 bond lengths
* X:		X(12), 12 Cartesian coordinates
* rM(0:111):    Array to store monomials
* P(0:305):     Array to store polynomials
* dVdX:		dVdX(12), derivatives of V w.r.t. Cartesian coordinates 
* dVdR:		dVdR(6), derivatives of V w.r.t.6 bond length
* dRdX:		dRdX(6,12), derivatives of R(6) w.r.t. 12  
*		Cartesian coordinates
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j
      double precision dVdX(12),X(12)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)

C Initialize dVdX
      do i=1,12
        dVdX(i)=0.0d0
      enddo

C Call EvdVdR to evaluate dVdR(6)
      Call evdvdr

C Call EvdRdX to evaluate dRdX(6,12)
      Call evdrdx(X)  

C Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi), j=1 to 6
      do i=1,12
        do j=1,6
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo

C      write(*,*) 'The 12 dVdX are:'
C      Write(*,9999) (dVdX(i),i=1,12) 
C 9999 Format(1x,3F15.8)

      return
      end 
      subroutine EvMorse
      
      implicit double precision (a-h,o-z)

      integer i

      common /coord/    R(6)
      common /msprmt/   a,re
      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      
C Morse term ms = exp(-(r-re)/a)
C re:	equlibrium bond length(1.098A)
C a: 	nonlinear paramter optimized for diatomic molecules, unit Anstrom
      do i=1,6
         rms(i)=Exp(-(r(i)-re)/a)
      enddo

      end 
      subroutine EvMono
***********************************************************************
*  The subroutine reads six Morse variables(X) and calculate the
*  monomials(M) that do not have usable decomposition.
*  For A4 with max. degree 9, the number of monomials is nom.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)

      rm(0) = 1.0d0
      rm(1) = rms(6)
      rm(2) = rms(5)
      rm(3) = rms(4)
      rm(4) = rms(3)
      rm(5) = rms(2)
      rm(6) = rms(1)
      rm(7) = rm(3)*rm(4)
      rm(8) = rm(2)*rm(5)
      rm(9) = rm(1)*rm(6)
      rm(10) = rm(1)*rm(2)
      rm(11) = rm(1)*rm(3)
      rm(12) = rm(2)*rm(3)
      rm(13) = rm(1)*rm(4)
      rm(14) = rm(2)*rm(4)
      rm(15) = rm(1)*rm(5)
      rm(16) = rm(3)*rm(5)
      rm(17) = rm(4)*rm(5)
      rm(18) = rm(2)*rm(6)
      rm(19) = rm(3)*rm(6)
      rm(20) = rm(4)*rm(6)
      rm(21) = rm(5)*rm(6)
      rm(22) = rm(1)*rm(7)
      rm(23) = rm(2)*rm(7)
      rm(24) = rm(1)*rm(8)
      rm(25) = rm(2)*rm(16)
      rm(26) = rm(2)*rm(17)
      rm(27) = rm(3)*rm(17)
      rm(28) = rm(1)*rm(18)
      rm(29) = rm(1)*rm(19)
      rm(30) = rm(1)*rm(20)
      rm(31) = rm(3)*rm(20)
      rm(32) = rm(1)*rm(21)
      rm(33) = rm(2)*rm(21)
      rm(34) = rm(1)*rm(12)
      rm(35) = rm(1)*rm(17)
      rm(36) = rm(2)*rm(20)
      rm(37) = rm(3)*rm(21)
      rm(38) = rm(1)*rm(14)
      rm(39) = rm(1)*rm(16)
      rm(40) = rm(2)*rm(19)
      rm(41) = rm(4)*rm(21)
      rm(42) = rm(2)*rm(27)
      rm(43) = rm(1)*rm(31)
      rm(44) = rm(1)*rm(33)
      rm(45) = rm(1)*rm(23)
      rm(46) = rm(1)*rm(25)
      rm(47) = rm(1)*rm(26)
      rm(48) = rm(1)*rm(27)
      rm(49) = rm(1)*rm(40)
      rm(50) = rm(1)*rm(36)
      rm(51) = rm(2)*rm(31)
      rm(52) = rm(1)*rm(37)
      rm(53) = rm(2)*rm(37)
      rm(54) = rm(1)*rm(41)
      rm(55) = rm(2)*rm(41)
      rm(56) = rm(3)*rm(41)
      rm(57) = rm(1)*rm(42)
      rm(58) = rm(1)*rm(51)
      rm(59) = rm(1)*rm(53)
      rm(60) = rm(1)*rm(55)
      rm(61) = rm(1)*rm(56)
      rm(62) = rm(2)*rm(56)
      rm(63) = rm(1)*rm(62)
      rm(64) = rm(2)*rm(57)
      rm(65) = rm(3)*rm(57)
      rm(66) = rm(4)*rm(57)
      rm(67) = rm(5)*rm(57)
      rm(68) = rm(1)*rm(58)
      rm(69) = rm(3)*rm(58)
      rm(70) = rm(4)*rm(58)
      rm(71) = rm(1)*rm(59)
      rm(72) = rm(2)*rm(59)
      rm(73) = rm(1)*rm(60)
      rm(74) = rm(2)*rm(60)
      rm(75) = rm(1)*rm(61)
      rm(76) = rm(2)*rm(62)
      rm(77) = rm(3)*rm(61)
      rm(78) = rm(3)*rm(62)
      rm(79) = rm(4)*rm(61)
      rm(80) = rm(4)*rm(62)
      rm(81) = rm(5)*rm(59)
      rm(82) = rm(5)*rm(60)
      rm(83) = rm(5)*rm(62)
      rm(84) = rm(6)*rm(58)
      rm(85) = rm(6)*rm(59)
      rm(86) = rm(6)*rm(60)
      rm(87) = rm(6)*rm(61)
      rm(88) = rm(2)*rm(64)
      rm(89) = rm(3)*rm(65)
      rm(90) = rm(4)*rm(66)
      rm(91) = rm(5)*rm(67)
      rm(92) = rm(1)*rm(68)
      rm(93) = rm(3)*rm(69)
      rm(94) = rm(4)*rm(70)
      rm(95) = rm(1)*rm(71)
      rm(96) = rm(2)*rm(72)
      rm(97) = rm(1)*rm(73)
      rm(98) = rm(2)*rm(74)
      rm(99) = rm(1)*rm(75)
      rm(100) = rm(2)*rm(76)
      rm(101) = rm(3)*rm(77)
      rm(102) = rm(3)*rm(78)
      rm(103) = rm(4)*rm(79)
      rm(104) = rm(4)*rm(80)
      rm(105) = rm(5)*rm(81)
      rm(106) = rm(5)*rm(82)
      rm(107) = rm(5)*rm(83)
      rm(108) = rm(6)*rm(84)
      rm(109) = rm(6)*rm(85)
      rm(110) = rm(6)*rm(86)
      rm(111) = rm(6)*rm(87)


      return

      end subroutine EvMono
      subroutine EvPoly
***********************************************************************
*  The subroutine reads monomials(m) and calculate the
*  permutation invariant polynomials(p)
*  For A4 with max. degree 9, the number of polynomials is nob.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)

      p(0) = rm(0)
      p(1) = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
      p(2) = rm(7) + rm(8) + rm(9)
      p(3) = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) +
     $       rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
      p(4) = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
      p(5) = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) +
     $       rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
      p(6) = rm(34) + rm(35) + rm(36) + rm(37)
      p(7) = rm(38) + rm(39) + rm(40) + rm(41)
      p(8) = p(1)*p(2) - p(5)
      p(9) = p(1)*p(3) - p(6) - p(7) - p(5) - p(6)
     $      - p(7) - p(5) - p(6) - p(7)
      p(10) = p(1)*p(4) - p(9) - p(8)
      p(11) = rm(42) + rm(43) + rm(44)
      p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) +
     $       rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
      p(13) = p(2)*p(3) - p(12)
      p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12)
     $      - p(11) - p(11) - p(11)
      p(15) = p(1)*p(6) - p(12)
      p(16) = p(1)*p(7) - p(12)
      p(17) = p(2)*p(2) - p(11) - p(11)
      p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16)
     $      - p(14) - p(12) - p(11) - p(15) - p(16) - p(14)
     $      - p(12) - p(11) - p(12) - p(11)
      p(19) = p(2)*p(4) - p(14)
      p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
      p(21) = p(1)*p(10) - p(20) - p(19)
      p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
      p(23) = p(1)*p(11) - p(22)
      p(24) = p(2)*p(6)
      p(25) = p(2)*p(7)
      p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22)
     $      - p(22) - p(22)
      p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
      p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25)
     $      - p(23) - p(22) - p(24) - p(25) - p(23) - p(22)
     $      - p(22)
      p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
      p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
      p(31) = p(2)*p(9) - p(26) - p(28)
      p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
      p(33) = p(4)*p(6) - p(25)
      p(34) = p(4)*p(7) - p(24)
      p(35) = p(1)*p(17) - p(27)
      p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
      p(37) = p(2)*p(10) - p(32)
      p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
      p(39) = p(1)*p(21) - p(38) - p(37)
      p(40) = rm(63)
      p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) +
     $       rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + rm(75) +
     $       rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) +
     $       rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + rm(87)
      p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40)
     $      - p(40) - p(40) - p(40)
      p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
      p(44) = p(2)*p(12) - p(41)
      p(45) = p(3)*p(11) - p(41)
      p(46) = p(5)*p(6) - p(41)
      p(47) = p(5)*p(7) - p(41)
      p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
      p(49) = p(4)*p(11) - p(42)
      p(50) = p(2)*p(15) - p(46)
      p(51) = p(2)*p(16) - p(47)
      p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
      p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
      p(54) = p(6)*p(6) - p(42) - p(42)
      p(55) = p(7)*p(7) - p(42) - p(42)
      p(56) = p(3)*p(17) - p(44)
      p(57) = p(2)*p(18) - p(48)
      p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47)
     $      - p(45) - p(45)
      p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
      p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
      p(61) = p(2)*p(20) - p(52) - p(58)
      p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
      p(63) = p(6)*p(10) - p(51)
      p(64) = p(7)*p(10) - p(50)
      p(65) = p(2)*p(17) - p(43)
      p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59)
     $      - p(60) - p(58)
      p(67) = p(2)*p(19) - p(49)
      p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57)
     $      - p(66) - p(66)
      p(69) = p(2)*p(21) - p(62)
      p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
      p(71) = p(1)*p(39) - p(70) - p(69)
      p(72) = p(40)*p(1)
      p(73) = p(2)*p(22) - p(72)
      p(74) = p(6)*p(11)
      p(75) = p(7)*p(11)
      p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72)
     $      - p(72) - p(72)
      p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) +
     $       rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + rm(99) +
     $       rm(100) + rm(101) + rm(102) + rm(103) + rm(104) + rm(105) +
     $       rm(106) + rm(107) + rm(108) + rm(109) + rm(110) + rm(111)
      p(78) = p(1)*p(42) - p(72) - p(76)
      p(79) = p(5)*p(11) - p(72) - p(73) - p(72)
      p(80) = p(2)*p(26) - p(76) - p(77)
      p(81) = p(6)*p(12) - p(72) - p(76) - p(72)
      p(82) = p(7)*p(12) - p(72) - p(76) - p(72)
      p(83) = p(8)*p(11) - p(72)
      p(84) = p(6)*p(17)
      p(85) = p(7)*p(17)
      p(86) = p(9)*p(11) - p(76) - p(77)
      p(87) = p(2)*p(29) - p(81)
      p(88) = p(6)*p(14) - p(75) - p(75)
      p(89) = p(2)*p(30) - p(82)
      p(90) = p(7)*p(14) - p(74) - p(74)
      p(91) = p(1)*p(48) - p(76) - p(81) - p(82)
      p(92) = p(10)*p(11) - p(78)
      p(93) = p(2)*p(33) - p(88)
      p(94) = p(2)*p(34) - p(90)
      p(95) = p(10)*p(12) - p(77) - p(93) - p(94)
      p(96) = p(2)*p(27) - p(73) - p(79)
      p(97) = p(2)*p(28) - p(76) - p(86)
      p(98) = p(1)*p(53) - p(80) - p(79) - p(97)
      p(99) = p(1)*p(54) - p(81)
      p(100) = p(1)*p(55) - p(82)
      p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89)
     $      - p(86)
      p(102) = p(6)*p(18) - p(74) - p(90)
      p(103) = p(7)*p(18) - p(75) - p(88)
      p(104) = p(2)*p(31) - p(77) - p(86)
      p(105) = p(2)*p(36) - p(91) - p(101)
      p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90)
     $      - p(86)
      p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
      p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
      p(109) = p(2)*p(38) - p(95) - p(106)
      p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
      p(111) = p(6)*p(21) - p(94)
      p(112) = p(7)*p(21) - p(93)
      p(113) = p(1)*p(65) - p(96)
      p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
      p(115) = p(2)*p(37) - p(92)
      p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
      p(117) = p(2)*p(39) - p(110)
      p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
      p(119) = p(1)*p(71) - p(118) - p(117)
      p(120) = p(40)*p(2)
      p(121) = p(40)*p(3)
      p(122) = p(40)*p(4)
      p(123) = p(11)*p(12) - p(121)
      p(124) = p(2)*p(42) - p(122)
      p(125) = p(6)*p(22) - p(121)
      p(126) = p(7)*p(22) - p(121)
      p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
      p(128) = p(6)*p(23) - p(123)
      p(129) = p(7)*p(23) - p(123)
      p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125)
     $      - p(128) - p(126) - p(129) - p(124) - p(123) - p(122)
     $      - p(121) - p(120) - p(125) - p(126) - p(124) - p(123)
     $      - p(122) - p(121) - p(120) - p(122) - p(121) - p(120)
     $      - p(120) - p(120) - p(120) - p(120)
      p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
      p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129)
     $      - p(127) - p(121)
      p(133) = p(1)*p(78) - p(122) - p(131)
      p(134) = p(11)*p(11) - p(120) - p(120)
      p(135) = p(2)*p(48) - p(130)
      p(136) = p(11)*p(17) - p(120)
      p(137) = p(2)*p(44) - p(123)
      p(138) = p(2)*p(45) - p(121)
      p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
      p(140) = p(6)*p(27) - p(127)
      p(141) = p(2)*p(54)
      p(142) = p(7)*p(27) - p(127)
      p(143) = p(2)*p(52) - p(131) - p(132)
      p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
      p(145) = p(2)*p(55)
      p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
      p(147) = p(11)*p(18) - p(130)
      p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
      p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
      p(150) = p(6)*p(30) - p(121) - p(146)
      p(151) = p(11)*p(19) - p(122)
      p(152) = p(2)*p(50) - p(128)
      p(153) = p(2)*p(51) - p(129)
      p(154) = p(11)*p(20) - p(131) - p(132)
      p(155) = p(2)*p(59) - p(144) - p(148)
      p(156) = p(6)*p(32) - p(129)
      p(157) = p(2)*p(60) - p(146) - p(149)
      p(158) = p(7)*p(32) - p(128)
      p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
      p(160) = p(11)*p(21) - p(133)
      p(161) = p(2)*p(63) - p(156)
      p(162) = p(2)*p(64) - p(158)
      p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
      p(164) = p(2)*p(53) - p(124) - p(139)
      p(165) = p(2)*p(58) - p(131) - p(154)
      p(166) = p(3)*p(54) - p(125) - p(144)
      p(167) = p(3)*p(55) - p(126) - p(146)
      p(168) = p(3)*p(65) - p(137)
      p(169) = p(17)*p(18) - p(135)
      p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
      p(171) = p(4)*p(54) - p(135)
      p(172) = p(4)*p(55) - p(135)
      p(173) = p(2)*p(66) - p(150)
      p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173)
     $      - p(165) - p(147)
      p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
      p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
      p(177) = p(2)*p(61) - p(132) - p(154)
      p(178) = p(2)*p(68) - p(159) - p(174)
      p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158)
     $      - p(154)
      p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
      p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
      p(182) = p(2)*p(70) - p(163) - p(179)
      p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
      p(184) = p(6)*p(39) - p(162)
      p(185) = p(7)*p(39) - p(161)
      p(186) = p(2)*p(65) - p(136)
      p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175)
     $      - p(176) - p(174)
      p(188) = p(2)*p(67) - p(151)
      p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
      p(190) = p(2)*p(69) - p(160)
      p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
      p(192) = p(2)*p(71) - p(183)
      p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
      p(194) = p(1)*p(119) - p(193) - p(192)
      p(195) = p(40)*p(5)
      p(196) = p(40)*p(6)
      p(197) = p(40)*p(7)
      p(198) = p(40)*p(8)
      p(199) = p(40)*p(9)
      p(200) = p(40)*p(10)
      p(201) = p(11)*p(22) - p(195)
      p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196)
     $      - p(197) - p(195) - p(196) - p(197)
      p(203) = p(17)*p(22) - p(198)
      p(204) = p(6)*p(43)
      p(205) = p(7)*p(43)
      p(206) = p(11)*p(26) - p(199) - p(202)
      p(207) = p(2)*p(76) - p(199) - p(202)
      p(208) = p(2)*p(78) - p(200)
      p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
      p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
      p(213) = p(11)*p(29) - p(209)
      p(214) = p(11)*p(30) - p(211)
      p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
      p(216) = p(2)*p(77) - p(199) - p(206)
      p(217) = p(6)*p(49) - p(205)
      p(218) = p(7)*p(49) - p(204)
      p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209)
     $      - p(217) - p(211) - p(218) - p(207) - p(204) - p(205)
     $      - p(200) - p(199) - p(198) - p(200) - p(198) - p(200)
     $      - p(198)
      p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
      p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218)
     $      - p(216)
      p(222) = p(1)*p(133) - p(200) - p(220)
      p(223) = p(11)*p(27) - p(195) - p(203)
      p(224) = p(1)*p(134) - p(201)
      p(225) = p(2)*p(81) - p(209)
      p(226) = p(2)*p(80) - p(202) - p(206)
      p(227) = p(2)*p(82) - p(211)
      p(228) = p(2)*p(91) - p(215) - p(219)
      p(229) = p(11)*p(28) - p(199) - p(207)
      p(230) = p(6)*p(53) - p(205)
      p(231) = p(5)*p(54) - p(209)
      p(232) = p(7)*p(53) - p(204)
      p(233) = p(7)*p(54) - p(196)
      p(234) = p(5)*p(55) - p(211)
      p(235) = p(6)*p(55) - p(197)
      p(236) = p(11)*p(35) - p(198)
      p(237) = p(6)*p(65)
      p(238) = p(7)*p(65)
      p(239) = p(2)*p(86) - p(199) - p(229)
      p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
      p(241) = p(17)*p(29) - p(225)
      p(242) = p(2)*p(99) - p(231)
      p(243) = p(17)*p(30) - p(227)
      p(244) = p(2)*p(95) - p(220) - p(221)
      p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
      p(246) = p(2)*p(100) - p(234)
      p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
      p(248) = p(11)*p(36) - p(215) - p(219)
      p(249) = p(2)*p(102) - p(233)
      p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
      p(251) = p(2)*p(103) - p(235)
      p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
      p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
      p(254) = p(11)*p(37) - p(200)
      p(255) = p(2)*p(93) - p(217)
      p(256) = p(2)*p(94) - p(218)
      p(257) = p(11)*p(38) - p(220) - p(221)
      p(258) = p(2)*p(107) - p(245) - p(250)
      p(259) = p(6)*p(62) - p(218)
      p(260) = p(2)*p(108) - p(247) - p(252)
      p(261) = p(7)*p(62) - p(217)
      p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
      p(263) = p(11)*p(39) - p(222)
      p(264) = p(2)*p(111) - p(259)
      p(265) = p(2)*p(112) - p(261)
      p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
      p(267) = p(2)*p(98) - p(208) - p(240)
      p(268) = p(6)*p(54) - p(210)
      p(269) = p(7)*p(55) - p(212)
      p(270) = p(2)*p(96) - p(203) - p(223)
      p(271) = p(2)*p(97) - p(207) - p(229)
      p(272) = p(2)*p(101) - p(215) - p(248)
      p(273) = p(2)*p(106) - p(220) - p(257)
      p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
      p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
      p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251)
     $      - p(248)
      p(277) = p(6)*p(66) - p(213) - p(252)
      p(278) = p(7)*p(66) - p(214) - p(250)
      p(279) = p(2)*p(104) - p(216) - p(239)
      p(280) = p(2)*p(105) - p(219) - p(248)
      p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
      p(282) = p(10)*p(54) - p(227)
      p(283) = p(10)*p(55) - p(225)
      p(284) = p(2)*p(114) - p(253) - p(276)
      p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276)
     $      - p(273)
      p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
      p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
      p(288) = p(2)*p(109) - p(221) - p(257)
      p(289) = p(2)*p(116) - p(262) - p(285)
      p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261)
     $      - p(257)
      p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
      p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
      p(293) = p(2)*p(118) - p(266) - p(290)
      p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
      p(295) = p(6)*p(71) - p(265)
      p(296) = p(7)*p(71) - p(264)
      p(297) = p(1)*p(186) - p(270)
      p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
      p(299) = p(2)*p(115) - p(254)
      p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284)
     $      - p(298)
      p(301) = p(2)*p(117) - p(263)
      p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
      p(303) = p(2)*p(119) - p(294)
      p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
      p(305) = p(1)*p(194) - p(304) - p(303)

      return

      end subroutine EvPoly
      subroutine evbas
***********************************************************************
*  The subroutine elminate the 2-body terms in Bowman's approach
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      
      integer i
      double precision b1(306) 

C Pass P(0:305) to BM1(1:306)
      do i=1,306
        b1(i)=p(i-1)
      enddo


C Remove unconnected terms and 2-body terms and pass to B(1:276)
      b(1)=b1(4)

      do i=2,4
        b(i)=b1(i+4)
      enddo

      b(5)=b1(10)

      do i=6,11
        b(i)=b1(i+6)
      enddo

      b(12)=b1(19)
      b(13)=b1(21)

      do i=14,26
        b(i)=b1(i+9)
      enddo

      b(27)=b1(37)
      b(28)=b1(39)

      do i=29,53
        b(i)=b1(i+12)
      enddo

      b(54)=b1(67)
      b(55)=b1(69)
      b(56)=b1(71)

      do i=57,97
        b(i)=b1(i+16)
      enddo

      b(98)=b1(115)
      b(99)=b1(117)
      b(100)=b1(119)

      do i=101,166
        b(i)=b1(i+20)
      enddo

      b(167)=b1(188)
      b(168)=b1(190)
      b(169)=b1(192)
      b(170)=b1(194)

      do i=171,272
        b(i)=b1(i+25)
      enddo

      b(273)=b1(299)
      b(274)=b1(301)
      b(275)=b1(303)
      b(276)=b1(305)

      return

      end 
      subroutine ev2gm2(r,v2,dv2dr,imol,igrad) 
***********************************************************************
*
* Subroutine to evaluate the 2-body potential energy and gradient
* for given r with generalized Morse potential with Scheme 2.
* V(r) = 0 for r -> inifity
* V(r) = De*(1-exp(-f(y)*(r-re)))^2 - De       
* f(y) = c0 + c1*y + c2*y^2 + c3*y^3 + c4*y^4 + c5*y^5 + c6*y^6
* y = (r^4 - re^4)/(r^4 + re^4)
* 
* imol .eq. 1   Paramemters for N2 dissociation without SEC 
*               (for Yuliya's N4 data)
* imol .eq. 2   Paramemters for N2 dissociation with SEC 
*               (for Zoltan's N2O2 data)
* imol .eq. 3   Paramemters for NO dissociation with SEC 
*               (for Zoltan's N2O2 data)
* imol .eq. 4   Paramemters for O2 dissociation with SEC 
*               (for Zoltan's N2O2 data and future O4 data)
***********************************************************************
      
      implicit double precision (a-h,o-z) 

      integer igrad,imol
      double precision c(0:10)
      double precision re,de
      double precision fy,dfdr,r,u,v2,dv2dr,y,dydr

      if (imol .eq. 1) then
C Parameter for N2 dissociation without SEC (for Yuliya's N4 data)
        re=1.098d0
        de=228.7d0
        c(0) = 2.70963254293d0
        c(1) = 1.32620177271d-1
        c(2) = 2.96757048793d-1
        c(3) = 1.97112432229d-1
        c(4) =-5.02002309588d-1
        c(5) = 3.80734244606d-1
        c(6) = 1.21001628750d0
      else if (imol .eq. 2) then  
C Parameter for N2 dissociation with SEC (for Zoltan's N2O2)
        re=1.098d0
        de=228.4d0
        c(0) = 2.71405774451d0
        c(1) = 1.32757649829d-1
        c(2) = 2.66756890408d-1
        c(3) = 1.95350725241d-1
        c(4) =-4.08663480982d-1
        c(5) = 3.92451705557d-1
        c(6) = 1.13006674877d0
      else if (imol .eq. 3) then
C Parameter for NO dissociation with SEC (for Zoltan's N2O2)
        re=1.1508d0
        de=149.9d0
        c(0) = 2.81134495569d0
        c(1) = 1.43241169611d-1
        c(2) = 1.35760038863d-2
        c(3) = 3.92892178507d-1
        c(4) = 9.29495534058d-1
        c(5) = 2.66966672332d-1
        c(6) =-3.68118714223d-1
      else if (imol .eq. 4) then
C Parameter for O2 dissociation with SEC (for Zoltan's N2O2 and KRY's O4)
        re=1.208d0
        de=120.243d0
        c(0) = 2.69132890094d0
        c(1) = 3.39550045614d-1
        c(2) = 3.46401777195d-1
        c(3) =-7.76983671636d-1
        c(4) =-3.29632972405d-1
        c(5) = 2.40883331247d0
        c(6) = 2.09264029009d0
      endif

      y=(r*r*r*r - re*re*re*re)/(r*r*r*r + re*re*re*re)

      fy = c(0) + c(1)*y + c(2)*y*y + c(3)*y*y*y + c(4)*y*y*y*y
     $   + c(5)*y*y*y*y*y + c(6)*y*y*y*y*y*y

      u=exp(-fy*(r-re))

      v2=de*(1.0d0-u)*(1.0d0-u)-de

      if (igrad .eq. 1) then
        dfdy = c(1) + 2.0d0*c(2)*y + 3.0d0*c(3)*y*y + 4.0d0*c(4)*y*y*y 
     $       + 5.0d0*c(5)*y*y*y*y + 6.0d0*c(6)*y*y*y*y*y
        v = r*r*r*r + re*re*re*re 
        dydr = 8.0d0*r*r*r*re*re*re*re/(v*v)
        dfdr = dfdy*dydr
        dv2dr = 2.0d0*de*(1-u)*u*(dfdr*(r-re)+fy)
      endif
  
      end 
      subroutine EvdVdR
***********************************************************************
* Subroutine to evaluate dVdR for giving R 
* dVdR = dV2dR + C*dBdR
* C:		Coefficients, stored in 'dim.inc' 
* P:		Basis functions evaluated for giving R
* M:		Monomials evluated for giving R
* dV2dR:        Gradient of 2-body interactions
* dMsdR:	dMsdR(6,6), 6 Morse terms w.r.t. 6 bond length
* dMdR:		dMdR(6,nom), nom monomials w.r.t.6 bond length
* dPdR:		dPdR(6,nob), nop polynomial basis functions 
*		w.r.t. 6 bond length
* nom:  	number of monomials
* nob:  	number of basis functions(polynomials)
* M(nom):	Array to store monomials
* P(nob):	Array to store polynomials
***********************************************************************
      
      implicit double precision (a-h,o-z)
      
      integer i,j
      double precision dist,v2,dv2dr

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)

C Initialize dVdR(6)
      do i=1,6
        dVdR(i)=0.0d0
      enddo

C Add dV2dR(i) to dVdR
      do i=1,6
        dist=R(i)
        call ev2gm2(dist,v2,dv2dr,1,1)
        dVdR(i)=dv2dr
      enddo

C Calculate dMorse/dr(6,6) for giving R(6)
      call evdmsdr

C Calculate the monomials for each point by using six Morse terms
      call evdmdr

C Calculate the polynomials by using monomials
      call evdpdr 

C Remove 2-body interactions and unconnected terms from polynomials
      call evdbdr

C Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
      do i=1,6      
        do j=1,276
         dVdR(i)=dVdR(i) + c(j)*dBdR(i,j)
        enddo
      enddo

      return
      end 
      subroutine EvdRdX(X)
***********************************************************************
* Subroutine to evaluate dRdX for giving R and X 
* R:		R(6), 6 bond lengths
* X:		X(12), 12 Cartesian coordinates
* 
* dMdR:		dMdR(6,nom), nom monomials w.r.t.6 bond length
* dPdR:		dPdR(6,nob), nop polynomial basis functions 
*		w.r.t. 6 bond length
* M(nom):	Array to store monomials
* P(nob):	Array to store polynomials
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j
      double precision X(12)

      common /coord/    R(6)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)

C Initialize dRdX(6,12)
      do i=1,6
        do j=1,12
          dRdX(i,j)=0.0d0
        enddo
      enddo

C Start to calculate the non-zero dRdX
C dr1dx
      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

C dr2dx
      dRdX(2,1)=(x(1)-x(7))/r(2)
      dRdX(2,2)=(x(2)-x(8))/r(2)
      dRdX(2,3)=(x(3)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)

C dr3dx
      dRdX(3,1)=(x(1)-x(10))/r(3)
      dRdX(3,2)=(x(2)-x(11))/r(3)
      dRdX(3,3)=(x(3)-x(12))/r(3)
      dRdX(3,10)=-dRdX(3,1)
      dRdX(3,11)=-dRdX(3,2)
      dRdX(3,12)=-dRdX(3,3)

C dr4dx
      dRdX(4,4)=(x(4)-x(7))/r(4)
      dRdX(4,5)=(x(5)-x(8))/r(4)
      dRdX(4,6)=(x(6)-x(9))/r(4)
      dRdX(4,7)=-dRdX(4,4)
      dRdX(4,8)=-dRdX(4,5)
      dRdX(4,9)=-dRdX(4,6)

C dr5dx
      dRdX(5,4)=(x(4)-x(10))/r(5)
      dRdX(5,5)=(x(5)-x(11))/r(5)
      dRdX(5,6)=(x(6)-x(12))/r(5)
      dRdX(5,10)=-dRdX(5,4)
      dRdX(5,11)=-dRdX(5,5)
      dRdX(5,12)=-dRdX(5,6)

C dr6dx
      dRdX(6,7)=(x(7)-x(10))/r(6)
      dRdX(6,8)=(x(8)-x(11))/r(6)
      dRdX(6,9)=(x(9)-x(12))/r(6)
      dRdX(6,10)=-dRdX(6,7)
      dRdX(6,11)=-dRdX(6,8)
      dRdX(6,12)=-dRdX(6,9)
C Finish the calculation of non-zero dRdX

      return

      end 
      subroutine EvdMsdR
***********************************************************************
* Subroutine to evalute the derivatives of Morse term X
* w.r.t. interatomic distance R(6)
* dmsdR:	Local variables, dirm(6,6)
* a:		Nonlinear pamameter(Angstrom)
* re:		equilibrium bond length(Angstrom)
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j

      common /coord/    R(6)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)
      common /msprmt/   a,re

C Initialize dmsdr
      do i=1,6
        do j=1,6
          dmsdr(i,j)=0.0d0
        enddo
      enddo
C
C Morse term dmsdr = exp(-(r-re)/a)
C dmsdr(i,j)=0	i!=j
C dmsdr(i,i)= -(1/a)*Exp(-r(i)-re)/a)
C
      do i=1,6
         dmsdr(i,i)=-(1/a)*Exp(-(r(i)-re)/a)
      enddo 

      return

      end 
      subroutine EvdMdR
***********************************************************************
*  The subroutine reads M(nom) and dMSdR(6,6) and calculate the
*  dMdR(6,nom) that do not have usable decomposition.
*  For A4 with max. degree 9, the number of monomials is nom.
***********************************************************************

      implicit double precision (a-h,o-z)

      integer i

      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)

      do i=1,6
        dmdr(i,0) = 0.0d0
        dmdr(i,1) = dmsdr(i,6)
        dmdr(i,2) = dmsdr(i,5)
        dmdr(i,3) = dmsdr(i,4)
        dmdr(i,4) = dmsdr(i,3)
        dmdr(i,5) = dmsdr(i,2)
        dmdr(i,6) = dmsdr(i,1)
        dmdr(i,7) = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
        dmdr(i,8) = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
        dmdr(i,9) = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
        dmdr(i,10) = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
        dmdr(i,11) = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
        dmdr(i,12) = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
        dmdr(i,13) = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
        dmdr(i,14) = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
        dmdr(i,15) = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
        dmdr(i,16) = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
        dmdr(i,17) = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
        dmdr(i,18) = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
        dmdr(i,19) = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
        dmdr(i,20) = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
        dmdr(i,21) = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
        dmdr(i,22) = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
        dmdr(i,23) = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
        dmdr(i,24) = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
        dmdr(i,25) = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
        dmdr(i,26) = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
        dmdr(i,27) = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
        dmdr(i,28) = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
        dmdr(i,29) = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
        dmdr(i,30) = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
        dmdr(i,31) = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
        dmdr(i,32) = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
        dmdr(i,33) = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
        dmdr(i,34) = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
        dmdr(i,35) = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
        dmdr(i,36) = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
        dmdr(i,37) = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
        dmdr(i,38) = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
        dmdr(i,39) = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
        dmdr(i,40) = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
        dmdr(i,41) = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
        dmdr(i,42) = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
        dmdr(i,43) = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
        dmdr(i,44) = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
        dmdr(i,45) = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
        dmdr(i,46) = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
        dmdr(i,47) = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
        dmdr(i,48) = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
        dmdr(i,49) = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
        dmdr(i,50) = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
        dmdr(i,51) = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
        dmdr(i,52) = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
        dmdr(i,53) = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
        dmdr(i,54) = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
        dmdr(i,55) = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
        dmdr(i,56) = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
        dmdr(i,57) = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
        dmdr(i,58) = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
        dmdr(i,59) = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
        dmdr(i,60) = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
        dmdr(i,61) = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
        dmdr(i,62) = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
        dmdr(i,63) = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
        dmdr(i,64) = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
        dmdr(i,65) = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
        dmdr(i,66) = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
        dmdr(i,67) = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
        dmdr(i,68) = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
        dmdr(i,69) = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
        dmdr(i,70) = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
        dmdr(i,71) = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
        dmdr(i,72) = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
        dmdr(i,73) = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
        dmdr(i,74) = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
        dmdr(i,75) = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
        dmdr(i,76) = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
        dmdr(i,77) = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
        dmdr(i,78) = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
        dmdr(i,79) = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
        dmdr(i,80) = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
        dmdr(i,81) = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
        dmdr(i,82) = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
        dmdr(i,83) = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
        dmdr(i,84) = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
        dmdr(i,85) = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
        dmdr(i,86) = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
        dmdr(i,87) = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
        dmdr(i,88) = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
        dmdr(i,89) = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
        dmdr(i,90) = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
        dmdr(i,91) = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
        dmdr(i,92) = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
        dmdr(i,93) = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
        dmdr(i,94) = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
        dmdr(i,95) = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
        dmdr(i,96) = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
        dmdr(i,97) = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
        dmdr(i,98) = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
        dmdr(i,99) = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
        dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
        dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
        dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
        dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
        dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
        dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
        dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
        dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
        dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
        dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
        dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
        dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
      enddo

      return

      end subroutine EvdMdR
      subroutine EvdPdr
***********************************************************************
*  The subroutine reads monomials(m) and calculate the
*  permutation invariant polynomials(p)
*  For A4 with max. degree 9, the number of polynomials is nob.
***********************************************************************

      implicit double precision (a-h,o-z)

      integer i

      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)

      do i=1,6
        dpdr(i,0) = dmdr(i,0)
        dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) +
     $         dmdr(i,5) + dmdr(i,6)
        dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
        dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) +
     $         dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i,17) +
     $         dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
        dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1)
     $        - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
        dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) +
     $         dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i,29) +
     $         dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
        dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
        dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
        dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2)
     $        - dpdr(i,5)
        dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3)
     $        - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6)
     $        - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
        dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4)
     $        - dpdr(i,9) - dpdr(i,8)
        dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
        dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) +
     $         dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + dmdr(i,52) +
     $         dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
        dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3)
     $        - dpdr(i,12)
        dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5)
     $        - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12)
     $        - dpdr(i,11) - dpdr(i,11) - dpdr(i,11)
        dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6)
     $        - dpdr(i,12)
        dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7)
     $        - dpdr(i,12)
        dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2)
     $        - dpdr(i,11) - dpdr(i,11)
        dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3)
     $        - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16)
     $        - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15)
     $        - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11)
     $        - dpdr(i,12) - dpdr(i,11)
        dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4)
     $        - dpdr(i,14)
        dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4)
     $        - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
        dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10)
     $        - dpdr(i,20) - dpdr(i,19)
        dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) +
     $         dmdr(i,61) + dmdr(i,62)
        dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11)
     $        - dpdr(i,22)
        dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
        dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
        dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12)
     $        - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22)
     $        - dpdr(i,22) - dpdr(i,22)
        dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5)
     $        - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
        dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5)
     $        - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25)
     $        - dpdr(i,23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25)
     $        - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
        dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6)
     $        - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7)
     $        - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
        dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9)
     $        - dpdr(i,26) - dpdr(i,28)
        dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14)
     $        - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
        dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6)
     $        - dpdr(i,25)
        dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7)
     $        - dpdr(i,24)
        dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17)
     $        - dpdr(i,27)
        dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18)
     $        - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
        dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10)
     $        - dpdr(i,32)
        dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10)
     $        - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
        dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21)
     $        - dpdr(i,38) - dpdr(i,37)
        dpdr(i,40) = dmdr(i,63)
        dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) +
     $         dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + dmdr(i,71) +
     $         dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) +
     $         dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + dmdr(i,79) +
     $         dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) +
     $         dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + dmdr(i,87)
        dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22)
     $        - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40)
     $        - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11)
     $        - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12)
     $        - dpdr(i,41)
        dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11)
     $        - dpdr(i,41)
        dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6)
     $        - dpdr(i,41)
        dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7)
     $        - dpdr(i,41)
        dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7)
     $        - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
        dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11)
     $        - dpdr(i,42)
        dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15)
     $        - dpdr(i,46)
        dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16)
     $        - dpdr(i,47)
        dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12)
     $        - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
        dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14)
     $        - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
        dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6)
     $        - dpdr(i,42) - dpdr(i,42)
        dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7)
     $        - dpdr(i,42) - dpdr(i,42)
        dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17)
     $        - dpdr(i,44)
        dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18)
     $        - dpdr(i,48)
        dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14)
     $        - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47)
     $        - dpdr(i,45) - dpdr(i,45)
        dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9)
     $        - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
        dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9)
     $        - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
        dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20)
     $        - dpdr(i,52) - dpdr(i,58)
        dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32)
     $        - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
        dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10)
     $        - dpdr(i,51)
        dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10)
     $        - dpdr(i,50)
        dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17)
     $        - dpdr(i,43)
        dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18)
     $        - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59)
     $        - dpdr(i,60) - dpdr(i,58)
        dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19)
     $        - dpdr(i,49)
        dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36)
     $        - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57)
     $        - dpdr(i,66) - dpdr(i,66)
        dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21)
     $        - dpdr(i,62)
        dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21)
     $        - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
        dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39)
     $        - dpdr(i,70) - dpdr(i,69)
        dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
        dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22)
     $        - dpdr(i,72)
        dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
        dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
        dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22)
     $        - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72)
     $        - dpdr(i,72) - dpdr(i,72)
        dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) +
     $         dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + dmdr(i,95) +
     $         dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) +
     $         dmdr(i,100) + dmdr(i,101) + dmdr(i,102) + dmdr(i,103) +
     $         dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) +
     $         dmdr(i,108) + dmdr(i,109) + dmdr(i,110) + dmdr(i,111)
        dpdr(i,78) = dpdr(i,1)*p(42) + p(1)*dpdr(i,42)
     $        - dpdr(i,72) - dpdr(i,76)
        dpdr(i,79) = dpdr(i,5)*p(11) + p(5)*dpdr(i,11)
     $        - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
        dpdr(i,80) = dpdr(i,2)*p(26) + p(2)*dpdr(i,26)
     $        - dpdr(i,76) - dpdr(i,77)
        dpdr(i,81) = dpdr(i,6)*p(12) + p(6)*dpdr(i,12)
     $        - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,82) = dpdr(i,7)*p(12) + p(7)*dpdr(i,12)
     $        - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
        dpdr(i,83) = dpdr(i,8)*p(11) + p(8)*dpdr(i,11)
     $        - dpdr(i,72)
        dpdr(i,84) = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
        dpdr(i,85) = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
        dpdr(i,86) = dpdr(i,9)*p(11) + p(9)*dpdr(i,11)
     $        - dpdr(i,76) - dpdr(i,77)
        dpdr(i,87) = dpdr(i,2)*p(29) + p(2)*dpdr(i,29)
     $        - dpdr(i,81)
        dpdr(i,88) = dpdr(i,6)*p(14) + p(6)*dpdr(i,14)
     $        - dpdr(i,75) - dpdr(i,75)
        dpdr(i,89) = dpdr(i,2)*p(30) + p(2)*dpdr(i,30)
     $        - dpdr(i,82)
        dpdr(i,90) = dpdr(i,7)*p(14) + p(7)*dpdr(i,14)
     $        - dpdr(i,74) - dpdr(i,74)
        dpdr(i,91) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48)
     $        - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
        dpdr(i,92) = dpdr(i,10)*p(11) + p(10)*dpdr(i,11)
     $        - dpdr(i,78)
        dpdr(i,93) = dpdr(i,2)*p(33) + p(2)*dpdr(i,33)
     $        - dpdr(i,88)
        dpdr(i,94) = dpdr(i,2)*p(34) + p(2)*dpdr(i,34)
     $        - dpdr(i,90)
        dpdr(i,95) = dpdr(i,10)*p(12) + p(10)*dpdr(i,12)
     $        - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
        dpdr(i,96) = dpdr(i,2)*p(27) + p(2)*dpdr(i,27)
     $        - dpdr(i,73) - dpdr(i,79)
        dpdr(i,97) = dpdr(i,2)*p(28) + p(2)*dpdr(i,28)
     $        - dpdr(i,76) - dpdr(i,86)
        dpdr(i,98) = dpdr(i,1)*p(53) + p(1)*dpdr(i,53)
     $        - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
        dpdr(i,99) = dpdr(i,1)*p(54) + p(1)*dpdr(i,54)
     $        - dpdr(i,81)
        dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55)
     $        - dpdr(i,82)
        dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18)
     $        - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89)
     $        - dpdr(i,86)
        dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18)
     $        - dpdr(i,74) - dpdr(i,90)
        dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18)
     $        - dpdr(i,75) - dpdr(i,88)
        dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31)
     $        - dpdr(i,77) - dpdr(i,86)
        dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36)
     $        - dpdr(i,91) - dpdr(i,101)
        dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32)
     $        - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90)
     $        - dpdr(i,86)
        dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29)
     $        - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
        dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30)
     $        - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
        dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38)
     $        - dpdr(i,95) - dpdr(i,106)
        dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62)
     $        - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
        dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21)
     $        - dpdr(i,94)
        dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21)
     $        - dpdr(i,93)
        dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65)
     $        - dpdr(i,96)
        dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66)
     $        - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
        dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37)
     $        - dpdr(i,92)
        dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18)
     $        - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
        dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39)
     $        - dpdr(i,110)
        dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39)
     $        - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
        dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71)
     $        - dpdr(i,118) - dpdr(i,117)
        dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
        dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
        dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
        dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12)
     $        - dpdr(i,121)
        dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42)
     $        - dpdr(i,122)
        dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22)
     $        - dpdr(i,121)
        dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22)
     $        - dpdr(i,121)
        dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41)
     $        - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
        dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23)
     $        - dpdr(i,123)
        dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23)
     $        - dpdr(i,123)
        dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41)
     $        - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125)
     $        - dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124)
     $        - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120)
     $        - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123)
     $        - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,122)
     $        - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120)
     $        - dpdr(i,120) - dpdr(i,120)
        dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42)
     $        - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
        dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41)
     $        - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129)
     $        - dpdr(i,127) - dpdr(i,121)
        dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78)
     $        - dpdr(i,122) - dpdr(i,131)
        dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11)
     $        - dpdr(i,120) - dpdr(i,120)
        dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48)
     $        - dpdr(i,130)
        dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17)
     $        - dpdr(i,120)
        dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44)
     $        - dpdr(i,123)
        dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45)
     $        - dpdr(i,121)
        dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14)
     $        - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
        dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27)
     $        - dpdr(i,127)
        dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
        dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27)
     $        - dpdr(i,127)
        dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52)
     $        - dpdr(i,131) - dpdr(i,132)
        dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81)
     $        - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
        dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
        dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82)
     $        - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
        dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18)
     $        - dpdr(i,130)
        dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28)
     $        - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28)
     $        - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
        dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30)
     $        - dpdr(i,121) - dpdr(i,146)
        dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19)
     $        - dpdr(i,122)
        dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50)
     $        - dpdr(i,128)
        dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51)
     $        - dpdr(i,129)
        dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20)
     $        - dpdr(i,131) - dpdr(i,132)
        dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59)
     $        - dpdr(i,144) - dpdr(i,148)
        dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32)
     $        - dpdr(i,129)
        dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60)
     $        - dpdr(i,146) - dpdr(i,149)
        dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32)
     $        - dpdr(i,128)
        dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34)
     $        - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
        dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21)
     $        - dpdr(i,133)
        dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63)
     $        - dpdr(i,156)
        dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64)
     $        - dpdr(i,158)
        dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21)
     $        - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
        dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53)
     $        - dpdr(i,124) - dpdr(i,139)
        dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58)
     $        - dpdr(i,131) - dpdr(i,154)
        dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54)
     $        - dpdr(i,125) - dpdr(i,144)
        dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55)
     $        - dpdr(i,126) - dpdr(i,146)
        dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65)
     $        - dpdr(i,137)
        dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18)
     $        - dpdr(i,135)
        dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98)
     $        - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
        dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54)
     $        - dpdr(i,135)
        dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55)
     $        - dpdr(i,135)
        dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66)
     $        - dpdr(i,150)
        dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101)
     $        - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173)
     $        - dpdr(i,165) - dpdr(i,147)
        dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102)
     $        - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
        dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103)
     $        - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
        dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61)
     $        - dpdr(i,132) - dpdr(i,154)
        dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68)
     $        - dpdr(i,159) - dpdr(i,174)
        dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62)
     $        - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158)
     $        - dpdr(i,154)
        dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38)
     $        - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
        dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38)
     $        - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
        dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70)
     $        - dpdr(i,163) - dpdr(i,179)
        dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110)
     $        - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
        dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39)
     $        - dpdr(i,162)
        dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39)
     $        - dpdr(i,161)
        dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65)
     $        - dpdr(i,136)
        dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66)
     $        - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175)
     $        - dpdr(i,176) - dpdr(i,174)
        dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67)
     $        - dpdr(i,151)
        dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66)
     $        - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
        dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69)
     $        - dpdr(i,160)
        dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21)
     $        - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
        dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71)
     $        - dpdr(i,183)
        dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71)
     $        - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
        dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119)
     $        - dpdr(i,193) - dpdr(i,192)
        dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
        dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
        dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
        dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
        dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
        dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
        dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22)
     $        - dpdr(i,195)
        dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22)
     $        - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196)
     $        - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
        dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22)
     $        - dpdr(i,198)
        dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
        dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
        dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26)
     $        - dpdr(i,199) - dpdr(i,202)
        dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76)
     $        - dpdr(i,199) - dpdr(i,202)
        dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78)
     $        - dpdr(i,200)
        dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41)
     $        - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42)
     $        - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
        dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41)
     $        - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
        dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42)
     $        - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
        dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29)
     $        - dpdr(i,209)
        dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30)
     $        - dpdr(i,211)
        dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22)
     $        - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
        dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77)
     $        - dpdr(i,199) - dpdr(i,206)
        dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49)
     $        - dpdr(i,205)
        dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49)
     $        - dpdr(i,204)
        dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77)
     $        - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209)
     $        - dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207)
     $        - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - dpdr(i,199)
     $        - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200)
     $        - dpdr(i,198)
        dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78)
     $        - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
        dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41)
     $        - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218)
     $        - dpdr(i,216)
        dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133)
     $        - dpdr(i,200) - dpdr(i,220)
        dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27)
     $        - dpdr(i,195) - dpdr(i,203)
        dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134)
     $        - dpdr(i,201)
        dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81)
     $        - dpdr(i,209)
        dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80)
     $        - dpdr(i,202) - dpdr(i,206)
        dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82)
     $        - dpdr(i,211)
        dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91)
     $        - dpdr(i,215) - dpdr(i,219)
        dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28)
     $        - dpdr(i,199) - dpdr(i,207)
        dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53)
     $        - dpdr(i,205)
        dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54)
     $        - dpdr(i,209)
        dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53)
     $        - dpdr(i,204)
        dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54)
     $        - dpdr(i,196)
        dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55)
     $        - dpdr(i,211)
        dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55)
     $        - dpdr(i,197)
        dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35)
     $        - dpdr(i,198)
        dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
        dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
        dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86)
     $        - dpdr(i,199) - dpdr(i,229)
        dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139)
     $        - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
        dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29)
     $        - dpdr(i,225)
        dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99)
     $        - dpdr(i,231)
        dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30)
     $        - dpdr(i,227)
        dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95)
     $        - dpdr(i,220) - dpdr(i,221)
        dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81)
     $        - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
        dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100)
     $        - dpdr(i,234)
        dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82)
     $        - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
        dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36)
     $        - dpdr(i,215) - dpdr(i,219)
        dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102)
     $        - dpdr(i,233)
        dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58)
     $        - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
        dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103)
     $        - dpdr(i,235)
        dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58)
     $        - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
        dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150)
     $        - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
        dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37)
     $        - dpdr(i,200)
        dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93)
     $        - dpdr(i,217)
        dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94)
     $        - dpdr(i,218)
        dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38)
     $        - dpdr(i,220) - dpdr(i,221)
        dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107)
     $        - dpdr(i,245) - dpdr(i,250)
        dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62)
     $        - dpdr(i,218)
        dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108)
     $        - dpdr(i,247) - dpdr(i,252)
        dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62)
     $        - dpdr(i,217)
        dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64)
     $        - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
        dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39)
     $        - dpdr(i,222)
        dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111)
     $        - dpdr(i,259)
        dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112)
     $        - dpdr(i,261)
        dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39)
     $        - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
        dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98)
     $        - dpdr(i,208) - dpdr(i,240)
        dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54)
     $        - dpdr(i,210)
        dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55)
     $        - dpdr(i,212)
        dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96)
     $        - dpdr(i,203) - dpdr(i,223)
        dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97)
     $        - dpdr(i,207) - dpdr(i,229)
        dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101)
     $        - dpdr(i,215) - dpdr(i,248)
        dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106)
     $        - dpdr(i,220) - dpdr(i,257)
        dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59)
     $        - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
        dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60)
     $        - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
        dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66)
     $        - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251)
     $        - dpdr(i,248)
        dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66)
     $        - dpdr(i,213) - dpdr(i,252)
        dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66)
     $        - dpdr(i,214) - dpdr(i,250)
        dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104)
     $        - dpdr(i,216) - dpdr(i,239)
        dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105)
     $        - dpdr(i,219) - dpdr(i,248)
        dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170)
     $        - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
        dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54)
     $        - dpdr(i,227)
        dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55)
     $        - dpdr(i,225)
        dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114)
     $        - dpdr(i,253) - dpdr(i,276)
        dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174)
     $        - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276)
     $        - dpdr(i,273)
        dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68)
     $        - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
        dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68)
     $        - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
        dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109)
     $        - dpdr(i,221) - dpdr(i,257)
        dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116)
     $        - dpdr(i,262) - dpdr(i,285)
        dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110)
     $        - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261)
     $        - dpdr(i,257)
        dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70)
     $        - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
        dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70)
     $        - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
        dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118)
     $        - dpdr(i,266) - dpdr(i,290)
        dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183)
     $        - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
        dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71)
     $        - dpdr(i,265)
        dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71)
     $        - dpdr(i,264)
        dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186)
     $        - dpdr(i,270)
        dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187)
     $        - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
        dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115)
     $        - dpdr(i,254)
        dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189)
     $        - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284)
     $        - dpdr(i,298)
        dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117)
     $        - dpdr(i,263)
        dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39)
     $        - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
        dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119)
     $        - dpdr(i,294)
        dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119)
     $        - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
        dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194)
     $        - dpdr(i,304) - dpdr(i,303)
      enddo

      return

      end subroutine EvdPdR
      subroutine evdbdr
***********************************************************************
*  The subroutine elminate the 2-body terms in Bowman's approach
***********************************************************************

      implicit double precision (a-h,o-z)

      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)
      
      integer i
      double precision db1dr(6,306) 

C Pass P(0:305) to BM1(1:306)
      do j=1,6
      do i=1,306
        db1dr(j,i)=dpdr(j,i-1)
      enddo
      enddo

C Remove unconnected terms and 2-body terms and pass to B(1:276)
      do j=1,6
      dbdr(j,1)=db1dr(j,4)

      do i=2,4
        dbdr(j,i)=db1dr(j,i+4)
      enddo

      dbdr(j,5)=db1dr(j,10)

      do i=6,11
        dbdr(j,i)=db1dr(j,i+6)
      enddo

      dbdr(j,12)=db1dr(j,19)
      dbdr(j,13)=db1dr(j,21)

      do i=14,26
        dbdr(j,i)=db1dr(j,i+9)
      enddo

      dbdr(j,27)=db1dr(j,37)
      dbdr(j,28)=db1dr(j,39)

      do i=29,53
        dbdr(j,i)=db1dr(j,i+12)
      enddo

      dbdr(j,54)=db1dr(j,67)
      dbdr(j,55)=db1dr(j,69)
      dbdr(j,56)=db1dr(j,71)

      do i=57,97
        dbdr(j,i)=db1dr(j,i+16)
      enddo

      dbdr(j,98)=db1dr(j,115)
      dbdr(j,99)=db1dr(j,117)
      dbdr(j,100)=db1dr(j,119)

      do i=101,166
        dbdr(j,i)=db1dr(j,i+20)
      enddo

      dbdr(j,167)=db1dr(j,188)
      dbdr(j,168)=db1dr(j,190)
      dbdr(j,169)=db1dr(j,192)
      dbdr(j,170)=db1dr(j,194)

      do i=171,272
        dbdr(j,i)=db1dr(j,i+25)
      enddo

      dbdr(j,273)=db1dr(j,299)
      dbdr(j,274)=db1dr(j,301)
      dbdr(j,275)=db1dr(j,303)
      dbdr(j,276)=db1dr(j,305)
      enddo

      return

      end 
C Begin
      block data prmt

      implicit double precision (a-h,o-z)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:305),C(276),B(276)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
     $                  dRdX(6,12),dBdR(6,276)
      common /msprmt/   a,re

C Nonlinear parameters of a(0.9A) and re (1.098A)
      data a    /0.9d0/
      data re   /1.098d0/

C Linear parameters optimized by the weighted-least square fitting
      data C(  1)   /  0.108232000847D+04/
      data C(  2)   / -0.383052557989D+04/
      data C(  3)   /  0.839494388810D+04/
      data C(  4)   / -0.838396243277D+04/
      data C(  5)   / -0.357056476708D+04/
      data C(  6)   /  0.370412755496D+04/
      data C(  7)   /  0.402655542449D+04/
      data C(  8)   /  0.459302228421D+04/
      data C(  9)   /  0.180315790229D+05/
      data C( 10)   / -0.981148332230D+04/
      data C( 11)   /  0.155594405528D+05/
      data C( 12)   /  0.382446322072D+04/
      data C( 13)   /  0.426778294068D+04/
      data C( 14)   /  0.239830807947D+05/
      data C( 15)   / -0.169709044789D+05/
      data C( 16)   / -0.178869989161D+05/
      data C( 17)   / -0.332020758916D+05/
      data C( 18)   / -0.194019326760D+04/
      data C( 19)   /  0.576837149814D+04/
      data C( 20)   / -0.173228888404D+05/
      data C( 21)   /  0.284285880672D+05/
      data C( 22)   / -0.177670404692D+05/
      data C( 23)   /  0.515290193571D+02/
      data C( 24)   / -0.291983524136D+05/
      data C( 25)   / -0.555729302469D+04/
      data C( 26)   / -0.148968347726D+05/
      data C( 27)   /  0.334284288765D+04/
      data C( 28)   / -0.348791539328D+04/
      data C( 29)   /  0.622603586957D+05/
      data C( 30)   / -0.535166008070D+05/
      data C( 31)   /  0.115395863503D+06/
      data C( 32)   /  0.723720465223D+05/
      data C( 33)   /  0.116927755758D+05/
      data C( 34)   /  0.166456259441D+05/
      data C( 35)   / -0.181795815292D+05/
      data C( 36)   /  0.238192165988D+05/
      data C( 37)   / -0.498809470818D+04/
      data C( 38)   /  0.247564729701D+05/
      data C( 39)   /  0.734392788468D+05/
      data C( 40)   /  0.411055945411D+05/
      data C( 41)   /  0.548126860613D+04/
      data C( 42)   /  0.147311445667D+05/
      data C( 43)   /  0.455993650858D+05/
      data C( 44)   /  0.472974769530D+05/
      data C( 45)   / -0.134931567320D+05/
      data C( 46)   / -0.919312454993D+03/
      data C( 47)   /  0.182371842748D+05/
      data C( 48)   / -0.219062324555D+05/
      data C( 49)   / -0.143134239177D+04/
      data C( 50)   / -0.612606279476D+04/
      data C( 51)   /  0.329026080805D+05/
      data C( 52)   / -0.317273136510D+05/
      data C( 53)   /  0.197206345486D+05/
      data C( 54)   / -0.146703151032D+04/
      data C( 55)   / -0.286804318308D+04/
      data C( 56)   /  0.482070970629D+04/
      data C( 57)   / -0.202923650697D+05/
      data C( 58)   /  0.120893278276D+06/
      data C( 59)   / -0.396500084405D+05/
      data C( 60)   /  0.210724336430D+05/
      data C( 61)   / -0.461690024947D+05/
      data C( 62)   /  0.344112211278D+05/
      data C( 63)   / -0.492964436067D+05/
      data C( 64)   / -0.102381662936D+05/
      data C( 65)   / -0.359160759193D+05/
      data C( 66)   / -0.198647773743D+05/
      data C( 67)   /  0.539501764423D+05/
      data C( 68)   / -0.912958065639D+05/
      data C( 69)   /  0.275238240217D+04/
      data C( 70)   / -0.146785486834D+05/
      data C( 71)   / -0.243200816860D+04/
      data C( 72)   /  0.275943650312D+05/
      data C( 73)   /  0.245524867413D+05/
      data C( 74)   / -0.197157661163D+05/
      data C( 75)   / -0.130926931129D+05/
      data C( 76)   /  0.110602206074D+05/
      data C( 77)   / -0.131259741275D+05/
      data C( 78)   / -0.826582073835D+05/
      data C( 79)   / -0.173526822600D+05/
      data C( 80)   / -0.142885462712D+05/
      data C( 81)   /  0.170396271144D+05/
      data C( 82)   /  0.215051573077D+04/
      data C( 83)   / -0.172147788865D+05/
      data C( 84)   / -0.506527993442D+04/
      data C( 85)   / -0.429683705379D+05/
      data C( 86)   / -0.179966097736D+05/
      data C( 87)   /  0.295786064666D+05/
      data C( 88)   / -0.651746790389D+03/
      data C( 89)   /  0.909987620115D+04/
      data C( 90)   /  0.881730414410D+04/
      data C( 91)   /  0.491469905908D+04/
      data C( 92)   / -0.238479435714D+05/
      data C( 93)   /  0.181819093903D+05/
      data C( 94)   /  0.618897436087D+04/
      data C( 95)   / -0.358396457607D+05/
      data C( 96)   /  0.983174669571D+05/
      data C( 97)   / -0.222757600410D+05/
      data C( 98)   /  0.941883986329D+04/
      data C( 99)   / -0.129931322597D+05/
      data C(100)   / -0.518957451256D+04/
      data C(101)   /  0.502150797894D+06/
      data C(102)   / -0.127189053433D+06/
      data C(103)   / -0.471858710777D+05/
      data C(104)   /  0.170093300456D+05/
      data C(105)   / -0.733635609633D+05/
      data C(106)   / -0.378019237891D+04/
      data C(107)   /  0.105960323206D+05/
      data C(108)   / -0.109725341209D+06/
      data C(109)   /  0.599571151662D+05/
      data C(110)   /  0.408194536570D+05/
      data C(111)   / -0.368634396864D+04/
      data C(112)   /  0.762014036496D+05/
      data C(113)   /  0.178190103614D+04/
      data C(114)   / -0.722949186981D+05/
      data C(115)   / -0.781807583040D+05/
      data C(116)   /  0.194483251924D+05/
      data C(117)   /  0.656474114607D+05/
      data C(118)   / -0.128036515935D+05/
      data C(119)   /  0.398287271535D+05/
      data C(120)   / -0.625092472673D+05/
      data C(121)   /  0.172090063917D+05/
      data C(122)   /  0.675522525128D+05/
      data C(123)   /  0.257466408403D+05/
      data C(124)   /  0.974708624146D+04/
      data C(125)   / -0.215479021597D+05/
      data C(126)   / -0.431236813229D+05/
      data C(127)   / -0.304026519580D+05/
      data C(128)   / -0.737078112603D+04/
      data C(129)   / -0.820114526299D+05/
      data C(130)   /  0.599194678229D+04/
      data C(131)   / -0.895063829603D+04/
      data C(132)   /  0.404871762568D+05/
      data C(133)   /  0.121937973858D+03/
      data C(134)   /  0.269158902167D+05/
      data C(135)   / -0.145174235284D+05/
      data C(136)   /  0.890472515492D+03/
      data C(137)   /  0.210324686348D+05/
      data C(138)   /  0.967791167291D+03/
      data C(139)   /  0.184339407545D+04/
      data C(140)   /  0.572045376723D+04/
      data C(141)   /  0.829376767965D+04/
      data C(142)   /  0.234398941256D+05/
      data C(143)   / -0.669617011142D+04/
      data C(144)   /  0.509990643324D+04/
      data C(145)   / -0.563004865548D+03/
      data C(146)   /  0.866048806511D+04/
      data C(147)   /  0.205542074504D+05/
      data C(148)   /  0.306364914983D+05/
      data C(149)   / -0.694109944383D+04/
      data C(150)   / -0.111742753767D+05/
      data C(151)   / -0.176041281355D+03/
      data C(152)   / -0.483046575615D+05/
      data C(153)   /  0.145296629902D+05/
      data C(154)   /  0.585931484562D+04/
      data C(155)   /  0.135429147648D+05/
      data C(156)   /  0.121292140229D+05/
      data C(157)   / -0.906761042822D+04/
      data C(158)   / -0.292627880245D+04/
      data C(159)   /  0.210808391881D+04/
      data C(160)   / -0.102070471820D+05/
      data C(161)   /  0.222598746489D+05/
      data C(162)   / -0.432823625479D+04/
      data C(163)   / -0.227093181806D+04/
      data C(164)   /  0.234458564401D+05/
      data C(165)   / -0.731425566413D+05/
      data C(166)   /  0.838292463406D+04/
      data C(167)   /  0.478249492913D+04/
      data C(168)   / -0.166527064993D+05/
      data C(169)   /  0.158230244351D+05/
      data C(170)   /  0.247039528120D+04/
      data C(171)   /  0.204703921213D+06/
      data C(172)   / -0.193631041360D+06/
      data C(173)   / -0.552819147588D+06/
      data C(174)   / -0.329399675737D+06/
      data C(175)   /  0.120444143730D+06/
      data C(176)   / -0.197630330093D+04/
      data C(177)   / -0.183906826997D+06/
      data C(178)   /  0.683725744752D+05/
      data C(179)   /  0.404594263788D+05/
      data C(180)   /  0.393637852485D+05/
      data C(181)   /  0.473359743314D+05/
      data C(182)   / -0.891076659528D+05/
      data C(183)   / -0.448680528184D+05/
      data C(184)   / -0.367933441170D+03/
      data C(185)   /  0.338358168891D+05/
      data C(186)   / -0.863952552119D+05/
      data C(187)   /  0.653622881890D+05/
      data C(188)   / -0.609133314712D+05/
      data C(189)   / -0.453094488263D+05/
      data C(190)   / -0.448525047759D+05/
      data C(191)   / -0.355816225174D+05/
      data C(192)   /  0.417393783369D+05/
      data C(193)   / -0.133867821370D+05/
      data C(194)   / -0.314938469107D+05/
      data C(195)   / -0.263453916428D+04/
      data C(196)   /  0.328296813300D+05/
      data C(197)   / -0.507351810098D+04/
      data C(198)   / -0.192121276506D+05/
      data C(199)   / -0.450290587850D+05/
      data C(200)   /  0.300953063670D+05/
      data C(201)   / -0.373528166249D+05/
      data C(202)   /  0.100479785157D+05/
      data C(203)   / -0.296954363580D+05/
      data C(204)   /  0.136109142122D+05/
      data C(205)   /  0.126486286466D+05/
      data C(206)   /  0.221929308562D+05/
      data C(207)   / -0.296851538640D+05/
      data C(208)   /  0.102897645572D+05/
      data C(209)   /  0.945556756944D+05/
      data C(210)   /  0.819802511124D+04/
      data C(211)   / -0.198440697793D+05/
      data C(212)   / -0.559611854489D+04/
      data C(213)   / -0.211708607842D+05/
      data C(214)   /  0.127805397273D+05/
      data C(215)   / -0.820494210221D+04/
      data C(216)   /  0.263847446371D+05/
      data C(217)   / -0.124201821079D+05/
      data C(218)   /  0.113594562612D+05/
      data C(219)   / -0.152698768710D+05/
      data C(220)   / -0.183741935900D+05/
      data C(221)   / -0.351149667011D+05/
      data C(222)   /  0.252731998280D+04/
      data C(223)   /  0.107632118740D+05/
      data C(224)   /  0.256967496163D+04/
      data C(225)   /  0.155593585914D+05/
      data C(226)   /  0.372031759011D+05/
      data C(227)   /  0.147888834018D+05/
      data C(228)   / -0.140840658024D+05/
      data C(229)   / -0.638873122608D+04/
      data C(230)   / -0.158719841223D+05/
      data C(231)   /  0.149476667169D+05/
      data C(232)   / -0.115486384117D+05/
      data C(233)   /  0.892787012460D+04/
      data C(234)   / -0.853673138571D+04/
      data C(235)   / -0.112698691434D+05/
      data C(236)   /  0.272741601222D+04/
      data C(237)   /  0.119543853963D+05/
      data C(238)   /  0.473314740470D+04/
      data C(239)   / -0.292984535559D+04/
      data C(240)   /  0.170134527007D+04/
      data C(241)   /  0.535130434607D+04/
      data C(242)   /  0.721489330540D+03/
      data C(243)   / -0.465342829798D+04/
      data C(244)   / -0.553783732812D+05/
      data C(245)   / -0.992508668200D+04/
      data C(246)   /  0.187020363744D+04/
      data C(247)   /  0.246466819967D+04/
      data C(248)   / -0.987652029997D+03/
      data C(249)   / -0.309558455425D+04/
      data C(250)   /  0.975083114532D+04/
      data C(251)   / -0.114148798777D+04/
      data C(252)   / -0.816808630241D+04/
      data C(253)   /  0.102317405015D+05/
      data C(254)   / -0.225335183555D+04/
      data C(255)   /  0.699682112656D+03/
      data C(256)   /  0.275854224075D+04/
      data C(257)   /  0.283424282325D+04/
      data C(258)   /  0.241129666745D+05/
      data C(259)   / -0.983285202228D+04/
      data C(260)   /  0.174060899439D+04/
      data C(261)   / -0.137237094271D+04/
      data C(262)   / -0.244032321366D+05/
      data C(263)   /  0.613665161967D+04/
      data C(264)   /  0.445702123602D+03/
      data C(265)   / -0.340473475708D+04/
      data C(266)   /  0.123141974300D+04/
      data C(267)   / -0.100305108824D+04/
      data C(268)   / -0.287149527130D+04/
      data C(269)   /  0.228487016478D+03/
      data C(270)   / -0.544749868439D+04/
      data C(271)   /  0.154812736178D+05/
      data C(272)   /  0.286973907526D+03/
      data C(273)   /  0.602858747186D+03/
      data C(274)   /  0.583595840631D+04/
      data C(275)   / -0.475565936385D+04/
      data C(276)   / -0.384111331120D+03/

      end
