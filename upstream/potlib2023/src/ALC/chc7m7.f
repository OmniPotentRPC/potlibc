C
      SUBROUTINE PREPOT
C
C   System:          ALC
C   Functional form: eLEPS + Va
C                    Extended LEPS plus an auxiliary barrier term
C   Common name:     CHC7M7
C   Reference:       M. M. Kreevoy, D. Ostovic, 
C                    D. G. Truhlar, and B. C. Garrett
C                    J. Phys. Chem. 90, 3766-3774 (1986)
C
C   Note:            In A-L-C, A = acceptor group, L = H (D), 
C                    and C = methyl donor
C
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Al
C      infinitely far from the LC diatomic and R(LC) set equal to the
C      LC equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Al-L)
C                            R(2) = R(L-C)
C                            R(3) = R(Al-C)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and
C        derivatives of the ground electronic state energy with respect
C        to the coordinates.
C      /PT2CM/ NSURF, NDER, NDUM(21)
C        passes the control flags where
C        NDER  = 0 => no derivatives are computed
C              = 1 => derivatives of the energy for the ground electronic
C                     state with respect to the coordinates are computed
C        NDUM  - not used
C      /PT4CM/ IPTPRT, IDUM(19)
C        IPTPRT passes the FORTRAN unit number used for potential output
C        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to
C        Re(AB), the equilibrium bond length for the AB diatomic.
C        Similarly, the terms EASYBC and EASYAC represent the energies
C        in the BC and the AC valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1
C      IPTPRT           6
C
C*****
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER (ROOT2=1.41421356D0)
         PARAMETER (CKCAL = 627.5095D0)
         PARAMETER (CANGS =   0.529177106D0)
         COMMON /PT1CM/ RR(3),E,DEDRR(3)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
         COMMON /VACOM/ AN1, FNS
C
      DIMENSION ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),DO4Z(3),B(3)
      DIMENSION X(3),COUL(3),EXCH(3)
      SAVE
C
C   Echo the potential parameters to the file linked to FORTRAN unit IPTPRT
      WRITE (IPTPRT, 600) D, RE, BETA, Z
      WRITE (IPTPRT, 610) AN1,FNS
C
      DO  10 I = 1,3
C   CONVERT TO ATOMIC UNITS
      D(I)=D(I)/CKCAL 
      RE(I) = RE(I)/CANGS
      BETA(I) = BETA(I)*CANGS 
C   COMPUTE USEFUL CONSTANTS
      ZPO(I) = 1.0D0 + Z(I)
      OP3Z(I) = 1.0D0 + 3.0D0*Z(I)
      TOP3Z(I) = 2.0D0*OP3Z(I)
      ZP3(I) = Z(I) + 3.0D0
      TZP3(I) = 2.0D0*ZP3(I)
      DO4Z(I) = D(I)/4.0D0/ZPO(I)
   10 B(I) = BETA(I)*DO4Z(I)*2.0D0
C
C   Initialize the energies in the asymptotic valleys
         EASYAB = D(1)
         EASYBC = D(2)
         EASYAC = D(3)
      RETURN
C
      ENTRY POT
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C    Compute the eLEPS part of the potential energy
      DO 21 I = 1,3
            X(I) = EXP(-BETA(I)*(RR(I)-RE(I)))
            COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
            EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
21    CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+
     1      (EXCH(3)-EXCH(1))**2)
      E = COUL(1) + COUL(2) + COUL(3) - (RAD/ROOT2) + D(2)
C
C   Compute the derivative of the eLEPS part of the potential, if NDER=1
      IF (NDER .NE. 1) GO TO 700
      S = EXCH(1) + EXCH(2) + EXCH(3)
      DO 22 I = 1,3
      DEDRR(I) = 0.D0     
      IF(X(I).LT.1.D-30) GO TO 22 
      T= (3.0D0*EXCH(I)-S)/ROOT2*(OP3Z(I)*X(I)-ZP3(I))
C
C   If RAD is very small compared to T print a warning to the file linked 
C   to FORTRAN IPTPRT and set T/RAD equal to T.
      IF(ABS(RAD).LT.1.D-32.AND.ABS(T).GT.1.D-12) THEN
        WRITE(IPTPRT,6000) T,RAD     
      ELSE IF(ABS(RAD).GT.1.D-32) THEN
        T = T/RAD 
      END IF
C
      DEDRR(I) = B(I)*X(I)*(T        
     1          -ZP3(I)*X(I)+OP3Z(I)) 
22    CONTINUE
700   CONTINUE
C
C   Compute the auxiliary barrier term, Va, where 
C   Va = a{[(R1-R2)*(R2-R3)*(R3-R1)]^2}*exp[-alpha*{(R1+R2+R3)^3}]
      R=RR(1)+RR(2)+RR(3)
      R2=R*R
      R3=R2*R
      EXNS=EXP(-FNS*R3)
      ENS = 0.D0
      WNT=(RR(1)-RR(2))*(RR(2)-RR(3))*(RR(3)-RR(1))
      WN=ABS(WNT)
      WN2=WNT*WNT
      ENS=(AN1*WN2)*EXNS
C    Add the energy terms together
      E=E+ENS
C
C    Compute the derivative of the Va term with respect to the 
C    coordinates, if NDER = 1
      IF (NDER .EQ. 1) THEN
      DELTA =-1.D0
      IF (WN.EQ.WNT) DELTA =1.D0
      ENSP1 = 0.D0
      ENSP2 = 0.D0
      ENSP3 = 0.D0
      ENSPWN=(AN1*WN*2.D0)*EXNS
      ENSPR=(-3.D0)*FNS*R2*ENS
C DENS/DXI=(DENS/DWN)(DWN/DXI)+(DENS/DR)(DR/DXI)
      X21=RR(1)*RR(1)
      X22=RR(2)*RR(2)
      X23=RR(3)*RR(3)
      WNP1=(2.D0*RR(1)*(RR(3)-RR(2))+X22-X23)*DELTA
      WNP2=(2.D0*RR(2)*(RR(1)-RR(3))+X23-X21)*DELTA
      WNP3=(2.D0*RR(3)*(RR(2)-RR(1))+X21-X22)*DELTA
      ENSP1=ENSPWN*WNP1+ENSPR
      ENSP2=ENSPWN*WNP2+ENSPR
      ENSP3=ENSPWN*WNP3+ENSPR
      DEDRR(1)=DEDRR(1)+ENSP1
      DEDRR(2)=DEDRR(2)+ENSP2
      DEDRR(3)=DEDRR(3)+ENSP3
      ENDIF
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'Al + LC CHC7M7 potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /,2X,T5,'Bond', T46, 'AI-L', T58, 'L-C', T69, 'C-AI',
     *        /,2X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT(/,2X,T5,'Parameters for the auxiliary barrier term ',
     *       '(a.u.)',/,2X,T20,'A = ',F15.8,T45,'alpha = ',T53,F15.8,
     *       //,1X,'*****')
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
6000  FORMAT(/,2X,'In the potential routine: T, RAD = ',1P,2E15.7,
     1       /,2X,'T/RAD has been set equal to T')  
C
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /VACOM/ AN1, FNS
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER, NSURF /6, 1, 0/
         DATA NDUM /21*0/ 
C
C   Initialize the eLEPS potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA DE/ 80.00000D0, 73.00000D0, 80.00000D0/
         DATA RE/ 1.10000D0, 1.10000D0, 1.54000D0/
         DATA BETA/ 2.07800D0, 2.17500D0, 2.07800D0/
         DATA Z/ 0.00000D0, 0.00000D0, -0.18900D0/
C   Initialize the potential parameters for the auxiliary barrier term; 
C   the parameters are in Hartree atomic units.
         DATA AN1, FNS / 0.026D0, 0.004D0/
         END
C*****
C
