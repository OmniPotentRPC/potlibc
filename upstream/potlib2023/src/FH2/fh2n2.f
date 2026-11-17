C
      SUBROUTINE PREPOT
C
C   System:          FH2
C   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
C   Common name:     no. 2
C   Reference:       D. G. Truhlar, B. C. Garrett, and N. C. Blais
C                    J. Chem. Phys. 80, 232 (1984).  
C
C   Notes:           In this version of the FH2 potential the AB and AC
C                    Sato parameters are functions of PHI, where PHI is
C                    the angle of A-to-BC vector with the BC axis.
C
C   THE DERIVATIVES OF THE SATO PARATMETERS WITH RESPECT TO PHI HAVE
C   NOT BEEN CODED FOR THIS POTENTIAL!
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in data statements and in the
C   block data subprogram PPARMS.
C   Coordinates, potential energy, and derivatives are passed 
C   through the common block POTCM:
C                  COMMON /POTCM/ R(3), ENERGY, DEDR(3)
C   The potential energy in the three asymptotic valleys are 
C   stored in the common block ASYCOM:
C                  COMMON /ASYCOM/ EASYAB, EASYBC, EASYAC
C   The potential energy in the AB valley, EASYAB, is equal to the potential 
C   energy of the H "infinitely" far from the FH diatomic, with the 
C   FH diatomic at its equilibrium configuration.  Similarly, the terms 
C   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys, 
C   respectively.
C   All the information passed through the common blocks POTCM and ASYCOM 
C   is in in Hartree atomic units.  
C
C   This potential is written such that:
C                  R(1) = R(F-H)
C                  R(2) = R(H-H)
C                  R(3) = R(H-F)
C   The zero of energy is defined at F "infinitely" far from the H2 diatomic.
C
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and
C        derivatives of the ground electronic state energy with respect
C        to the coordinates.
C      /PT2CM/ NSURF, NDER, NDUM(21)
C        passes the control flags where
C        NSURF = 0    ground electronic surface
C        NDER  = 0 => no derivatives are computed
C              = 1 => derivatives of the energy for the ground electronic
C                     state with respect to the coordinates are computed
C        LCHI  - Flag that controls printing in the POT and FPHI subroutines
C              = 0, do not print extra info
C              = 1, print extra info to the file linked to FORTRAN unit IPRT
C              - LCHI should be nonzero for debugging
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
C      LCHI             0
C      IPTPRT           6
C
C*****
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER (CKCAU  = 627.5095D0)
         PARAMETER (CANGAU = 0.52917706D0)
         PARAMETER (R2     = 1.41421356D0)
         COMMON /PT1CM/ R(3), E, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, LCHI, NDUM(20)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
C
         COMMON /COSCOM/ COSPHI
         COMMON /PHICOM/ DEG1,Z0,FDER,DEG
         COMMON /Z3COM/ Z(3)
         DIMENSION D(3), RE(3), BETA(3), ZPO(3), OP3Z(3), ZP3(3),
     1             TZP3(3), TOP3Z(3), DO4Z(3), B(3)
         DIMENSION X(3),COUL(3),EXCH(3)
         DATA D /141.196D0, 109.449D0, 141.196D0/
         DATA RE /0.917D0, 0.7419D0, 0.917D0/
         DATA BETA /2.21870D0, 1.942D0, 2.2187D0/
         SAVE
C
C   Echo potential parameters to FORTRAN unit IPRT
         WRITE (IPRT, 600) D, RE, BETA, Z
         WRITE (IPRT ,610) DEG1
C
      Z0=Z(1)
      DO  10 I = 1,3
C   CONVERT TO ATOMIC UNITS
      D(I)=D(I)/CKCAU
      RE(I) = RE(I)/CANGAU
      BETA(I) = BETA(I)*CANGAU
C   COMPUTE USEFUL CONSTANTS
      ZPO(I) = 1.0D0 + Z(I)
      OP3Z(I) = 1.0D0 + 3.0D0*Z(I)
      TOP3Z(I) = 2.0D0*OP3Z(I)
      ZP3(I) = Z(I) + 3.0D0
      TZP3(I) = 2.0D0*ZP3(I)
      DO4Z(I) = D(I)/4.0D0/ZPO(I)
   10 B(I) = BETA(I)*DO4Z(I)*2.0D0
C
C    Set the values of the classical energy in the three asymptotic valleys
             EASYAB = D(1)
             EASYBC = D(2)
             EASYAC = D(3)
C
600   FORMAT (/, 2X, T5, 'PREPOT has been called for the FH2 ',
     *                   'extended-LEPS potential energy surface',
     *        /, 2X, T5, 'called no. 2',
     *        /, 2X, T5, 'Potential energy surface parameters:',
     *        /, 2X, T5, 'Bond', T47, 'F-H', T58, 'H-H', T69, 'H-F',
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT (/,2X,T5,'The above Sato parameters apply only for PHI ',
     *                'less than',F10.3,1X,'degrees',/)
C
      RETURN
C
      ENTRY POT
C
C   Check the value of NSURF
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
      CALL PHICAL (PHI)
      IF (LCHI.NE.0) WRITE (IPRT,620) R(1),R(2),R(3),PHI
620   FORMAT(2X,3F8.3,T35,F10.4)
      Z(1)=FPHI(PHI)
      Z(3)=Z(1)
      DO 11 I=1,3,2
      ZPO(I) = 1.0D0 + Z(I)
      OP3Z(I) = 1.0D0 + 3.0D0*Z(I)
      TOP3Z(I) = 2.0D0*OP3Z(I)
      ZP3(I) = Z(I) + 3.0D0
      TZP3(I) = 2.0D0*ZP3(I)
   11 DO4Z(I) = D(I)/4.0D0/ZPO(I)
      S = 0.0D0
      DO 21 I = 1,3
      X(I) = EXP(-BETA(I)*(R(I)-RE(I)))
      COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
      EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
   21 S = S+EXCH(I)
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+
     1      (EXCH(3)-EXCH(1))**2)
      E = -RAD/R2
         IF (NDER .EQ. 1) THEN
             DO 22 I = 1,3
                   DEDR(I) = 0.D0 
                   IF(X(I).LT.1.D-30) GO TO 22
                   T= (3.0D0*EXCH(I)-S)/R2*(OP3Z(I)*X(I)-ZP3(I))
                   IF(RAD.NE.0.D0) GO TO 2222                         
                   IF(T.EQ.0.D0) GO TO 2223
                   WRITE(IPRT,6000) T,RAD 
                   GO TO 2223 
 2222              T = T/RAD 
 2223              CONTINUE  
                   DEDR(I) = B(I)*X(I)*(T-ZP3(I)*X(I)+OP3Z(I))
22           CONTINUE
         ENDIF
         E = E + COUL(1) + COUL(2) + COUL(3) + D(2)
         IF (DEG.LE.DEG1) RETURN
C   DERIVATIVE CHANGES FROM LEPS TO LEPS2 ARE TO BE ADDED HERE.
C   THEY ARE COMPUTED USING FDER FROM FPHI AND COSPHI FROM COSCOM.
C   FDER IS DERIVATIVE OF SATO PARAMETER WITH RESPECT TO PHI.
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
6000  FORMAT(/,2X,'In the eLEPS potential T, RAD = ',1P,2E15.7,
     1       /,2X,'T/RAD has been set equal to T')  
C
      RETURN
      END
      FUNCTION FPHI(PHI)
C   FH2 NEW SURFACE NO. 2   DON TRUHLAR
C   AB,AC SATO PARAMETERS ARE Z(1)=Z(3) FOR DEG.LE.DEG1.
C   DEG1 MUST BE SET IN PREPOT.
C   Z0 IS SET EQUAL TO Z(1)=Z(3) IN PREPOT.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT2CM/ NSURF, NDER, LCHI, NDUM(20)
      COMMON /PT4CM/ IPRT, IDUM(19)
      COMMON /PHICOM/ DEG1,Z0,FDER,DEG
      COMMON /Z3COM/ Z(3)
      DEG=57.29577951D0*PHI
      IF (LCHI.NE.0) WRITE (IPRT,100) DEG
  100 FORMAT (2X,T35,F10.4)
      IF (DEG.GT.180.0D0)DEG=DEG-180.0D0
      IF (DEG.LT.0.0D0)  DEG=-DEG
      IF (DEG.GT.90.0D0) DEG=180.0D0-DEG
      IF (DEG.LT.0.0D0) STOP 89
      IF (DEG.GT.DEG1) GO TO 200
      IF (LCHI.NE.0) WRITE (IPRT,100) DEG
      FPHI=Z0
      FDER=0.0D0
      RETURN
  200 IF (DEG.GT.20.0D0) GO TO 300
      DEGDIF=DEG-10.0D0
      DD2=DEGDIF*DEGDIF
      DD3 = DEGDIF*DD2
      FPHI=0.175D0 +1.15D-4*DD2    -6.0D-6*DD3
      FDER=         2.30D-4*DEGDIF -1.8D-5*DD2
      RETURN
  300 IF (DEG.GT.60.0D0) GO TO 400
      DEGDIF=DEG-41.0D0
      FPHI=0.191D0 +0.0005D0*DEGDIF
      FDER=         0.0005D0
      RETURN
  400 IF (DEG.GT.90.0D0) GO TO 500
      SINDEG=SIN(PHI)
      S2=SINDEG*SINDEG
      FPHI=0.17569D0 +3.308D-2*S2
      COSDEG=COS(PHI)
      FDER=           6.616D-2*SINDEG*COSDEG
      RETURN
  500 STOP 90
      END
      SUBROUTINE PHICAL (PHI)
C   FH2 NEW SURFACE NO. 2   DON TRUHLAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /COSCOM/ COSPHI
      COMMON /PT1CM/ R(3),E,DEDR(3)
      HR2=0.5D0*R(2)
      R12=R(1)*R(1)
      R22=R(2)*R(2)
      R32=R(3)*R(3)
      HR22=0.25D0*R22
      COS23=(R12-R22-R32)/(-2.0D0*R(2)*R(3))
      FHH2=HR22+R32-2.0D0*HR2*R(3)*COS23
      FHH=SQRT(FHH2)
      COSPHI=(R32-HR22-FHH2)/(-2.0D0*HR2*FHH)
      IF (COSPHI.GT.1.0D0) COSPHI=1.0D0
      IF (COSPHI.LT.-1.0D0) COSPHI=-1.0D0
      PHI=ACOS(COSPHI)
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, LCHI, NDUM(20)
         COMMON /PHICOM/ DEG1,Z0,FDER,DEG
         COMMON /Z3COM/ Z(3)
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NDER /1/
         DATA NSURF /0/
         DATA LCHI /0/
         DATA NDUM /20*0/
C
C   Initialize potential variables
         DATA DEG1 /10.0D0/
         DATA Z /0.175D0, 0.104D0, 0.175D0/
C
         END
C*****
