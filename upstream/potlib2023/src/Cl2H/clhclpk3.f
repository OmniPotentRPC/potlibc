C
      SUBROUTINE PREPOT
C
C   System:          ClHCl
C   Functional form: LEPS (London-Eyring-Polanyi-Sato)
C   Common name:     PK3
C   Reference:       A. Persky and H. Kornweitz
C                    J. Phys. Chem. 91, 5496 (1987)
C
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 1
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
C      The classical potential energy is set equal to zero for the Cl
C      infinitely far from the HCl diatomic and R(HCl) set equal to the
C      HCl equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Cl-H)
C                            R(2) = R(H-Cl)
C                            R(3) = R(Cl-Cl)
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
         PARAMETER (R2     = 1.41421356D0)
         PARAMETER (CKCAU  = 627.5095D0)
         PARAMETER (CANGAU = 0.52917706D0)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
C
      DIMENSION ZPO(3), OP3Z(3), ZP3(3),TZP3(3), TOP3Z(3), DO4Z(3), B(3)
      DIMENSION X(3),COUL(3),EXCH(3)
      SAVE
C
         WRITE (IPTPRT, 600) D, RE, BETA, Z
C
C   CONVERT TO ATOMIC UNITS
      DO  10 I = 1,3
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
C   Initialize the energy in the three asymptotic valleys
      EASYAB = D(1)
      EASYBC = D(2)
      EASYAC = D(3)
C
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
C   Compute the LEPS potential - this is equivalent to the JW OH2 potential.
      DO 21 I = 1,3
            X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
            COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
            EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
21    CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+
     1      (EXCH(3)-EXCH(1))**2)
      ENERGY = COUL(1) + COUL(2) + COUL(3) - (RAD/R2) + D(2)
C
C   Compute the derivatives of the LEPS energy w.r.t.the coordinates
      IF (NDER .EQ. 1) THEN
          SVAL = EXCH(1) + EXCH(2) + EXCH(3)
          DO 22 I = 1,3
                DEDR(I)=0.0D0 
                IF(X(I).LT.1.0D-30) GO TO 22
                TVAL = (3.0D0*EXCH(I)-SVAL)/R2*(OP3Z(I)*X(I)-ZP3(I)) 
                IF(ABS(RAD) .LT. 1.0D-32 .AND. 
     1             ABS(TVAL) .GT. 1.0D-12) THEN 
                   WRITE(IPTPRT, 6000) TVAL, RAD
                ELSE IF(ABS(RAD).GT.1.0D-32) THEN 
                   TVAL = TVAL/RAD 
                END IF 
                DEDR(I)=B(I)*X(I)*(TVAL-ZP3(I)*X(I)+OP3Z(I))  
22       CONTINUE
      ENDIF
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClHCl PK3 potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /,1X,T5,'Bond', T46, 'Cl-H', T58, 'H-Cl', T68, 'Cl-Cl',
     *        /,1X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
6000  FORMAT(/,2X,'In the ClHCl PK3 potential, T, RAD = ',1P,2E15.7,
     1       /,2X,'T/RAD has been set equal to T')  
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT /6/
         DATA NDER /1/
         DATA NDUM /21*0/
         DATA NSURF /0/
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA D    / 106.477D0, 106.477D0, 57.983D0/
         DATA RE   / 1.275D0, 1.275D0, 1.988D0/    
         DATA BETA / 1.868D0, 1.868D0, 2.002D0/
         DATA Z    / 0.220D0, 0.220D0, -0.241D0/
C
         END
C*****
