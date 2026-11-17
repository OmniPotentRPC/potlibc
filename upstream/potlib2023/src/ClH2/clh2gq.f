C
      SUBROUTINE PREPOT
C
C   System:          ClH2
C   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)
C                    plus E3C (a three-center term).
C   Common name:     GQ
C   Reference:       D. W. Schwenke, S. C. Tucker, R. Steckler, F. B. Brown,
C                    G. C. Lynch, D. G. Truhlar, and B. C. Garrett
C                    J. Chem. Phys. 90, 3110-3120 (1989)
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
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Cl-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Cl-H)
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
         PARAMETER (CKCAU = 627.5095D0)
         PARAMETER (CANGAU = 0.52917706D0)
         PARAMETER (R2     = 1.41421356D0)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
         COMMON /V3CCOM/ CON, ALP, APRIME, BPRIME, Q2, Q4
C
      DIMENSION ZPO(3),OP3Z(3),ZP3(3),TZP3(3), TOP3Z(3),DO4Z(3),B(3)
      DIMENSION X(3),COUL(3),EXCH(3)
      SAVE
C
      WRITE (IPTPRT, 600) D, RE, BETA, Z
      WRITE (IPTPRT, 610) CON, APRIME, BPRIME, ALP, Q2, Q4
C
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
C   Set the values of the variables EASYBC, EASYAC, and EASYBA
         EASYAB = D(1)
         EASYBC = D(2)
         EASYAC = D(3)
C
      RETURN
C
      ENTRY POT
C   Check the value of NDER
         IF (NDER .GT. 1) THEN 
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C   Compute the eLEPS part of the potential
C
      DO 20 I = 1,3
            X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
            COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
            EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
20    CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+
     1      (EXCH(3)-EXCH(1))**2)
      ENERGY = COUL(1) + COUL(2) + COUL(3) - (RAD/R2) + D(2)
C
C   Compute the derivatives of the eLEPS energy w.r.t. the bond lengths
C
         IF (NDER .EQ. 1) THEN
             SVAL = EXCH(1) + EXCH(2) + EXCH(3)
             DO 30 I = 1, 3
                   DEDR(I) = 0.D0
                   IF (X(I) .LT. 1.D-30) GO TO 30
                   TVAL = (3.0D0*EXCH(I)-SVAL)/R2*(OP3Z(I)*X(I)-ZP3(I)) 
                   IF (ABS(RAD) .LT. 1.D-32 .AND. 
     1                 ABS(TVAL) .GT. 1.D-12) THEN
                       WRITE (IPTPRT, 6000) TVAL, RAD
                   ELSE
                       TVAL = TVAL/RAD
                   ENDIF
                   DEDR(I) = B(I)*X(I)*(TVAL-ZP3(I)*X(I)+OP3Z(I)) 
30           CONTINUE
         ENDIF
C
C   Compute the three-center term
C
C   F1VAL is equal to f1 in the literature reference, 
C   similarly F2VAL = f2, and F3VAL = f3
      F1VAL = EXP(-ALP*(R(1)-R(3))**2)
      F2VAL = EXP(-APRIME*(R(1)+R(3))**4) 
      F3VAL = EXP(-BPRIME*(R(1)+R(3)-R(2))**2)
      E3C   = CON * F1VAL * F2VAL * F3VAL
C
C   Compute the angular dependent term, f4, which is called F4VAL
      DIST = R(1)**2 + R(3)**2 - R(2)**2
      CTHETA = DIST/(2.0D0*R(1)*R(3))
      IF (CTHETA .LT. -1.0D0)THEN
         THETA = ACOS(-1.0D0)
      ELSE IF (CTHETA .GE. 1.0D0)THEN
         THETA = 0.0D0
      ELSE
         THETA = ACOS(CTHETA)
      END IF
      STHETA = SIN(THETA)
C
      F4VAL  = (Q2*(STHETA**2))+(Q4 * (STHETA**4)) + 1.0D0 
      E3C    = E3C * F4VAL
      ENERGY = E3C + ENERGY
C
C   Compute the derivatives of the three-center term
C
         IF (NDER .EQ. 1) THEN
C   Recall DIST = R(1)**2 + R(3)**2 - R(2)**2
             RSUM3 = (R(1) + R(3))**3
             RDIF  = R(1) - R(3)
             RCO   = R(1) + R(3) - R(2)
             ABCVAL = E3C/F4VAL
             R1P5 = R(1)**5
             R1P4 = R(1)**4
             R1P3 = R(1)**3
             R1P2 = R(1)**2
             R3P5 = R(3)**5
             R3P4 = R(3)**4
             R3P3 = R(3)**3
             R3P2 = R(3)**2
             TRM1 = 4.0D0*APRIME*RSUM3
             TRM2 = 2.0D0*BPRIME*RCO
             TRM3 = 2.0D0*ALP*RDIF
             ANG11 = -(DIST**4)/(4.0D0*R1P5*R3P4)
             ANG12 =  (DIST**3)/(2.0D0*R1P3*R3P4)
             ANG13 =  (DIST**2)/(R1P3*R3P2)
             ANG14 = (-2.0D0*DIST)/(R(1)*R3P2)
             ANG1 = Q4 * (ANG11 + ANG12 + ANG13 + ANG14) +  
     1              Q2 * ((ANG13/2.0D0) + (ANG14/2.0D0))
             ANG31 = -(DIST**4)/(4.0D0*R1P4*R3P5)
             ANG32 =  (DIST**3)/(2.0D0*R1P4*R3P3)
             ANG33 =  (DIST**2)/(R1P2*R3P3)
             ANG34 = (-2.0D0*DIST)/(R1P2*R(3))
             ANG3 = Q4 * (ANG31 + ANG32 + ANG33 + ANG34) + 
     1              Q2 * (ANG33/2.0D0 + ANG34/2.0D0)
             ANG21 = -R(2)*(DIST**3)/(2.0D0*R1P4*R3P4)
             ANG22 =  R(2)*2.0D0*DIST/(R1P2*R3P2)
             ANG2 = Q4 * (ANG21 + ANG22) +  
     1              Q2 * (ANG22/2.0D0)
             D3CDR1 = E3C * (-TRM1 - TRM2 - TRM3) + ABCVAL * ANG1
             D3CDR2 = E3C * TRM2 + ABCVAL * ANG2
             D3CDR3 = E3C * (-TRM1 - TRM2 + TRM3) + ABCVAL * ANG3
             DEDR(1) = DEDR(1) + D3CDR1                                         
             DEDR(2) = DEDR(2) + D3CDR2
             DEDR(3) = DEDR(3) + D3CDR3
         ENDIF
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClH2 GQ potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /,1X,T5,'Bond', T46, 'Cl-H', T58, 'H-H', T69, 'H-Cl',
     *        /,1X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Sato parameters:',
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT(/, 1X, T5, 'Parameters for the three-center term:',
     * /,1X,T5,'J', T10,'=',T12,E13.5,T40,'ALFP', T46,'=',T48,E13.5,
     * /,1X,T5,'BETP',T10,'=',T12,E13.5,T40,'ALFT',T46,'=',T48,E13.5,
     * /,1X,T5,'Q2', T10,'=',T12,E13.5,T40,'Q4', T46,'=',T48,E13.5,
     * //,1X,'*****')
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
6000  FORMAT(/,1X,'In the ClH2 potential GQ, T, RAD = ',1P,2E15.7,
     *       /,1X,'T/RAD has been set equal to T')
C
      RETURN
      END                                                 
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
         COMMON /V3CCOM/ CON, ALP, APRIME, BPRIME, Q2, Q4
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER, NSURF /6, 1, 0/
         DATA NDUM /21*0/ 
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA D    / 106.447D0, 109.458D0, 106.447D0 /
         DATA RE   / 1.2732D0, 0.74127D0, 1.2732D0 /
         DATA BETA / 1.8674D0, 1.94130D0, 1.8674D0 /
         DATA Z    / 0.187D0, 0.167D0, 0.187D0 /
C
C   Initialize the parameters for the three-center term
         DATA CON, ALP / 0.07672D0, 0.3D0/
         DATA APRIME, BPRIME / 0.000863D0, 0.15D0/
         DATA Q2, Q4 / 0.0D0, 1.725D0/
C
         END
C*****
