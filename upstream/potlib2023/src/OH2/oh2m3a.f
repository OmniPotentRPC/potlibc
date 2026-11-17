C
      SUBROUTINE PREPOT
C
C   System:          OH2
C   Common name:     M3na, symmetrized M2 potential (A'' surface)
C   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential
C                    for collinear geometeries plus anti-Morse 
C                    bending potential. Gamma is a function of the
C                    collinear OH distance.
C
C   Data file:       potoh2m3a.dat
C   Reference:       unpublished
C   Cross reference: T. Joseph, D. G. Truhlar, and B. C. Garrett
C                    J. Chem. Phys. 88, 6982-6990 (1988).
C
C   Notes:           This is the A'' surface for OH2.
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PPARMS.
C   The potetinal parameters for the collinear part of the potential are 
C   read in from the file 'potoh2m3n.dat' which is linked to FORTRAN 
C   unit IREAD.
C
C        This potential is written such that:
C                       R(1) = R(O-H)
C                       R(2) = R(H-H)
C                       R(3) = R(H-O)
C   The classical potential energy is set equal to zero for the O
C   infinitely far from the equilibrium H2 diatomic.
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
C        NDUM  - not used
C      /PT4CM/ IPTPRT, IREAD, IDUM(18)
C        IPTPRT passes the FORTRAN unit number used for potential output
C        IREAD  passed the FORTRAN unit number used for data file input
C        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to
C        Re(AB), the equilibrium bond length for the AB diatomic.
C        Similarly, the terms EASYBC and EASYAC represent the energies 
C        in the BC and the AC valleys, respectively.
C
C   Potential parameters' default settings
C                  Variable            Default value
C                  NSURF               0
C                  NDER                1
C                  IPRT                6
C                  IREAD               4
C
C
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ R1,R2,R3,V,DER(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         DIMENSION DE1(3), DE2(3), DCOSDR(3), DFDR(3)
         DATA ONE/1.0D0/
         SAVE
C
      HPI = 2.0D0*ATAN(ONE)
      WRITE (IPRT, 600)
C   Open the file containing the potential parameters.
      OPEN (UNIT=IREAD, FILE='potoh2m3a.dat', FORM='FORMATTED',
     *      STATUS='OLD', ERR=1000)
      READ (IREAD, *) IOPSW, A
      IOPSW = MAX(0, IOPSW)
      A = MAX(AMN, A)
      WRITE (IPRT, 601) IOPSW, A
      CALL PREM2
C   Close the potential data file
      CLOSE (UNIT=IREAD)
      ASQ = A*A
      RETURN
1000  WRITE (IPRT, 6000)
      STOP 'PREPOT 1'
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the A'''' OH2 ',
     *                'potential energy surface M3na',
     *       //,2X,T5,'Potential energy surface parameters:')
601   FORMAT (/,2X,T5,'Switching function parameter, n=', T50,I5,
     *        /,2X,T5,'Switching function, Fn(chi) = ',
     *                'sin**2((pi/2)*Fn-1(chi))',
     *        /,2X,T5,'where F0(chi) = sin**2(chi/2)',
     *        /,2X,T5,'cos(chi) = 0.5*(r1**2 - r3**2)/r2*rho',
     *        /,2X,T5,'rho**2 = R**2 + a**2,  a=',T50,1PE15.7)
6000  FORMAT (/,2X,T5,'Error opening potential input data file')
C
      ENTRY POT
C
C   Check the values of NSURF and NDER for validity.
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
C
C  CHECK FOR SMALL VALUES OF R
      IF (R1 .LT. 0.001D0 .OR. R2. LT. 0.001D0 .OR. 
     *    R3 .LT. 0.001D0) THEN
         V = 1.D10
         IF (NDER .NE. 1) RETURN
         DER(1) = -V
         DER(2) = -V
         DER(3) = -V
         RETURN
      END IF
C
C    EVALUATE SWITCHING FUNCTION AND ITS DERIVATIVE
      R12 = R1*R1
      R22 = R2*R2
      R32 = R3*R3
      BIGR2 = 0.5D0*(R12 - 0.5D0*R22 + R32)
      BIGR2 = MAX(ZERO, BIGR2)
      RHO2 = BIGR2 + ASQ
      FRHO2 = 4.0D0*RHO2
      RHO = SQRT(RHO2)
      DR13 = R12 - R32
      COSCHI = DR13/(2.0D0*R2*RHO)
C    DERIVATIVES OF COS
      IF (NDER .EQ. 1) THEN
          T = 0.25D0/(R2*RHO2)
          DCOSDR(2) = -COSCHI*(FRHO2 - R22)*T
          T = T/RHO
          DCOSDR(1) = R1*(FRHO2 - DR13)*T
          DCOSDR(3) = -R3*(FRHO2 + DR13)*T
      ENDIF
C    EVALUATE SWITCHING FUNCTION
      G = COSCHI
      DGDX = 1.0D0
      IF (IOPSW .GT. 0) THEN
         DO 20 I = 1,IOPSW
            T = HPI*G
            G = SIN(T)
            DGDX = HPI*COS(T)*DGDX
   20    CONTINUE
      END IF
      F = 0.5D0*(1.0D0 - G)
      IF (NDER .EQ. 1) THEN
          DFDX = -0.5D0*DGDX
          DO 30 I = 1,3
                DFDR(I) = DFDX*DCOSDR(I)
30        CONTINUE
      ENDIF 
      OMF = 1.0D0 - F
      IF (F .GT. 1.D-10) THEN
C    EVALUATE POTENTIAL FOR O + H2 DIRECTION
         CALL POTM2
         V1 = V
         IF (NDER .EQ. 1) THEN
             DO 50 I = 1,3
50                 DE1(I) = DER(I)
         ENDIF
      ELSE
         V1 = 0.0D0
         IF (NDER .EQ. 1) THEN
             DO 55 I = 1,3
55                 DE1(I) = 0.0D0
         ENDIF
      END IF
      IF (ABS(F-0.5D0) .LT. 1.D-10) THEN
         V2 = V
         IF (NDER .EQ. 1) THEN
             DE2(1) = DER(3)
             DE2(2) = DER(2)
             DE2(3) = DER(1)
         ENDIF 
      ELSE IF (OMF .GT. 1.D-10) THEN
C    INTERCHANGE R1 AND R3 AND RECOMPUTE POTENTIAL
         T = R1
         R1 = R3
         R3 = T
         CALL POTM2
         V2 = V
         IF (NDER .EQ. 1) THEN
             DE2(1) = DER(3)
             DE2(2) = DER(2)
             DE2(3) = DER(1)
         ENDIF
         T = R1
         R1 = R3
         R3 = T
      ELSE
         V2 = 0.0D0
         IF (NDER .EQ. 1) THEN
             DO 60 I = 1,3
60                 DE2(I) = 0.0D0
         ENDIF
      END IF
C
      V = F*V1 + OMF*V2
      IF (NDER .NE. 1) RETURN
      DELV = V1 - V2
      DO 300 I = 1,3
             DER(I) = F*DE1(I) + OMF*DE2(I) + DELV*DFDR(I)
300   CONTINUE
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,
     *       /,2X,T5,'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C
      SUBROUTINE PREM2
C
C   The subprogram PREM2 determines the bending correction.
C
C   The potential parameters for the collinear potential are 
C   initialized in the block data subprogram PPARMS.
C   The potential parameters for the bending potential are read
C   in in the subprogram PREM22.  The subprograms PREM2 and PREM22
C   do not need to be called in the MAIN program.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1
      DIMENSION DBCOOR(2)
      COMMON /PT1CM/ RAB,RBC,RAC,ENERGY,DER(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(21)
      COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      COMMON /AMBPCM/ D(3), BETA, REQ(3), A, GAMP(3), RBB
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)
      PARAMETER (BOHR = .52917706D0)
C   Conversion for kcal to Hartree (kcal/Hartree)
      PARAMETER (HARTRE = 627.5095D0)
      DATA TOL /1.D-5/
      SAVE
C
C   Initialize the potentials for the bending potential
      CALL PREM22
C
C   Echo the potential parameters for the bending correction.
      WRITE (IPRT,600) REQ, D
      WRITE (IPRT,601) BETA, A, RBB, GAMP
C
C       CONVERT TO ATOMIC UNITS
      DO 50 IT = 1,3
         REQ(IT) = REQ(IT)/BOHR
         D(IT) = D(IT)/HARTRE
   50 CONTINUE
      BETA = BETA * BOHR
      DE = D(1)
      A = A/BOHR
      RHH = 0.74144D0 / BOHR
      RBB = RBB /BOHR
C
C   Initialize the energy in the asymptotic valleys
      EASYAB = D(1)
      EASYBC = D(2)
      EASYAC = D(3)
C
600   FORMAT(/,2X,T5,'Parameters for the bending correction',
     *       /,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',
     *       /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *       T44, F10.5, T55, F10.5, T66, F10.5,
     *       /,2X,T5,'Dissociation energy (kcal/mol):',
     *       T44, F10.5, T55, F10.5, T66, F10.5)
601   FORMAT(2X,T5,'Morse beta parameter (Angstroms**-1):', T44,1PE15.8,
     *     /,2X,T5,'Pauling parameter (Angstroms):',T44,1PE15.8,
     *     /,2X,T5,'RBB (Angstroms):',T44,1PE15.8,
     *     /,2X,T5,'Gamma for the bending correction:',T44,1PE15.8,
     *          T60,1PE15.8,/,2X,T44,1PE15.8)
C
      RETURN
C
      ENTRY POTM2
C
      IF (RAB .GT. TOL .AND. RBC .GT. TOL .AND. RAC .GT. TOL) THEN
         RACCOL = RAB + RBC
         BCOOR = 0.0D0
         IF (NDER .EQ. 1) THEN
             DBCOOR(1) = 0.0D0
             DBCOOR(2) = 0.0D0
         ENDIF
         DAC = 0.0D0
         ROB = 1.0D0
         IF(ABS(RHH-RBB).GT.1.D-10) THEN
C
C     FIRST CALCULATE SOME USEFUL NUMBERS
C        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
C        AB AND BC RESPECTIVELY. (EQUATION 2.2)
            NAB = EXP((REQ(1) - RAB)/A)
            IF (NAB.LT.1.D-15) NAB = 1.D-15
            NBC = EXP((REQ(2) - RBC)/A)
            IF (NBC.LT.1.D-15) NBC = 1.D-15
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
C       SEE EQUATION 2.7A, 2.7B
            C = 1.D0 - NAB/NBC
            IF (ABS(C) .LT. 1.D-14) THEN
               NAB1=.5D0
               NBC1=.5D0
            ELSE
               C = C/NAB
               NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)
               NBC1 = 1.D0 - NAB1
            END IF
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
C          FIRST TRAP OUT ANY ZERO ARGUEMENTS
            IF (NAB1*NBC1 .LE. 0.0D0) THEN
               ROB = 1.
            ELSE
               STUFF = A * LOG(NAB1*NBC1)
               T = RHH - STUFF
               ROB = 1.0D0
               IF(ABS(T).GT.1.D-15) ROB = (RBB - STUFF)/T
            END IF
         END IF
C
C        CALCULATE BENDING CORRECTION
C           EVALUATE GAMMA AND DERIVATIVE
         EX = EXP(RACCOL*(GAMP(2) + GAMP(3)*RACCOL))
         GAMMA = GAMP(1) + EX
         IF (NDER .EQ. 1) DGAM = (GAMP(2) + 2.D0*GAMP(3)*RACCOL)*EX
         GAMMA = GAMMA*DE
         IF (NDER .EQ. 1) DGAM = DGAM*DE
C           CALCULATE V(RAC) TERM
         EX = EXP(-BETA*(RAC-REQ(3))/ROB)
         EX1 = EX*(1.D0 + 0.5D0*EX)
         IF (NDER .EQ. 1) DEX1 = EX*(1.D0 + EX)
C           CALCULATE V(RAB+RBC) TERM
         EX = EXP(BETA*(REQ(3)-RBC-RAB)/ROB)
         EX2 = EX*(1.D0 + 0.5D0*EX)
         IF (NDER .EQ. 1) DEX2 = EX*(1.D0 + EX)
C        HERE IS THE BENDING CORRECTION
         BCOOR = GAMMA*(EX1 - EX2)
C        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
C        OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
C          CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO RAC)
         GAMMA = GAMMA*BETA/ROB
         IF (NDER .EQ. 1) THEN
             DAC = -GAMMA*DEX1
C          CALCULATE CORRECTIONS TO DAB, DBC
             DBCOOR(1) = GAMMA*DEX2 + DGAM*(EX1-EX2)
             DBCOOR(2)= DBCOOR(1)
         ENDIF
C       NOW CALCULATE THE COLLINEAR ENERGY
C          STORE RAC IN RACTMP THEN SET RAC TO THE COLLINEAR
C          GEOMETRY AND CALL POTM22
         RACTMP = RAC
         RAC = RBC + RAB
         CALL POTM22
C
         RAC = RACTMP
C
         ENERGY = ENERGY + BCOOR
C
         IF (NDER .EQ. 1) THEN
             DO 300 IT = 1,2
                    DER(IT) = DBCOOR(IT) + DER(IT) + DER(3)
300          CONTINUE
             DER(3) = DAC
         ENDIF
      ELSE
        ENERGY = 1.D30
        IF (NDER .NE. 1) RETURN
        DER(1) = 0.0D0
        IF (RAB .LE. TOL) DER(1) = -ENERGY
        DER(2) = 0.0D0
        IF (RAB .LE. TOL) DER(2) = -ENERGY
        DER(3) = 0.0D0
        IF (RAB .LE. TOL) DER(3) = -ENERGY
      END IF
C
      RETURN
      END
C
      SUBROUTINE PREM22
C  COLLINEAR A + BC POTENTIAL FUNCTION
C     RMO-SPLINE TYPE FIT
C     COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PHI(20),BSP(20,3),CSP(22,3)
      COMMON /PT1CM/ R(3),E,D(3)
      COMMON /PT2CM/ NDER, NDUM(21)
      COMMON /PT4CM/ IPRT, IREAD
      CHARACTER*20 TITLE
      DATA DETRAD/0.0174532D0/, CKCAL/627.5095D0/, BOHR/0.529177D0/
      SAVE
C
C  READ IN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)
      WRITE (IPRT, 599)
C     TITLE
      READ (IREAD, 800) TITLE
      A1INF = 0.0D0
      A2INF = 0.0D0
      IAYES = 0
C  READ IN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS
      READ (IREAD, 801) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,A2INF
      IF(A2INF .EQ.0.0D0) A2INF = D1INF-D2INF
C  READ IN THE PIVOT POINT, POTENTIAL AT THE PIVOT, FLAG, AND
C     ASYMPTOTIC EXPONENTIAL PARAMETERS
      READ (IREAD, 802) R1S,R2S,VSP,IAYES,ALPH1,ALPH2
      WRITE (IPRT, 600) TITLE
      WRITE (IPRT, 601) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,A2INF
      WRITE (IPRT, 602) R1S,R2S,VSP,ALPH1,ALPH2
      IF(IAYES.NE.0) WRITE (IPRT, 603)
      IF(IAYES.EQ.0) WRITE (IPRT, 604)
C  READ IN NUMBER OF ANGLES IN ROTATED MORSE OSCILLATOR FIT
      READ (IREAD, 803) NPHI
C  READ IN ANGLES
      READ (IREAD, 804) (PHI(I),I = 1,NPHI)
C  READ IN FUNCTION VALUES AT THESE ANGLE
      DO 10 J = 1,3
         READ (IREAD, 804) (BSP(I,J),I = 1,NPHI)
   10 CONTINUE
      NPHIP2 = NPHI + 2
C  READ IN SPLINE FIT COEFFICIENTS
      DO 20 J = 1,3
         READ (IREAD, 804)(CSP(I,J),I = 1,NPHIP2)
   20 CONTINUE
      WRITE (IPRT, 605) NPHI
      WRITE (IPRT, 610) 
      DO 30 I = 1,NPHI
         WRITE (IPRT, 606) I,PHI(I),(BSP(I,J),CSP(I,J),J=1,3)
   30 CONTINUE
      WRITE (IPRT, 611) 
      DO 40 I=1,2
         K = I + NPHI
         WRITE (IPRT, 607) K,(CSP(K,J),J=1,3)
   40 CONTINUE
C  CONVERT TO ATOMIC UNITS
      D1INF = D1INF/CKCAL
      A1INF = A1INF/CKCAL
      R1INF = R1INF/BOHR
      B1INF = B1INF*BOHR
      D2INF = D2INF/CKCAL
      A2INF = A2INF/CKCAL
      R2INF = R2INF/BOHR
      B2INF = B2INF*BOHR
      R1S = R1S/BOHR
      R2S = R2S/BOHR
      VSP = VSP/CKCAL
      ALPH1 = ALPH1*BOHR
      ALPH2 = ALPH2*BOHR
      DO 50 I = 1,NPHI
         BSP(I,1) = BSP(I,1)/CKCAL
         BSP(I,2) = BSP(I,2)/BOHR
         BSP(I,3) = BSP(I,3)*BOHR
   50 CONTINUE
      NNN = NPHI + 2
      DO 60 I = 1,NNN
         CSP(I,1) = CSP(I,1)/CKCAL
         CSP(I,2) = CSP(I,2)/BOHR
         CSP(I,3) = CSP(I,3)*BOHR
   60 CONTINUE
C  CONVERT FROM DEGREES TO RADIANS
      T = 1.0D0/DETRAD
      DO 70 J = 1,3
        CSP(2,J) = CSP(2,J)*T
   70 CONTINUE
      T = T**3
      DO 80 II = 1,NPHI
         PHI(II) = PHI(II)*DETRAD
         I = II + 2
         DO 80 J = 1,3
            CSP(I,J) = CSP(I,J)*T
   80 CONTINUE
C
C  SET UP CONSTANTS
C     REACTANT ASYMPTOTE
      DIS = BSP(1,1)
      REQ = BSP(1,2)
      BET = BSP(1,3)
      T1 = 1.0D0 - EXP(-BET*REQ)
      AER = VSP - DIS*T1*T1
      REQ = R2S - REQ
      D1DIF = DIS - D1INF
      R1DIF = REQ - R1INF
      B1DIF = BET - B1INF
      A1DIF = AER - A1INF
C     PRODUCT ASYMPTOTE
      DIS = BSP(NPHI,1)
      REQ = BSP(NPHI,2)
      BET = BSP(NPHI,3)
      T1 = 1.0D0 - EXP(-BET*REQ)
      AER = VSP - DIS*T1*T1
      REQ = R1S - REQ
      D2DIF = DIS - D2INF
      R2DIF = REQ - R2INF
      B2DIF = BET - B2INF
      A2DIF = AER - A2INF
C  DEFINE ASYMPTOTIC ENERGY AND ZERO OF ENERGY
      DASY = D1INF
C
599   FORMAT(/,2X,T5,'Parameters for the collinear part of the ',
     *               'potential')
600   FORMAT(2X,T5,'Title in potential data file: ',A20)
601   FORMAT(2X,T5,'D1INF, R1INF, B1INF, D2INF, R2INF, B2INF, ',
     *               'A1INF, A2INF (kcal/mol):',
     *       /,(2X,T15,F10.5,T30,F10.5,T45,F10.5,T60,F10.5))
602   FORMAT(2X,T5,'R1S, R2S, VSP, ALPH1, ALPH2 (ANG/KCAL):',
     *       /,(2X,T15,F10.5,T30,F10.5,T45,F10.5))
603   FORMAT(2X,T5,'Note: the pivot point is fixed at VSP')
604   FORMAT(2X,T5,'Note: the pivot point is not fixed; VSP is ',
     *               'the 3-body break-up energy')
605   FORMAT(2X,T5,'NPHI = ',I5)
610   FORMAT(2X,T13,'I',T25,'PHI',T39,'BSP',T56,'CSP')
606   FORMAT(2X,T11,I5,T20,F10.5,(T35,F10.5,T50,F15.5))
611   FORMAT(2X,T11,'K',T20,'CSP')
607   FORMAT(2X,T11,I5,(T20,E15.5,T40,E15.5))
  800 FORMAT (A20)
  801 FORMAT (8F9.5)
  802 FORMAT (3F10.5,I5,2F10.5)
  803 FORMAT (16I5)
  804 FORMAT (4G20.10)
      RETURN
C
      ENTRY POTM22
C
      R1 = R(1)
      R2 = R(2)
      R3 = R(3)
      IF (R1 .GT. R1S) THEN
         IF (R2 .GT. R2S) THEN
C   3-ATOM BREAK UP REGION
            E = DASY
            IF (NDER .EQ. 1) THEN
                D(1) = 0.0D0
                D(2) = 0.0D0
            ENDIF 
         ELSE
C   REACTANT ASYMPTOTIC REGION
            DR1 = R1 - R1S
            IF (DR1*ALPH1 .LT. -70.0D0) WRITE (IPRT, 6601) R1, R2, R3,
     *         DR1, ALPH1
            T1 = EXP(-DR1*ALPH1)
            T2 = -ALPH1*T1
C.....INTERPOLATE FOR MORSE PARAMETERS
            DIS = D1INF + D1DIF*T1
            DDIS = D1DIF*T2
            REQ = R1INF + R1DIF*T1
            DREQ = R1DIF*T2
            BET = B1INF + B1DIF*T1
            DBET = B1DIF*T2
            IF(IAYES.EQ.0) THEN
               AER = VSP - DIS
               DAER = -DDIS
            ELSE
               AER = A1INF + A1DIF*T1
               DAER = A1DIF*T2
            END IF
            DELR = R2 - REQ
            IF (DELR*BET .LT. -70.0D0) WRITE (IPRT, 6602) R1, R2, R3, 
     *                                        DELR, BET
            T1 = EXP(-BET*DELR)
            T2 = 1.0D0 - T1
            E = AER + DIS*T2*T2
            IF (NDER .EQ. 1) THEN
                T3 = 2.0D0*DIS*T1
                D(1) = DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))
                D(2) = T3*BET*T2
            ENDIF
         END IF
      ELSE IF (R2 .GT. R2S) THEN
C   PRODUCT ASYMPTOTIC REGION
         DR2 = R2 - R2S
         IF (DR2*ALPH2 .LT. -70.D0) WRITE (IPRT, 6603) R1, R2, R3, DR2,
     *      ALPH2
C.....INTERPOLATE FOR MORSE PARAMETERS
         T1 = EXP(-DR2*ALPH2)
         T2 = -ALPH2*T1
         DIS = D2INF + D2DIF*T1
         DDIS = D2DIF*T2
         REQ = R2INF + R2DIF*T1
         DREQ = R2DIF*T2
         BET = B2INF + B2DIF*T1
         DBET = B2DIF*T2
         IF(IAYES.EQ.0) THEN
            AER = VSP - DIS
            DAER = -DDIS
         ELSE
            AER = A2INF + A2DIF*T1
            DAER = A2DIF*T2
         END IF
         DELR = R1 - REQ
         IF (DELR*BET .LT. -70.0D0) WRITE (IPRT, 6604) R1, R2, R3, 
     *                                     DELR, BET
         T1 = EXP(-BET*DELR)
         T2 = 1.0D0 - T1
         E = AER + DIS*T2*T2
         IF (NDER .EQ. 1) THEN
             T3 = 2.0D0*DIS*T1
             D(2) = DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))
             D(1) = T3*BET*T2
         ENDIF
      ELSE
C   SPLINE INTERPOLATION REGION
         DR1 = R1S - R1
         DR2 = R2S - R2
         T1 = DR1*DR1 + DR2*DR2
         RR = SQRT(T1)
         DRRD1 = -DR1/RR
         DRRD2 = -DR2/RR
         ANGLE = ATAN(DR1/DR2)
         DPHD1 = -DR2/T1
         DPHD2 = DR1/T1
         DIS = POTSPL(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
         DDIS = POTSPL(4,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
         REQ = POTSPL(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
         DREQ = POTSPL(4,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
         BET =  POTSPL(3,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)
         DBET = POTSPL(4,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)
         IF(IAYES.EQ.0) THEN
            AER = VSP - DIS
            DAER = -DDIS
         ELSE
            T1 = EXP(-BET*REQ)
            T2 = 1.0D0 - T1
            AER = VSP - DIS*T2*T2
            DAER = -T2*(DDIS*T2 + 2.0D0*DIS*T1*(BET*DREQ + REQ*DBET))
         END IF
         IF (BET*(RR-REQ) .GT.  70) WRITE (IPRT, 6600) R1, R2, R3, RR,
     *      ANGLE, DIS, REQ, BET
         DELR = RR - REQ
         T1 = EXP(BET*DELR)
         T2 = 1.0D0 - T1
         E = AER + DIS*T2*T2
         IF (NDER .EQ. 1) THEN
             T3 = 2.0D0*DIS*T1
             DVDP = DAER + T2*(DDIS*T2 - T3*(DBET*DELR - BET*DREQ))
             T3 = T3*T2*BET
             D(1) =  DVDP*DPHD1 - T3*DRRD1
             D(2) =  DVDP*DPHD2 - T3*DRRD2
         ENDIF
      END IF
      E = E 
      IF (NDER .NE. 1) RETURN
      D(1) = D(1)
      D(2) = D(2)
      D(3) = 0.0D0
C
6600  FORMAT(/,2X,T5,'R1, R2, R3, RR, ANGLE, DIS, REQ, BET:',
     *       /,(2X,T5,1PE20.10,T30,1PE20.10,T55,1PE20.10))
6601  FORMAT(/,2X,T5,'R1, R2, R3, DELR, ALPH1:',/,2X,T10,5(1PE13.5,1X))
6602  FORMAT(/,2X,T5,'R1, R2, R3, DR2, BET:',/,2X,T10,5(1PE13.5,1X))
6603  FORMAT(/,2X,T5,'R1, R2, R3, DELR, ALPH2:',/,2X,T10,5(1PE13.5,1X))
6604  FORMAT(/,2X,T5,'R1, R2, R3, DR2, BET:',/,2X,T10,5(1PE13.5,1X))
C
      RETURN
      END
C
      FUNCTION POTSPL(ISW,NN,X,Y,C,D,XPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),C(*),D(*)
C
C  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
C  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
C  OR ITS DERIVATIVE AT A GIVEN POINT
C
C  CALLING SEQUENCE .......
C     ISW ... CONTROL OPTION
C         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS
C                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS
C                ARE STORED IN THE ARRAY C.
C         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS
C                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY
C                THE PARAMETER 'XPT'.
C         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATE
C         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WIT
C                IS RETURNED.
C
C     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO
C            BE FITTED
C
C     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN
C             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST
C             NN.
C
C     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.
C           MUST BE DIMENSIONED AT LEAST NN+2 .
C
C     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .
C
C     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS
C              SET TO 2).  THE VALUE OF POTSPL IS SET TO THE
C              INTERPOLATED VALUE.
C
C   Output from this subprogram is written to FORTRAN unit IUNIT
C
C *****  USER NOTES  *****
C
C     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......
C
C       A.  CALL POTSPL WITH THE KNOTS.  THIS SETS UP THE
C           COEFFICIENT ARRAY C.
C           EG.  DUMY=POTSPL(1,NN,X,Y,C,D,XPT)
C
C       B.  CALL POTSPL WITH THE ARRAY C WHICH WAS DEFINED BY THE
C           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE
C           POINT 'XPT' .
C           EG.   VALUE=POTSPL(2,NN,X,Y,C,D,XPT)
C
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.
C
      N = NN
      NP1 = N+1
      NP2 = N+2
      Z = XPT
      IF (ISW .EQ. 1) THEN
         C(1) = Y(1)
         D(1) = 1.0D0
         C(NP1) = 0.0D0
         D(NP1) = 0.0D0
         C(NP2) = 0.0D0
         D(NP2) = 0.0D0
         DO 41 I = 2,N
            C(I) = Y(I)-Y(1)
41          D(I) = X(I)-X(1)
         DO 410 I = 3,NP2
            IF(D(I-1).EQ.0.0D0) THEN
               WRITE(IUNIT,1001)
               STOP 'POTSPL 1'
            END IF
            PIVOT = 1.0D0/D(I-1)
            IF(I.LT.NP2) THEN
               SUPD = X(I-1)-X(I-2)
               IF(SUPD.LT.0.0D0) THEN
                  WRITE(IUNIT,1000)
                  STOP 'POTSPL 2'
               END IF
               SUPD = SUPD*SUPD*SUPD
            ELSE
               SUPD = 1.0D0
            END IF
            DFACT = SUPD*PIVOT
            CFACT = C(I-1)*PIVOT
            IF(I.LE.N) THEN
               DO 47 J = I,N
                  V = X(J)-X(I-2)
                  C(J) = C(J)-D(J)*CFACT
47                D(J) = V*V*V-D(J)*DFACT
            END IF
            IF(I.LT.NP2) THEN
               C(NP1) = C(NP1)-D(NP1)*CFACT
               D(NP1) = 1.0D0-D(NP1)*DFACT
            END IF
            C(NP2) = C(NP2)-D(NP2)*CFACT
            D(NP2) = X(I-2)-D(NP2)*DFACT
410      CONTINUE
         DO 411 I = 1,N
            J = NP2-I
            IF(J.EQ.NP1) THEN
               V = 1.0D0
            ELSE
               V = X(J)-X(J-1)
               V = V*V*V
            END IF
            IF(D(J+1).EQ.0.0D0) THEN
               WRITE(IUNIT,1001)
               STOP 'POTSPL 3'
            END IF
            C(J+1) = C(J+1)/D(J+1)
            C(J) = C(J)-C(J+1)*V
411      CONTINUE
         IF(D(2).EQ.0.0D0) THEN
            WRITE(IUNIT,1001)
            STOP 'POTSPL 4'
         END IF
         C(2) = C(2)/D(2)
      ELSE IF (ISW .EQ. 2) THEN
         POTSPL = C(1) + C(2)*(Z-X(1))
         DO 51 I = 1,N
            V = Z-X(I)
            IF(V.LE.0) RETURN
            POTSPL = POTSPL + C(I+2)*V**3
51       CONTINUE
      ELSE IF (ISW .EQ. 3) THEN
         POTSPL = C(1) + C(2)*(Z-X(1))
         DERIV = C(2)/3.D0
         DO 53 I = 1,N
            V = Z-X(I)
            IF(V.LE.0) RETURN
            V2 = V*V
            POTSPL = POTSPL + C(I+2)*V2*V
            DERIV = DERIV + C(I+2)*V2
   53    CONTINUE
      ELSE IF (ISW .EQ. 4) THEN
         POTSPL = 3.D0*DERIV
      END IF
      RETURN
1000  FORMAT(1X,5X,'***** ERROR IN POTSPL ... UNORDERED X-VALUES
     1*****')
1001  FORMAT(1X,5X,' ***** ERROR IN POTSPL ... DIVIDE FAULT *****')
      END
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
         COMMON /AMBPCM/ D(3), BETA, REQ(3), A, GAMP(3), RBB
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT, IREAD /6, 4/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters for the collinear part of the 
C   potential; the energy parameters are in kcal/mol, and the lengths 
C   are in Angstroms.  These are the parameters for the A'' surface.
         DATA D / 106.56D0, 109.472D0, 106.56D0/
         DATA BETA / 2.07942D0/
         DATA REQ / 0.96966D0, 0.74144D0, 0.96966D0/
         DATA A / 0.26D0/
         DATA GAMP / 0.210820D0, 5.42027D0, -1.44235D0/
         DATA RBB / 0.74144D0/
C
         END
C*****
