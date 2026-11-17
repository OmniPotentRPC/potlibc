C
      SUBROUTINE PREPOT
C
C   System:          H3
C   Common name:     TK
C   Reference:       D. G. Truhlar and A. Kuppermann
C                    J. Chem. Phys. 56, 2232 (1984).  
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements.
C   and in the block data subprogram PPARMS.
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
         COMMON /PT1CM/ RAB,RBC,RAC,ENERGY,DER(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
C
C   Conversion constants
C   BOHR (Angstrom/bohr) converts Angstroms to Bohr
C   HARTRE (kcal/Hartree) converts from Hartree to kcal/mol
          PARAMETER (BOHR = 0.52917706D0)
          PARAMETER (HARTRE = 627.5095D0)
C
         DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK
         DIMENSION DBCOOR(2),JUNK(4),REQ(3)
         DATA DE /109.47D0/
         DATA BETA /1.97354/
         DATA REQ /0.741287D0, 0.741287D0, 0.741287D0/
         DATA A /0.26D0/
         DATA GAMMA /0.5D0/
         DATA RBB /0.741270D0/
         SAVE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SET COLLINEAR POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL PREPT2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Echo the bending potential parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   DE    = dissociation energy in kcal/mol
C   BETA  = Morse beta in reciprocal Angstroms
C   A     = Pauling parameter in Angstroms
C   GAMMA = gamma for the bend correction
C           gamma is dimensionless and ranges from 0.5 to 0.65
C   RBB   = parameter used in the scaling calculation for the 
C           bending correction
C
C
      WRITE (IPRT,600) DE, DE, DE, (REQ(I),I=1,3),BETA, BETA, BETA
      WRITE (IPRT, 610) A, GAMMA, RBB
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CONVERT TO ATOMIC UNITS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 50 IT = 1,3
      REQ(IT) = REQ(IT)/BOHR
   50 CONTINUE
      BETA = BETA * BOHR
      DE = DE/HARTRE
      A = A/BOHR
      RHH = 1.4007977D0
      RBB = RBB /BOHR
C
CCCCCCCCCCCC
      RETURN
CCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        MOVE ONTO THE MAIN PART OF THE ROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCC
      ENTRY POT 
CCCCCCCCCCCCCCCC
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        FIRST CALCULATE SOME USEFUL NUMBERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
C        AB AND BC RESPECTIVELY. (EQUATION 2.2)
C
      NAB = EXP((REQ(1) - RAB)/A)
C
      IF (NAB.LT.1.D-15) NAB = 1.D-15
C
      NBC = EXP((REQ(2) - RBC)/A)
C
      IF (NBC.LT.1.D-15) NBC = 1.D-15
C
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
C       SEE EQUATION 2.7A, 2.7B
C
      C = 1.D0 - NAB/NBC
      IF (ABS(C).GT.1.D-14) GOTO 90
      NAB1=.5D0
      NBC1=.5D0
      GO TO 95
   90 C = C/NAB
      NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)
      NBC1 = 1.D0 - NAB1
C
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
C
C          FIRST TRAP OUT ANY ZERO ARGUEMENTS
C
   95 IF (NAB1*NBC1.GT.0.0D0) GO TO 103
           ROB = 1.D0
      GO TO 107
C
  103 STUFF = A * LOG(NAB1*NBC1)
      ROB = (RBB - STUFF)/(RHH - STUFF)
C
  107 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        CALCULATE BENDING CORRECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        CALCULATE V(RAC)
C
      JUNK(3) = EXP(-BETA*(RAC-REQ(3))/ROB)
      VAC = GAMMA * DE * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))
C
C        CALCULATE V(RAB+RBC)
C
      JUNK(4) = EXP(BETA*(REQ(3)-RBC-RAB)/ROB)
      VABBC= GAMMA * DE * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))
C
C        HERE IS THE BENDING CORRECTION
C
      BCOOR = VAC - VABBC
C
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
C          OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
C
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO RAC)
C
      DAC = -GAMMA*DE*BETA*JUNK(3)*(1.D0+JUNK(3))/ROB
      DBCOOR(1) = GAMMA * DE * BETA * JUNK(4) * (1.D0 + JUNK(4))/ROB
      DBCOOR(2)= DBCOOR(1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        NOW CALCULATE THE COLLINEAR ENERGY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          STORE RAC IN RACTMP THEN SET RAC TO THE COLLINEAR
C          GEOMETRY AND CALL POT
C
      RACTMP = RAC
      RAC = RBC + RAB
C
      CALL POT2 
C
      DER(1) = DER(1) + DER(3)
      DER(2) = DER(2) + DER(3)
      DER(3) = DAC
      RAC = RACTMP
C
      ENERGY = ENERGY + BCOOR
C
      DO 300 IT = 1,2
           DER(IT) = DBCOOR(IT) + DER(IT)
  300 CONTINUE
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',
     *                'potential energy surface TK',
     *       //,2X,T5,'Potential energy surface parameters:',
     *        /,2X,T5,'Bond', T47, 'H-H', T58, 'H-H', T69, 'H-H',
     *        /,2X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT(/,2X,T5,'Pauling parameter (Angstroms):',T44,F10.5,
     *       /,2X,T5,'Gamma for the bending correction:',T44,F10.5,
     *       /,2X,T5,'RBB (Angstroms):',T44,F10.5)
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not supported in this ',
     *              'version of the potential.')
C
      RETURN
      END
C
      SUBROUTINE PREPT2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RABL,RBCL
C         COORDINATES,POTENTIAL ENERGY, AND DERIVATIVES ARE
C         PASSED THROUGH COMMON /POTCM/ RAB,RBC,RAC,V,DAB,DBC,DAC
      COMMON /PT1CM/ RAB,RBC,R3,V,D1,D2,D3
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      PARAMETER (HARTRE = 627.5095D0)
C        DATA HAS ENERGY IN UNITS OF KCAL/MOLE  IN DATA STATEMENT
C        THIS IN CONVERTED TO HARTREES.  LENGTH IS IN BOHR.
      DATA ALFD,A,SN,D,ALFN,XN/0.70152D0,0.08939D0,1.765D0,109.47D0,
     1                         1.04435D0,1.40083D0/
      DATA L,RNMSN/3,1.85630D0/
      SAVE
C  THIS SUBPROGRAM SHOULD BE CALLED ONCE BEFORE THE FIRST CALL
C TO POT
      RN=RNMSN+SN
      RNMXN=RN-XN
      ALFDIF=ALFD-ALFN
      R2T=1.4142136D0*RNMSN
      D = D /HARTRE
C
C    Set the values of the classical energy in the three asymptotic valleys
             EASYAB = D
             EASYBC = EASYAB
             EASYAC = EASYAB
C
      RETURN
C
      ENTRY POT2
C THIS VERSION WORKS ONLY FOR L .LT. 4
      RABL=RAB.GT.RN
      RBCL=RBC.GT.RN
      IF (RABL.OR.RBCL) GO TO 12
      RNMX=RN-RAB
      RNMY = RN - RBC
      THETA= ATAN(RNMY/RNMX)
      ST = SIN(THETA)
      CT = COS(THETA)
      STT = 2.D0*ST*CT
      CTT = CT*CT - ST*ST
      LM4=L-4
    3 STTL=STT**L
      STT4=STTL*STT**(-LM4)
    6 GAM=D*(1.D0-A*STTL)
      ALF=ALFN+ALFDIF*STT4
      BRACK=R2T-RNMXN
      RCRC=2.D0*RN*RN+RAB*RAB+RBC*RBC-2.D0*RN*(RAB+RBC)
      SQUIG= BRACK*STT4+RNMXN-SQRT(RCRC)
C      DERIVATIVE INFORMATION
      T1 = STTL*CTT*3.D0*A*D
       DGAMX = -T1/RNMX
       DGAMY = T1/RNMY
       T2 = STT4*CTT*4.D0
       T1 = T2*ALFDIF
       DALFX = T1/RNMX
       DALFY = -T1/RNMY
       T1 = T2*BRACK
      DSQIGX = CT + T1/RNMX
       DSQIGY = ST - T1/RNMY
      GO TO 13
   12 GAM=D
      ALF=ALFN
      IF(RABL) SQUIG = RBC-XN
      IF(RBCL) SQUIG = RAB - XN
      IF(RABL.AND.RBCL) SQUIG = SQUIG + RBC - RN
   13 CONTINUE
      T1 = EXP(-ALF*SQUIG)
      QUAN = 1.D0-T1
      T2 = QUAN*QUAN - 1.D0
      V = GAM*T2 + D
      IF(RABL.OR.RBCL) GO TO 20
      D1 = DGAMX*T2
       D2 = DGAMY * T2
      T2 = T1*QUAN*2.D0*GAM
      D1 = D1 + T2*(ALF*DSQIGX + SQUIG*DALFX)
      D2 = D2 + T2*(ALF*DSQIGY + SQUIG*DALFY)
      D3 = 0.D0
      RETURN
   20 T2 = 2.D0*ALF*GAM*T1*QUAN
      D1 = 0.D0
      D2 = 0.D0
      IF(RABL) D2 = T2
      IF(RBCL) D1 = T2
      D3 = 0.D0
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NDER /1/
         DATA NDUM /21*0/
C
         END
C*****
