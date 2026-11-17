C
      SUBROUTINE PREPOT
C
C   System:          IH2
C   Functional form: 
C   Common name:     RMC
C   References:      J. W. Duff and D. G. Truhlar, J. Chem. Phys., Vol. 52, 2477 (1975).
C
C   Note: This is the three-dimensional version of the IH2 
C         RMC potential energy surface; the potential energy is 
C         defined as E = Ecollinear + Ebend.
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PPARMS.
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
         COMMON /RMCRCM/ RR(2), S(2), D(2), ALPHA(2), RE(2), ALFD, CLBH
         COMMON /RMCICM/ IRXN, L(2)
         COMMON /BENDCM/ DE, BETA, REQ(3), A, GAMMA, RBB
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
         PARAMETER (CKCAL = 627.5095D0)
         PARAMETER (CANGS = 0.52917706D0)
         SAVE
C
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1,JUNK
      DIMENSION DBCOOR(2),JUNK(4)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      COLLINEAR POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Initialize the parameters for the collinear part of the potential
      CALL PREPT2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      BENDING POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Echo the parameters for the bend potential to the file linked to 
C   FORTRAN unit IPRT
         WRITE (IPRT, 650) REQ, DE, BETA, A, GAMMA, RBB
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CONVERT TO ATOMIC UNITS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 50 IT = 1,3
      REQ(IT) = REQ(IT)/CANGS 
   50 CONTINUE
      BETA = BETA * CANGS 
      DE = DE/CKCAL 
      A = A/CANGS
      RHH = 1.4007977D0
      RBB = RBB /CANGS
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
      IF (NDER .EQ. 1) THEN
          DAC       = -GAMMA*DE*BETA*JUNK(3)*(1.D0+JUNK(3))/ROB
          DBCOOR(1) = GAMMA * DE * BETA * JUNK(4) * 
     *                (1.D0 + JUNK(4))/ROB
          DBCOOR(2) = DBCOOR(1)
      ENDIF
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
C   Reset RAC the its correct value
      RAC = RACTMP
C
C   Add the energy terms together, i.e., E = Ecollinear + Ebend
      ENERGY = ENERGY + BCOOR
C
      IF (NDER .NE. 1) RETURN
C
C   If NDER = 1, then compute the derivatives of the energy with respect to
C   the internal coordinates, i.e., dE/dRi = dEcollinear/dRi + dEbend/dRi
      DER(1) = DER(1) + DER(3)
      DER(2) = DER(2) + DER(3)
      DER(3) = DAC
C
      DO 300 IT = 1,2
           DER(IT) = DBCOOR(IT) + DER(IT)
  300 CONTINUE
C
650   FORMAT (/,2X,T5,'Parameters for the bend potential',
     *        /,2X,T5,'Bond', T47, 'I-H', T58, 'H-H', T69, 'H-I',
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Dissociation energy (kcal/mol):', 
     *        T44, F15.8, 
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F15.8, 
     *        /,2X,T5,'Pauling parameter (Angstroms):', 
     *        T44, F15.8, 
     *        /,2X,T5,'Gamma (unitless):', 
     *        T44, F15.8, 
     *        /,2X,T5,'RBB (Angstroms):', 
     *        T44, F15.8)
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,
     *       /,2X,T5,'This value of NDER is not allowed in this ',
     *               'version of the potential.')
C
      RETURN
      END
C*****
      SUBROUTINE PREPT2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT1CM/ R(3),V,DVDR(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(21)
      COMMON /PT4CM/ IPRT, IDUM(19)
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      COMMON /RMCRCM/ RR(2), S(2), D(2), ALPHA(2), RE(2), ALFD, CLBH
      COMMON /RMCICM/ IRXN, L(2)
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr
      PARAMETER (CKCAL = 627.5095D0)
      PARAMETER (CANGS = 0.52917706D0)
C
C   This subprogram initializes the constants for the collinear part of 
C   the potential.
C   PREPT2 must be called once before any calls to POT2.
C
      LOGICAL LRR1,LRR2
      DIMENSION C(2,2),PRE(2),RRMRE(2)
      DIMENSION A(2),DADP(2),DPDR(2),RRMR(2)
C
      PID4 = ATAN(1.D0)
      PID2 = 2.D0*PID4
      M = IRXN
      N = 1
      IF (IRXN.EQ.1) N = 2
      C(M,1) = CLBH/D(M)
      IF (IRXN.EQ.1) C(M,1) = C(M,1) + (D(M)-D(N))/D(M)
      C(M,2) = 0.0D0
      DMDN = D(M)/D(N)
      C(N,2) = 2.0D0*DMDN*C(M,1)*DBLE(L(N)-L(M)) +
     1   2.0D0*DBLE(L(N))*(1.0D0-DMDN)
      C(N,1) = 1.0D0-DMDN*(1.0D0-C(M,1))+PID4*PID4*C(N,2)
C
C   Echo the potential parameters to the file linked to FORTRAN unit IPRT.
         WRITE (IPRT, 600) D(1), D(2), D(1), RE(1), RE(2), RE(1),
     *                     ALPHA(1), ALPHA(2), ALPHA(1)
         WRITE (IPRT, 610) RR, S, L, ALFD
         WRITE (IPRT, 620) ((C(I,J),I=1,2),J=1,2)
C   CONVERT TO ATOMIC UNITS
      DO   2 I = 1,2
      RR(I) = RR(I)/CANGS
      S(I) = S(I)/CANGS
      RE(I) = RE(I)/CANGS
      ALPHA(I) = ALPHA(I)*CANGS
   2   D(I)=D(I)/CKCAL
      ALFD = ALFD*CANGS
C   CALCULATE USEFUL CONSTANTS
      RR1MS1=RR(1)-S(1)
      RR2MS2=RR(2)-S(2)
      RRS1D2=RR1MS1/RR2MS2
      SR1R2=SQRT(RR1MS1**2+RR2MS2**2)
      DO   5 N=1,2
      RRMRE(N)=RR(N)-RE(N)
    5 PRE(N)=SR1R2-RRMRE(N)
C   Initialize the asymptotic energy terms
         EASYAB = D(1)
         EASYBC = D(2)
         EASYAC = D(1)
      RETURN
C
      ENTRY POT2
C
C   This subprogram computes the collinear part of the potential energy
C
      RRMR(1) = RR(1)-R(1)
      RRMR(2) = RR(2)-R(2)
      LRR1 = RRMR(1).LE.0.0D0
      LRR2 = RRMR(2).LE.0.0D0
      IF (NDER .EQ. 1) THEN
          DO 10 J = 1,3
10              DVDR(J) = 0.0D0
      ENDIF
      V = D(2)
      IF (LRR1.AND.LRR2) RETURN
      IF (LRR1.OR.LRR2) GO TO 50
      PHIP = ATAN (RRS1D2*RRMR(2)/RRMR(1))
      SN2PP = SIN (2.0D0*PHIP)
      CS2PP = COS (2.0D0*PHIP)
      TN2PP = SN2PP/CS2PP
      SN2PP3 = SN2PP**3
      SN2PP4 = SN2PP3*SN2PP
      I = 1
      IF (PHIP.GT.PID4) I = 2
      A(I) = C(I,1) + C(I,2)*(PHIP-PID2)*PHIP
      DADP(I) = C(I,2)*(2.0D0*PHIP-PID2)
      RR12 = SQRT (RRMR(1)**2+RRMR(2)**2)
      LM3 = L(I) - 3
      SN2PPL = SN2PP3*SN2PP**LM3
      ALF = ALPHA(I) +(ALFD-ALPHA(I))*SN2PP4
      GAM = D(I)*(1.0D0-A(I)*SN2PPL)
      XI = PRE(I)*SN2PP4+RRMRE(I)-RR12
      QUAN = 1.0D0 - EXP(-ALF*XI)
      QUAN2 = QUAN*QUAN
      V = GAM*(QUAN2-1.0D0)+D(2)
C
      IF (NDER .NE. 1) RETURN
C   If NDER=1, compute the derivatives of the energy with respect to the 
C   internal coordinates.
      FAC = RRS1D2/(RRMR(1)**2+(RRS1D2*RRMR(2))**2)
      DPDR(1) = FAC*RRMR(2)
      DPDR(2) = -FAC*RRMR(1)
      DO 40 J = 1,2
      DS4DR = 8.0D0*SN2PP4/TN2PP*DPDR(J)
      DSLDR = 2.0D0*DBLE(L(I))*SN2PPL/TN2PP*DPDR(J)
      DALDR = (ALFD-ALPHA(I))*DS4DR
      DGADR = -D(I)*A(I)*DSLDR - D(I)*DADP(I)*SN2PPL*DPDR(J)
      DXIDR = PRE(I)*DS4DR + RRMR(J)/RR12
   40 DVDR(J) = DGADR*(QUAN2-1.0D0)+2.0D0*GAM*QUAN*(1.0D0-QUAN)*
     1          (ALF*DXIDR + XI*DALDR)
      RETURN
   50 J = 1
      IF (LRR1) J = 2
      QUAN = 1.0D0-EXP(-ALPHA(J)*(R(J)-RE(J)))
      V = D(J)*(QUAN**2-1.0D0)+D(2)
C   If NDER =1, compute the derivatives of the energy with respect to the 
C   internal coordinates.
      IF (NDER .EQ. 1) DVDR(J) = 2.0D0*D(J)*ALPHA(J)*QUAN*(1.0D0-QUAN)
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the IH2 ',
     *                '3D potential energy surface RMC',
     *       //,2X,T5,'Parameters for the collinear part of the ',
     *                'potential',
     *        /,2X,T5,'Bond', T47, 'I-H', T58, 'H-H', T69, 'H-I',
     *        /,2X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT (/,2X,T5,'Point of rotation (Angstroms):',
     *        T44, F10.5, T55, F10.5,
     *        /,2X,T5,'Saddle point position (Angstroms):', 
     *        T44, F10.5, T55, F10.5,
     *        /,2X,T5,'Power in Gamma:', 
     *        T44, I4, T55, I4, 
     *        /,2X,T5,'Alpha at saddle point (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5) 
620   FORMAT (2X,T5,'Coefficients in Gamma:', 
     *        (T44, F10.5, T55, F10.5)) 
C
      RETURN
      END
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /RMCRCM/ RR(2), S(2), D(2), ALPHA(2), RE(2), ALFD, CLBH
         COMMON /RMCICM/ IRXN, L(2)
         COMMON /BENDCM/ DE, BETA, REQ(3), A, GAMMA, RBB
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA ALPHA / 1.789D0, 1.974D0/
         DATA D     / 73.682D0, 109.499D0/
         DATA RE    / 1.604D0, 0.742D0/
         DATA S     / 1.6783D0, 1.2801D0/
         DATA RR    / 2.750D0, 1.750D0/
         DATA ALFD, CLBH /1.500D0, 35.879D0/
         DATA L     / 12, 6/
C
C   If IRXN = 1, then the reaction is in the endothermic direction
C   If IRXN = 2, then the reaction is in the exothermic direction
         DATA  IRXN / 1/
C
C   Initialize the paramters for the bend potential
         DATA DE    / 73.682D0/
         DATA BETA  / 1.789D0/
         DATA REQ   / 1.604D0, 0.742D0, 1.604D0/
         DATA A     / 0.26D0/
         DATA GAMMA / 0.50D0/
         DATA RBB   / 0.7412871D0/
C
         END
C*****
C
