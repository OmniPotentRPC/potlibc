C
      SUBROUTINE PREPOT
C
C   System:          ClH2
C   Functional form: Empirical triatomic potential surface defined
C                    over orthogonal bond order coordinates
C   Common name:     AL/AB
C   References:      B. C. Garrett, D. G. Truhlar, and A. W. Magnuson
C                    J. Chem. Phys. 76, 2321 (1982)
C                    N. Agmon and R. D. Levine
C                    J. Chem. Phys. 71, 3034 (1979)
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
C      coordinates - bohrs
C      derivatives - hartrees/bohr
C
C   Surface(s): 
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
C      Internal, Definition: R(1) = R(Cl-first H)
C                            R(2) = R(first H-second H)
C                            R(3) = R(Cl-second H)
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
         DOUBLE PRECISION LAMDA, MORSE, M, NAB, NAB1, NBC, NBC1, 
     *                    NUM, JUNK
         PARAMETER (BOHR = .52917706D0)
         PARAMETER (HRTREE = 627.5095D0)
         COMMON /PT1CM/ RAB,RBC,RAC,ENERGY,DER(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PPARCM/ DE(3), REQ(3), BETA(3), LAMDA, A, GAMMA, RBB
         DIMENSION OMEGA(2),JUNK(4),DVSTAR(2),MORSE(2),COEF(2),DBCOOR(2)

         SAVE RHH
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Echo the potential parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE (IPTPRT,600) DE, REQ, BETA
      WRITE (IPTPRT,610) A, GAMMA, RBB, LAMDA
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CONVERT TO ATOMIC UNITS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 50 IT = 1,3
      REQ(IT) = REQ(IT)/BOHR
      BETA(IT) = BETA(IT) * BOHR
      DE(IT) = DE(IT)/HRTREE
   50 CONTINUE
      A = A/BOHR
      LAMDA = LAMDA/HRTREE
      RHH = 1.4007977D0
      RBB = RBB /BOHR
C
C   Set the values for the potenetial in the three asymptotic valleys
         EASYAB = DE(1)
         EASYBC = DE(2)
         EASYAC = DE(3)
C
CCCCCCCCCCCC
      RETURN
CCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        MOVE ONTO THE MAIN PART OF THE ROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ENTRY POT
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C        FIRST CALCULATE SOME USEFUL NUMBERS
C
      NAB = EXP((REQ(1) - RAB)/A)
C
      IF (NAB.LT.1.D-15) NAB = 1.D-15
C
      NBC = EXP((REQ(2) - RBC)/A)
C
      IF (NBC.LT.1.D-15) NBC = 1.D-15
C
C        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS
C        AB AND BC RESPECTIVELY. (EQUATION 2.2)
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
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE
C       SEE EQUATION 2.7A, 2.7B
C
   95 M = NAB + NBC
C
C        M IS THE SUM OF THE BOND ORDERS (EQUATION 2.3)
C
      OMEGA(1) = NAB / M
      OMEGA(2) = NBC / M
C
C       THE OMEGAS ARE USED IN EQUATION 3.15A
C
      SIGMA = -A * LOG(M)
C
C        SIGMA IS THE MYSTERY NUMBER REFERENCED IN EQUATION 3.15B
C        SE EQUATION 2.13
C
      DO 100 IT = 1,2
      JUNK(IT) = EXP(BETA(IT)*(-SIGMA))
      MORSE(IT) = DE(IT)*(JUNK(IT)*JUNK(IT)-2.D0*JUNK(IT))
  100 CONTINUE
C
C        MORSE FUNCTION--EQUATION 4-16, GAS PHASE REACTION
C        RATE THEORY, H. S. JOHNSTON.
C
C        ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS
C
C          FIRST TRAP OUT ANY ZERO ARGUEMENTS
C
      IF (NAB1*NBC1.GT.0.D0) GO TO 103
           ROB = 1.D0
      GO TO 107
C
  103 STUFF = A * LOG(NAB1*NBC1)
      ROB = (RBB - STUFF)/(RHH - STUFF)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CALCULATE DE(OMEGA)
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        THIS IS FORMULA 3.17
C
  107 DEW = 0.D0
      DO 110 IT = 1,2
      DEW = OMEGA(IT)*DE(IT) + DEW
  110 CONTINUE
      DEW = DEW - LAMDA * BIGM(OMEGA(1))
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        CALCULATE V*(SIGMA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      VSTAR=0.D0
      NUM = 0.D0
      DENOM = 0.D0
      DO 210 IT = 1,2
      NUM = OMEGA(IT) * MORSE(IT) + NUM
      DENOM = OMEGA(IT)*DE(IT) + DENOM
  210 CONTINUE
      VSTAR = NUM/DENOM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        CALCULATE BENDING CORRECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        CALCULATE V(RAC)
C
      JUNK(3) = EXP(-BETA(3)*(RAC-REQ(3))/ROB)
      VAC = GAMMA * DE(3) * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))
C
C        CALCULATE V(RAB+RBC)
C
      JUNK(4) = EXP(BETA(3)*(REQ(3)-RBC-RAB)/ROB)
      VABBC= GAMMA * DE(3) * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))
C
C        HERE IS THE BENDING CORRECTION
C
      BCOOR = VAC - VABBC
C
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
C          OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
C
C        CALCULATE DER(3) (THE DERIVATIVE WITH RESPECT TO RAC)
C
      IF (NDER .EQ. 1) THEN
          DER(3) = -GAMMA*DE(3)*BETA(3)*JUNK(3)*(1.D0+JUNK(3))/ROB
          DBCOOR(1) = GAMMA * DE(3) * BETA(3) * JUNK(4) * 
     *                (1.D0+ JUNK(4))/ROB
          DBCOOR(2)= DBCOOR(1)
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        NOW CALCULATE THE ENERGY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ENERGY = DEW * VSTAR + BCOOR + DE(2)
C
C           DE(2) ADJUSTS THE POTENTIAL SO THAT IT IS ZERO
C           IN THE ASYMPTOTIC REACTANTS REGION.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        NOW CALCULATE THE DERIVATIVES, IF NDER = 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF (NDER .NE. 1) GO TO 920
C
C        NEXT CALCULATE DER(1) AND DER(2)
C
      DO 300 IT = 1,2
C
C        FIRST CALCULATE DVSTAR(1) AND DVSTAR(2) SEE EQUATION B.5
C
      V1=OMEGA(IT)/(DE(1)*OMEGA(1) + DE(2) * OMEGA(2))
C
      IF (IT.EQ.1) GO TO 280
      V2 = OMEGA(1)*(MORSE(1) - MORSE(2) - VSTAR * (DE(1)-DE(2)))/A
      GO TO 285
  280 V2 = OMEGA(2)*(MORSE(2) - MORSE(1) - VSTAR*(DE(2)-DE(1)))/A
  285 CONTINUE
C
      V3 = -OMEGA(IT)*2.D0*DE(IT)*BETA(IT)*(JUNK(IT)*JUNK(IT)-JUNK(IT))
C
      IF (IT.EQ.1) GO TO 290
      V4 = -OMEGA(1)*2.D0*DE(1)*BETA(1)*(JUNK(1)*JUNK(1) - JUNK(1))
      GO TO 295
  290 V4 = -OMEGA(2)*2*DE(2)*BETA(2)*(JUNK(2)*JUNK(2) - JUNK(2))
  295 CONTINUE
C
      DVSTAR(IT) = V1*(V2+V3+V4)
C
C       NOW CALCULATE THE COEFFICIENT OF VSTAR IN EQUATION B.1.
C
      COEF(IT) = (NAB*NBC/(M*M*A))
     1    * ((DE(1)-DE(2)+LAMDA * LOG(NAB/NBC)) * ((-1.D0)**IT))
C
C       FINALLY PUT IT ALL TOGETHER TO CALCULATE THE DERIVITIVE
C       SEE EQUATION B.1
C
      DER(IT) = COEF(IT) * VSTAR + DEW * DVSTAR(IT) + DBCOOR(IT)
C
  300 CONTINUE
  310 CONTINUE
920   CONTINUE
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClH2 AL/AB potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /,1X,T5,'Bond', T46, 'Cl-H', T58, 'H-H', T69, 'H-Cl',
     *        /,1X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
610   FORMAT(/,1X,T5,'Pauling parameter (Angstroms):',T44,F10.5,
     *       /,1X,T5,'Gamma for the bending correction:',T44,F10.5,
     *       /,1X,T5,'RBB (Angstroms):',T44,F10.5,
     *       /,1X,T5,'LAMDA (kcal/mol):',T44,F10.5,//,1X,'*****')
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C 
      RETURN
      END
C
      FUNCTION BIGM(DUMMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C        THIS FUNCTION RETURNS THE M THAT IS REFERENCED
C        IN EQUATION 3.12.
C
      IF (DUMMY.LT.1.D0) GO TO 10
      BIGM = 0.D0
      RETURN
C
   10 IF (DUMMY.GT.0.D0) GO TO 20
      BIGM = 0.D0
      RETURN
C
   20 BIGM = -DUMMY*LOG(DUMMY) - (1.D0-DUMMY)*LOG(1.D0-DUMMY)
C
      RETURN
      END
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         DOUBLE PRECISION LAMDA
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PPARCM/ DE(3), REQ(3), BETA(3), LAMDA, A, GAMMA, RBB
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER, NSURF /6, 1, 0/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters
C   DE in kcal/mol, REQ in Angstroms, and BETA in reciprocal Angstroms
         DATA DE/ 106.447D0, 109.458D0, 106.447D0/
         DATA REQ/ 1.2732D0, 0.74127D0, 1.2732D0/
         DATA BETA/ 1.8674D0, 1.9413D0, 1.8674D0/
C
C   LAMDA in kcal/mol, Pauling parameter, A, in Angstroms, 
C   GAMMA for the bend correction, unitless, and RBB in Angstroms.
         DATA LAMDA / 8.2444D0/
         DATA A / 0.272D0/
         DATA GAMMA / 0.5D0/
         DATA RBB / 0.74127D0/
C
         END
C*****
