C
      SUBROUTINE PREPOT
C
C   System:          OH2
C   Functional form: Anti-Morse Bends plus Bowman's collinear spline
C   Common name:     OH2DIM
C   Data file: potoh2dim.dat
C   Reference:       A. Whitlock, J. T. Muckerman, E. R. Fisher, \"Theoretical Invertigations of 
C                    the Energetics and Dynamics of the Reaction of O(3P, 1D) + H2 and C(1D) + H2.\" 
C                    Report from Research Institute for Engineering Sciences, Wayne State University, 1976;
C                    A. F. Wagner, G. C. Schatz, J. M. Bowman, J. Chem. Phys. 74, 4960 (1981); 
C                    G. C. Schatz, A. F. Wagner, S. P. Walch, J. M. Bowman, J. Chem. Phys. 74, 4984 (1981); 
C                    J. M. Bowman and K. T. Lee in \"Potential Energy Surfaces and Dynamics Calculations,\" 
C                    D. G. Truhlar, ed. (Plenum, 1981), p. 359. 
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are read in from a the file "potoh2dim.dat" 
C   which is linked to FORTRAN unit IREAD.
C   The control flags for the potential, such as IREAD, are
C   initialized in the block data subprogram PPARMS.
C
C   This potential is written such that:
C                  R(1) = R(O-H)
C                  R(2) = R(H-H)
C                  R(3) = R(H-O)
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          GENERAL INFORMATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   POTA.FOR HAS ANTI-MORSE BEND PLUS BOWMAN'S COLLINEAR SPLINE ROUTINE.
C   ANTI-MORSE BENDING POTENTIAL IS ADDED TO COLLINEAR POTENTIAL.
C   THIS SUBROUTINE CAN BE USED TO ADD A BENDING CORRECTION TO 
C   ANY COLLINEAR POTENTIAL SUBROUTINE.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C          DEFINE VARIABLES
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ RAB,RBC,RAC,ENERGY,DER(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         DOUBLE PRECISION NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR
         DIMENSION DBCOOR(2),JUNK(4),REQ(3)
C Conversion constant from Angstroms to Bohr (Angstrom/Bohr)
         PARAMETER (BOHR = .52917706D0)
C Conversion constant from Hartree to kcal (kcal/Hartree)
         PARAMETER (HRTREE = 627.5095D0)
      save
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      READ COLLINEAR POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      OPEN (UNIT=IREAD, FILE='potoh2dim.dat', FORM='FORMATTED', 
     *      STATUS='OLD', ERR=100)
      WRITE(IPRT, 600)
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',
     *                'potential energy surface DIM',
     *        /,2X,T5,'Parameters for the potential energy surface',/)
C
      CALL PRPOT2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      READ BENDING POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE(IPRT,3)
3     FORMAT(/,2X,T5,'Parameters for the bending potential')
C
C      DISSOCIATION ENERGY IN KCAL/MOLE
C
      READ(IREAD,*) DE
C
C       READ IN MORSE BETAS--RECIPROCAL ANGSTROMS
C
      READ(IREAD,*) BETA
C
C        READ IN THE EQUILIBRIUM DISTANCES IN ANGSTROMS
C
      READ(IREAD,*) (REQ(IT), IT = 1,3)
C
C        READ IN THE PAULING PARAMETER IN ANGSTROMS
C
      READ(IREAD,*) A
C
C       READ IN GAMMA FOR THE BEND CORRECTION
C       GAMMA IS DIMENSIONLESS AND RANGES FROM .5 TO .65
C
      READ(IREAD,*) GAMMA
C
C       READ IN RBB IN ANGSTROMS--THIS NUMBER IS USED
C       IN A SCALING CALCULATION FOR THE BENDING CORRECTION
C
      READ(IREAD,*) RBB
C
C   Close the input data file
      CLOSE (UNIT = IREAD)
C
C   Echo the potential parameters to the file linked to FORTRAN unit IPRT
      WRITE (IPRT,25) (REQ(IT), IT = 1,3)
      WRITE(IPRT,10) DE
      WRITE(IPRT,21) BETA
      WRITE (IPRT,30) A
      WRITE (IPRT,31) GAMMA
      WRITE(IPRT,35) RBB
25    FORMAT(2X,T10,'Equilibrium bond distances (Angstroms):',
     *       /,2X,T15,3(1PE20.10,1X))
10    FORMAT(2X,T10,'Dissociation energies (kcal/mol):',
     *          T50,1PE20.10)
21    FORMAT(2X,T10,'Morse Beta parameters (Angstroms**-1):',
     *          T50,1PE20.10)
35    FORMAT(2X,T10,'RBB (Angstroms):', 
     *          T50,1PE20.10)
30    FORMAT(2X,T10,'Pauling parameter (Angstroms):',
     *          T50,1PE20.10)
31    FORMAT(2X,T10,'Gamma (unitless):',
     *          T50,1PE20.10)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       CONVERT TO ATOMIC UNITS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO 50 IT = 1,3
      REQ(IT) = REQ(IT)/BOHR
   50 CONTINUE
      BETA = BETA * BOHR
      DE = DE/HRTREE
      A = A/BOHR
      RHH = 1.4007977D0
      RBB = RBB /BOHR
C
CCCCCCCCCCCC
      RETURN
100   WRITE (IPRT, 1000)
      STOP 'PREPOT 1'
1000  FORMAT (/,2X,T5,'Error opening potential input data file')
CCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        MOVE ONTO THE MAIN PART OF THE ROUTINE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCC
      ENTRY POT                                                         8/17R79
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
      VAC = GAMMA *1.0D-14 * DE * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))
C
C        CALCULATE V(RAB+RBC)
C
      JUNK(4) = EXP(BETA*(REQ(3)-RBC-RAB)/ROB)
      VABBC= GAMMA *1.0D-14 * DE * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))
C
C        HERE IS THE BENDING CORRECTION
C
      BCOOR = (VAC - VABBC)*1.0D14
C
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE
C          OF BCOOR WITH RESPECT TO RAB  RBC AND RAC.
C
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO RAC)
C
      IF (NDER .EQ. 1) THEN
          DAC = -GAMMA*DE*BETA*JUNK(3)*(1.D0+JUNK(3))/ROB
          DBCOOR(1) = GAMMA * DE * BETA * JUNK(4) * 
     *                (1.D0 + JUNK(4))/ROB
          DBCOOR(2)= DBCOOR(1)
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
      CALL POT2                                                         8/17R79
C
      IF (NDER .EQ. 1) THEN
          DER(1) = DER(1) + DER(3)
          DER(2) = DER(2) + DER(3)
          DER(3) = DAC
      ENDIF 
      RAC = RACTMP
C
      ENERGY = ENERGY + BCOOR
C
      IF (NDER .EQ. 1) THEN
          DO 300 IT = 1,2
                 DER(IT) = DBCOOR(IT) + DER(IT)
300       CONTINUE
      ENDIF
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C
C POTENTIAL FUNCTION FOR RATE1D SERIES
C RMO-SPLINE TYPE FIT WITH BEBO TYPE CONTINUATION INTO THE ASYMBTOTIC REGION
C POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM WITH UNITS OF EV/BOHR
C POTENTIAL ROUTINE ASSUMES SPLINE INFO IS IN UNITS OF KCAL/ANGSTROM
C
      SUBROUTINE PRPOT2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (HRTREE = 627.5095D0)
      COMMON /PT1CM/ R(3),E,D(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(21)
      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
      COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
      DIMENSION V(4)
      COMMON/NSTEP/STEP
      DIMENSION PHI(40),BSP(40,3),XBC(3),XAB(3),CBC(3),ABC(3)
     1,CAB(3),AAB(3)
      DIMENSION CSP(48,3),SCRAP(48)
      PARAMETER (DETRAD=.0174532D0)
      CHARACTER*80 TITLE
      save
C
C    Read in the potential parameters from the file linked to FORTRAN IREAD
C
      READ (IREAD,843)TITLE
      READ (IREAD,2) STEP
      READ (IREAD,401)NPHI
      READ (IREAD,403)R1S,R2S,VSP,XBC,XAB,CBC,ABC,CAB,AAB
      READ (IREAD,407)(PHI(I),I = 1,NPHI)
C
      DO 402 J = 1,3
             READ(IREAD,407)(BSP(I,J),I = 1,NPHI)
402   CONTINUE
C
      WRITE(IPRT,845)TITLE
      WRITE(IPRT,125)STEP
      WRITE(IPRT,111)XAB(1), XBC(1), XAB(1),
     *               XAB(2), XBC(2), XAB(2),
     *               XAB(3), XBC(3), XAB(3)
      WRITE(IPRT,105)R1S,R2S,VSP
      WRITE(IPRT,113)(PHI(I),(BSP(I,J),J = 1,3),I = 1,NPHI)
      WRITE(IPRT,117)(CBC(I),ABC(I),I = 1,3),(CAB(I),AAB(I),I = 1,3)
C
C   Initialize the energy in the asymptotic valleys
      EASYAB = XAB(1)/HRTREE
      EASYBC = XBC(1)/HRTREE 
      EASYAC = EASYAB
C
2     FORMAT(F12.10)
105   FORMAT(/,2X,T5,'R1S, R2S, VSP = ',3(F10.5,1X))
111   FORMAT (/,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',
     *        /,2X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5)
113   FORMAT(/,2X,T24,'Phi',T39,'De',T59,'Re',T73,'Beta',/,
     *      (2X,T20,F10.5,T35,F10.5,T55,F10.5,T70,F10.5))
117   FORMAT(/,2X,T5,'For the asymptote: MO parameter value = ',
     *               '(asymp. value)+C*EXP(-A*(RI-RIS))',
     *       /,2X,T10,'CBC(I); I=1,3: ',3(F15.5, 1X),
     *       /,2X,T10,'ABC(I); I=1,3: ',3(F15.5, 1X),
     *       /,2X,T10,'CAB(I); I=1,3: ',3(F15.5, 1X),
     *       /,2X,T10,'AAB(I); I=1,3: ',3(F15.5, 1X))
125   FORMAT(2X,T5,'Step size for the numerical derivative ',
     *               'calculation = ',1PE20.10)
401   FORMAT(16I5)
403   FORMAT(3E20.7)
407   FORMAT(4E20.7)
843   FORMAT(A80)
845   FORMAT(2X,T5,'Title card in potential input data file: ',/,A80)
C
      RETURN
C
C COCK SPLINES
C
      ENTRY POT2
      IC=1
   10 GO TO (81,82,83,84),IC
  81    R1=R(1)
      R2=R(2)
      R3=R(3)
      GO TO 85
   82 R1=R(1) + STEP
      R2=R(2)
      R3=R(3)
      GO TO 85
   83 R1=R(1)
      R2=R(2) + STEP
      R3=R(3)
      GO TO 85
   84 R1=R(1)
      R2=R(2)
      R3=R(3) + STEP
 85     R1=R1*0.529177D0
      R2=R2*0.529177D0
      R3=R3*0.529177D0
      DR1=R1S - R1
      DR2=R2S - R2
      ANGLE=ATAN(DR1/DR2)
      ANGLE=ANGLE/DETRAD
      DO 500 K=1,3
      PPP=SPLINE(1,NPHI,PHI,BSP(1,K),CSP(1,K),SCRAP,ANGLE)
 500  CONTINUE
C
C SETUP ANY CONSTANTS
C
      DBC = XBC(1)
      DERG = XBC(1)/23.061D0
      RERG = (R2S-XBC(2))/.529177D0
      BETRG = XBC(3)*.529177D0
      DEPR = XAB(1)/23.061D0
      REPR = (R1S-XAB(2))/.529177D0
      BETPR = XAB(3)*.529177D0
C
C ENTRY INTO POTENTIAL CALCULATOR
C
      IF(R1.GT.R1S) GO TO 100
      IF(R2.GT.R2S) GO TO 102
C
C INTERACTION REGION OF POLAR COORDINATES
C
      DR1=R1S-R1
      DR2=R2S-R2
      ANGLE=ATAN(DR1/DR2)
      ANGLE=ANGLE/DETRAD
      RR=SQRT(DR1*DR1+DR2*DR2)
      DIS=SPLINE(2,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
      REQ=SPLINE(2,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
      BET=SPLINE(2,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)
C      D(1)=SPLINE(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
C      D(2)=SPLINE(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
C30      D(1)=D(1)/627.5095D0*0.529177D0
C      D(2)=D(2)/627.5095D0*0.529177D0
30      X = DIS*((1.D0-EXP(BET*(RR-REQ)))**2.0D0 -1.D0) + VSP          12AUG84
      V(IC) = X
      IC=IC+1
      IF(IC.LT.5) GO TO 10
      DO 333 I=1,3
 333  V(I)=V(I)/627.5095D0
      IF (NDER .EQ. 1) THEN
          D(1)=(V(2)-V(1))/STEP
          D(2)=(V(3)-V(1))/STEP
          D(3)=0.0D0
      ENDIF
      E=V(1)
   93 RETURN
C
C LARGE R1 ASYMPTOTIC REGION
C
  100 IF(R2.GT.R2S) GO TO 104
      RR = R2S-R2
      SEP = R1-R1S
      DIS = XBC(1) + CBC(1)*EXP(-MIN(ABC(1)*SEP,25.D0))
      REQ = XBC(2) + CBC(2)*EXP(-MIN(ABC(2)*SEP,25.D0))
      BET = XBC(3) + CBC(3)*EXP(-MIN(ABC(3)*SEP,25.D0))
C      IF(R1.GT.10.D0)D(1)=1.0D1
C      IF(R1.GT.10.D0)D(2)=1.0D-8
CI     IF(R1.GT.10.D0)GO TO 30
C      D(1)=-CBC(1)*ABC(1)*((1.D0-EXP(-BET*(RR-REQ)))**2-1.D0) -DIS*
C     12.D0*(1.D0-EXP(-BET*(RR-REQ)))*EXP(-BET*(RR-REQ))*((-CBC(3)*ABC(
C     13)*EXP(-ABC(3)*SEP)*(RR-REQ) - CBC(2)*ABC(2)*EXP(ABC(2)*SEP)*
C     1BET))*BET*(RR-REQ)
C      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ)))*BET*(RR-REQ)*EXP(-BET*
C     1(RR-REQ))*BET
      GO TO 30
C
C LARGE R2 ASYMPTOTIC REGION
C
  102 CONTINUE
      RR = R1S-R1
      SEP = R2-R2S
      DIS = XAB(1) + CAB(1)*EXP(-MIN(AAB(1)*SEP,25.D0))
      REQ = XAB(2) + CAB(2)*EXP(-MIN(AAB(2)*SEP,25.D0))
      BET = XAB(3) + CAB(3)*EXP(-MIN(AAB(3)*SEP,25.D0))
C      IF(R2.GT.10.D0)D(2)=1.0D1
C      IF(R2.GT.10.D0)D(1)=1.0D-8
C      IF(R2.GT.10.D0)GO TO 30
C      D(1)=-CAB(1)*AAB(1)*((1.D0-EXP(-BET*(RR-REQ)))**2-1.D0) -DIS*
C     12.D0*(1.D0-EXP(-BET*(RR-REQ)))*EXP(-BET*(RR-REQ))*((-CAB(3)*AAB(
C     13)*EXP(-AAC(3)*SEP)*(RR-REQ) - CAB(2)*AAB(2)*EXP(AAB(2)*SEP)*
C     1BET))*BET*(RR-REQ)
C      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ)))*BET*(RR-REQ)*EXP(-BET*
C     1(RR-REQ))*BET
      GO TO 30
C
C LARGE R1,R2 3 BODY BREAK-UP REGION
C
  104 CONTINUE
      J=J+1
      E=XBC(1)/627.5095D0
      IF (NDER .EQ. 1) THEN 
          D(1)=1.0D-9
          D(2)=1.0D-9
          D(3)=0.0D0
      ENDIF
      GO TO 93
      END
      FUNCTION SPLINE(ISW,NN,X,Y,C,D,XPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),C(*),D(*)
C
C
C  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE
C  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT
C  OR ITS DERIVATIVE AT A GIVEN POINT
C
C  CALLING SEQUENCE .......
C
C     ISW ... CONTROL OPTION
C         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS
C                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS
C                ARE STORED IN THE ARRAY C.
C         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS
C                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY
C                THE PARAMETER 'XPT'.
C         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATED AT XPT
C         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WITH ISW=3
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
C              SET TO 2).  THE VALUE OF SPLINE IS SET TO THE
C              INTERPOLATED VALUE.
C
C *****  USER NOTES  *****
C
C     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......
C
C       A.  CALL SPLINE WITH THE KNOTS.  THIS SETS UP THE
C           COEFFICIENT ARRAY C.
C           EG.  DUMY=SPLINE(1,NN,X,Y,C,D,XPT)
C
C       B.  CALL SPLINE WITH THE ARRAY C WHICH WAS DEFINED BY THE
C           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE
C           POINT 'XPT' .
C           EG.   VALUE=SPLINE(2,NN,X,Y,C,D,XPT)
C
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.
C
C     Output from this subprogram is written to unit 6.
C
2     N=NN
      NP1=N+1
      NP2=N+2
      Z=XPT
24    GO TO (4,5,6,7),ISW
4     C(1)=Y(1)
      D(1)=1.0D0
      C(NP1)=0.0D0
      D(NP1)=0.0D0
      C(NP2)=0.0D0
      D(NP2)=0.0D0
      DO 41 I=2,N
      C(I)=Y(I)-Y(1)
41    D(I)=X(I)-X(1)
      DO 410 I=3,NP2
      IF(D(I-1).NE.0)GO TO 43
      WRITE(6,1001)
      STOP 'SPLINE 1'
43    PIVOT=1.0D0/D(I-1)
      IF(I.GE.NP2)GO TO 45
      SUPD=X(I-1)-X(I-2)
      IF(SUPD.GE.0)GO TO 44
      WRITE(6,1000)
      STOP 'SPLINE 2'
44    SUPD=SUPD*SUPD*SUPD
      GO TO 46
45    SUPD=1.0D0
46    DFACT=SUPD*PIVOT
      CFACT=C(I-1)*PIVOT
      IF(I.GT.N)GO TO 48
      DO 47 J=I,N
      V=X(J)-X(I-2)
      C(J)=C(J)-D(J)*CFACT
47    D(J)=V*V*V-D(J)*DFACT
48    CONTINUE
      IF(I.GE.NP2)GO TO 49
      C(NP1)=C(NP1)-D(NP1)*CFACT
      D(NP1)=1.0D0-D(NP1)*DFACT
49    C(NP2)=C(NP2)-D(NP2)*CFACT
410   D(NP2)=X(I-2)-D(NP2)*DFACT
      DO 411 I=1,N
      J=NP2-I
      IF(J.NE.NP1)GO TO 413
      V=1.0D0
      GO TO 414
413   V=X(J)-X(J-1)
      V=V*V*V
414   IF(D(J+1).NE.0)GO TO 415
      WRITE(6,1001)
      STOP 'SPLINE 3'
415   C(J+1)=C(J+1)/D(J+1)
411   C(J)=C(J)-C(J+1)*V
      IF(D(2).NE.0)GO TO 416
      WRITE(6,1001)
      STOP 'SPLINE 4'
416   C(2)=C(2)/D(2)
      RETURN
5     SPLINE=C(1)+C(2)*(Z-X(1))
      DO 51 I=1,N
      V=Z-X(I)
      IF(V.GT.0)GO TO 51
      RETURN
51    SPLINE=SPLINE+C(I+2)*V*V*V
      RETURN
    6 CONTINUE
      SPLINE=C(1)+C(2)*(Z-X(1))
      DERIV = C(2)/3.D0
      DO 53 I = 1,N
      V=Z-X(I)
      IF(V.LE.0) RETURN
      V2 = V*V
      SPLINE = SPLINE + C(I+2)*V2*V
      DERIV = DERIV + C(I+2)*V2
   53 CONTINUE
      RETURN
    7 CONTINUE
      SPLINE = 3.D0*DERIV
      RETURN
1000  FORMAT(1X,5X,'***** ERROR IN SPLINE ... UNORDERED X-VALUES
     1*****')
1001  FORMAT(1X,5X,' ***** ERROR IN SPLINE ... DIVIDE FAULT *****')
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA IREAD /4/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/
C
         END
C*****
