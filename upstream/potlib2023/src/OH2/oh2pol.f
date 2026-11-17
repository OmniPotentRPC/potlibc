C
      SUBROUTINE PREPOT
C
C   System:          OH2
C   Functional form: RMO-Spline fit
C   Common name:     POL
C   Data file:       potoh2pol.dat
C   Reference:       S. P. Walch, T. H. Dunning, Jr., F. W. Bobrowicz, R. Raffenetti, 
C                    J. Chem. Phys. 72, 406 (1980); S. P. Walch, A. F. Wagner, G. C. Schatz, 
C                    and T. H. Dunning, Jr., J. Chem. Phys. 72, 2894 (1980); 
C                    A. F. Wagner, G. C. Schatz, J. M. Bowman, J. Chem. Phys. 74, 4960 (1981); 
C                    G. C. Schatz, A. F. Wagner, S. P. Walch, J. M. Bowman, J. Chem. Phys. 74, 4984 (1981); 
C                    J. M. Bowman and K. T. Lee in \"Potential Energy Surfaces and Dynamics Calculations,\"
C                    D. G. Truhlar, ed. (Plenum, 1981), p. 359.
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are read in from the file "potoh2pol.dat" which 
C   is linked to FORTRAN unit number IREAD.
C   The control flags for the potential are initialized in the block 
C   data subprogram PPARMS.
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
C          ANTI-MORSE BENDING POTENTIAL
C          SUBROUTINE BEND WILL ADD A BENDING CORRECTION TO ANY
C          COLLINEAR POTENTIAL SUBROUTINE.
C
C          PREBEND MUST BE CALLED BEFORE ANY CALLS TO BEND
C          CAN BE MADE.  THE POTENTIAL PARAMETERS FOR THE COL-
C          LINEAR POTENTIAL ARE READ IN BY A CALL TO PREPOT
C          (CONTAINED IN PREBEND).  THEN THE BENDING PARAMETERS ARE
C          COORDINATES, POTENTIAL ENERGY, AND DERIVATIVES ARE
C          PASSED THROUGH COMMON /POTCM/ RAB,RBC,RAC,ENERGY,DER(1--3).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C          DEFINE VARIABLES
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PT1CM/ RAB,RBC,RAC,ENERGY,DER(3)
      COMMON /PT2CM/ NSURF, NDER, NFLAG(20)
      COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
      DOUBLE PRECISION NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR
      DIMENSION DBCOOR(2),JUNK(4),REQ(3)
C   Conversion for Angstroms to Bohr (Angstrom/Bohr)
      PARAMETER (BOHR = .52917706D0)
C   Conversion for kcal to Hartree (kcal/Hartree)
      PARAMETER (HRTREE = 627.5095D0)
      SAVE
C
C   Echo the name of the potential routine to FORTRAN unit IPRT
      WRITE (IPRT, 600)
C
C   Open the file containing the parameters for the potential
      OPEN (UNIT=IREAD, FILE = 'potoh2pol.dat', FORM='FORMATTED', 
     *      STATUS='OLD', ERR = 100)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      READ COLLINEAR POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE(IPRT, 610)
      CALL PREPT2       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      READ BENDING POTENTIAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE(IPRT, 620)
C
C      DISSOCIATION ENERGY IN KCAL/MOLE
      READ(IREAD,*) DE
C
C       READ IN MORSE BETAS--RECIPROCAL ANGSTROMS
      READ (IREAD,*) BETA
C
C        READ IN THE EQUILIBRIUM DISTANCES IN ANGSTROMS
      READ (IREAD,*) (REQ(IT), IT = 1,3)
C
C        READ IN THE PAULING PARAMETER IN ANGSTROMS
      READ (IREAD,*) A
C
C       READ IN GAMMA FOR THE BEND CORRECTION
C       GAMMA IS DIMENSIONLESS AND RANGES FROM .5 TO .65
      READ (IREAD,*) GAMMA
C
C       READ IN RBB IN ANGSTROMS--THIS NUMBER IS USED
C       IN A SCALING CALCULATION FOR THE BENDING CORRECTION
      READ(IREAD,*) RBB
C
C   Close the input data file
      CLOSE (UNIT = IREAD)
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
C    Echo the potential parameters to the file linked to FORTRAN unit IREAD
      WRITE (IPRT,630) (REQ(IT), IT = 1,3)
      WRITE(IPRT,631) DE
      WRITE(IPRT,632) BETA
      WRITE(IPRT,635) RBB
      WRITE (IPRT,633) A
      WRITE (IPRT,634) GAMMA
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',
     *                'potential energy surface POL')
610   FORMAT (/,2X,T5,'Parameters for the collinear part of the ',
     *                'potential:')
620   FORMAT (/,2X,T5,'Parameters for the bending potential:')
630   FORMAT(2X,T10,'Equilibrium bond distances (Angstroms):',
     *       /,2X,T15,3(1PE20.10,1X))
631   FORMAT(2X,T10,'Dissociation energy (kcal/mol):',
     *          T50,1PE20.10)
632   FORMAT(2X,T10,'Morse Beta parameter (Angstroms**-1):',
     *          T50,1PE20.10)
635   FORMAT(2X,T10,'RBB (Angstroms):', 
     *          T50,1PE20.10)
634   FORMAT(2X,T10,'Gamma (unitless):',
     *          T50,1PE20.10)
633   FORMAT(2X,T10,'Pauling parameter (Angstroms):',
     *          T50,1PE20.10)
6000  FORMAT (/,2X,T5,'Error opening potential input data file')
CCCCCCCCCCCC
      RETURN
CCCCCCCCCCCC
100   WRITE (IPRT, 6000)
      STOP 'PREPOT 1'
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
      CALL POT2
C   Add the energy terms together
      ENERGY = ENERGY + BCOOR
C
C   If NDER = 1, compute the derivatives
      IF (NDER .NE. 1) RETURN
      DER(1) = DER(1) + DER(3)
      DER(2) = DER(2) + DER(3)
      DER(3) = DAC
      RAC = RACTMP
C
      DO 300 IT = 1,2
           DER(IT) = DBCOOR(IT) + DER(IT)
  300 CONTINUE
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
C POTENTIAL FUNCTION FOR RATEID SERIES
C RMO-SPLINE TYPE FIT
C POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM IN UNITS OF HARTREE/
C BOHR
C POTENTIAL ROUTINE ASSUMES SPLINE INFO IN UNITS OF KCAL/ANGSTROM
      SUBROUTINE PREPT2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(4),PHI(20),BSP(20,5),CSP(22,5)
      COMMON /PT1CM/ R(3),E,D(3)
      COMMON /PT2CM/ NSURF, NDER, NDUM(21)
      COMMON /PT4CM/ IPRT, IREAD, IDUM(18)
      CHARACTER*4 TITLE
      SAVE
      READ(IREAD,843)TITLE
 843  FORMAT(5A4)
      A1INF=0.0D0
      A2INF=0.0D0
      IAYES= 0
      READ(IREAD,4)D1INF,R1EINF,B1INF,D2INF,R2INF,B2INF,A1INF,A2INF
   4  FORMAT(8F9.5)
      IF(A2INF .EQ.0.0D0) A2INF= D1INF-D2INF
      READ(IREAD,6)R1S,R2S,VSP,IAYES,ALPH1,ALPH2
   6  FORMAT(3F10.5,I5,2F10.5)
      READ(IREAD,8)STEP
   8  FORMAT(F15.10)
      WRITE(IPRT,845)TITLE
845   FORMAT(2X,T5,'Title card from input data file:',5A4)
      WRITE(IPRT,121)D1INF,R1EINF,B1INF,D2INF,R2INF,B2INF,A1INF,A2INF
121   FORMAT(2X,T5,'D1INF,R1EINF,B1INF,D2INF,R2EINF,B2INF,A1INF,',
     *             'A2INF (kcal/mol):',/,(2X,T5,4(F10.5,1X)))
      WRITE(IPRT,123)R1S,R2S,VSP,ALPH1,ALPH2
123   FORMAT(2X,T5,'R1S,R2S,VSP,ALPH1,ALPH2 (ang/kcal):',
     *       /,2X,T10,5(F10.5,1X))
      IF(IAYES.NE.0) WRITE(IPRT,133)
133   FORMAT(2X,T5,'Warning: Pivot point is fixed at VSP')
      IF(IAYES.EQ.0) WRITE(IPRT,131)
131   FORMAT(2X,T5,'Warning: Pivot point is not fixed; VSP is ',
     *               'is 3-body break-up energy')
      WRITE(IPRT,125)STEP
125   FORMAT(2X,T5,'Step for numerical derivative = ', F15.10)
      READ(IREAD,401)NPHI,NV
      READ(IREAD,403)(PHI(I),I = 1,NPHI)
      DO 402 J = 1,NV
      READ(IREAD,403)(BSP(I,J),I = 1,NPHI)
 402  CONTINUE
      NPHIP2 = NPHI + 2
      DO 404 J = 1,NV
      READ(IREAD,403)(CSP(I,J),I = 1,NPHIP2)
 404  CONTINUE
 401  FORMAT(16I5)
 403  FORMAT(4E20.7)
      WRITE(IPRT,111)NPHI,NV
111   FORMAT(2X,T5,'NPHI, NV = ',2(I5,1X))
      WRITE (IPRT, 700)
700   FORMAT(/,2X,T5,'I, PHI(I), (BSP(I,J),CSP(I,J);J=1,NV):')
      DO 20 I = 1,NPHI
      WRITE(IPRT,103)I,PHI(I),(BSP(I,J),CSP(I,J),J=1,NV)
103   FORMAT(2X,T5,I5,1X,F10.5,(T25,F10.5,1X,E15.5))
  20  CONTINUE
      WRITE (IPRT, 701)
701   FORMAT(/,2X,T5,'K,CSP(K,J);J=1,NV')
      DO 22 I=1,2
      K = I + NPHI
      WRITE(IPRT,105)K,(CSP(K,J),J=1,NV)
105   FORMAT(2X,T5,I5,5X,4(1X,E15.5))
  22  CONTINUE
      D1 = BSP(1,1)
      D2 = BSP(NPHI,1)
      R1E = BSP(1,2)
      R2E = BSP(NPHI,2)
      B1 = BSP(1,3)
      B2 = BSP(NPHI,3)
      A1 = VSP-D1*(1.0D0-EXP(-B1*R1E))**2.0D0
      A2 = VSP-D2*(1.0D0-EXP(-B2*R2E))**2.0D0
      R1E = R2S -R1E
      R2E = R1S -R2E
      DEPR = D2INF/23.061D0
      BETPR = B2INF*0.529177D0
      REPR = R2INF/0.529177D0
      DERG = D1INF/23.061D0
      BETRG = B1INF*0.529177D0
      RERG = R1EINF/0.529177D0
      RETURN
      ENTRY POT2
      IC=1
  86  GO TO (81,82,83,84),IC
  81  R1=R(1)
      R2=R(2)
      R3=R(3)
      GO TO 85
  82  R1=R(1) + STEP
      R2=R(2)
      R3=R(3)
      GO TO 85
  83  R1=R(1)
      R2=R(2) + STEP
      R3=R(3)
      GO TO 85
  84  R1=R(1)
      R2=R(2)
      R3=R(3) + STEP
  85  R1=R1*0.529177D0
      R2=R2*0.529177D0
      R3=R3*0.529177D0
      DATA DETRAD/0.0174532D0/
      IF(R1.GT.R1S)GO TO 100
      IF(R2.GT.R2S)GO TO 102
C.....DO NORMAL SPLINE INTERPOLATION
      DR1=R1S-R1
      DR2=R2S-R2
      ANGLE=ATAN(DR1/DR2)
      ANGLE=ANGLE/DETRAD
      RR=SQRT(DR1*DR1 + DR2*DR2)
      DIS=SPLINE(2,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)
      REQ=SPLINE(2,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)
      BET=SPLINE(2,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)
      IF(IAYES.EQ.0) A=VSP - DIS
      IF(IAYES.NE.0)
     1A=VSP-DIS*(1.0D0-EXP(-BET*REQ))**2.0D0
      POTSP=DIS*(1.0D0-EXP(BET*(RR-REQ)))**2.0D0 + A
      GO TO 50
 100  IF(R2.LT.R2S) GO TO 101
      E=VSP/627.5095D0
      D(1)=1.0D-9
      D(2)=1.0D-9
      D(3)=0.0D0
      RETURN
  101 CONTINUE
      DR=D1-D1INF
      DRE=R1E-R1EINF
      DBE=B1-B1INF
      DA=A1-A1INF
      DELR=R1-R1S
      FACTR=EXP(-DELR*ALPH1)
C.....INTERPOLATE FOR MORSE PARAMETERS
      DER=D1INF + DR*FACTR
      RER=R1EINF+DRE*FACTR
      BER=B1INF+DBE*FACTR
      IF(IAYES.EQ.0) A = VSP-DER
      IF(IAYES.NE.0) A=A1INF+DA*FACTR
      DR2=R2-RER
      POTSP=DER*(1.0D0-EXP(-BER*DR2))**2.0D0 + A                        13OCT83
      GO TO 50
 102  CONTINUE
      DR=D2-D2INF
      DRE=R2E-R2INF
      DBE=B2-B2INF
      DA=A2-A2INF
      DELR=R2-R2S
      FACTR=EXP(-DELR*ALPH2)
      DER=D2INF+DR*FACTR
      RER=R2INF + DRE*FACTR
      BER=B2INF + DBE*FACTR
      IF(IAYES.EQ.0) A = VSP-DER
      IF(IAYES.NE.0) A=A2INF+DA*FACTR
      DR1=R1-RER
      POTSP=DER*(1.0D0-EXP(-BER*DR1))**2.0D0+A
  50  V(IC)=POTSP
      IC=IC+1
      IF(IC.LT.4) GO TO 86
      DO 33 I=1,3
  33  V(I)=V(I)/627.5095D0
      E=V(1)
      IF (NDER .NE. 1) RETURN
      D(1)=(V(2)-V(1))/STEP
      D(2)=(V(3)-V(1))/STEP
      D(3)=0.0D0
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SPLINE(ISW,NN,X,Y,C,D,XPT) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION X(*),Y(*),C(*),D(*)                                             
      PARAMETER (IUNIT = 6)
C                                                                               
C                                                                               
C  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE                    
C  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT                   
C  OR ITS DERIVATIVE AT A GIVEN POINT                                           
C                                                                               
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
C                                                                               
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.             
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.                        
C                                                                               
C     Output from this subprogram is written to the file linked to 
C     FORTRAN unit IUNIT
C                                                                               
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
      WRITE(IUNIT,1001) 
      STOP 'SPLINE 1'
43    PIVOT=1.0D0/D(I-1)                                                        
      IF(I.GE.NP2)GO TO 45                                                      
      SUPD=X(I-1)-X(I-2)                                                        
      IF(SUPD.GE.0)GO TO 44                                                     
      WRITE(IUNIT,1000)
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
      WRITE(IUNIT,1001)  
      STOP 'SPLINE 3'
415   C(J+1)=C(J+1)/D(J+1)                                                      
411   C(J)=C(J)-C(J+1)*V                                                        
      IF(D(2).NE.0)GO TO 416                                                    
      WRITE(IUNIT,1001)
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
         DATA NDER /1/
         DATA NSURF /0/
         DATA NDUM /21*0/
C
         END
C*****
