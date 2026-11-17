C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          OH2                                                        
C   Functional form: RMO-Spline fit                                             
C   Common name:     POL                                                        
C   Interface:       potlib2001
C   Reference:
C
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are assigned in the Block Data Subprogram          
C   PTPACM.  The control flags for the potential also are initialized           
C   in the Block Data Subprogram PPTARM.                                        
C   Coordinates, potential energy, and derivatives are passed                   
C   All the information passed through the common block PT1CM                   
C   is in Hartree atomic units.                                                 
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(O-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-O)                                           
C   The classical potential energy is set equal to zero for the O               
C   infinitely far from the equilibrium H2 diatomic.                            
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states availalble                            
C                 (1,1) = 1 as only gs state available                          
C        NDER  = 0 => no derivatives should be calculated                       
C        NDER  = 1 => calculate first derivatives                               
C        NFLAG - these integer values can be used to flag options               
C                within the potential;                                          
C                                                                               
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
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
C          PASSED THROUGH COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3).               
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCC                                                     
C          DEFINE VARIABLES                                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCC                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 JUNK
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
C   Conversion for Angstroms to Bohr (Angstrom/Bohr)                            
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
      PARAMETER (HRTREE = 627.5095D0)                                           
C                                                                               
      COMMON /AMBPCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                            
      COMMON /PRECM/  DBCOOR(2),JUNK(4),BCOOR                                   
      COMMON /PRE2CM/ D1INF,R1EINF,B1INF,D2INF,                                 
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                V(4),PHI(20),BSP(20,5),CSP(22,5),                         
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                STEP,D1,D2,R1E,R2E,B1,B2,A1,A2,                           
     +                IAYES,NPHI,NV                                             
C                                                                               
C   Echo the name of the potential routine to UNIT NFLAG(18)                    
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      WRITE (NFLAG(18), 600)                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C      ASSIGN COLLINEAR POTENTIAL PARAMETERS                                    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C                                                                               
      WRITE(NFLAG(18), 610)                                                     
C                                                                               
      CALL PREPOT2                                                              
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C      ASSIGN BENDING POTENTIAL PARAMETERS                                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C                                                                               
      WRITE(NFLAG(18), 620)                                                     
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
C    Echo the potential parameters to the file linked to UNIT NFLAG(18)         
C                                                                               
      WRITE(NFLAG(18),630) (REQ(IT), IT = 1,3)                                  
      WRITE(NFLAG(18),631) DE                                                   
      WRITE(NFLAG(18),632) BETA                                                 
      WRITE(NFLAG(18),635) RBB                                                  
      WRITE(NFLAG(18),633) A                                                    
      WRITE(NFLAG(18),634) GAMMA                                                
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
CCCCCCCCCCCC                                                                    
C                                                                               
      EZERO(1)=0.0d0                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
C                                                                               
      INDEXES(1) = 8                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
CCCCCCCCCCCC                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C        MOVE TO THE POT SUBROUTINE                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(O-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-O)                                           
C   The classical potential energy is set equal to zero for the O               
C   infinitely far from the equilibrium H2 diatomic.                            
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C      ENTRY POT                                                                
CCCCCCCCCCCCCCCC                                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK                                             
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
C   Conversion for Angstroms to Bohr (Angstrom/Bohr)                            
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
      PARAMETER (HRTREE = 627.5095D0)                                           
C                                                                               
      COMMON /AMBPCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                            
      COMMON /PRECM/  DBCOOR(2),JUNK(4),BCOOR                                   
      COMMON /PRE2CM/ D1INF,R1EINF,B1INF,D2INF,                                 
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                V(4),PHI(20),BSP(20,5),CSP(22,5),                         
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                STEP,D1,D2,R1E,R2E,B1,B2,A1,A2,                           
     +                IAYES,NPHI,NV                                             
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Check the values of NASURF and NDER for validity.                           
C                                                                               
      IF (NASURF(1,1) .EQ. 0) THEN                                              
         WRITE(NFLAG(18), 900) NASURF(1,1)                                      
         STOP                                                                   
      ENDIF                                                                     
         IF (NDER .GT. 1) THEN                                                  
             WRITE (NFLAG(18), 910) NDER                                        
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
      NAB = EXP((REQ(1) - R(1))/A)                                              
C                                                                               
      IF (NAB.LT.1.D-15) NAB = 1.D-15                                           
C                                                                               
      NBC = EXP((REQ(2) - R(2))/A)                                              
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
C        CALCULATE V(R(3))                                                      
C                                                                               
      JUNK(3) = EXP(-BETA*(R(3)-REQ(3))/ROB)                                    
      VAC = GAMMA * DE * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))                        
C                                                                               
C        CALCULATE V(R(1)+R(2))                                                 
C                                                                               
      JUNK(4) = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)                                
      VABBC= GAMMA * DE * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))                      
C                                                                               
C        HERE IS THE BENDING CORRECTION                                         
C                                                                               
      BCOOR = VAC - VABBC                                                       
C                                                                               
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE                      
C          OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).                        
C                                                                               
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
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
C          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR                  
C          GEOMETRY AND CALL POT                                                
C                                                                               
      RACTMP = R(3)                                                             
      R(3) = R(2) + R(1)                                                        
C                                                                               
      CALL POT2                                                                 
C                                                                               
C   Add the energy terms together                                               
C                                                                               
      ENGYGS = ENGYGS + BCOOR                                                   
C                                                                               
C   If NDER = 1, compute the derivatives                                        
C                                                                               
      IF (NDER .NE. 1) THEN                                                     
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
          RETURN                                                                
      END IF                                                                    
      DEGSDR(1) = DEGSDR(1) + DEGSDR(3)                                         
      DEGSDR(2) = DEGSDR(2) + DEGSDR(3)                                         
      DEGSDR(3) = DAC                                                           
      R(3) = RACTMP                                                             
C                                                                               
      DO 300 IT = 1,2                                                           
           DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT)                                 
  300 CONTINUE                                                                  
C                                                                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT2                                                        
C                                                                               
C   POTENTIAL FUNCTION FOR RATEID SERIES                                        
C   RMO-SPLINE TYPE FIT                                                         
C   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM                            
C   IN UNITS OF HARTREE/BOHR                                                    
C   POTENTIAL ROUTINE ASSUMES SPLINE INFO IN UNITS OF KCAL/ANGSTROM             
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK                                             
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
      PARAMETER (DETRAD = 0.0174532D0)                                          
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (BOHR = 0.529177D0)                                             
C                                                                               
      COMMON /AMBPCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                            
      COMMON /PRECM/  DBCOOR(2),JUNK(4),BCOOR                                   
      COMMON /PRE2CM/ D1INF,R1EINF,B1INF,D2INF,                                 
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                V(4),PHI(20),BSP(20,5),CSP(22,5),                         
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                STEP,D1,D2,R1E,R2E,B1,B2,A1,A2,                           
     +                IAYES,NPHI,NV                                             
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
C  TITLE                                                                        
C                                                                               
      TITLE='   O+H2 AB 9-12-79 '                                               
C                                                                               
C  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)               
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
C      A1INF=0.0D0                                                              
C      A2INF=0.0D0                                                              
C      IAYES= 0                                                                 
C                                                                               
C  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS            
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
      IF(A2INF.EQ.0.0D0) A2INF= D1INF-D2INF                                     
C                                                                               
C  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,                              
C  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS                                  
C  IN BLOCK DATA SUBROUTINE.  ECHO.                                             
C                                                                               
      WRITE(NFLAG(18),845)TITLE                                                 
      WRITE(NFLAG(18),121)D1INF,R1EINF,B1INF,D2INF,R2INF,B2INF,A1INF,A2I        
C     NF                                                                        
      WRITE(NFLAG(18),123)R1S,R2S,VSP,ALPH1,ALPH2                               
C                                                                               
      IF(IAYES.NE.0) WRITE(NFLAG(18),133)                                       
      IF(IAYES.EQ.0) WRITE(NFLAG(18),131)                                       
C                                                                               
      WRITE(NFLAG(18),125)STEP                                                  
C                                                                               
      NPHIP2 = NPHI + 2                                                         
C                                                                               
C  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE                   
C                                                                               
      WRITE(NFLAG(18),111)NPHI,NV                                               
      WRITE (NFLAG(18), 700)                                                    
C                                                                               
      DO 20 I = 1,NPHI                                                          
          WRITE(NFLAG(18),103)I,PHI(I),(BSP(I,J),CSP(I,J),J=1,NV)               
  20  CONTINUE                                                                  
                                                                                
      WRITE (NFLAG(18), 701)                                                    
                                                                                
      DO 22 I=1,2                                                               
         K = I + NPHI                                                           
         WRITE(NFLAG(18),105)K,(CSP(K,J),J=1,NV)                                
  22  CONTINUE                                                                  
C                                                                               
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
C                                                                               
      DEPR = D2INF/23.061D0                                                     
      BETPR = B2INF*0.529177D0                                                  
      REPR = R2INF/0.529177D0                                                   
      DERG = D1INF/23.061D0                                                     
      BETRG = B1INF*0.529177D0                                                  
      RERG = R1EINF/0.529177D0                                                  
C                                                                               
103   FORMAT(2X,T5,I5,1X,F10.5,(T25,F10.5,1X,E15.5))                            
105   FORMAT(2X,T5,I5,5X,4(1X,E15.5))                                           
111   FORMAT(2X,T5,'NPHI, NV = ',2(I5,1X))                                      
121   FORMAT(2X,T5,'D1INF,R1EINF,B1INF,D2INF,R2EINF,B2INF,A1INF,',              
     *             'A2INF (kcal/mol):',/,(2X,T5,4(F10.5,1X)))                   
123   FORMAT(2X,T5,'R1S,R2S,VSP,ALPH1,ALPH2 (ang/kcal):',                       
     *       /,2X,T10,5(F10.5,1X))                                              
125   FORMAT(2X,T5,'Step for numerical derivative = ', F15.10)                  
133   FORMAT(2X,T5,'Warning: Pivot point is fixed at VSP')                      
131   FORMAT(2X,T5,'Warning: Pivot point is not fixed; VSP is ',                
     *               'is 3-body break-up energy')                               
700   FORMAT(/,2X,T5,'I, PHI(I), (BSP(I,J),CSP(I,J);J=1,NV):')                  
701   FORMAT(/,2X,T5,'K,CSP(K,J);J=1,NV')                                       
845   FORMAT(2X,T5,'Title card for potential: ',A20)                            
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT2                                                           
C                                                                               
C      ENTRY POT2                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK                                             
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
      PARAMETER (DETRAD = 0.0174532D0)                                          
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (BOHR = 0.529177D0)                                             
C                                                                               
      COMMON /AMBPCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                            
      COMMON /PRECM/  DBCOOR(2),JUNK(4),BCOOR                                   
      COMMON /PRE2CM/ D1INF,R1EINF,B1INF,D2INF,                                 
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                V(4),PHI(20),BSP(20,5),CSP(22,5),                         
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                STEP,D1,D2,R1E,R2E,B1,B2,A1,A2,                           
     +                IAYES,NPHI,NV                                             
      DIMENSION SCRAP(20)                                                       
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
      IC=1                                                                      
  86  GO TO (81,82,83,84),IC                                                    
C                                                                               
  81  R1=R(1)                                                                   
      R2=R(2)                                                                   
      R3=R(3)                                                                   
      GO TO 85                                                                  
C                                                                               
  82  R1=R(1) + STEP                                                            
      R2=R(2)                                                                   
      R3=R(3)                                                                   
      GO TO 85                                                                  
C                                                                               
  83  R1=R(1)                                                                   
      R2=R(2) + STEP                                                            
      R3=R(3)                                                                   
      GO TO 85                                                                  
C                                                                               
  84  R1=R(1)                                                                   
      R2=R(2)                                                                   
      R3=R(3) + STEP                                                            
C                                                                               
  85  R1=R1*0.529177D0                                                          
      R2=R2*0.529177D0                                                          
      R3=R3*0.529177D0                                                          
C                                                                               
      IF(R1.GT.R1S)GO TO 100                                                    
      IF(R2.GT.R2S)GO TO 102                                                    
C                                                                               
C.....DO NORMAL SPLINE INTERPOLATION                                            
C                                                                               
      DR1=R1S-R1                                                                
      DR2=R2S-R2                                                                
      ANGLE=ATAN(DR1/DR2)                                                       
      ANGLE=ANGLE/DETRAD                                                        
      RR=SQRT(DR1*DR1 + DR2*DR2)                                                
      DIS= SPLINE(2,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                     
      REQ2=SPLINE(2,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                     
      BET= SPLINE(2,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                     
      IF(IAYES.EQ.0) ANU=VSP - DIS                                              
      IF(IAYES.NE.0) ANU=VSP-DIS*(1.0D0-EXP(-BET*REQ2))**2.0D0                  
      POTSP=DIS*(1.0D0-EXP(BET*(RR-REQ2)))**2.0D0 + ANU                         
      GO TO 50                                                                  
 100  IF(R2.LT.R2S) GO TO 101                                                   
      ENGYGS=VSP/627.5095D0                                                     
      DEGSDR(1)=1.0D-9                                                          
      DEGSDR(2)=1.0D-9                                                          
      DEGSDR(3)=0.0D0                                                           
      RETURN                                                                    
  101 CONTINUE                                                                  
      DR=D1-D1INF                                                               
      DRE=R1E-R1EINF                                                            
      DBE=B1-B1INF                                                              
      DA=A1-A1INF                                                               
      DELR=R1-R1S                                                               
      FACTR=EXP(-DELR*ALPH1)                                                    
C                                                                               
C.....INTERPOLATE FOR MORSE PARAMETERS                                          
C                                                                               
      DER=D1INF + DR*FACTR                                                      
      RER=R1EINF+DRE*FACTR                                                      
      BER=B1INF+DBE*FACTR                                                       
      IF(IAYES.EQ.0) ANU = VSP-DER                                              
      IF(IAYES.NE.0) ANU=A1INF+DA*FACTR                                         
      DR2=R2-RER                                                                
      POTSP=DER*(1.0D0-EXP(-BER*DR2))**2.0D0 + ANU                        13OCT8
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
      IF(IAYES.EQ.0) ANU = VSP-DER                                              
      IF(IAYES.NE.0) ANU=A2INF+DA*FACTR                                         
      DR1=R1-RER                                                                
      POTSP=DER*(1.0D0-EXP(-BER*DR1))**2.0D0+ANU                                
  50  V(IC)=POTSP                                                               
      IC=IC+1                                                                   
      IF(IC.LT.4) GO TO 86                                                      
      DO 33 I=1,3                                                               
  33  V(I)=V(I)/627.5095D0                                                      
      ENGYGS=V(1)                                                               
      IF (NDER .NE. 1) RETURN                                                   
      DEGSDR(1)=(V(2)-V(1))/STEP                                                
      DEGSDR(2)=(V(3)-V(1))/STEP                                                
      DEGSDR(3)=0.0D0                                                           
      RETURN                                                                    
      END                                                                       
C                                                                               
      DOUBLE PRECISION FUNCTION SPLINE(ISW,NN,X,Y,C,D,XPT)                      
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
      DIMENSION X(NN),Y(NN),C(NN+2),D(NN+2)                                     
C     PARAMETER (IPRT = 6)                                                      
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
C     UNIT NFLAG(18)                                                            
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
      WRITE(NFLAG(18),1001)                                                     
      STOP 'SPLINE 1'                                                           
43    PIVOT=1.0D0/D(I-1)                                                        
      IF(I.GE.NP2)GO TO 45                                                      
      SUPD=X(I-1)-X(I-2)                                                        
      IF(SUPD.GE.0)GO TO 44                                                     
      WRITE(NFLAG(18),1000)                                                     
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
      WRITE(NFLAG(18),1001)                                                     
      STOP 'SPLINE 3'                                                           
415   C(J+1)=C(J+1)/D(J+1)                                                      
411   C(J)=C(J)-C(J+1)*V                                                        
      IF(D(2).NE.0)GO TO 416                                                    
      WRITE(NFLAG(18),1001)                                                     
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
     1       *****')                                                            
1001  FORMAT(1X,5X,' ***** ERROR IN SPLINE ... DIVIDE FAULT *****')             
      END                                                                       
C                                                                               
C*****                                                                          
C                                                                               
      BLOCK DATA PTPACM                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                              
C                                                                               
C         REAL*8 NAB,NAB1,NBC,NBC1,JUNK                                         
      REAL*8 JUNK
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
         COMMON /AMBPCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                         
         COMMON /PRECM/  DBCOOR(2),JUNK(4),BCOOR                                
         COMMON /PRE2CM/ D1INF,R1EINF,B1INF,D2INF,                              
     +                   R2INF,B2INF,A1INF,A2INF,                               
     +                   R1S,R2S,VSP,ALPH1,ALPH2,                               
     +                   V(4),PHI(20),BSP(20,5),CSP(22,5),                      
     +                   D1DIF,R1DIF,B1DIF,A1DIF,                               
     +                   D2DIF,R2DIF,B2DIF,A2DIF,                               
     +                   STEP,D1,D2,R1E,R2E,B1,B2,A1,A2,                        
     +                   IAYES,NPHI,NV                                          
C                                                                               
C   Initialize the flags and the I/O unit numbers for the potential             
C                                                                               
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
         DATA NFLAG /1,1,15*0,6,0,0/                                            
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,0,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
C   Dissociation energy in kcal/mol                                             
C                                                                               
         DATA DE / 106.56D0 /                                                   
C                                                                               
C   Morse beta in reciprocal Angstroms                                          
C                                                                               
         DATA BETA / 2.07942D0/                                                 
C                                                                               
C   Equilibrium bond lengths in Angstroms                                       
C                                                                               
         DATA REQ / 0.96966D0, 0.74144D0, 0.96966D0/                            
C                                                                               
C   Pauling parameters in Angstroms                                             
C                                                                               
         DATA A / 0.26D0/                                                       
C                                                                               
C   Gamma parameters for the bend correction                                    
C                                                                               
         DATA GAMMA /  0.4131D0 /                                               
C                                                                               
C   RBB in Angstroms - this number is used in                                   
C   a scaling calculation for the bending correction                            
C                                                                               
         DATA RBB / 0.74144D0 /                                                 
C                                                                               
         DATA D1INF / 97.23777d0 /                                              
         DATA R1EINF /  .75933d0 /                                              
         DATA B1INF /  2.01627d0 /                                              
         DATA D2INF / 98.81423d0 /                                              
         DATA R2INF /   .98712d0 /                                              
         DATA B2INF /  2.30153d0 /                                              
         DATA A1INF /   0.0d0    /                                              
         DATA A2INF /  -0.1d0    /                                              
         DATA R1S   / 2.21208d0  /                                              
         DATA R2S   / 2.21531d0  /                                              
         DATA VSP   / 87.65999d0 /                                              
         DATA IAYES /1/                                                         
         DATA ALPH1 /  1.0d0     /                                              
         DATA ALPH2 /  1.0d0     /                                              
         DATA STEP  /0.000000010d0 /                                            
         DATA NPHI /11/                                                         
         DATA NV   /3/                                                          
         DATA (PHI(I), I=1,11)                                                  
     +        / 0.0000000d+00,  0.1000000d+02,                                  
     +          0.2000000d+02,  0.3000000d+02,                                  
     +          0.3732999d+02,  0.4000000d+02,                                  
     +          0.5000000d+02,  0.6000000d+02,                                  
     +          0.7000000d+02,  0.8000000d+02,                                  
     +          0.9000000d+02/                                                  
         DATA (BSP(I,1), I=1,11)                                                
     +        / 0.9723777d+02,  0.9597664d+02,                                  
     +          0.9292151d+02,  0.9047624d+02,                                  
     +          0.9036362d+02,  0.9132085d+02,                                  
     +          0.9458319d+02,  0.9728180d+02,                                  
     +          0.9809436d+02,  0.9854700d+02,                                  
     +          0.9881423d+02/                                                  
         DATA (BSP(I,2), I=1,11)                                                
     +        / 0.1455982d+01,  0.1474749d+01,                                  
     +          0.1529910d+01,  0.1605620d+01,                                  
     +          0.1633813d+01,  0.1631247d+01,                                  
     +          0.1533492d+01,  0.1395812d+01,                                  
     +          0.1298821d+01,  0.1242818d+01,                                  
     +          0.1224964d+01/                                                  
         DATA (BSP(I,3), I=1,11)                                                
     +        / 0.2016270d+01,  0.1955010d+01,                                  
     +          0.1879873d+01,  0.1619415d+01,                                  
     +          0.1484336d+01,  0.1467517d+01,                                  
     +          0.1689022d+01,  0.1935001d+01,                                  
     +          0.2141747d+01,  0.2261585d+01,                                  
     +          0.2301529d+01/                                                  
         DATA (CSP(I,1), I=1,13)                                                
     +        / 0.9723777d+02, -0.7558417d-01,                                  
     +         -0.5052923d-03,  0.1237764d-02,                                  
     +         -0.7532287d-03,  0.2249910d-02,                                  
     +         -0.1043168d-01,  0.8624236d-02,                                  
     +         -0.9957934d-03,  0.1126430d-02,                                  
     +         -0.6615264d-03,  0.1681010d-03,                                  
     +         -0.5891778d-04/                                                  
         DATA (CSP(I,2), I=1,13)                                                
     +        / 0.1455982d+01,  0.1101632d-02,                                  
     +          0.7749966d-05, -0.1010505d-04,                                  
     +         -0.1182017d-04,  0.3430826d-05,                                  
     +         -0.4369002d-04,  0.8096568d-04,                                  
     +         -0.1197233d-04, -0.1871181d-04,                                  
     +          0.6506310d-05, -0.1045387d-04,                                  
     +          0.8100479d-05/                                                  
         DATA (CSP(I,3), I=1,13)                                                
     +        / 0.2016270d+01, -0.7200198d-02,                                  
     +          0.1074172d-04, -0.7832711d-04,                                  
     +          0.1557412d-03, -0.8529179d-04,                                  
     +          0.3140934d-03, -0.4585322d-03,                                  
     +          0.1656866d-03, -0.4269112d-04,                                  
     +          0.2110969d-04,  0.1294235d-04,                                  
     +         -0.1547269d-04/                                                  
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
