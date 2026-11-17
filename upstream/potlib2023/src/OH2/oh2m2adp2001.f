C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          OH2  A'' surface                                           
C   Common name:     M2                                                         
C   Functional form: Rotated-Morse-oscillator-spline (RMOS) potential           
C                    for collinear geometeries plus anti-Morse                  
C                    bending potential.                                         
C   References:      B. C. Garrett and D. G. Truhlar                            
C                    Int. J. Quantum Chem. 29, 1463 (1986)                      
C                                                                               
C   Interface:       potlib2001
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   PREPOT must be called once before any calls to POT.                         
C   The control flags for the potential are contained in the block data         
C   subprogram PTPACM.                                                          
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in hartree atomic units.                                                 
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
C                                                                               
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
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
         DOUBLE PRECISION NAB,NAB1,NBC,NBC1                                     
C                                                                               
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)                           
C                                                                               
         PARAMETER (BOHR = .52917706D0)                                         
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
         PARAMETER (HARTRE = 627.5095D0)                                        
C                                                                               
C  Vaule of tolerance                                                           
C                                                                               
         PARAMETER (TOL = 1.D-05)                                               
C                                                                               
         COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                  
         COMMON /PRECM/  DBCOOR(2)                                              
C                                                                               
C   Echo the name of the potential energy surface                               
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 600)                                                 
C                                                                               
C   Assign the parameters for the collinear part of the potential               
C                                                                               
      CALL PREPOT2                                                              
C                                                                               
C   Echo the parameters for the bending potential                               
C                                                                               
      WRITE (NFLAG(18), 620) REQ,D                                              
      WRITE (NFLAG(18), 630) BETA,A,RBB,GAMP                                    
C                                                                               
C       CONVERT TO ATOMIC UNITS                                                 
C                                                                               
      DO 50 IT = 1,3                                                            
         REQ(IT) = REQ(IT)/BOHR                                                 
         D(IT) = D(IT)/HARTRE                                                   
   50 CONTINUE                                                                  
C                                                                               
      BETA = BETA * BOHR                                                        
      DE = D(1)                                                                 
      A = A/BOHR                                                                
      RHH = 0.74144D0 / BOHR                                                    
      RBB = RBB /BOHR                                                           
C                                                                               
C   Initialize the energy in the asymptotic valleys                             
C                                                                               
      EASYAB = D(1)                                                             
      EASYBC = D(2)                                                             
      EASYAC = D(3)                                                             
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the A '''' OH2 ',             
     *                'potential energy surface M2',                            
     *       //,2X,T5,'Potential parameters for the collinear ',                
     *                'part of the potential')                                  
620   FORMAT(/,2X,T5,'Parameters for the bending correction',                   
     *       /,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',                
     *       /,2X,T5,'Equilibrium bond lengths (Angstroms):',                   
     *       T44, F10.5, T55, F10.5, T66, F10.5,                                
     *       /,2X,T5,'Dissociation energy (kcal/mol):',                         
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
630   FORMAT(2X,T5,'Morse beta parameter (Angstroms**-1):', T44,1PE15.8,        
     *     /,2X,T5,'Pauling parameter (Angstroms):',T44,1PE15.8,                
     *     /,2X,T5,'RBB (Angstroms):',T44,1PE15.8,                              
     *     /,2X,T5,'Gamma for the bending correction:',T44,1PE15.8,             
     *          T60,1PE15.8,/,2X,T44,1PE15.8)                                   
C                                                                               
C                                                                               
      EZERO(1)=D(2)                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='B. C. Garrett, D. G. Truhlar'                                    
       REF(2)='Int. J. Quantum Chem. 29, 1463(1986)'                            
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
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(O-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-O)                                           
C   The classical potential energy is set equal to zero for the O               
C   infinitely far from the equilibrium H2 diatomic.                            
C                                                                               
C      ENTRY POT                                                                
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
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
         DOUBLE PRECISION NAB,NAB1,NBC,NBC1                                     
C                                                                               
C   Conversion for Angstroms to Bohr (Angstroms/Bohr)                           
C                                                                               
         PARAMETER (BOHR = .52917706D0)                                         
C                                                                               
C   Conversion for kcal to Hartree (kcal/Hartree)                               
C                                                                               
         PARAMETER (HARTRE = 627.5095D0)                                        
C                                                                               
C  Vaule of tolerance                                                           
C                                                                               
         PARAMETER (TOL = 1.D-05)                                               
C                                                                               
         COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                  
         COMMON /PRECM/  DBCOOR(2)                                              
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
      IF (R(1) .GT. TOL .AND. R(2) .GT. TOL .AND. R(3) .GT. TOL) THEN           
         RACCOL = R(1) + R(2)                                                   
         BCOOR = 0.0D0                                                          
         IF (NDER .EQ. 1) THEN                                                  
             DBCOOR(1) = 0.0D0                                                  
             DBCOOR(2) = 0.0D0                                                  
             DAC = 0.0D0                                                        
         ENDIF                                                                  
         ROB = 1.0D0                                                            
         IF(ABS(RHH-RBB).GT.1.D-10) THEN                                        
C                                                                               
C     FIRST CALCULATE SOME USEFUL NUMBERS                                       
C     NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS                          
C     AB AND BC RESPECTIVELY. (EQUATION 2.2)                                    
C                                                                               
            NAB = EXP((REQ(1) - R(1))/A)                                        
            IF (NAB.LT.1.D-15) NAB = 1.D-15                                     
            NBC = EXP((REQ(2) - R(2))/A)                                        
            IF (NBC.LT.1.D-15) NBC = 1.D-15                                     
C                                                                               
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE                     
C       SEE EQUATION 2.7A, 2.7B                                                 
C                                                                               
            C = 1.D0 - NAB/NBC                                                  
            IF (ABS(C) .LT. 1.D-14) THEN                                        
               NAB1=.5D0                                                        
               NBC1=.5D0                                                        
            ELSE                                                                
               C = C/NAB                                                        
               NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)                      
               NBC1 = 1.D0 - NAB1                                               
            END IF                                                              
C                                                                               
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS                    
C          FIRST TRAP OUT ANY ZERO ARGUMENTS                                    
C                                                                               
            IF (NAB1*NBC1 .LE. 0.0D0) THEN                                      
               ROB = 1.D0                                                       
            ELSE                                                                
               STUFF = A * LOG(NAB1*NBC1)                                       
               T = RHH - STUFF                                                  
               ROB = 1.0D0                                                      
               IF(ABS(T).GT.1.D-15) ROB = (RBB - STUFF)/T                       
            END IF                                                              
         END IF                                                                 
C                                                                               
C        CALCULATE BENDING CORRECTION                                           
C        EVALUATE GAMMA AND DERIVATIVE                                          
C                                                                               
         EX = EXP(RACCOL*(GAMP(2) + GAMP(3)*RACCOL))                            
         GAMMA = GAMP(1) + EX                                                   
         IF (NDER .EQ. 1) DGAM = (GAMP(2) + 2.D0*GAMP(3)*RACCOL)*EX             
         GAMMA = GAMMA*DE                                                       
         IF (NDER .EQ. 1) DGAM = DGAM*DE                                        
C                                                                               
C           CALCULATE V(R(3)) TERM                                              
C                                                                               
         EX = EXP(-BETA*(R(3)-REQ(3))/ROB)                                      
         EX1 = EX*(1.D0 + 0.5D0*EX)                                             
         IF (NDER .EQ. 1) DEX1 = EX*(1.D0 + EX)                                 
C                                                                               
C           CALCULATE V(R(1)+R(2)) TERM                                         
C                                                                               
         EX = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)                                  
         EX2 = EX*(1.D0 + 0.5D0*EX)                                             
         IF (NDER .EQ. 1) DEX2 = EX*(1.D0 + EX)                                 
C                                                                               
C        HERE IS THE BENDING CORRECTION                                         
C                                                                               
         BCOOR = GAMMA*(EX1 - EX2)                                              
C                                                                               
C        WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE                        
C        OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).                          
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             GAMMA = GAMMA*BETA/ROB                                             
             DAC = -GAMMA*DEX1                                                  
C                                                                               
C          CALCULATE CORRECTIONS TO DAB, DBC                                    
C                                                                               
             DBCOOR(1) = GAMMA*DEX2 + DGAM*(EX1-EX2)                            
             DBCOOR(2)= DBCOOR(1)                                               
         ENDIF                                                                  
C                                                                               
C       NOW CALCULATE THE COLLINEAR ENERGY                                      
C       STORE R(3) IN RACTMP THEN SET R(3) TO                                   
C       THE COLLINEAR GEOMETRY AND CALL POT2                                    
C                                                                               
         RACTMP = R(3)                                                          
         R(3) = R(2) + R(1)                                                     
C                                                                               
         CALL POT2                                                              
C                                                                               
         R(3) = RACTMP                                                          
C                                                                               
         ENGYGS = ENGYGS + BCOOR                                                
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             DO 300 IT = 1,2                                                    
                    DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT) + DEGSDR(3)            
300          CONTINUE                                                           
             DEGSDR(3) = DAC                                                    
         ENDIF                                                                  
      ELSE                                                                      
        ENGYGS = 1.D30                                                          
        IF (NDER .EQ. 1) THEN                                                   
           DEGSDR(1) = 0.0D0                                                    
           IF (R(1) .LE. TOL) DEGSDR(1) = -ENGYGS                               
           DEGSDR(2) = 0.0D0                                                    
           IF (R(1) .LE. TOL) DEGSDR(2) = -ENGYGS                               
           DEGSDR(3) = 0.0D0                                                    
           IF (R(1) .LE. TOL) DEGSDR(3) = -ENGYGS                               
        END IF                                                                  
      END IF                                                                    
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,                     
     *       /,2X,T5,'This value of NDER is not allowed in this ',              
     *               'version of the potential.')                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT2                                                        
C                                                                               
C  COLLINEAR A + BC POTENTIAL FUNCTION                                          
C  RMO-SPLINE TYPE FIT                                                          
C  COORDINATES AND ENERGIES PASSED TO POTENTIAL IN ATOMIC UNITS                 
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
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
      COMMON /PRE2CM/ D1INF,R1INF,B1INF,D2INF,                                  
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                PHI(20),BSP(20,3),CSP(22,3),                              
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                DASY,IAYES,NPHI                                           
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
C  TITLE                                                                        
C                                                                               
      TITLE=' O+H2 MAB 9-12-79 '                                                
C                                                                               
C  ASSIGN COLLINEAR PARAMETERS (UNITS OF KCAL, ANGSTROM, DEGREES)               
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
      A1INF = 0.0D0                                                             
      A2INF = 0.0D0                                                             
      IAYES = 0                                                                 
C                                                                               
C  ASSIGN PARAMETERS FOR THE ASYMPTOTIC REACTANT AND PRODUCT REGIONS            
C  IN BLOCK DATA SUBROUTINE                                                     
C                                                                               
      IF(A2INF .EQ.0.0D0) A2INF = D1INF-D2INF                                   
C                                                                               
C  ASSIGN THE PIVOT POINT, POTENTIAL AT THE PIVOT,                              
C  FLAG, AND ASYMPTOTIC EXPONENTIAL PARAMETERS                                  
C  IN BLOCK DATA SUBROUTINE.  ECHO.                                             
C                                                                               
      WRITE (NFLAG(18), 600) TITLE                                              
      WRITE (NFLAG(18), 601) D1INF,R1INF,B1INF,D2INF,R2INF,B2INF,A1INF,A        
C     2INF                                                                      
      WRITE (NFLAG(18), 602) R1S,R2S,VSP,ALPH1,ALPH2                            
C                                                                               
      IF(IAYES.NE.0) WRITE (NFLAG(18), 603)                                     
      IF(IAYES.EQ.0) WRITE (NFLAG(18), 604)                                     
C                                                                               
      NPHIP2 = NPHI + 2                                                         
C                                                                               
C  ECHO POTENTIAL PARAMTERS ASSIGNED BY BLOCK DATA SUBROUTINE                   
C                                                                               
      WRITE (NFLAG(18), 605) NPHI                                               
      WRITE (NFLAG(18), 610)                                                    
C                                                                               
      DO 30 I = 1,NPHI                                                          
         WRITE (NFLAG(18), 606) I,PHI(I),(BSP(I,J),CSP(I,J),J=1,3)              
   30 CONTINUE                                                                  
C                                                                               
      WRITE (NFLAG(18), 611)                                                    
      DO 40 I=1,2                                                               
         K = I + NPHI                                                           
         WRITE (NFLAG(18), 607) K,(CSP(K,J),J=1,3)                              
   40 CONTINUE                                                                  
C                                                                               
C  CONVERT TO ATOMIC UNITS                                                      
C                                                                               
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
C                                                                               
C  CONVERT FROM DEGREES TO RADIANS                                              
C                                                                               
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
C                                                                               
C     REACTANT ASYMPTOTE                                                        
C                                                                               
      DIS = BSP(1,1)                                                            
      REQ2 = BSP(1,2)                                                           
      BET = BSP(1,3)                                                            
      T1 = 1.0D0 - EXP(-BET*REQ2)                                               
      AER = VSP - DIS*T1*T1                                                     
      REQ2 = R2S - REQ2                                                         
      D1DIF = DIS - D1INF                                                       
      R1DIF = REQ2 - R1INF                                                      
      B1DIF = BET - B1INF                                                       
      A1DIF = AER - A1INF                                                       
C                                                                               
C     PRODUCT ASYMPTOTE                                                         
C                                                                               
      DIS = BSP(NPHI,1)                                                         
      REQ2 = BSP(NPHI,2)                                                        
      BET = BSP(NPHI,3)                                                         
      T1 = 1.0D0 - EXP(-BET*REQ2)                                               
      AER = VSP - DIS*T1*T1                                                     
      REQ2 = R1S - REQ2                                                         
      D2DIF = DIS - D2INF                                                       
      R2DIF = REQ2 - R2INF                                                      
      B2DIF = BET - B2INF                                                       
      A2DIF = AER - A2INF                                                       
C                                                                               
C  DEFINE ZERO OF ENERGY                                                        
C                                                                               
      DASY = D1INF                                                              
C                                                                               
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
  800 FORMAT (2x,//,A20,//)                                                     
  801 FORMAT (8F9.5)                                                            
  802 FORMAT (3F10.5,I5,2F10.5)                                                 
  803 FORMAT (16I5)                                                             
  804 FORMAT (4G20.10)                                                          
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT2                                                           
C                                                                               
C      ENTRY POT2                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
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
      COMMON /PRE2CM/ D1INF,R1INF,B1INF,D2INF,                                  
     +                R2INF,B2INF,A1INF,A2INF,                                  
     +                R1S,R2S,VSP,ALPH1,ALPH2,                                  
     +                PHI(20),BSP(20,3),CSP(22,3),                              
     +                D1DIF,R1DIF,B1DIF,A1DIF,                                  
     +                D2DIF,R2DIF,B2DIF,A2DIF,                                  
     +                DASY,IAYES,NPHI                                           
C                                                                               
      CHARACTER*20 TITLE                                                        
C                                                                               
      R1 = R(1)                                                                 
      R2 = R(2)                                                                 
      R3 = R(3)                                                                 
      IF (R1 .GT. R1S) THEN                                                     
         IF (R2 .GT. R2S) THEN                                                  
C                                                                               
C   3-ATOM BREAK UP REGION                                                      
C                                                                               
           ENGYGS = DASY                                                        
            IF (NDER .EQ. 1) THEN                                               
                DEGSDR(1) = 0.0D0                                               
                DEGSDR(2) = 0.0D0                                               
            ENDIF                                                               
         ELSE                                                                   
C                                                                               
C   REACTANT ASYMPTOTIC REGION                                                  
C                                                                               
            DR1 = R1 - R1S                                                      
            IF (DR1*ALPH1 .LT. -70.0D0) WRITE (NFLAG(18), 6601) R1, R2,         
C     R3,                                                                       
     *         DR1, ALPH1                                                       
            T1 = EXP(-DR1*ALPH1)                                                
            T2 = -ALPH1*T1                                                      
C                                                                               
C.....INTERPOLATE FOR MORSE PARAMETERS                                          
C                                                                               
            DIS = D1INF + D1DIF*T1                                              
            DDIS = D1DIF*T2                                                     
            REQ2 = R1INF + R1DIF*T1                                             
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
            DELR = R2 - REQ2                                                    
            IF (DELR*BET .LT.  -70.0D0)                                         
     +          WRITE (NFLAG(18), 6602) R1, R2, R3, DELR, BET                   
            T1 = EXP(-BET*DELR)                                                 
            T2 = 1.0D0 - T1                                                     
           ENGYGS = AER + DIS*T2*T2                                             
            IF (NDER .EQ. 1) THEN                                               
               T3 = 2.0D0*DIS*T1                                                
               DEGSDR(1) =                                                      
     +                 DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))          
               DEGSDR(2) = T3*BET*T2                                            
            ENDIF                                                               
         END IF                                                                 
      ELSE IF (R2 .GT. R2S) THEN                                                
C                                                                               
C   PRODUCT ASYMPTOTIC REGION                                                   
C                                                                               
         DR2 = R2 - R2S                                                         
         IF (DR2*ALPH2 .LT.  -70.D0)                                            
     *      WRITE (NFLAG(18), 6603) R1, R2, R3, DR2,ALPH2                       
C                                                                               
C.....INTERPOLATE FOR MORSE PARAMETERS                                          
C                                                                               
         T1 = EXP(-DR2*ALPH2)                                                   
         T2 = -ALPH2*T1                                                         
         DIS = D2INF + D2DIF*T1                                                 
         DDIS = D2DIF*T2                                                        
         REQ2 = R2INF + R2DIF*T1                                                
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
         DELR = R1 - REQ2                                                       
         IF (DELR*BET .LT.  -70.0D0)                                            
     *       WRITE (NFLAG(18), 6604) R1, R2, R3, DELR, BET                      
         T1 = EXP(-BET*DELR)                                                    
         T2 = 1.0D0 - T1                                                        
        ENGYGS = AER + DIS*T2*T2                                                
         IF (NDER .EQ. 1) THEN                                                  
             T3 = 2.0D0*DIS*T1                                                  
             DEGSDR(2) = DAER + T2*(DDIS*T2 - T3*(BET*DREQ - DBET*DELR))        
C                                                                               
             DEGSDR(1) = T3*BET*T2                                              
         ENDIF                                                                  
      ELSE                                                                      
C                                                                               
C   SPLINE INTERPOLATION REGION                                                 
C                                                                               
         DR1 = R1S - R1                                                         
         DR2 = R2S - R2                                                         
         T1 = DR1*DR1 + DR2*DR2                                                 
         RR = SQRT(T1)                                                          
         DRRD1 = -DR1/RR                                                        
         DRRD2 = -DR2/RR                                                        
         ANGLE = ATAN(DR1/DR2)                                                  
         DPHD1 = -DR2/T1                                                        
         DPHD2 = DR1/T1                                                         
         DIS = SPLINE(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                 
         DDIS = SPLINE(4,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                
         REQ2 = SPLINE(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                
         DREQ = SPLINE(4,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                
         BET =  SPLINE(3,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                
         DBET = SPLINE(4,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                
         IF(IAYES.EQ.0) THEN                                                    
            AER = VSP - DIS                                                     
            DAER = -DDIS                                                        
         ELSE                                                                   
            T1 = EXP(-BET*REQ2)                                                 
            T2 = 1.0D0 - T1                                                     
            AER = VSP - DIS*T2*T2                                               
            DAER = -T2*(DDIS*T2 + 2.0D0*DIS*T1*(BET*DREQ + REQ2*DBET))          
         END IF                                                                 
         IF (BET*(RR-REQ2) .GT.  70)                                            
     *       WRITE (NFLAG(18), 6600) R1,R2,R3,RR,ANGLE,DIS,REQ2,BET             
         DELR = RR - REQ2                                                       
         T1 = EXP(BET*DELR)                                                     
         T2 = 1.0D0 - T1                                                        
        ENGYGS = AER + DIS*T2*T2                                                
        IF(NDER.EQ.1) THEN                                                      
            T3 = 2.0D0*DIS*T1                                                   
            DVDP = DAER + T2*(DDIS*T2 - T3*(DBET*DELR - BET*DREQ))              
            T3 = T3*T2*BET                                                      
            DEGSDR(1) =  DVDP*DPHD1 - T3*DRRD1                                  
            DEGSDR(2) =  DVDP*DPHD2 - T3*DRRD2                                  
        ENDIF                                                                   
      END IF                                                                    
      IF(NDER.EQ.1) THEN                                                        
         DEGSDR(1) = DEGSDR(1)                                                  
         DEGSDR(2) = DEGSDR(2)                                                  
         DEGSDR(3) = 0.0D0                                                      
      ENDIF                                                                     
C                                                                               
6600  FORMAT(/,2X,T5,'R1, R2, R3, RR, ANGLE, DIS, REQ2, BET:',                  
     *       /,(2X,T5,1PE20.10,T30,1PE20.10,T55,1PE20.10))                      
6601  FORMAT(/,2X,T5,'R1, R2, R3, DELR, ALPH1:',/,2X,T10,5(1PE13.5,1X))         
6602  FORMAT(/,2X,T5,'R1, R2, R3, DR2, BET:',/,2X,T10,5(1PE13.5,1X))            
6603  FORMAT(/,2X,T5,'R1, R2, R3, DELR, ALPH2:',/,2X,T10,5(1PE13.5,1X))         
6604  FORMAT(/,2X,T5,'R1, R2, R3, DR2, BET:',/,2X,T10,5(1PE13.5,1X))            
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      FUNCTION SPLINE(ISW,NN,X,Y,C,D,XPT)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
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
C     PARAMETER (IPRT = 6)                                                      
      DIMENSION X(NN),Y(NN),C(NN+2),D(NN+2)                                     
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
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.             
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.                        
C                                                                               
C     Output from this file is written to UNIT NFLAG(18) which is               
C     defined in a parameter statement.                                         
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
               WRITE(NFLAG(18),1001)                                            
               STOP 'SPLINE 1'                                                  
            END IF                                                              
            PIVOT = 1.0D0/D(I-1)                                                
            IF(I.LT.NP2) THEN                                                   
               SUPD = X(I-1)-X(I-2)                                             
               IF(SUPD.LT.0.0D0) THEN                                           
                  WRITE(NFLAG(18),1000)                                         
                  STOP 'SPLINE 2'                                               
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
               WRITE(NFLAG(18),1001)                                            
               STOP 'SPLINE 3'                                                  
            END IF                                                              
            C(J+1) = C(J+1)/D(J+1)                                              
            C(J) = C(J)-C(J+1)*V                                                
411      CONTINUE                                                               
         IF(D(2).EQ.0.0D0) THEN                                                 
            WRITE(NFLAG(18),1001)                                               
            STOP 'SPLINE 4'                                                     
         END IF                                                                 
         C(2) = C(2)/D(2)                                                       
      ELSE IF (ISW .EQ. 2) THEN                                                 
         SPLINE = C(1) + C(2)*(Z-X(1))                                          
         DO 51 I = 1,N                                                          
            V = Z-X(I)                                                          
            IF(V.LE.0) RETURN                                                   
            SPLINE = SPLINE + C(I+2)*V**3                                       
51       CONTINUE                                                               
      ELSE IF (ISW .EQ. 3) THEN                                                 
         SPLINE = C(1) + C(2)*(Z-X(1))                                          
         DERIV = C(2)/3.D0                                                      
         DO 53 I = 1,N                                                          
            V = Z-X(I)                                                          
            IF(V.LE.0) RETURN                                                   
            V2 = V*V                                                            
            SPLINE = SPLINE + C(I+2)*V2*V                                       
            DERIV = DERIV + C(I+2)*V2                                           
   53    CONTINUE                                                               
      ELSE IF (ISW .EQ. 4) THEN                                                 
         SPLINE = 3.D0*DERIV                                                    
      END IF                                                                    
      RETURN                                                                    
1000  FORMAT(1X,5X,'***** ERROR IN SPLINE ... UNORDERED X-VALUES                
     1*****')                                                                   
1001  FORMAT(1X,5X,' ***** ERROR IN SPLINE ... DIVIDE FAULT *****')             
      END                                                                       
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
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
         COMMON /AMBPCM/ D(3),BETA,REQ(3),A,GAMP(3),RBB,RHH,DE                  
         COMMON /PRE2CM/ D1INF,R1INF,B1INF,D2INF,                               
     +                   R2INF,B2INF,A1INF,A2INF,                               
     +                   R1S,R2S,VSP,ALPH1,ALPH2,                               
     +                   PHI(20),BSP(20,3),CSP(22,3),                           
     +                   D1DIF,R1DIF,B1DIF,A1DIF,                               
     +                   D2DIF,R2DIF,B2DIF,A2DIF,                               
     +                   DASY,IAYES,NPHI                                        
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
C--- bend potential parameters for A'' surface of O + H2                        
C   Assign the parameters for the non-collinear part of the potential           
C                                                                               
C   Dissociation energy in kcal/mol                                             
C                                                                               
         DATA D /106.56D0, 109.472D0, 106.56D0/                                 
C                                                                               
C   Morse beta in reciprocal Angstroms                                          
C                                                                               
         DATA BETA /2.07942D0/                                                  
C                                                                               
C   Equilibrium bond lengths in Angstroms                                       
C                                                                               
         DATA REQ /0.96966D0, 0.74144D0, 0.96966D0/                             
C                                                                               
C   Pauling parameters in Angstroms                                             
C                                                                               
         DATA A /0.26D0/                                                        
C                                                                               
C   Gamma parameters for the bend correction                                    
C                                                                               
         DATA GAMP /0.210820D0, 5.42027D0, -1.44235D0/                          
C                                                                               
C   RBB in Angstroms - this number is used in                                   
C   a scaling calculation for the bending correction                            
C                                                                               
         DATA RBB /0.74144D0/                                                   
C                                                                               
         DATA D1INF /  95.3d0    /                                              
         DATA R1INF /    .76d0   /                                              
         DATA B1INF /   2.0202d0 /                                              
         DATA D2INF /  92.4d0    /                                              
         DATA R2INF /    .98d0   /                                              
         DATA B2INF /   2.41105d0/                                              
         DATA A1INF /   0.0d0    /                                              
         DATA A2INF /   0.0d0    /                                              
         DATA R1S   /  2.21208d0 /                                              
         DATA R2S   /  2.21531d0 /                                              
         DATA VSP   / 95.3d0     /                                              
         DATA IAYES /  0         /                                              
         DATA ALPH1 /  1.6d0     /                                              
         DATA ALPH2 /  1.6d0     /                                              
         DATA NPHI / 11 /                                                       
         DATA (PHI(I), I=1,11)                                                  
     +        / 0.0000000d+00,  0.1000000d+02,                                  
     +          0.2000000d+02,  0.3000000d+02,                                  
     +          0.3733000d+02,  0.4000000d+02,                                  
     +          0.5000000d+02,  0.6000000d+02,                                  
     +          0.7000000d+02,  0.8000000d+02,                                  
     +          0.9000000d+02/                                                  
         DATA (BSP(I,1), I=1,11)                                                
     +        / 0.9477920d+02,  0.9316200d+02,                                  
     +          0.9035320d+02,  0.8517960d+02,                                  
     +          0.8272720d+02,  0.8304320d+02,                                  
     +          0.8855580d+02,  0.9003100d+02,                                  
     +          0.9129800d+02,  0.9184900d+02,                                  
     +          0.9209100d+02/                                                  
         DATA (BSP(I,2), I=1,11)                                                
     +        / 0.1456500d+01,  0.1475500d+01,                                  
     +          0.1530500d+01,  0.1606400d+01,                                  
     +          0.1634900d+01,  0.1632300d+01,                                  
     +          0.1534600d+01,  0.1396600d+01,                                  
     +          0.1299800d+01,  0.1243800d+01,                                  
     +          0.1225600d+01/                                                  
         DATA (BSP(I,3), I=1,11)                                                
     +        / 0.2041900d+01,  0.1987500d+01,                                  
     +          0.1905400d+01,  0.1672700d+01,                                  
     +          0.1554000d+01,  0.1539400d+01,                                  
     +          0.1749300d+01,  0.2014500d+01,                                  
     +          0.2227100d+01,  0.2348800d+01,                                  
     +          0.2383800d+01/                                                  
         DATA (CSP(I,1), I=1,13)                                                
     +        / 0.9477920d+02,  -0.1478277d+00,                                 
     +         -0.1389235d-03,  -0.3580591d-03,                                 
     +          0.1450636d-02,   0.1148717d-02,                                 
     +         -0.9458414d-03,  -0.4983968d-02,                                 
     +          0.5921762d-02,  -0.2814975d-02,                                 
     +          0.1001140d-02,  -0.2747850d-03,                                 
     +         -0.5702495d-05/                                                  
         DATA (CSP(I,2), I=1,13)                                                
     +        / 0.1456500d+01,   0.1137130d-02,                                 
     +          0.7628700d-05,  -0.9772201d-05,                                 
     +         -0.1201120d-04,   0.2800151d-05,                                 
     +         -0.4114549d-04,   0.7866376d-04,                                 
     +         -0.1122307d-04,  -0.1936703d-04,                                 
     +          0.6791209d-05,  -0.1039780d-04,                                 
     +          0.8032967d-05/                                                  
         DATA (CSP(I,3), I=1,13)                                                
     +        / 0.2041900d+01,  -0.5887558d-02,                                 
     +          0.4475579d-05,  -0.5455348d-04,                                 
     +          0.1230139d-03,  -0.6997501d-04,                                 
     +          0.2822753d-03,  -0.4019424d-03,                                 
     +          0.1215636d-03,  -0.1548133d-04,                                 
     +          0.9961687d-05,   0.1813458d-04,                                 
     +         -0.1747243d-04/                                                  
         END                                                                    
C                                                                               
C*****                                                                          
                                                                                
