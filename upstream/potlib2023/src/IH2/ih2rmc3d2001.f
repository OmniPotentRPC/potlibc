      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          IH2                                                        
C   Functional form:                                                            
C   Common name:     RMC                                                        
C   References:                                                                 
C                                                                               
C   Note: This is the three-dimensional version of the IH2                      
C         RMC potential energy surface; the potential energy is                 
C         defined as E = Ecollinear + Ebend.                                    
C                                                                               
C
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 1
c
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in the block data subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the IH diatomic, with the             
C   IH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(I-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-I)                                           
C   The classical potential energy is set equal to zero for the I               
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
         DOUBLE PRECISION NAB,NAB1,NBC,NBC1,JUNK                                
C                                                                               
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr                   
C                                                                               
         PARAMETER (CKCAL = 627.5095D0)                                         
         PARAMETER (CANGS = 0.52917706D0)                                       
C                                                                               
         COMMON /RMCRCM/ RR(2),S(2),D(2),ALPHA(2),RE(2),ALFD,CLBH               
         COMMON /RMCICM/ IRXN,L(2)                                              
         COMMON /BENDCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                         
C                                                                               
         COMMON /PRECM/  DBCOOR(2),JUNK(4)                                      
         COMMON /PRE2CM/ C(2,2),PRE(2),RRMRE(2),                                
     +                   A2(2),DADP(2),DPDR(2),RRMR(2),                         
     +                   PID4,PID2,DMDN,RR1MS1,RR2MS2,                          
     +                   RRS1D2,SR1R2                                           
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C      COLLINEAR POTENTIAL PARAMETERS                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C                                                                               
C   Initialize the parameters for the collinear part of the potential           
C                                                                               
      CALL PREPOT2                                                              
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C      BENDING POTENTIAL PARAMETERS                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C                                                                               
C   Echo the parameters for the bend potential to the                           
C   file linked to UNIT NFLAG(18)                                               
C                                                                               
         WRITE (NFLAG(18), 650) REQ,DE,BETA,A,GAMMA,RBB                         
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
C                                                                               
      EZERO(1)=D(2)                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
C                                                                               
      INDEXES(1) = 53                                                           
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
C        MOVE POTENTIAL SUBROIUTINE                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the IH diatomic, with the             
C   IH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(I-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-I)                                           
C   The classical potential energy is set equal to zero for the I               
C   infinitely far from the equilibrium H2 diatomic.                            
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C      ENTRY POT                                                                
CCCCCCCCCCCCCCCC                                                                
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
         DOUBLE PRECISION NAB,NAB1,NBC,NBC1,JUNK                                
                                                                                
C                                                                               
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr                   
C                                                                               
         PARAMETER (CKCAL = 627.5095D0)                                         
         PARAMETER (CANGS = 0.52917706D0)                                       
C                                                                               
         COMMON /RMCRCM/ RR(2),S(2),D(2),ALPHA(2),RE(2),ALFD,CLBH               
         COMMON /RMCICM/ IRXN, L(2)                                             
         COMMON /BENDCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                         
C                                                                               
         COMMON /PRECM/  DBCOOR(2),JUNK(4)                                      
         COMMON /PRE2CM/ C(2,2),PRE(2),RRMRE(2),                                
     +                   A2(2),DADP(2),DPDR(2),RRMR(2),                         
     +                   PID4,PID2,DMDN,RR1MS1,RR2MS2,                          
     +                   RRS1D2,SR1R2                                           
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
      C2 = 1.D0 - NAB/NBC                                                       
      IF (ABS(C2).GT.1.D-14) GOTO 90                                            
      NAB1=.5D0                                                                 
      NBC1=.5D0                                                                 
      GO TO 95                                                                  
   90 C2 = C2/NAB                                                               
      NAB1 = (2.D0 + C2 - SQRT(4.D0+C2*C2))/(2.D0*C2)                           
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
C          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).                       
C                                                                               
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
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
C          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR                  
C          GEOMETRY AND CALL POT                                                
C                                                                               
      RACTMP = R(3)                                                             
      R(3) = R(2) + R(1)                                                        
C                                                                               
      CALL POT2                                                                 
C                                                                               
C   Reset R(3) the its correct value                                            
C                                                                               
      R(3) = RACTMP                                                             
C                                                                               
C   Add the energy terms together, i.e., E = Ecollinear + Ebend                 
C                                                                               
      ENGYGS = ENGYGS + BCOOR                                                   
C                                                                               
      IF (NDER .NE. 1) THEN                                                     
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
          RETURN                                                                
      END IF                                                                    
C                                                                               
C   If NDER = 1, then compute the derivatives of the energy with respect to     
C   the internal coordinates, i.e., dE/dRi = dEcollinear/dRi + dEbend/dRi       
C                                                                               
      DEGSDR(1) = DEGSDR(1) + DEGSDR(3)                                         
      DEGSDR(2) = DEGSDR(2) + DEGSDR(3)                                         
      DEGSDR(3) = DAC                                                           
C                                                                               
      DO 300 IT = 1,2                                                           
           DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT)                                 
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
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/,2X,T5,'POT has been called with NDER = ',I5,                     
     *       /,2X,T5,'This value of NDER is not allowed in this ',              
     *               'version of the potential.')                               
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
C*****                                                                          
C                                                                               
      SUBROUTINE PREPOT2                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      DOUBLE PRECISION JUNK
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
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr                   
C                                                                               
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (CANGS = 0.52917706D0)                                          
C                                                                               
      COMMON /RMCRCM/ RR(2),S(2),D(2),ALPHA(2),RE(2),ALFD,CLBH                  
      COMMON /RMCICM/ IRXN, L(2)                                                
C                                                                               
C   This subprogram initializes the constants                                   
C   for the collinear part of the potential.                                    
C   PREPOT2 must be called once before any calls to POT2.                       
C                                                                               
      COMMON /PRECM/  DBCOOR(2),JUNK(4)                                         
      COMMON /PRE2CM/ C(2,2),PRE(2),RRMRE(2),                                   
     +                A2(2),DADP(2),DPDR(2),RRMR(2),                            
     +                PID4,PID2,DMDN,RR1MS1,RR2MS2,                             
     +                RRS1D2,SR1R2                                              
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
     +         2.0D0*DBLE(L(N))*(1.0D0-DMDN)                                    
      C(N,1) = 1.0D0-DMDN*(1.0D0-C(M,1))+PID4*PID4*C(N,2)                       
C                                                                               
C   Echo the potential parameters to the file linked to UNIT NFLAG(18).         
C                                                                               
         WRITE (NFLAG(18), 600) D(1), D(2), D(1), RE(1), RE(2), RE(1),          
     *                     ALPHA(1), ALPHA(2), ALPHA(1)                         
         WRITE (NFLAG(18), 610) RR, S, L, ALFD                                  
         WRITE (NFLAG(18), 620) ((C(I,J),I=1,2),J=1,2)                          
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
C   CONVERT TO ATOMIC UNITS                                                     
C                                                                               
      DO   2 I = 1,2                                                            
           RR(I) = RR(I)/CANGS                                                  
           S(I) = S(I)/CANGS                                                    
           RE(I) = RE(I)/CANGS                                                  
           ALPHA(I) = ALPHA(I)*CANGS                                            
           D(I)=D(I)/CKCAL                                                      
 2    CONTINUE                                                                  
      ALFD = ALFD*CANGS                                                         
C                                                                               
C   CALCULATE USEFUL CONSTANTS                                                  
C                                                                               
      RR1MS1=RR(1)-S(1)                                                         
      RR2MS2=RR(2)-S(2)                                                         
      RRS1D2=RR1MS1/RR2MS2                                                      
      SR1R2=SQRT(RR1MS1**2+RR2MS2**2)                                           
      DO   5 N=1,2                                                              
           RRMRE(N)=RR(N)-RE(N)                                                 
           PRE(N)=SR1R2-RRMRE(N)                                                
 5    CONTINUE                                                                  
C                                                                               
C   Initialize the asymptotic energy terms                                      
C                                                                               
         EASYAB = D(1)                                                          
         EASYBC = D(2)                                                          
         EASYAC = D(1)                                                          
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT2                                                           
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      DOUBLE PRECISION JUNK
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
C      ENTRY POT2                                                               
C                                                                               
C   Conversion constants: kcal/mol/Hartree and Angstroms/Bohr                   
C                                                                               
      PARAMETER (CKCAL = 627.5095D0)                                            
      PARAMETER (CANGS = 0.52917706D0)                                          
C                                                                               
      COMMON /RMCRCM/ RR(2),S(2),D(2),ALPHA(2),RE(2),ALFD,CLBH                  
      COMMON /RMCICM/ IRXN, L(2)                                                
      COMMON /PRECM/  DBCOOR(2),JUNK(4)                                         
      COMMON /PRE2CM/ C(2,2),PRE(2),RRMRE(2),                                   
     +                A2(2),DADP(2),DPDR(2),RRMR(2),                            
     +                PID4,PID2,DMDN,RR1MS1,RR2MS2,                             
     +                RRS1D2,SR1R2                                              
C                                                                               
C   This subprogram initializes the constants                                   
C   for the collinear part of the potential.                                    
C   PREPT2 must be called once before any calls to POT2.                        
C                                                                               
      LOGICAL LRR1,LRR2                                                         
C                                                                               
C   This subprogram computes the collinear part of the potential energy         
C                                                                               
      RRMR(1) = RR(1)-R(1)                                                      
      RRMR(2) = RR(2)-R(2)                                                      
      LRR1 = RRMR(1).LE.0.0D0                                                   
      LRR2 = RRMR(2).LE.0.0D0                                                   
      IF (NDER .EQ. 1) THEN                                                     
          DO 10 J = 1,3                                                         
10              DEGSDR(J) = 0.0D0                                               
      ENDIF                                                                     
      ENGYGS = EZERO(1)                                                         
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
      A2(I) = C(I,1) + C(I,2)*(PHIP-PID2)*PHIP                                  
      DADP(I) = C(I,2)*(2.0D0*PHIP-PID2)                                        
      RR12 = SQRT (RRMR(1)**2+RRMR(2)**2)                                       
      LM3 = L(I) - 3                                                            
      SN2PPL = SN2PP3*SN2PP**LM3                                                
      ALF = ALPHA(I) +(ALFD-ALPHA(I))*SN2PP4                                    
      GAM = D(I)*(1.0D0-A2(I)*SN2PPL)                                           
      XI = PRE(I)*SN2PP4+RRMRE(I)-RR12                                          
      QUAN = 1.0D0 - EXP(-ALF*XI)                                               
      QUAN2 = QUAN*QUAN                                                         
      ENGYGS = GAM*(QUAN2-1.0D0)+ EZERO(1)                                      
C                                                                               
      IF (NDER .NE. 1) RETURN                                                   
C                                                                               
C   If NDER=1, compute the derivatives of the energy with respect to the        
C   with respect to the internal coordinates.                                   
C                                                                               
      FAC = RRS1D2/(RRMR(1)**2+(RRS1D2*RRMR(2))**2)                             
      DPDR(1) = FAC*RRMR(2)                                                     
      DPDR(2) = -FAC*RRMR(1)                                                    
      DO 40 J = 1,2                                                             
      DS4DR = 8.0D0*SN2PP4/TN2PP*DPDR(J)                                        
      DSLDR = 2.0D0*DBLE(L(I))*SN2PPL/TN2PP*DPDR(J)                             
      DALDR = (ALFD-ALPHA(I))*DS4DR                                             
      DGADR = -D(I)*A2(I)*DSLDR - D(I)*DADP(I)*SN2PPL*DPDR(J)                   
      DXIDR = PRE(I)*DS4DR + RRMR(J)/RR12                                       
   40 DEGSDR(J) = DGADR*(QUAN2-1.0D0)+2.0D0*GAM*QUAN*(1.0D0-QUAN)*              
     1          (ALF*DXIDR + XI*DALDR)                                          
      RETURN                                                                    
   50 J = 1                                                                     
      IF (LRR1) J = 2                                                           
      QUAN = 1.0D0-EXP(-ALPHA(J)*(R(J)-RE(J)))                                  
      ENGYGS = D(J)*(QUAN**2-1.0D0)+ EZERO(1)                                   
C                                                                               
C   If NDER =1, compute the derivatives of the energy                           
C   with respect to the internal coordinates.                                   
C                                                                               
      IF (NDER .EQ. 1) DEGSDR(J) = 2.0D0*D(J)*ALPHA(J)*QUAN*(1.0D0-QUAN)        
C                                                                               
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
         COMMON /RMCRCM/ RR(2), S(2), D(2), ALPHA(2), RE(2), ALFD, CLBH         
         COMMON /RMCICM/ IRXN, L(2)                                             
         COMMON /BENDCM/ DE,BETA,REQ(3),A,GAMMA,RBB,RHH                         
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
C   Initialize the potential parameters; the energy parameters are in           
C   kcal/mol, and the lengths are in Angstroms.                                 
C                                                                               
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
C                                                                               
         DATA  IRXN / 1/                                                        
C                                                                               
C   Initialize the parameters for the bend potential                            
C                                                                               
         DATA DE    / 73.682D0/                                                 
         DATA BETA  / 1.789D0/                                                  
         DATA REQ   / 1.604D0, 0.742D0, 1.604D0/                                
         DATA A     / 0.26D0/                                                   
         DATA GAMMA / 0.50D0/                                                   
         DATA RBB   / 0.7412871D0/                                              
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
