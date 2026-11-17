C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          OH2                                                        
C   Functional form: LEPS (London-Eyring-Polanyi-Sato) plus a correction term   
C   Common name:     JWS                                                        
C   References:      G. C. Schatz                                               
C                    J. Chem. Phys. 83, 5677 (1985)                             
C                    B. R. Johnson and N. W. Winter                             
C                    J. Chem. Phys. 66, 4116 (1977).                            
C                                                                               
C   Interface:       potlib2001
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
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
C      coordinates  - bohrs                                                     
C      derivatives - hartrees/bohr                                              
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the O            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in the BLOCK DATA subprogram PTPACM                                  
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(O-first H)                                
C                            R(2) = R(first H-second H)                         
C                            R(3) = R(second H-O)                               
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, energy, and derivatives of                     
C        the energy with respect to the coordinates.                            
C        passes the control flags where                                         
C        NDER  = 0 => no derivatives are calculated                             
C        NDER  = 1 => calculate first derivatives of the energy for the         
C                     ground electronic state with respect to the coordinates   
C        NFLAG  Control flags                                                   
C      NFLAG(18-20)                                                             
C        passes the FORTRAN unit number used for potential output               
C      /ASYCM/ EASYAB, EASYBC, EASYAC                                           
C        passes the energy in the three asymptotic valleys for an A+BC system.  
C        The energy in the AB valley, EASYAB, is equal to the energy            
C        of the C atom "infinitely" far from the AB diatomic and R(AB) set      
C        equal to Re(AB), the equilibrium bond length for the AB diatomic.      
C        In this potential the AB valley represents H infinitely far from       
C        the OH diatomic and R(OH) equal to Re(OH).                             
C        Similarly, the terms EASYBC and EASYAC represent the energies in the   
C        H2 and the other OH asymptotic valleys, respectively.                  
C                                                                               
C   Default Parameter Values:                                                   
C      Variable      Default value                                              
C      NDER             1                                                       
C      NFLAG(18)        6                                                       
C                                                                               
C*****                                                                          
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
         PARAMETER (R2     = 1.41421356D0)                                      
         PARAMETER (CKCAU  = 627.5095D0)                                        
         PARAMETER (CANGAU = 0.52917706D0)                                      
C                                                                               
         COMMON /SATOCM/ D(3),RE(3),BETA(3),Z(3)                                
C                                                                               
         COMMON /PRECM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),                
     +                   DO4Z(3),B(3),X(3),COUL(3),EXCH(3)                      
C                                                                               
C   Echo the potential parameters to the file linked to UNIT NFLAG(18)          
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 600) D, RE, BETA, Z                                  
C                                                                               
600   FORMAT(/,1X,'*****',1X,'Potential Energy Surface',1X,'*****',/,           
     *       /,1X,T5,'OH2 JWS potential energy surface',                        
     *       //,1X,T5,'Parameters:',                                            
     *        /,1X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',               
     *        /,1X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Sato parameters:',                                       
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')                 
C                                                                               
C   CONVERT TO ATOMIC UNITS                                                     
C                                                                               
      DO  10 I = 1,3                                                            
          D(I)=D(I)/CKCAU                                                       
          RE(I) = RE(I)/CANGAU                                                  
          BETA(I) = BETA(I)*CANGAU                                              
C                                                                               
C   COMPUTE USEFUL CONSTANTS                                                    
C                                                                               
          ZPO(I) = 1.0D0 + Z(I)                                                 
          OP3Z(I) = 1.0D0 + 3.0D0*Z(I)                                          
          TOP3Z(I) = 2.0D0*OP3Z(I)                                              
          ZP3(I) = Z(I) + 3.0D0                                                 
          TZP3(I) = 2.0D0*ZP3(I)                                                
          DO4Z(I) = D(I)/4.0D0/ZPO(I)                                           
          B(I) = BETA(I)*DO4Z(I)*2.0D0                                          
 10   CONTINUE                                                                  
C                                                                               
C    Initialize the asymptotic energy values                                    
C                                                                               
         EASYAB = D(1)                                                          
         EASYBC = D(2)                                                          
         EASYAC = D(3)                                                          
C                                                                               
C                                                                               
      EZERO(1)=D(2)                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='B. R. Johnson, N. W. Winter'                                     
       REF(2)='J. Chem. Phys. 66, 4116(1977)'                                   
       REF(3)='G. C. Schatz'                                                    
       REF(4)='J. Chem. Phys. 83, 5677(1985)'                                   
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
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the O            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(O-first H)                                
C                            R(2) = R(first H-second H)                         
C                            R(3) = R(second H-O)                               
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
         PARAMETER (R2     = 1.41421356D0)                                      
         PARAMETER (CKCAU  = 627.5095D0)                                        
         PARAMETER (CANGAU = 0.52917706D0)                                      
C                                                                               
         COMMON /SATOCM/ D(3),RE(3),BETA(3),Z(3)                                
C                                                                               
         COMMON /PRECM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),                
     +                   DO4Z(3),B(3),X(3),COUL(3),EXCH(3)                      
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Check the value of NDER                                                     
C                                                                               
         IF (NDER .GT. 1) THEN                                                  
             WRITE (NFLAG(18), 900) NDER                                        
             STOP 'POT 1'                                                       
         ENDIF                                                                  
C                                                                               
C   Compute the LEPS potential - this is equivalent to the JW OH2 potential.    
C                                                                               
      DO 21 I = 1,3                                                             
            X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))                                
            COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)                       
            EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)                       
21    CONTINUE                                                                  
      RAD = SQRT((EXCH(1)-EXCH(2))**2 + (EXCH(2)-EXCH(3))**2 +                  
     +      (EXCH(3)-EXCH(1))**2)                                               
      ENGYGS = COUL(1) + COUL(2) + COUL(3) - (RAD/R2) +                         
     +         EZERO(1) - ANUZERO                                               
C                                                                               
C   Compute the derivatives of the LEPS energy w.r.t.the coordinates            
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          SVAL = EXCH(1) + EXCH(2) + EXCH(3)                                    
          DO 22 I = 1,3                                                         
                DEGSDR(I)=0.0D0                                                 
                IF(X(I).LT.1.0D-30) GO TO 22                                    
                TVAL = (3.0D0*EXCH(I)-SVAL)/R2*(OP3Z(I)*X(I)-ZP3(I))            
                IF(ABS(RAD) .LT. 1.0D-32 .AND.                                  
     +             ABS(TVAL) .GT. 1.0D-12) THEN                                 
                   WRITE(NFLAG(18), 6000) TVAL, RAD                             
                ELSE IF(ABS(RAD).GT.1.0D-32) THEN                               
                   TVAL = TVAL/RAD                                              
                END IF                                                          
                DEGSDR(I)=B(I)*X(I)*(TVAL-ZP3(I)*X(I)+OP3Z(I))                  
22       CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
C  George Schatz's correction term to JW LEPS surface                           
C                                                                               
      T1 = R(1) - 2.105D0                                                       
      T2 = R(3) - 2.105D0                                                       
      T3 = R(2) - 4.21D0                                                        
      VCOR = 0.04382D0*EXP(-T1*T1 - T2*T2 - 0.5D0*T3*T3)                        
      ENGYGS  = ENGYGS + VCOR                                                   
         IF (NDER .EQ. 1) THEN                                                  
             DEGSDR(1) = DEGSDR(1) - 2.D0*T1*VCOR                               
             DEGSDR(3) = DEGSDR(3) - 2.D0*T2*VCOR                               
             DEGSDR(2) = DEGSDR(2) - T3*VCOR                                    
         ENDIF                                                                  
C                                                                               
600   FORMAT(/,1X,'*****',1X,'Potential Energy Surface',1X,'*****',/,           
     *       /,1X,T5,'OH2 JWS potential energy surface',                        
     *       //,1X,T5,'Parameters:',                                            
     *        /,1X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',               
     *        /,1X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,1X,T5,'Sato parameters:',                                       
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')                 
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,             
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',             
     *                'coded in this potential')                                
6000  FORMAT(/,2X,'In the OH2 potential JWS, T, RAD = ',1P,2E15.7,              
     1       /,2X,'T/RAD has been set equal to T')                              
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
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3)                             
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
C   Initialize the potential parameters; the energy parameters                  
C   are in kcal/mol, and the lengths are in Angstroms.                          
C                                                                               
         DATA D    / 106.6D0, 109.4D0, 106.6D0/                                 
         DATA RE   / 0.9706D0, 0.7417D0, 0.9706D0/                              
         DATA BETA / 2.294D0, 1.942D0, 2.294D0/                                 
         DATA Z    / 0.0885D0, 0.0885D0, 0.0885D0/                              
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
