C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          FH2                                                        
C   Functional form: eLEPS (extended-London-Eyring-Polanyi-Sato)                
C   Common name:     no. 3                                                      
C   Reference:       R. Steckler, D. G. Truhlar, B. C. Garrett,                 
C                    N. C. Blais, and R. B. Walker,                             
C                    J. Chem. Phys. 81, 5700 (1984).                            
C                                                                               
C   Notes:           In this version of the FH2 potential the AB and AC         
C                    Sato parameters are functions of PHI, where PHI is         
C                    the angle of A-to-BC vector with the BC axis.              
C                                                                               
C
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 0
C
C   THE DERIVATIVES OF THE SATO PARATMETERS WITH RESPECT TO PHI HAVE            
C   NOT BEEN CODED FOR THIS POTENTIAL!                                          
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in DATA statements                    
C   and in the block data subprogram PTPACM.                                    
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys is                     
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the FH diatomic, with the             
C   FH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(F-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-F)                                                
C   The zero of energy is defined at F "infinitely" far from the H2 diatomic.   
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states availalble                            
C                 (1,1) = 1 as only gs state available                          
C        NDER  = 0 => no derivatives should be calculated                       
C        NDER  = 1 => calculate first derivatives                               
C        NFLAG - these integer values can be used to flag options               
C                within the potential.                                          
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
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
         PARAMETER (CKCAU  = 627.5095D0)                                        
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (R2     = 1.41421356D0)                                      
C                                                                               
         COMMON /COSCM/ COSPHI                                                  
         COMMON /PHICM/ DEG1,Z0,FDER,DEG                                        
         COMMON /Z3CM/  Z(3)                                                    
         COMMON /PRECM/  D(3),RE(3),BETA(3),ZPO(3),OP3Z(3),                     
     +                   ZP3(3),TZP3(3),TOP3Z(3),DO4Z(3),B(3),                  
     +                   X(3),COUL(3),EXCH(3)                                   
C                                                                               
C   Echo the potential parameters to UNIT NFLAG(18)                             
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 600) D, RE, BETA, Z                                  
         WRITE (NFLAG(18), 610) DEG1                                            
C                                                                               
      Z0=Z(1)                                                                   
      DO  10 I = 1,3                                                            
C                                                                               
C   CONVERT TO ATOMIC UNITS                                                     
C                                                                               
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
C   10 B(I) = BETA(I)*DO4Z(I)*2.0D0                                             
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = D(1)                                                      
             EASYBC = D(2)                                                      
             EASYAC = D(3)                                                      
C                                                                               
600   FORMAT (/, 2X, T5, 'PREPOT has been called for the FH2 ',                 
     *                   'extended-LEPS potential energy surface',              
     *        /, 2X, T5, 'called no. 3',                                        
     *        /, 2X, T5, 'Potential energy surface parameters:',                
     *        /, 2X, T5, 'Bond', T47, 'F-H', T58, 'H-H', T69, 'H-F',            
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):',                   
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):',               
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):',              
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Sato parameters:',                                    
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
610   FORMAT (/,2X,T5,'The above Sato parameters apply only for PHI ',          
     *                'less than',F10.3,1X,'degrees',/)                         
C                                                                               
C                                                                               
      EZERO(1)=D(2)                                                             
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='R. Steckler, D. G. Truhlar, B. C. Garrett,'                      
       REF(2)='N. C. Blais, R. B. Walker,'                                      
       REF(3)='J. Chem. Phys. 81, 5700(1984)'                                   
C                                                                               
      INDEXES(1) = 9                                                            
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
C   energy of the H "infinitely" far from the FH diatomic, with the             
C   FH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the FH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(F-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-F)                                                
C   The zero of energy is defined at F "infinitely" far from the H2 diatomic.   
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
         PARAMETER (CKCAU  = 627.5095D0)                                        
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (R2     = 1.41421356D0)                                      
C                                                                               
         COMMON /COSCM/ COSPHI                                                  
         COMMON /PHICM/ DEG1,Z0,FDER,DEG                                        
         COMMON /Z3CM/  Z(3)                                                    
         COMMON /PRECM/  D(3),RE(3),BETA(3),ZPO(3),OP3Z(3),                     
     +                   ZP3(3),TZP3(3),TOP3Z(3),DO4Z(3),B(3),                  
     +                   X(3),COUL(3),EXCH(3)                                   
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
C                                                                               
C   Check the value of NDER                                                     
C                                                                               
         IF (NDER .GT. 1) THEN                                                  
             WRITE (NFLAG(18), 910) NDER                                        
             STOP 'POT 2'                                                       
         ENDIF                                                                  
C                                                                               
      CALL PHICAL (PHI)                                                         
C      IF (LCHI.NE.0) WRITE (NFLAG(18),620) R(1),R(2),R(3),PHI                  
620   FORMAT(2X,3F8.3,T35,F10.4)                                                
      Z(1)=FPHI(PHI)                                                            
      Z(3)=Z(1)                                                                 
      DO 11 I=1,3,2                                                             
      ZPO(I) = 1.0D0 + Z(I)                                                     
      OP3Z(I) = 1.0D0 + 3.0D0*Z(I)                                              
      TOP3Z(I) = 2.0D0*OP3Z(I)                                                  
      ZP3(I) = Z(I) + 3.0D0                                                     
      TZP3(I) = 2.0D0*ZP3(I)                                                    
   11 DO4Z(I) = D(I)/4.0D0/ZPO(I)                                               
      S = 0.0D0                                                                 
      DO 21 I = 1,3                                                             
      X(I) = EXP(-BETA(I)*(R(I)-RE(I)))                                         
      COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)                             
      EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)                             
   21 S = S+EXCH(I)                                                             
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+                     
     1      (EXCH(3)-EXCH(1))**2)                                               
      ENGYGS = -RAD/R2                                                          
      IF (NDER .EQ. 1) THEN                                                     
          DO 22 I = 1,3                                                         
                DEGSDR(I) = 0.D0                                                
                IF(X(I).LT.1.D-30) GO TO 22                                     
                T= (3.0D0*EXCH(I)-S)/R2*(OP3Z(I)*X(I)-ZP3(I))                   
                IF(RAD.NE.0.D0) GO TO 2222                                      
                IF(T.EQ.0.D0) GO TO 2223                                        
                WRITE(NFLAG(18),6000) T,RAD                                     
                GO TO 2223                                                      
 2222           T = T/RAD                                                       
 2223           CONTINUE                                                        
                DEGSDR(I) = B(I)*X(I)*(T-ZP3(I)*X(I)+OP3Z(I))                   
22        CONTINUE                                                              
      ENDIF                                                                     
      ENGYGS = ENGYGS + COUL(1) + COUL(2) + COUL(3) + EZERO(1) - ANUZERO        
C                                                                               
      IF (DEG.LE.DEG1) THEN                                                     
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
          RETURN                                                                
      END IF                                                                    
C                                                                               
C   DERIVATIVE CHANGES FROM LEPS TO LEPS2 ARE TO BE ADDED HERE.                 
C   THEY ARE COMPUTED USING FDER FROM FPHI AND COSPHI FROM COSCOM.              
C   FDER IS DERIVATIVE OF SATO PARAMETER WITH RESPECT TO PHI.                   
C                                                                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
6000  FORMAT(/,2X,'In the eLEPS potential T, RAD = ',1P,2E15.7,                 
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
      FUNCTION FPHI(PHI)                                                        
C                                                                               
C   FH2 NEW SURFACE NO. 3   DON TRUHLAR                                         
C   AB,AC SATO PARAMETERS ARE Z(1)=Z(3) FOR DEG.LE.DEG1.                        
C   DEG1 MUST BE SET IN PREPOT.                                                 
C   Z0 IS SET EQUAL TO Z(1)=Z(3) IN PREPOT.                                     
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
      COMMON /PHICM/ DEG1,Z0,FDER,DEG                                           
      COMMON /Z3CM/ Z(3)                                                        
      DEG=57.29577951D0*PHI                                                     
C      IF (LCHI.NE.0) WRITE (NFLAG(18),100) DEG                                 
  100 FORMAT (2X,T35,F10.4)                                                     
      IF (DEG.GT.180.0D0)DEG=DEG-180.0D0                                        
      IF (DEG.LT.0.0D0)  DEG=-DEG                                               
      IF (DEG.GT.90.0D0) DEG=180.0D0-DEG                                        
      IF (DEG.LT.0.0D0) STOP 89                                                 
      IF (DEG.GT.DEG1) GO TO 200                                                
C      IF (LCHI.NE.0) WRITE (NFLAG(18),100) DEG                                 
      FPHI=Z0                                                                   
      FDER=0.0D0                                                                
      RETURN                                                                    
  200 IF (DEG.GT.20.0D0) GO TO 300                                              
      DEGDIF=DEG-10.0D0                                                         
      DD2=DEGDIF*DEGDIF                                                         
      DD3 = DEGDIF*DD2                                                          
      FPHI=0.175D0 +1.15D-4*DD2    -6.0D-6*DD3                                  
      FDER=         2.30D-4*DEGDIF -1.8D-5*DD2                                  
      RETURN                                                                    
  300 IF (DEG.GT.60.0D0) GO TO 400                                              
      DEGDIF=DEG-41.0D0                                                         
      FPHI=0.191D0 +0.0005D0*DEGDIF                                             
      FDER=         0.0005D0                                                    
      RETURN                                                                    
  400 IF (DEG.GT.90.0D0) GO TO 500                                              
      SINDEG=SIN(PHI)                                                           
      S2=SINDEG*SINDEG                                                          
      FPHI=0.17569D0 +3.308D-2*S2                                               
      COSDEG=COS(PHI)                                                           
      FDER=           6.616D-2*SINDEG*COSDEG                                    
      RETURN                                                                    
  500 STOP 90                                                                   
      END                                                                       
C                                                                               
      SUBROUTINE PHICAL (PHI)                                                   
C                                                                               
C   FH2 NEW SURFACE NO. 2   DON TRUHLAR                                         
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
      COMMON /COSCM/ COSPHI                                                     
      HR2=0.5D0*R(2)                                                            
      R12=R(1)*R(1)                                                             
      R22=R(2)*R(2)                                                             
      R32=R(3)*R(3)                                                             
      HR22=0.25D0*R22                                                           
      COS23=(R12-R22-R32)/(-2.0D0*R(2)*R(3))                                    
      FHH2=HR22+R32-2.0D0*HR2*R(3)*COS23                                        
      FHH=SQRT(FHH2)                                                            
      COSPHI=(R32-HR22-FHH2)/(-2.0D0*HR2*FHH)                                   
      IF (COSPHI.GT.1.0D0) COSPHI=1.0D0                                         
      IF (COSPHI.LT.-1.0D0) COSPHI=-1.0D0                                       
      PHI=ACOS(COSPHI)                                                          
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
         COMMON /PHICM/ DEG1,Z0,FDER,DEG                                        
         COMMON /Z3CM/  Z(3)                                                    
         COMMON /PRECM/  D(3),RE(3),BETA(3),ZPO(3),OP3Z(3),                     
     +                   ZP3(3),TZP3(3),TOP3Z(3),DO4Z(3),B(3),                  
     +                   X(3),COUL(3),EXCH(3)                                   
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
C   Initialize potential variables                                              
C                                                                               
         DATA DEG1 /10.0D0/                                                     
         DATA Z /0.175D0, 0.168D0, 0.175D0/                                     
C                                                                               
         DATA D /141.196D0,109.449D0,141.196D0/                                 
         DATA RE /0.917D0, 0.7419D0, 0.917D0/                                   
         DATA BETA /2.21870D0, 1.942D0, 2.2187D0/                               
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
