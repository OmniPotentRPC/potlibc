      Subroutine prepot                                                         
C                                                                               
C   System:          OH2                                                        
C   Functional form:                                                            
C   Common name:     OH2SL                                                      
C   Reference:                                                                  
C                                                                               
C   Interface:       potlib2001
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   Notes: This is a collinear potential energy surface, i.e., the potential    
C          energy is not defined for non-linear geometries.                     
C          The derivatives of the potential energy with respect to the          
C          coordinates are determined numerically using a step size STEP        
C          if R(1) and R(2) are less than 50 bohr.                              
C          The variable STEP is initialized in the block data subprogram        
C          PTPACM.                                                              
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are assigned in the Block Data Subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
C   through the common block PT1CM:                                             
C                  COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)                       
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C                                                                               
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
C                                                                               
C   The flags that indicate which calculations should be carried out in         
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
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C        IMPLICIT REAL*16(L)                                                    
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
C                                                                               
C      DIMENSION C(50),AL(3),X(3),S(3),VM(3),R(3),CHI(3)                        
C                                                                               
C   Echo potential data to the file linked to UNIT NFLAG(18)                    
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      WRITE (NFLAG(18),600)                                                     
      WRITE (NFLAG(18),601) STEP                                                
      WRITE (NFLAG(18),602)                                                     
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',                    
     *                'Schinke-Lester potential energy surface')                
601   FORMAT(/,2X,T5,'Step size for the numerical derivatives = ',E10.5)        
602   FORMAT(2X,T5,'The Schinke-Lester surface has analytical ',                
     *             'derivatives if R(1) and R(2) are ',                         
     *       /,2X,T5,'less than 50 bohr, otherwise the derivatives ',           
     *               'are computed numerically.')                               
C                                                                               
      J=1                                                                       
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = DE(1)/HART                                                
             EASYBC = DE(2)/HART                                                
             EASYAC = DE(3)/HART                                                
C                                                                               
C                                                                               
      EZERO(1)=DE(2)                                                            
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
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C     CALCULATE ENERGY AND DERIVATIVE                                           
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C                                                                               
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
C                                                                               
C      ENTRY POT                                                                
C                                                                               
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C                                                                               
C        IMPLICIT REAL*16(L)                                                    
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
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
      IC=4                                                                      
      IF(R(1).GT.30.D0)IC=1                                                     
      IF(R(2).GT.30.D0)IC=1                                                     
C                                                                               
   86 GO TO (84,82,83,81),IC                                                    
C                                                                               
 81   SOH=R(1)                                                                  
      SR2=R(2)                                                                  
      GO TO 85                                                                  
C                                                                               
 82   SR2=R(2)+STEP                                                             
      SOH=R(1)                                                                  
      GO TO 85                                                                  
C                                                                               
 83   SR2=R(2)                                                                  
      SOH=R(1)+STEP                                                             
      GO TO 85                                                                  
C                                                                               
   84 SR2=R(2)                                                                  
      SOH=R(1)                                                                  
      R(3)=R(3) + STEP                                                          
C                                                                               
   85 SR2=SR2*0.529177D0                                                        
      SOH=SOH*0.529177D0                                                        
      R(3)= R(3)*0.529177D0                                                     
      R(1) = SOH                                                                
      R(2) = SR2                                                                
C                                                                               
      Q = 1.0D0                                                                 
C                                                                               
      DO 10 I=1,3                                                               
         X(I)=R(I)                                                              
         CHI(I)=(1.0D0 - EXP(-ALF*X(I)**2))                                     
         S(I)=X(I)-1.D0                                                         
         VM(I)=DE(I)*(1.0D0-EXP(-BET(I)*(X(I)-XRE(I))))**2-DE(I)                
         Q=Q*(1.0D0-TANH(AL(I)*S(I)))                                           
  10  continue                                                                  
C                                                                               
      U = VM(1)*CHI(2)*CHI(3)+VM(2)*CHI(1)*CHI(3)                               
     *	+ VM(3)*CHI(1)*CHI(2)                                                    
C                                                                               
      BP = B*(1.0D0/X(1)**4+1.0D0/X(3)**4)                                      
C                                                                               
      CHECK = (X(1)*X(3)*(X(2)-0.2D0))                                          
C                                                                               
      IF(CHECK.GT.1000.D0) AP = 0.0D0                                           
      IF(CHECK.GT.1000.D0) GO TO 112                                            
C                                                                               
      AP = A/(X(1)*X(3)*(X(2)-0.2D0))**12                                       
C                                                                               
  112 U = U + AP + BP                                                           
C                                                                               
      T = C(1) + C(2)*(S(1)+S(3)) + C(3)*S(2) +                                 
     *    C(4)*(S(1)**2 + S(3)**2) +                                            
     *	  C(5)*S(2)**2 + C(6)*S(2)*(S(1) + S(3)) +                               
     *    C(7)*S(1)*S(3) + C(8)*(S(1)**3 + S(3)**3) +                           
     *    C(9)*S(2)**3 + C(10)*S(1)*S(3)*(S(1) + S(3))                          
C                                                                               
      TT = C(11)*S(2)**2*(S(1) + S(3)) +                                        
     *     C(12)*S(2)*(S(1)**2 + S(3)**2) +                                     
     *     C(13)*S(1)*S(2)*S(3) + C(14)*(S(1)**4 + S(3)**4) +                   
     *     C(15)*S(2)**4 + C(16)*(S(1)**3*S(3) + S(1)*S(3)**3)                  
C                                                                               
      WW = C(17)*S(1)**2*S(3)**2 + C(18)*S(2)*(S(1)**3 + S(3)**3) +             
     *     C(19)*S(2)**2*(S(1)**2 + S(3)**2) +                                  
     *     C(20)*S(2)**3*(S(1) + S(3)) + C(21)*S(2)**2*S(1)*S(3) +              
     *     C(22)*S(2)*(S(1)**2*S(3) + S(1)*S(3)**2)                             
C                                                                               
      W = C(23)*(S(1)**5 + S(3)**5) +                                           
     *    C(24)*(S(1)**4*S(3) + S(1)*S(3)**4) +                                 
     *    C(25)*(S(1)**3*S(3)**2 + S(1)**2*S(3)**3) +                           
     *    C(26)*S(2)*(S(1)**4 + S(3)**4) +                                      
     *    C(27)*S(2)*(S(1)**3*S(3) + S(1)*S(3)**3)                              
C                                                                               
      GG = C(28)*S(2)*S(1)**2*S(3)**2 +                                         
     *     C(29)*S(2)**2*(S(1)**3 + S(3)**3) +                                  
     *     C(30)*S(2)**2*(S(1)**2*S(3) + S(1)*S(3)**2) +                        
     *     C(31)*S(2)**3*(S(1)**2 + S(3)**2)                                    
C                                                                               
      Y = C(32)*S(2)**3*S(1)*S(3) + C(33)*S(2)**4*(S(1) + S(3)) +               
     *    C(34)*S(2)**5 + C(35)*(S(1)**6 + S(3)**6) +                           
     *    C(36)*(S(1)**5*S(3) + S(1)*S(3)**5) +                                 
     *    C(37)*(S(1)**4*S(3)**2 + S(1)**2*S(3)**4) +                           
     *    C(38)*S(1)**3*S(3)**3                                                 
C                                                                               
      Z = C(39)*S(2)*(S(1)**5 + S(3)**5) +                                      
     *    C(40)*S(2)*(S(1)**4*S(3) + S(1)*S(3)**4) +                            
     *    C(41)*S(2)*(S(1)**3*S(3)**2 + S(1)**2*S(3)**3) +                      
     *    C(42)*S(2)**2*(S(1)**4 + S(3)**4) +                                   
     *    C(43)*S(2)**2*(S(1)**3*S(3) + S(1)*S(3)**3) +                         
     *    C(44)*S(1)**2*S(2)**2*S(3)**2                                         
C                                                                               
      ZZ = C(45)*S(2)**3*(S(1)**3 + S(3)**3) +                                  
     *     C(46)*S(2)**3*(S(1)**2*S(3) + S(1)*S(3)**2) +                        
     *     C(47)*S(2)**4*(S(1)**2 + S(3)**2) +                                  
     *     C(48)*S(2)**4*S(1)*S(3) +                                            
     *     C(49)*S(2)**5*(S(1) + S(3)) + C(50)*S(2)**6                          
C                                                                               
      P = T + TT + WW + W + GG + Y + Z + ZZ                                     
C                                                                               
      TEMP = EXP(-GAM*(X(1) + X(3))**2)                                         
C                                                                               
      U = U + P*Q*EXP(-GAM*(X(1) + X(3))**2)                                    
C                                                                               
      XTER=U + EZERO(1)                                                         
      V(IC) = XTER/HART                                                         
C                                                                               
      IF(IC .LT. 4)GO TO 97                                                     
      IF(R(1).GT.30.D0)GO TO 77                                                 
      IF(R(2).LT.30.D0)GO TO 66                                                 
C                                                                               
 77   CONTINUE                                                                  
      IF(NDER.EQ.1) THEN                                                        
         DEGSDR(1)=(V(3) - V(4))/STEP                                           
         DEGSDR(2)=(V(2) - V(4))/STEP                                           
         DEGSDR(3)=(V(1) - V(4))/STEP                                           
      ENDIF                                                                     
C                                                                               
 98   IF(R(2).LT.20)GO TO 97                                                    
C                                                                               
      V(4)=V(4)*1.0D8                                                           
      V(4)=INT(V(4))                                                            
      V(4)=V(4)*1.0D-8                                                          
C                                                                               
   97 IC=IC+1                                                                   
C                                                                               
      IF(IC.LT.5) GO TO 86                                                      
C                                                                               
      ENGYGS=V(4)                                                               
C                                                                               
      GO TO 113                                                                 
C                                                                               
 66   CONTINUE                                                                  
C                                                                               
C BEGINNING OF NDER SWITCH AROUND DERIVATIVES                                   
C                                                                               
      IF(NDER.EQ.1) THEN                                                        
C                                                                               
      DA = C(2) + 2.0D0*S(1)*C(4) + S(2)*C(6) + S(3)*C(7) +                     
     *     3.0D0*S(1)**2*C(8) + (S(3)**2 + 2.0D0*S(1)*S(3))*C(10) +             
     *     S(2)**2*C(11) + 2.0D0*S(1)*S(2)*C(12) + S(2)*S(3)*C(13) +            
     *     4.0D0*S(1)**3*C(14)                                                  
C                                                                               
      DB = (3.0D0*S(1)**2*S(3) + S(3)**3)*C(16) +                               
     *      2.0D0*S(1)*S(3)**2*C(17) + 3.0D0*S(1)**2*S(2)*C(18) +               
     *      2.0D0*S(1)*S(2)**2*C(19) + S(2)**3*C(20) +                          
     *      S(2)**2*S(3)*C(21) + (2.0D0*S(1)*S(2)*S(3) +                        
     *      S(2)*S(3)**2)*C(22)                                                 
C                                                                               
      DC = 5.0D0*S(1)**4*C(23) + (4.0D0*S(1)**3*S(3) +                          
     *     S(3)**4)*C(24) + (3.0D0*S(1)**2*S(3)**2 +                            
     *     2.0D0*S(1)*S(3)**3)*C(25) + 4.0D0*S(1)**3*S(2)*C(26) +               
     *     (3.0D0*S(2)*S(1)**2*S(3) + S(2)*S(3)**3)*C(27)                       
C                                                                               
      DD = 2.0D0*S(1)*S(2)*S(3)**2*C(28) +                                      
     *     3.0D0*S(1)**2*S(2)**2*C(29) + (2.0D0*S(1)*S(2)**2*S(3) +             
     *     S(2)**2*S(3)**2)*C(30) + 2.0D0*S(1)*S(2)**3*C(31) +                  
     *     S(2)**3*S(3)*C(32) + S(2)**4*C(33) +                                 
     *     6.0D0*S(1)**5*C(35) + (5.0D0*S(1)**4*S(3) +                          
     *     S(3)**5)*C(36)                                                       
C                                                                               
      DF = (4.0D0*S(1)**3*S(3)**2 + 2.0D0*S(1)*S(3)**4)*C(37) +                 
     *      3.0D0*S(1)**2*S(3)**3*C(38) + 5.0D0*S(1)**4*S(2)*C(39) +            
     *     (4.0D0*S(2)*S(1)**3*S(3) + S(2)*S(3)**4)*C(40) +                     
     *     (3.0D0*S(1)**2*S(2)*S(3)**2 +                                        
     *      2.0D0*S(1)*S(2)*S(3)**3)*C(41) +                                    
     *      4.0D0*S(1)**3*S(2)**2*C(42)                                         
C                                                                               
      DG = (3.0D0*S(1)**2*S(2)**2*S(3) + S(2)**2*S(3)**3)*C(43) +               
     *      2.0D0*S(1)*S(2)**2*S(3)**2*C(44) +                                  
     *      3.0D0*S(1)**2*S(2)**3*C(45) + (2.0D0*S(1)*S(2)**3*S(3) +            
     *      S(2)**3*S(3)**2)*C(46) + 2.0D0*S(1)*S(2)**4*C(47) +                 
     *      S(2)**4*S(3)*C(48) + S(2)**5*C(49)                                  
C                                                                               
      D1D = DA + DB + DC + DD + DF + DG                                         
C                                                                               
      DA2 = C(3) + 2.0D0*S(2)*C(5) + (S(1) + S(3))*C(6) +                       
     *      3.0D0*S(2)**2*C(9) + 2.0D0*S(2)*(S(1) + S(3))*C(11) +               
     *     (S(1)**2 + S(3)**2)*C(12) + S(1)*S(3)*C(13) +                        
     *      4.0D0*S(2)**3*C(15) + (S(1)**3 + S(3)**3)*C(18) +                   
     *      2.0D0*S(2)*(S(1)**2 + S(3)**2)*C(19)                                
C                                                                               
      DB2 = 3.0D0*S(2)**2*(S(1) + S(3))*C(20) +                                 
     *      2.0D0*S(1)*S(2)*S(3)*C(21) + (S(1)**2*S(3) +                        
     *      S(1)*S(3)**2)*C(22) + (S(1)**4 + S(3)**4)*C(26) +                   
     *     (S(1)**3*S(3) + S(1)*S(3)**3)*C(27) +                                
     *      S(1)**2*S(3)**2*C(28) + 2.0D0*S(2)*(S(1)**3 +                       
     *      S(3)**3)*C(29)                                                      
C                                                                               
      DC2 = 2.0D0*S(2)*(S(1)**2*S(3) + S(1)*S(3)**2)*C(30) +                    
     *      3.0D0*S(2)**2*(S(1)**2 + S(3)**2)*C(31) +                           
     *      3.0D0*S(2)**2*S(1)*S(3)*C(32) + 4.0D0*S(2)**3*(S(1) +               
     *      S(3))*C(33) + 5.0D0*S(2)**4*C(34)                                   
C                                                                               
      DD2 = (S(1)**5 + S(3)**5)*C(39) + (S(1)**4*S(3) +                         
     *       S(1)*S(3)**4)*C(40) + (S(1)**3*S(3)**2 +                           
     *       S(1)**2*S(3)**3)*C(41) + 2.0D0*S(2)*(S(1)**4 +                     
     *       S(3)**4)*C(42) + 2.0D0*S(2)*(S(1)**3*S(3) +                        
     *       S(1)*S(3)**3)*C(43) + 2.0D0*S(1)**2*S(2)*S(3)**2*C(44)             
C                                                                               
      DF2 = 3.0D0*S(2)**2*(S(1)**3 + S(3)**3)*C(45) +                           
     *      3.0D0*S(2)**2*(S(1)**2*S(3) + S(1)*S(3)**2)*C(46) +                 
     *      4.0D0*S(2)**3*(S(1)**2 + S(3)**2)*C(47) +                           
     *      4.0D0*S(2)**3*S(1)*S(3)*C(48) +                                     
     *      5.0D0*S(2)**4*(S(1) + S(3))*C(49) +                                 
     *     6.0D0*S(2)**5*C(50)                                                  
C                                                                               
      D2D = DA2 + DB2 + DC2 + DD2 + DF2                                         
C                                                                               
      DA3 = C(2) + 2.0D0*S(3)*C(4) + S(2)*C(6) + S(1)*C(7) +                    
     *      3.0D0*S(3)**2*C(8) + (2.0D0*S(1)*S(3) + S(1)**2)*C(10) +            
     *      S(2)**2*C(11) + 2.0D0*S(2)*S(3)*C(12) + S(1)*S(2)*C(13) +           
     *      4.0D0*S(3)**3*C(14)                                                 
C                                                                               
      DB3 = (S(1)**3 + 3.0D0*S(1)*S(3)**2)*C(16) +                              
     *       2.0D0*S(1)**2*S(3)*C(17) + 3.0D0*S(2)*S(3)**2*C(18) +              
     *       2.0D0*S(2)**2*S(3)*C(19) + S(2)**3*C(20) +                         
     *       S(2)**2*S(1)*C(21) + S(2)*(S(1)**2 +                               
     *       2.0D0*S(1)*S(3))*C(22) + 5.0D0*S(3)**4*C(23) + (S(1)**4 +          
     *       4.0D0*S(1)*S(3)**3)*C(24)                                          
C                                                                               
      DC3 = (2.0D0*S(1)**3*S(3) + 3.0D0*S(1)**2*S(3)**2)*C(25) +                
     *       4.0D0*S(2)*S(3)**3*C(26) + S(2)*(S(1)**3 +                         
     *       3.0D0*S(1)*S(3)**2)*C(27) +                                        
     *       2.0D0*S(2)*S(1)**2*S(3)*C(28) + 3.0D0*S(2)**2*S(3)**2*C(29)        
C                                                                               
      DD3 = S(2)**2*(S(1)**2 + 2.0D0*S(1)*S(3))*C(30) +                         
     *      2.0D0*S(2)**3*S(3)*C(31) + S(2)**3*S(1)*C(32) +                     
     *      S(2)**4*C(33) + 6.0D0*S(3)**5*C(35) + (S(1)**5 +                    
     *      5.0D0*S(1)*S(3)**4)*C(36) + (2.0D0*S(1)**4*S(3) +                   
     *      4.0D0*S(1)**2*S(3)**3)*C(37) + 3.0D0*S(1)**3*S(3)**2*C(38)          
C                                                                               
      DF3 = 5.0D0*S(2)*S(3)**4*C(39) + S(2)*(S(1)**4 +                          
     *      4.0D0*S(1)*S(3)**3)*C(40) + S(2)*(2.0D0*S(1)**3*S(3) +              
     *      3.0D0*S(1)**2*S(3)**2)*C(41) +                                      
     *      4.0D0*S(2)**2*S(3)**3*C(42) + S(2)**2*(S(1)**3 +                    
     *      3.0D0*S(1)*S(3)**2)*C(43) +                                         
     *      2.0D0*S(1)**2*S(2)**2*S(3)*C(44) +                                  
     *      3.0D0*S(2)**3*S(3)**2*C(45) + S(2)**3*(S(1)**2 +                    
     *      2.0D0*S(1)*S(3))*C(46) + 2.0D0*S(2)**4*S(3)*C(47) +                 
     *      S(2)**4*S(1)*C(48) + S(2)**5*C(49)                                  
C                                                                               
      D3D = DA3 + DB3 + DC3 + DD3 + DF3                                         
C                                                                               
      DQ1=(1.0D0 - TANH(AL(2)*S(2)))*(1.0D0-TANH(AL(3)*S(3)))                   
      DQ2=(1.0D0 - TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(3)*S(3)))                   
      DQ3=(1.0D0 - TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(2)*S(2)))                   
C                                                                               
      DAP1=-12.0D0*AP/R(1)                                                      
      DAP2=-12.0D0*AP/(R(2) - 0.2D0)                                            
      DAP3=-12.0D0*AP/R(3)                                                      
C                                                                               
      DB1=-4.0D0*B/R(1)**5                                                      
      DB3=-4.0D0*B/R(3)**5                                                      
C                                                                               
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)DCRS1=1.0D30                             
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)GO TO 144                                
C                                                                               
C                                                                               
       DCSR1=-AL(1)/(COSH(AL(1)*S(1))**2)                                       
C                                                                               
C144   IF(COSH(AL(2)*S(2)).LT.1.0D-15)DCRS2=1.0D30                              
C      IF(COSH(AL(2)*S(2)).LT.1.0D-15)GO TO 145                                 
C                                                                               
      DCSR2=-AL(2)/(COSH(AL(2)*S(2))**2)                                        
C                                                                               
C145   IF(COSH(AL(3)*S(3)).LT.1.0D-15)DCRS3=1.0D30                              
C      IF(COSH(AL(3)*S(3)).LT.1.0D-15)GO TO 146                                 
C                                                                               
      DCSR3=-AL(3)/(COSH(AL(3)*S(3))**2)                                        
                                                                                
 146  PSI13=2.0D0*(-GAM)*TEMP*(R(1) + R(3))                                     
C                                                                               
      DCHI1=2.0D0*ALF*R(1)*EXP(-ALF*(R(1)**2))                                  
      DCHI2=2.0D0*ALF*R(2)*EXP(-ALF*(R(2)**2))                                  
      DCHI3=2.0D0*ALF*R(3)*EXP(-ALF*(R(3)**2))                                  
C                                                                               
C      IF(VM(1).GT.1.0D20)DVR1=1.0D30                                           
C      IF(VM(1).GT.1.0D20)GO TO 153                                             
      DVR1 = 2.0D0*DE(1)*(1.0D0 - EXP(-BET(1)*(R(1)-XRE(1))))*                  
     *       BET(1)*EXP(-BET(1)*(R(1)-XRE(1)))                                  
C                                                                               
C 153  IF(VM(2).GT.1.0D20)DVR2=1.0D30                                           
C      IF(VM(2).GT.1.0D20)GO TO 154                                             
C                                                                               
      DVR2 = 2.0D0*DE(2)*(1.0D0 - EXP(-BET(2)*(R(2)-XRE(2))))*                  
     *       BET(2)*EXP(-BET(2)*(R(2)-XRE(2)))                                  
C                                                                               
C 154  IF(VM(3).GT.1.0D20)DVR3=1.0D30                                           
C      IF(VM(3).GT.1.0D20)GO TO 155                                             
C                                                                               
      DVR3 = 2.0D0*DE(3)*(1.0D0 - EXP(-BET(3)*(R(3)-XRE(3))))*                  
     *       BET(3)*EXP(-BET(3)*(R(3)-XRE(3)))                                  
C                                                                               
 155  DEGSDR(1) = VM(2)*CHI(3)*DCHI1 + CHI(2)*CHI(3)*DVR1 +                     
     *          VM(3)*CHI(2)*DCHI1 + D1D*Q*TEMP + P*DQ1*TEMP*DCSR1 +            
     *          DAP1 + DB1 + P*Q*PSI13                                          
      DEGSDR(2) = DVR2*CHI(1)*CHI(3) + VM(1)*CHI(3)*DCHI2 +                     
     *          VM(3)*CHI(1)*DCHI2 + D2D*Q*TEMP + P*DQ2*TEMP*DCSR2 +            
     *          DAP2                                                            
      DEGSDR(3) = VM(2)*CHI(1)*DCHI3 + VM(1)*CHI(2)*DCHI3 +                     
     *          DVR3*CHI(1)*CHI(2) + D3D*Q*TEMP + P*DQ3*TEMP*DCSR3 +            
     *          P*Q*PSI13 + DAP3 + DB3                                          
C                                                                               
      DO 222 I=1,3                                                              
         DEGSDR(I) = DEGSDR(I)/HART*BOHR                                        
 222  continue                                                                  
C                                                                               
C END OF NDER SWITCH AROUND DERIVATIVES                                         
C                                                                               
      ENDIF                                                                     
C                                                                               
      J=J+1                                                                     
      GO TO 98                                                                  
C                                                                               
  113 CONTINUE                                                                  
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',                    
     *                'Schinke-Lester potential energy surface')                
601   FORMAT(/,2X,T5,'Step size for the numerical derivatives = ',E10.5)        
602   FORMAT(2X,T5,'The Schinke-Lester surface has analytical ',                
     *             'derivatives if R(1) and R(2) are ',                         
     *       /,2X,T5,'less than 50 bohr, otherwise the derivatives ',           
     *               'are computed numerically.')                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
       stop                                                                     
       END                                                                      
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   
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
         PARAMETER (BOHR = 0.529177D0)                                          
         PARAMETER (HART = 27.21161D0)                                          
C                                                                               
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)                                  
         COMMON /NDERCM/ STEP                                                   
         COMMON /PRECM/  C(50),A,B,ALF,GAM,AL(3),                               
     +                   X(3),S(3),VM(3),CHI(3),V(4)                            
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
C   are in eV and the lengths are in Angstroms.                                 
C                                                                               
         DATA DE  /4.035D0,4.1662D0,4.035D0   /                                 
         DATA XRE /0.9923D0,0.7590D0,0.9923D0 /                                 
         DATA BET /2.4113D0,2.0682D0, 2.4113D0/                                 
C                                                                               
C   Initialize the step size used for the numerical derivatives                 
C   of the energy with respect to the coordinates.                              
C                                                                               
      DATA STEP/5.0D-7/                                                         
C                                                                               
      DATA A   / 0.157002D-03/                                                  
      DATA B   /-0.162301D-01/                                                  
      DATA ALF / 0.9D0       /                                                  
      DATA GAM / 0.045D0     /                                                  
      DATA AL  / 1.0D0,0.8D0,1.0D0/                                             
C                                                                               
      DATA C/0.354195D+01, 0.302780D+01, 0.577312D+00, 0.539877D+01,            
     +       0.233027D+01, 0.157276D+02, 0.701851D+01,-0.478535D+02,            
     +      -0.106915D+01, 0.268137D+02,-0.104806D+02,-0.182660D+02,            
     +       0.833941D+01, 0.872497D+02, 0.917870D+01,-0.697766D+02,            
     +       0.434862D+02,-0.156633D+02, 0.380162D+02,-0.796542D+01,            
     +       0.602854D+02, 0.500800D+02,-0.545761D+02, 0.410210D+02,            
     +      -0.890739D+01, 0.191135D+02,-0.512653D+02,                          
     +       0.365096D+02,-0.235972D+02,-0.102322D+03, 0.837770D+00,            
     +       0.258200D+02,-0.198105D+01,-0.484189D+01, 0.105377D+02,            
     +       0.372352D+01,-0.207237D+02, 0.223342D+02, 0.244787D+01,            
     +       0.979572D+01,-0.135612D+01,-0.724496D+01, 0.234984D+02,            
     +       0.136364D+03, 0.253799D+01,-0.203821D+02, 0.699846D+01,            
     +       0.198529D+02,-0.475902D+01, 0.236328D+01/                          
C                                                                               
         END                                                                    
                                                                                
                                                                                
                                                                                
