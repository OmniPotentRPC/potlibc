C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          IH2                                                        
C   Functional form:                                                            
C   Common name:                                                                
C   Reference:       L. Raff et al.                                             
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
C   This potential is written such that:                                        
C                  R(1) = R(I-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-I)                                                
C   The zero of energy is defined at I "infinitely" far from the H2 diatomic.   
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
         COMMON /PARMCM/ D1(3),D3(3),ALPH(3),BETA(3),RE(3),                     
     +                   A(3),C(3),SIG(3),RC(3)                                 
C         DIMENSION E1(3),E3(3),FJ(3),DE1(3),DE3(3),Q(3),DJ(3),DQ(3)            
C                                                                               
C   Echo the potential name and the potential parameters to the                 
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 600) D1, D3, ALPH, BETA, RE                          
         WRITE (NFLAG(18), 610) A, C, SIG, RC                                   
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the IH2 ',                    
     *                'potential energy surface of Raff et al.',                
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T48, 'I-H', T60, 'H-H', T72, 'H-H',               
     *        /,2X,T5,'Dissociation energies (eV):',                            
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'D3 energy values (eV):',                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Alpha:',                                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Beta:',                                                  
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Bohr):',                       
     *        T44, F11.5, T56, F11.5, T68, F11.5)                               
610   FORMAT (/,2X,T5,'A:',                                                     
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'C:',                                                     
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Sigma:',                                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'RC (Bohr):',                                             
     *        T44, F11.5, T56, F11.5, T68, F11.5)                               
C                                                                               
C   Convert the potential parameters to atomic units.                           
C                                                                               
      DO 5 I = 1,3                                                              
      D1(I) = D1(I)/27.21161D0                                                  
      D3(I) = D3(I)/27.21161D0                                                  
      C(I) = C(I)/27.21161D0                                                    
 5    CONTINUE                                                                  
C                                                                               
C    5 C(I) = C(I)/27.21161D0                                                   
C                                                                               
C   Initialize the asymptotic energy values                                     
C                                                                               
         EASYAB = D1(1)                                                         
         EASYBC = D1(2)                                                         
         EASYAC = D1(3)                                                         
C                                                                               
C                                                                               
      EZERO(1)=D1(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='L. Raff et al.'                                                  
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
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the IH diatomic, with the             
C   IH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(I-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-I)                                                
C   The zero of energy is defined at I "infinitely" far from the H2 diatomic.   
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
         COMMON /PARMCM/ D1(3),D3(3),ALPH(3),BETA(3),RE(3),                     
     +                   A(3),C(3),SIG(3),RC(3)                                 
         DIMENSION E1(3),E3(3),FJ(3),DE1(3),DE3(3),Q(3),DJ(3),DQ(3)             
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
      QT=0.0D0                                                                  
      SUMJ=0.0D0                                                                
      DO 2 I=1,3                                                                
      REL=RE(I)-R(I)                                                            
      EXA=EXP(ALPH(I)*REL)                                                      
      E1(I)=D1(I)*(1.0D0-EXA)**2-D1(I)                                          
      IF (NDER .EQ. 1) DE1(I)=2.0D0*D1(I)*ALPH(I)*EXA*(1.0D0-EXA)               
      IF (R(I).GT.RC(I)) GO TO 10                                               
      EXB=EXP(BETA(I)*REL)                                                      
      E3(I)=D3(I)*(1.0D0+EXB)**2-D3(I)                                          
      IF (NDER .EQ. 1) DE3(I)=-2.0D0*D3(I)*BETA(I)*EXB*(1.0D0+EXB)              
      GO TO 1                                                                   
   10 EXS=EXP(-SIG(I)*R(I))                                                     
      E3(I)=C(I)*(R(I)+A(I))*EXS                                                
      IF (NDER .EQ. 1) DE3(I)=C(I)*EXS*(1.0D0-SIG(I)*(R(I)+A(I)))               
    1 CONTINUE                                                                  
      FJ(I)=0.5D0*(E1(I)-E3(I))                                                 
      Q(I)=0.5D0*(E1(I)+E3(I))                                                  
      IF (NDER .EQ. 1) THEN                                                     
          DJ(I)=0.5D0*(DE1(I)-DE3(I))                                           
          DQ(I)=0.5D0*(DE1(I)+DE3(I))                                           
      ENDIF                                                                     
      QT=QT+Q(I)                                                                
      SUMJ=SUMJ+FJ(I)                                                           
    2 CONTINUE                                                                  
      EX1=(FJ(1)-FJ(2))**2                                                      
      EX2=(FJ(2)-FJ(3))**2                                                      
      EX3=(FJ(3)-FJ(1))**2                                                      
      EXCH=SQRT(0.5D0*(EX1+EX2+EX3))                                            
      ENGYGS=QT-EXCH                                                            
      ENGYGS = ENGYGS + EZERO(1)                                                
      IF (NDER .EQ. 1) THEN                                                     
          DO 3 I=1,3                                                            
3              DEGSDR(I)=DQ(I)-0.5D0*DJ(I)*(3.0D0*FJ(I)-SUMJ)/EXCH              
      ENDIF                                                                     
C                                                                               
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the IH2 ',                    
     *                'potential energy surface of Raff et al.',                
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T48, 'I-H', T60, 'H-H', T72, 'H-H',               
     *        /,2X,T5,'Dissociation energies (eV):',                            
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'D3 energy values (eV):',                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Alpha:',                                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Beta:',                                                  
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Bohr):',                       
     *        T44, F11.5, T56, F11.5, T68, F11.5)                               
610   FORMAT (/,2X,T5,'A:',                                                     
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'C:',                                                     
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'Sigma:',                                                 
     *        T44, F11.5, T56, F11.5, T68, F11.5,                               
     *        /,2X,T5,'RC (Bohr):',                                             
     *        T44, F11.5, T56, F11.5, T68, F11.5)                               
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
C*****                                                                          
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)                                   
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
         COMMON /PARMCM/ D1(3),D3(3),ALPH(3),BETA(3),RE(3),                     
     +                   A(3),C(3),SIG(3),RC(3)                                 
C                                                                               
C   Initialize the flags for the potential                                      
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
C   The energy values are in eV.                                                
C                                                                               
         DATA D1   / 3.194D0,   4.7466D0,  3.194D0/                             
         DATA D3   / 1.44399D0, 1.9668D0,  1.44399D0/                           
         DATA ALPH / 0.9468D0,  1.04435D0, 0.9468D0/                            
         DATA BETA / 0.794D0,   1.0001D0,  0.794D0/                             
C                                                                               
C   The lengths are in bohr.                                                    
C                                                                               
         DATA RE  /  3.032D0,    1.402D0,    3.032D0/                           
         DATA A   / -3.60481D0,  1.0D0,     -3.60481D0/                         
         DATA C   / 73599.62D0, 25.5530D0,  73599.62D0/                         
         DATA SIG /  2.47075D0,  1.67564D0,  2.47075D0/                         
         DATA RC  /  4.25D0,     1.6D0,      4.25D0/                            
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
