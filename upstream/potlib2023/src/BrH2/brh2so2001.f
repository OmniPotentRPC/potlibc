      subroutine prepot                                                         
C
C   System:           BrH2
C   Functional form:  extended LEPS plus three-center term
C   Common name:      BrH2SO
C   Interface:        potlib2001
C   Reference:        G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao
C                     J. Phys. Chem. 99, 207 (1995)
C
C   Number of bodies: 3
C   Number of derivatives: 0
C   Number of electronic surfaces: 2
C
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         COMMON /PRE22CM/ D22(3),DE22(3),RE22(3),BETAA22(3),R222,               
     +                    CCON,BETA22(3)                                        
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      call prepot11                                                     4s19m94 
      call prepot12                                                     4s19m94 
      call prepot22                                                     4s19m94 
C                                                                               
      EZERO(1)=DE11(2)                                                          
      EZERO(2)=DE22(2)                                                          
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
C                                                                               
      INDEXES(1) = 35                                                           
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      return                                                            4s19m94 
      END                                                                       
C                                                                               
      SUBROUTINE POT                                                            
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
      if (NASURF(1,1).NE.0) call pot11                                          
      if (NASURF(2,2).NE.0) call pot22                                          
      if (NASURF(1,2)+NASURF(2,1).NE.0) call pot12                              
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                            4s19m94 
      end                                                               4s19m94 
C                                                                               
      SUBROUTINE PREPOT11                                                       
C                                                                               
C                                                                               
C   System:          BrH2                                                       
C   Functional form: extended-LEPS plus three-center term                       
C   Common name:     BrH2SEC                                                    
C   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, et. al.           
C                    to be published                                            
C                                                                               
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in DATA statements.                   
C   Coordinates, potential energy, and derivatives are passed                   
C   through the common block POTCM:                                             
C                                                                               
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCOM:                                          
C                  COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                        
C   For the this potential energy surface the potential energy in the           
C   BC valley is equal to the energy of the HH diatomic.                        
C   The energy in the AC valley is equal to the energy of the HBr diatomic,     
C   which is also the energy value in the BA valley.                            
C   All the information passed through the common blocks PV1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(Br-H)                                               
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-Br)                                               
C   The zero of energy is defined at Br "infinitely" far from the H2 diatomic.  
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states are available                         
C        NDER  = 0 => no derivatives should be calculated                       
C        NDER  = 1 => calculate first derivatives                               
C        NFLAG  - these 20 integer values can be used to flag options           
C                within the potential;                                          
C                                                                               
C                                                                               
C   IPRT - FORTRAN unit number to which all output will be written              
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NASURF(1,1)         1                                        
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
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
C         A1  =  V3CP1                                                          
C         A2  =  V3CP2                                                          
C         A3  =  V3CP3                                                          
C         A4  =  V3CP4                                                          
C         A5  =  V3CP5                                                          
C         A6  =  V3CP6                                                          
C                                                                               
         WRITE (NFLAG(18),100) D11, RE11, BETAA11,                              
     +                    ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,                         
     +                    V3CP1, V3CP2, V3CP3, V3CP4, V3CP5                     
C                                                                               
C   Convert the constants to atomic units                                       
C                                                                               
         DO 10 I = 1, 3                                                         
               DE11(I)   = D11(I) / CHARKC                                      
               RE11(I)   = RE11(I) / CANGAU                                     
               BETA11(I) = BETAA11(I) * CANGAU                                  
 10      CONTINUE                                                               
C                                                                               
C    Initialize the asyptotic energy terms                                      
C                                                                               
               EASYBC = DE11(2)                                                 
               EASYAC = DE11(3)                                                 
               EASYBA = DE11(1)                                                 
C                                                                               
 100  FORMAT(/,2X,T5,'PREPOT has been called for the BrH2 potential',           
     *               ' SEC',                                                    
     *       /,2X,T5,'Potential parameters:',                                   
     *       /,2X,T5,'Bond:',T38,'BrH',T49,'HH',T58,'HBr',                      
     *       /,2X,T5,'Dissociation energies:',                                  
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Equilibrium bond lengths:',                               
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Morse betas:',                                            
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Parameters for ZHH: ',T45,F10.5,T55,F10.5,                
     *       /,2X,T5,'Parameters for ZHBr: ',T45,F10.5,T55,F10.5,               
     *       /,2X,T5,'Parameters for V3c: ',5(T45,F10.5,T55,F10.5,/))           
C                                                                               
         RETURN                                                                 
         END                                                                    
C                                                                               
      SUBROUTINE POT11                                                          
C                                                                               
C         ENTRY POT11(rr,e,n)                                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
         PARAMETER(xlamb11 =  0.4568d0/27.21161d0)                              
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         DIMENSION RR(3),B(3)                                                   
C                                                                               
        rr(1)=r(1)                                                              
        rr(2)=r(2)                                                              
        rr(3)=r(3)                                                              
C                                                                               
C   Compute bond dependent Sato parameters                                      
C                                                                               
      R1P2 = RR(1) * RR(1)                                                      
      R2P2 = RR(2) * RR(2)                                                      
      R3P2 = RR(3) * RR(3)                                                      
      R2I  = 1.0D0/RR(2)                                                        
      ZCHI = 0.0D0                                                              
C                                                                               
C   CAPR IS THE Br TO H2 DISTANCE                                               
C                                                                               
      CAPR2 = 0.5D0 * (R1P2 + R3P2 - 0.5D0 * R2P2)                              
      IF (CAPR2 .LE. 0.D0) THEN                                                 
          CAPR2 = 0.0D0                                                         
          CAPR  = 0.0D0                                                         
      ELSE                                                                      
          CAPR = SQRT(CAPR2)                                                    
C                                                                               
C   Compute Chi, the angle between CAPR and RHH                                 
C                                                                               
          COSCHI = 0.5D0*(R3P2 - R1P2)/(RR(2)*CAPR)                             
          SCHIP2 = 1.0D0 - COSCHI*COSCHI                                        
          ZCHI = CAPR2*(ZHBRP1*SCHIP2 + ZHBRP2*SCHIP2*SCHIP2)                   
      ENDIF                                                                     
      Z(1) = 0.1675D0 + ZCHI                                                    
      Z(3) = Z(1)                                                               
      Z(2) = ZHHP1*R2P2/(ZHHP2+R2P2)                                            
C                                                                               
C   Compute the LEPS constants                                                  
C                                                                               
         DO 20 I = 1, 3                                                         
               ZPO(I)   = 1.0D0+Z(I)                                            
               OP3Z(I)  = 1.0D0+3.0D0*Z(I)                                      
               TOP3Z(I) = 2.0D0*OP3Z(I)                                         
               ZP3(I)   = Z(I)+3.0D0                                            
               TZP3(I)  = 2.0D0*ZP3(I)                                          
               DO4Z(I)  = DE11(I)/4.0D0/ZPO(I)                                  
               B(I) = BETA11(I)*DO4Z(I)*2.0D0                                   
 20      CONTINUE                                                               
         PE = 0.D0                                                              
         DO 30 I = 1, 3                                                         
               X(I)    = EXP(-BETA11(I)*(RR(I)-RE11(I)))                        
               COUL(I) = DO4Z(I) * (ZP3(I)*X(I) - TOP3Z(I)) * X(I)              
               EXCH(I) = DO4Z(I) * (OP3Z(I)*X(I) - TZP3(I)) * X(I)              
               PE      = PE + COUL(I)                                           
 30      CONTINUE                                                               
         RAD = SQRT((EXCH(1) - EXCH(2)) ** 2 +                                  
     &              (EXCH(2) - EXCH(3)) ** 2 +                                  
     &              (EXCH(3) - EXCH(1)) ** 2 )                                  
         PE  = PE - (RAD / R211) + EZERO(1)                                     
C                                                                               
C   Compute the three-center term                                               
C                                                                               
         R1P3   = RR(1) + RR(3)                                                 
         A1  =  V3CP1                                                           
         A2  =  V3CP2                                                           
         A3  =  V3CP3                                                           
         A4  =  V3CP4                                                           
         A5  =  V3CP5                                                           
         A6  =  V3CP6                                                           
C                                                                               
         R1M3   = RR(1) - RR(3)                                                 
         R1P3M2 = R1P3 - RR(2)                                                  
         R1M3P2 = R1M3*R1M3                                                     
         R132P2 = R1P3M2*R1P3M2                                                 
C                                                                               
C   Compute the bond angle from the three internal coordinates                  
C                                                                               
         DIST  = RR(1)*RR(1) + RR(3)*RR(3) - RR(2)*RR(2)                        
         COSTH = DIST/(2.0D0*RR(1)*RR(3))                                       
         SINTP2 = 1.0D0 - COSTH*COSTH                                           
         SINTP6 = SINTP2*SINTP2*SINTP2                                          
         SINTP8 = SINTP6*SINTP2                                                 
         FCN4  = 1.D0 + A5*SINTP2 + A6*SINTP8                                   
C                                                                               
         V3C = A1*FCN4*exp(-(a2*r1p3+a3*r1m3p2+a4*r132p2))                      
C                                                                               
C   Sum all the terms that make up the energy                                   
C                                                                               
         PE = PE + V3C                                                          
         ENGYGS = PE+xlamb11/3.d0                                               
c                                                                               
         RETURN                                                                 
         END                                                                    
C                                                                               
      SUBROUTINE PREPOT12                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         RETURN                                                                 
         END                                                                    
C                                                                               
      SUBROUTINE POT12                                                          
C                                                                               
C      ENTRY POT12(rr,e,n)                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
         PARAMETER(xlamb12 = 0.4568d0/27.21161d0)                               
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         ENGYIJ(1)=-dsqrt(2.d0)*xlamb12/3.d0                                    
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT22                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         COMMON /PRE22CM/ D22(3),DE22(3),RE22(3),BETAA22(3),R222,               
     +                    CCON,BETA22(3)                                        
C                                                                               
         WRITE (NFLAG(18),100) D22, RE22, BETAA22,                              
     +                 ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,                            
     +                 V3CP1, V3CP2, V3CP3, V3CP4, V3CP5                        
C                                                                               
C   Convert the constants to atomic units                                       
C                                                                               
         DO 10 I = 1, 3                                                         
               DE22(I)   = D22(I) / CHARKC                                      
               RE22(I)   = RE22(I) / CANGAU                                     
               BETA22(I) = BETAA22(I) * CANGAU                                  
 10      CONTINUE                                                               
C                                                                               
C    Initialize the asyptotic energy terms                                      
C                                                                               
               EASYBC = DE22(2)                                                 
               EASYAC = DE22(3)                                                 
               EASYBA = DE22(1)                                                 
C                                                                               
100   FORMAT(/,2X,T5,'PREPOT22 has been called for the BrH2 potential',         
     *               ' SEC',                                                    
     *       /,2X,T5,'Potential parameters:',                                   
     *       /,2X,T5,'Bond:',T38,'BrH',T49,'HH',T58,'HBr',                      
     *       /,2X,T5,'Dissociation energies:',                                  
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Equilibrium bond lengths:',                               
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Morse betas:',                                            
     *            T35,F10.5,T45,F10.5,T55,F10.5,                                
     *       /,2X,T5,'Parameters for ZHH: ',T45,F10.5,T55,F10.5,                
     *       /,2X,T5,'Parameters for ZHBr: ',T45,F10.5,T55,F10.5,               
     *       /,2X,T5,'Parameters for V3c: ',5(T45,F10.5,T55,F10.5,/))           
C                                                                               
         RETURN                                                                 
         END                                                                    
C                                                                               
      SUBROUTINE POT22                                                          
C                                                                               
C      ENTRY POT22(rr,e,n)                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (xlamb22 = 0.4568d0/27.21161d0)                              
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         COMMON /PRE22CM/ D22(3),DE22(3),RE22(3),BETAA22(3),R222,               
     +                    CCON,BETA22(3)                                        
C                                                                               
         DIMENSION RR(3)                                                        
C                                                                               
         rr(1)=r(1)                                                             
         rr(2)=r(2)                                                             
         rr(3)=r(3)                                                             
         y2=exp(-beta22(2)*(rr(2)-re22(2)))                                     
         y1=exp(-beta22(1)*(rr(1)-re22(1)))                                     
         y3=exp(-beta22(3)*(rr(3)-re22(3)))                                     
         vmorse=de22(2)*(y2**2-2.d0*y2)                                         
         vanti1=de22(1)*(y1**2+2.d0*y1)                                         
         vanti3=de22(3)*(y3**2+2.d0*y3)                                         
         pe=EZERO(2)  +                                                         
     +      vmorse+(vanti1+vanti3)*ccon+2.d0*xlamb22/3.d0                       
         ENGYES(1)=pe                                                           
c                                                                               
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
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER(CANGAU = 0.52917706D0)                                       
         PARAMETER(CHARKC = 627.5095D0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /SATOCM/  ZHHP1,ZHHP2,ZHBRP1,ZHBRP2,Z(3),DZ(3,3)                
         COMMON /V3CCM/   V3CP1,V3CP2,V3CP3,V3CP4,V3CP5,V3CP6,                  
     +                    V3C,DV3C(3)                                           
         COMMON /LENGCM/  X(3),COUL(3),EXCH(3),RAD                              
         COMMON /LSATCM/  ZPO(3),OP3Z(3),ZP3(3),TZP3(3),TOP3Z(3),               
     +                    DO4Z(3)                                               
         COMMON /LDERCM/  DXDR(3,3),DEQ(3),DERAD(3)                             
         COMMON /PRE11CM/ D11(3),DE11(3),RE11(3),BETAA11(3),R211,               
     +                    BETA11(3)                                             
C                                                                               
         COMMON /PRE22CM/ D22(3),DE22(3),RE22(3),BETAA22(3),R222,               
     +                    CCON,BETA22(3)                                        
C                                                                               
C   Initialize the potential parameters NASURF and NDER                         
C                                                                               
         DATA NASURF /1,35*0/                                                   
      DATA NDER /0/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                               
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,1,0/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
         DATA D11    /90.4473D0, 109.559D0, 90.4473D0/                          
         DATA RE11    /1.41443D0, 0.74144D0, 1.41443D0/                         
         DATA BETAA11 /1.81088D0, 1.94198D0, 1.81088D0/                         
         DATA R211    /1.41421356D0/                                            
C                                                                               
         DATA D22    /90.4473D0, 109.559D0, 90.4473D0/                          
         DATA RE22    /1.41443D0, 0.74144D0, 1.41443D0/                         
         DATA BETAA22 /1.81088D0, 1.94198D0, 1.81088D0/                         
         DATA R222    /1.41421356D0/                                            
C                                                                               
         DATA CCON /0.5d0/                                                      
C                                                                               
         DATA ZHHP1  / 2.91230499D-1 /                                          
         DATA ZHHP2  / 4.99160974D0  /                                          
         DATA ZHBRP1 / 3.84262775D-3 /                                          
         DATA ZHBRP2 / 2.87642885D-3 /                                          
         DATA V3CP1  / 1.67457884D0  /                                          
         DATA V3CP2  / 7.09518875D-1 /                                          
         DATA V3CP3  / 1.87582679D-1 /                                          
         DATA V3CP4  / 2.43517074D-1 /                                          
         DATA V3CP5  / 2.04145206D-1 /                                          
         DATA V3CP6  / 2.57833022D0  /                                          
C                                                                               
         END                                                                    
