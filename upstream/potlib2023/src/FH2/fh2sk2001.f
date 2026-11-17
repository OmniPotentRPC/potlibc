C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:           FH2                                                       
C   Functional form:  Extended LEPS (London-Erying-Polyani-Sato)                
C   Common name:      SK                                                        
C   Reference:        G. C. Schatz and A. Kuppermann                            
C                     J. Chem. Phys. 72, 2737 (1980).                           
C
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 1
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
C      coordinates - bohrs                                                      
C      derivatives - hartrees/bohr                                              
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the F            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in the BLOCK DATA subprogram PTPACM                                  
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R(F-first H)                                
C                            R(2) = R(first H-second H)                         
C                            R(3) = R(second H-F)                               
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, energy, and derivatives of                     
C        the energy with respect to the coordinates.                            
C        passes the control flags where                                         
C        NDER  = 0 => no derivatives are calculated                             
C        NDER  = 1 => calculate first derivatives of the energy for the         
C                     ground electronic state with respect to the coordinates   
C        NFLAG  - Control flags                                                 
C      NFLAG(18-20)                                                             
C        passes the FORTRAN unit number used for potential output               
C      /ASYCM/ EASYAB, EASYBC, EASYAC                                           
C        passes the energy in the three asymptotic valleys for an A+BC system.  
C        The energy in the AB valley, EASYAB, is equal to the energy            
C        of the C atom "infinitely" far from the AB diatomic and R(AB) set      
C        equal to Re(AB), the equilibrium bond length for the AB diatomic.      
C        In this potential the AB valley represents H infinitely far from       
C        the FH diatomic and R(FH) equal to Re(FH).                             
C        Similarly, the terms EASYBC and EASYAC represent the energies in the   
C        H2 and the other HF asymptotic valleys, respectively.                  
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),            
     *                   DO4Z(3), B(3)                                          
         PARAMETER (CKCAL = 627.5095D0)                                         
         PARAMETER (CANGS =   0.529177106D0)                                    
C                                                                               
C   Echo the potential parameters                                               
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 100) DE, RE, BETA, Z                                 
C                                                                               
100   FORMAT(/,1X,'*****',1X,'Potential Energy Surface',1X,'*****',/,           
     *       /,1X,T5,'FH2 SK extended LEPS potential energy surface',           
     *        //, 1X, T5, 'Parameters:',                                        
     *        /, 1X, T5, 'Bond', T47, 'F-H', T58, 'H-H', T69, 'H-F',            
     *        /, 1X, T5, 'Dissociation energies (kcal/mol):',                   
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 1X, T5, 'Equilibrium bond lengths (Angstroms):',               
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 1X, T5, 'Morse beta parameters (Angstroms**-1):',              
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 1X, T5, 'Sato parameters:',                                    
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')                 
C                                                                               
      DO  10 I = 1,3                                                            
C                                                                               
C   Convert to atomic units                                                     
C                                                                               
             DE(I)  = DE(I)/CKCAL                                               
             RE(I)   = RE(I)/CANGS                                              
             BETA(I) = BETA(I)*CANGS                                            
C                                                                               
C   Compute useful constants                                                    
C                                                                               
             ZPO(I)   = 1.0D0 + Z(I)                                            
             OP3Z(I)  = 1.0D0 + 3.0D0*Z(I)                                      
             TOP3Z(I) = 2.0D0*OP3Z(I)                                           
             ZP3(I)   = Z(I) + 3.0D0                                            
             TZP3(I)  = 2.0D0*ZP3(I)                                            
             DO4Z(I)  = DE(I)/4.0D0/ZPO(I)                                      
             B(I)     = BETA(I)*DO4Z(I)*2.0D0                                   
10    CONTINUE                                                                  
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = DE(1)                                                     
             EASYBC = DE(2)                                                     
             EASYAC = DE(3)                                                     
C                                                                               
C                                                                               
      EZERO(1)=DE(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='G. C. Schatz and A. Kuppermann'                                  
       REF(2)='J. Chem. Phys. 72, 2737(1980)'                                   
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
C*****                                                                          
C                                                                               
         SUBROUTINE POT                                                         
C                                                                               
C   System:          ABC                                                        
C   Functional form: Extended LEPS (London-Erying-Polyani-Sato)                 
C   References:      P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner,     
C                    and C. E. Young                                            
C                    J. Chem. Phys. 44, 1168 (1966)                             
C                                                                               
C   The potential parameters must be passed through the common blocks           
C   PT1CM, ASYCM, PT2CM, SATOCM, and LEPSCM.                                    
C   All information passed through the common blocks PT1CM, ASYCM,              
C   SATOCM, and LEPSCM must be in Hartree atomic units.                         
C                                                                               
C        For the reaction: A + BC -> AB + C we write:                           
C                          R(1) = R(A-B)                                        
C                          R(2) = R(B-C)                                        
C                          R(3) = R(C-A)                                        
C                                                                               
C   NOTE: The potential energy at the reactant asymptote, that is at            
C         A infinitely far from the BC diatomic, BC diatomic at its             
C         equilibrium configuration, is set equal to zero.                      
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),            
     *                   DO4Z(3), B(3)                                          
         DIMENSION X(3), COUL(3), EXCH(3)                                       
         PARAMETER (R2 = 1.41421356D0)                                          
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Initialize variables                                                        
C                                                                               
         ENGYGS = 0.D0                                                          
         DEGSDR(1) = 0.D0                                                       
         DEGSDR(2) = 0.D0                                                       
         DEGSDR(3) = 0.D0                                                       
C                                                                               
C   Check the value of NDER                                                     
C                                                                               
         IF (NDER .GT. 1) THEN                                                  
             WRITE (NFLAG(18), 900) NDER                                        
             STOP 'POT 1'                                                       
         ENDIF                                                                  
C                                                                               
C   Compute the energy.                                                         
C                                                                               
         DO 10 I = 1,3                                                          
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))                             
               COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)                    
               EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)                    
               ENGYGS  = ENGYGS + COUL(I)                                       
10       CONTINUE                                                               
C                                                                               
         RAD = SQRT((EXCH(1)-EXCH(2))**2 + (EXCH(2)-EXCH(3))**2 +               
     *              (EXCH(3)-EXCH(1))**2)                                       
C                                                                               
         ENGYGS = ENGYGS - RAD/R2 + EZERO(1)                                    
C                                                                               
C   Compute the derivatives of the energy with respect to the internal          
C   coordinates.                                                                
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             S = EXCH(1) + EXCH(2) + EXCH(3)                                    
             DO 20 I = 1,3                                                      
                   DEGSDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/R2*                 
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-                         
     *                       ZP3(I)*X(I)+OP3Z(I))                               
20           CONTINUE                                                           
         ENDIF                                                                  
C                                                                               
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,             
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',             
     *                'coded in this potential')                                
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
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
         DATA DE/ 141.196D0, 109.449D0, 141.196D0/                              
         DATA RE/ 0.917D0, 0.7419D0, 0.917D0/                                   
         DATA BETA/ 2.2087D0, 1.942D0, 2.2087D0/                                
         DATA Z/ 0.167D0, 0.106D0, 0.167D0/                                     
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
