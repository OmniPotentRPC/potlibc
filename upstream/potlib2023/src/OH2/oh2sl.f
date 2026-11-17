C
      SUBROUTINE PREPOT
C
C   System:          OH2
C   Functional form: 
C   Common name:     OH2SL
C   Reference:       R. Schinke and W. A. Lester, Jr., J. Chem. Phys., Vol. 70, 4893 (1979).
C
C   Notes: This is a collinear potential energy surface, i.e., the potential
C          energy is not defined for non-linear geometries.
C          The derivatives of the potential energy with respect to the 
C          coordinates are determined numerically using a step size STEP
C          if R1 and R2 are less than 50 bohr.
C          The variable STEP is initialized in the block data subprogram
C          PPARMS.
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in DATA statements
C   and in the block data subprogram PPARMS.
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
C      /PT4CM/ IPTPRT, IDUM(19)
C        IPTPRT passes the FORTRAN unit number used for potential output
C        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to
C        Re(AB), the equilibrium bond length for the AB diatomic.
C        Similarly, the terms EASYBC and EASYAC represent the energies
C        in the BC and the AC valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1
C      IPTPRT           6
C
C*****
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        IMPLICIT REAL*16(L)
         COMMON /PT1CM/R1,R2,R3,E,D(3)
         COMMON /PT2CM/ NSURF, NDER, NFLAG(21)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)
         COMMON /NDERCM/ STEP
C
      DIMENSION V(4)
      DIMENSION C(50),AL(3),X(3),S(3),VM(3),R(3),CHI(3)
      DATA A,B,ALF,GAM,AL/0.157002D-03,-0.162301D-01,0.9D0,0.045D0,             
     1   1.D0,0.8D0,1.D0/                                                       
      DATA BOHR/0.529177D0/,HART/27.21161D0/
      DATA C/0.354195D+01,0.302780D+01,0.577312D+00,0.539877D+01                
     1,0.233027D+01,0.157276D+02, 0.701851D+01,-0.478535D+02,                   
     1-0.106915D+01,0.268137D+02,-0.104806D+02,-0.182660D+02,                   
     1 0.833941D+01, 0.872497D+02, 0.917870D+01,-0.697766D+02,                  
     1 0.434862D+02,-0.156633D+02, 0.380162D+02,-0.796542D+01,                  
     1 0.602854D+02, 0.500800D+02,-0.545761D+02, 0.410210D+02,                  
     1-0.890739D+01, 0.191135D+02,-0.512653D+02,                                
     1 0.365096D+02,-0.235972D+02,-0.102322D+03, 0.837770D+00,                  
     1 0.258200D+02,-0.198105D+01,-0.484189D+01, 0.105377D+02,                  
     1 0.372352D+01,-0.207237D+02, 0.223342D+02, 0.244787D+01,                  
     1 0.979572D+01,-0.135612D+01,-0.724496D+01, 0.234984D+02,                  
     1 0.136364D+03, 0.253799D+01,-0.203821D+02, 0.699846D+01,                  
     1 0.198529D+02,-0.475902D+01, 0.236328D+01/                                
C
C   Echo potential data to the file linked to FORTRAN unit IPRT
      WRITE (IPRT, 600)
      WRITE (IPRT,601)STEP 
      WRITE (IPRT,602)
      J=1
C
C    Set the values of the classical energy in the three asymptotic valleys
             EASYAB = DE(1)/HART
             EASYBC = DE(2)/HART
             EASYAC = DE(3)/HART
C
      RETURN                                                                    
C                                                                               
C CALCULATE ENERGY AND DERIVATIVE                                               
C                                                                               
      ENTRY POT
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
      IC=4
      IF(R1.GT.30.D0)IC=1
      IF(R2.GT.30.D0)IC=1
   86 GO TO (84,82,83,81),IC
   81 SOH=R1
      SR2=R2
      R(3) = R3
      GO TO 85
82      SR2=R2+STEP                                                            
      SOH=R1                                                                 
      R(3)= R3
      GO TO 85
   83 SR2=R2                                                                 
      SOH=R1+STEP                                                            
      R(3)= R3
      GO TO 85
   84 SR2=R2
      SOH=R1
      R(3)=R3 + STEP
   85 SR2=SR2*0.529177D0                                                        
      SOH=SOH*0.529177D0                                                        
      R(3)= R(3) *0.529177D0
      R(1) = SOH
      R(2) = SR2
      Q = 1.D0                                                                  
      DO 10 I=1,3                                                               
      X(I)=R(I)                                                                 
      CHI(I)=(1.D0 - EXP(-ALF*X(I)**2))
      S(I)=X(I)-1.D0                                                            
      VM(I)=DE(I)*(1.D0-EXP(-BET(I)*(X(I)-XRE(I))))**2-DE(I)                    
   10 Q=Q*(1.D0-TANH(AL(I)*S(I)))                                              
      U=VM(1)*CHI(2)*CHI(3)+VM(2)*CHI(1)*CHI(3)                     
     1+VM(3)*CHI(1)*CHI(2)                                                
      BP=B*(1.D0/X(1)**4+1.0D0/X(3)**4)
      CHECK=(X(1)*X(3)*(X(2)-0.2D0))
      IF(CHECK.GT.1000.D0)AP=0
      IF(CHECK.GT.1000.D0)GO TO 112
      AP=A/(X(1)*X(3)*(X(2)-0.2D0))**12
  112 U=U + AP + BP
      T=C(1)+C(2)*(S(1)+S(3))+C(3)*S(2)+C(4)*(S(1)**2+S(3)**2)                  
     1+C(5)*S(2)**2+C(6)*S(2)*(S(1)+S(3))+C(7)*S(1)*S(3)+                       
     1C(8)*(S(1)**3+S(3)**3)+C(9)*S(2)**3+C(10)*S(1)*S(3)*(S(1)+S(3))
C
      TT=C(11)*S(2)**2*(S(1)+S(3))+C(12)*S(2)*(S(1)**2+S(3)**2)         
     1+C(13)*S(1)*S(2)*S(3)+C(14)*(S(1)**4+S(3)**4)                             
     1+C(15)*S(2)**4+C(16)*(S(1)**3*S(3)+S(1)*S(3)**3)                          
C
      WW=C(17)*S(1)**2*S(3)**2+C(18)*S(2)*(S(1)**3+S(3)**3)
     1+C(19)*S(2)**2*(S(1)**2+S(3)**2)+C(20)*S(2)**3*(S(1)+                     
     1S(3))+C(21)*S(2)**2*S(1)*S(3)+C(22)*S(2)*(S(1)**2*S(3)                    
     1+S(1)*S(3)**2)
C
      W=C(23)*(S(1)**5+S(3)**5)+C(24)*(S(1)**4*S(3)                
     1+S(1)*S(3)**4)+C(25)*(S(1)**3*S(3)**2+S(1)**2*S(3)**3)                    
     1+C(26)*S(2)*(S(1)**4+S(3)**4)+C(27)*S(2)*(S(1)**3*S(3)                    
     1+S(1)*S(3)**3)
C
      GG=C(28)*S(2)*S(1)**2*S(3)**2+C(29)*                          
     1S(2)**2*(S(1)**3+S(3)**3)+C(30)*S(2)**2*(S(1)**2*S(3)                     
     1+S(1)*S(3)**2)+C(31)*S(2)**3*(S(1)**2+S(3)**2)                            
C
      Y=C(32)*S(2)**3*S(1)*S(3)+C(33)*S(2)**4*(S(1)+S(3)) 
     1+C(34)*S(2)**5+C(35)*(S(1)**6+S(3)**6)+C(36)*                             
     1(S(1)**5*S(3)+S(1)*S(3)**5)+C(37)*(S(1)**4                                
     1*S(3)**2+S(1)**2*S(3)**4)+C(38)*S(1)**3*S(3)**3
C
      Z=C(39)*S(2)*(S(1)**5+S(3)**5)+C(40)*S(2)*(S(1)**4*S(3)    
     1+S(1)*S(3)**4)+C(41)*S(2)*(S(1)**3*S(3)**2+S(1)**2*S(3)**3)               
     1+C(42)*S(2)**2*(S(1)**4+S(3)**4)+C(43)*S(2)**2*                           
     1(S(1)**3*S(3)+S(1)*S(3)**3)+C(44)*S(1)**2*S(2)**2*S(3)**2                 
C
      ZZ=C(45)*S(2)**3*(S(1)**3+S(3)**3)+C(46)*S(2)**3*(S(1)**2*S(3) 
     1+S(1)*S(3)**2)+C(47)*S(2)**4*(S(1)**2+S(3)**2)+C(48)*S(2)**4*S(1)*        
     1   S(3)                                                                   
     1+C(49)*S(2)**5*(S(1)+S(3))+C(50)*S(2)**6                                  
      P=T+TT+WW+W+GG+Y+Z+ZZ
      TEMP = EXP(-GAM*(X(1)+X(3))**2)                                          
      U=U+P*Q*EXP(-GAM*(X(1)+X(3))**2)                                         
      XTER=U + DE(2)
      V(IC) = XTER/HART
      IF(IC .LT. 4)GO TO 97
      IF(R1.GT.30.D0)GO TO 77
      IF(R2.LT.30.D0)GO TO 66
77      D(1)=(V(3)-V(4))/STEP
       D(2)=(V(2)-V(4))/STEP
      D(3)=(V(1)-V(4))/STEP
98      IF(R(2).LT.20)GO TO 97
      V(4)=V(4)*1.0D8
      V(4)=INT(V(4))
      V(4)=V(4)*1.0D-8
   97 IC=IC+1                                                                   
      IF(IC.LT.5)GO TO 86
      E=V(4)
      GO TO 113
66    DA=C(2) + 2.0D0*S(1)*C(4) + S(2)*C(6) +S(3)*C(7) +3.0D0*S(1)**2*
     1C(8) + (S(3)**2 + 2.0D0*S(1)*S(3))*C(10) + S(2)**2*C(11) + 2.0D0
     1*S(1)*S(2)*C(12) + S(2)*S(3)*C(13) + 4.0D0*S(1)**3*C(14)
      DB=(3.0D0*S(1)**2*S(3) + S(3)**3)*C(16) + 2.0D0*S(1)*S(3)**2*
     1C(17) + 3.0D0*S(1)**2*S(2)*C(18) + 2.0D0*S(1)*S(2)**2*C(19) +
     1S(2)**3*C(20) + S(2)**2*S(3)*C(21) + (2.0D0*S(1)*S(2)*S(3) +
     1S(2)*S(3)**2)*C(22)
      DC=5.0D0*S(1)**4*C(23) + (4.0D0*S(1)**3*S(3) + S(3)**4)*C(24) +
     1(3.0D0*S(1)**2*S(3)**2 + 2.0D0*S(1)*S(3)**3)*C(25) + 4.0D0*S(1)**
     13*S(2)*C(26) + (3.0D0*S(2)*S(1)**2*S(3) + S(2)*S(3)**3)*C(27)
      DD=2.0D0*S(1)*S(2)*S(3)**2*C(28) + 3.0D0*S(1)**2*S(2)**2*C(29) +
     1(2.0D0*S(1)*S(2)**2*S(3) + S(2)**2*S(3)**2)*C(30) + 2.0D0*S(1)* 
     1S(2)**3*C(31) + S(2)**3*S(3)*C(32) +S(2)**4*C(33) + 6.0D0*S(1)**
     15*C(35) + (5.0D0*S(1)**4*S(3) + S(3)**5)*C(36)
      DF=(4.0D0*S(1)**3*S(3)**2 + 2.0D0*S(1)*S(3)**4)*C(37) + 3.0D0*
     1S(1)**2*S(3)**3*C(38) + 5.0D0*S(1)**4*S(2)*C(39) + (4.0D0*S(2)*
     1S(1)**3*S(3) + S(2)*S(3)**4)*C(40) + (3.0D0*S(1)**2*S(2)*S(3)**2
     1+ 2.0D0*S(1)*S(2)*S(3)**3)*C(41) + 4.0D0*S(1)**3*S(2)**2*C(42)
      DG=(3.0D0*S(1)**2*S(2)**2*S(3) + S(2)**2*S(3)**3)*C(43) +
     12.0D0*S(1)*S(2)**2*S(3)**2*C(44) + 3.0D0*S(1)**2*S(2)**3*C(45) +
     1(2.0D0*S(1)*S(2)**3*S(3) + S(2)**3*S(3)**2)*C(46) + 2.0D0*S(1)*
     1S(2)**4*C(47) + S(2)**4*S(3)*C(48) + S(2)**5*C(49)
      D1D=DA+DB+DC+DD+DF+DG
      DA2=C(3) + 2.0D0*S(2)*C(5) + (S(1) + S(3))*C(6) + 3.0D0*S(2)**2*
     1C(9) + 2.0D0*S(2)*(S(1) + S(3))*C(11) + (S(1)**2 + S(3)**2)*C(12) 
     1+ S(1)*S(3)*C(13) + 4.0D0*S(2)**3*C(15) + (S(1)**3 + S(3)**3)*
     1C(18) + 2.0D0*S(2)*(S(1)**2 + S(3)**2)*C(19)
      DB2=3.0D0*S(2)**2*(S(1) + S(3))*C(20) + 2.0D0*S(1)*S(2)*S(3)*
     1C(21) + (S(1)**2*S(3) + S(1)*S(3)**2)*C(22) + (S(1)**4 + S(3)**
     14)*C(26) + (S(1)**3*S(3) + S(1)*S(3)**3)*C(27) + S(1)**2*S(3)**
     12*C(28) + 2.0D0*S(2)*(S(1)**3 + S(3)**3)*C(29)
      DC2=2.0D0*S(2)*(S(1)**2*S(3) + S(1)*S(3)**2)*C(30) + 3.0D0*S(2)
     1**2*(S(1)**2 + S(3)**2)*C(31) + 3.0D0*S(2)**2*S(1)*S(3)*C(32) +
     14.0D0*S(2)**3*(S(1) + S(3))*C(33) + 5.0D0*S(2)**4*C(34)
      DD2=(S(1)**5 + S(3)**5)*C(39) + (S(1)**4*S(3) + S(1)*S(3)**4)*
     1C(40) + (S(1)**3*S(3)**2 + S(1)**2*S(3)**3)*C(41) + 2.0D0*S(2)*
     1(S(1)**4 + S(3)**4)*C(42) + 2.0D0*S(2)*(S(1)**3*S(3) + S(1)*S(3)
     1**3)*C(43) + 2.0D0*S(1)**2*S(2)*S(3)**2*C(44)
      DF2=3.0D0*S(2)**2*(S(1)**3 + S(3)**3)*C(45) + 3.0D0*S(2)**2*
     1(S(1)**2*S(3) + S(1)*S(3)**2)*C(46) + 4.0D0*S(2)**3*(S(1)**2 +
     1S(3)**2)*C(47) + 4.0D0*S(2)**3*S(1)*S(3)*C(48) + 5.0D0*S(2)**4
     1*(S(1) + S(3))*C(49) + 6.0D0*S(2)**5*C(50)
      D2D=DA2 +DB2+DC2+DD2+DF2
      DA3=C(2) + 2.0D0*S(3)*C(4) + S(2)*C(6) + S(1)*C(7) + 3.0D0*S(3)**
     12*C(8) + (2.0D0*S(1)*S(3) + S(1)**2)*C(10) + S(2)**2*C(11) + 2.0D0
     1*S(2)*S(3)*C(12) + S(1)*S(2)*C(13) + 4.0D0*S(3)**3*C(14)
      DB3=(S(1)**3 + 3.0D0*S(1)*S(3)**2)*C(16) +2.0D0*S(1)**2*S(3)*
     1C(17) + 3.0D0*S(2)*S(3)**2*C(18) + 2.0D0*S(2)**2*S(3)*C(19) +
     1S(2)**3*C(20) + S(2)**2*S(1)*C(21) + S(2)*(S(1)**2 + 2.0D0*S(1)*
     1S(3))*C(22) + 5.0D0*S(3)**4*C(23) + (S(1)**4 + 4.0D0*S(1)*S(3)**
     13)*C(24)
      DC3=(2.0D0*S(1)**3*S(3) + 3.0D0*S(1)**2*S(3)**2)*C(25) + 4.0D0*
     1S(2)*S(3)**3*C(26) + S(2)*(S(1)**3 + 3.0D0*S(1)*S(3)**2)*C(27) +
     12.0D0*S(2)*S(1)**2*S(3)*C(28) + 3.0D0*S(2)**2*S(3)**2*C(29)
      DD3=S(2)**2*(S(1)**2 + 2.0D0*S(1)*S(3))*C(30) + 2.0D0*S(2)**3*S(3)
     1*C(31) + S(2)**3*S(1)*C(32) + S(2)**4*C(33) + 6.0D0*S(3)**5*C(35)
     1+(S(1)**5 + 5.0D0*S(1)*S(3)**4)*C(36) + (2.0D0*S(1)**4*S(3) + 
     14.0D0*S(1)**2*S(3)**3)*C(37) + 3.0D0*S(1)**3*S(3)**2*C(38) 
      DF3=5.0D0*S(2)*S(3)**4*C(39) + S(2)*(S(1)**4 + 4.0D0*S(1)*S(3)**3)
     1*C(40) + S(2)*(2.0D0*S(1)**3*S(3) + 3.0D0*S(1)**2*S(3)**2)*C(41) +
     14.0D0*S(2)**2*S(3)**3*C(42) + S(2)**2*(S(1)**3 + 3.0D0*S(1)*S(3)
     1**2)*C(43) + 2.0D0*S(1)**2*S(2)**2*S(3)*C(44) + 3.0D0*S(2)**3*
     1S(3)**2*C(45) + S(2)**3*(S(1)
     1**2 + 2.0D0*S(1)*S(3))*C(46) + 2.0D0*S(2)**4*S(3)*C(47) +S(2)**4
     1*S(1)*C(48) + S(2)**5*C(49)
      D3D=DA3+DB3+DC3+DD3+DF3
      DQ1=(1.0D0-TANH(AL(2)*S(2)))*(1.0D0-TANH(AL(3)*S(3)))
      DQ2=(1.0D0-TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(3)*S(3)))
      DQ3=(1.0D0-TANH(AL(1)*S(1)))*(1.0D0-TANH(AL(2)*S(2)))
      DAP1=-12.0D0*AP/R(1)
      DAP2=-12.0D0*AP/(R(2)-0.2D0)
      DAP3=-12.0D0*AP/R(3)
      DB1=-4.0D0*B/R(1)**5
      DB3=-4.0D0*B/R(3)**5
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)DCRS1=1.0D30
C      IF(DCOSH(AL(1)*S(1)).LT.1.0D-15)GO TO 144
      DCSR1=-AL(1)/(COSH(AL(1)*S(1))**2)
C144      IF(COSH(AL(2)*S(2)).LT.1.0D-15)DCRS2=1.0D30
C      IF(COSH(AL(2)*S(2)).LT.1.0D-15)GO TO 145
      DCSR2=-AL(2)/(COSH(AL(2)*S(2))**2)
C145      IF(COSH(AL(3)*S(3)).LT.1.0D-15)DCRS3=1.030
C      IF(COSH(AL(3)*S(3)).LT.1.0D-15)GO TO 146
      DCSR3=-AL(3)/(COSH(AL(3)*S(3))**2)
146      PSI13=2.0D0*(-GAM)*TEMP*(R(1) + R(3))
      DCHI1=2.0D0*ALF*R(1)*EXP(-ALF*(R(1)**2))
      DCHI2=2.0D0*ALF*R(2)*EXP(-ALF*(R(2)**2))
      DCHI3=2.0D0*ALF*R(3)*EXP(-ALF*(R(3)**2))
C      IF(VM(1).GT.1.0D20)DVR1=1.0D30
C      IF(VM(1).GT.1.0D20)GO TO 153
      DVR1=2.0D0*DE(1)*(1.0D0-EXP(-BET(1)*(R(1)-XRE(1))))*BET(1)*EXP(
     1-BET(1)*(R(1)-XRE(1)))
C 153  IF(VM(2).GT.1.0D20)DVR2=1.0D30
C      IF(VM(2).GT.1.0D20)GO TO 154
      DVR2=2.0D0*DE(2)*(1.0D0-EXP(-BET(2)*(R(2)-XRE(2))))*BET(2)*EXP(
     1-BET(2)*(R(2)-XRE(2)))
C 154  IF(VM(3).GT.1.0D20)DVR3=1.0D30
C      IF(VM(3).GT.1.0D20)GO TO 155
      DVR3=2.0D0*DE(3)*(1.0D0-EXP(-BET(3)*(R(3)-XRE(3))))*BET(3)*EXP(
     1-BET(3)*(R(3)-XRE(3)))
155   D(1)=VM(2)*CHI(3)*DCHI1 + CHI(2)*CHI(3)*DVR1 + VM(3)*CHI(2)*DCHI1
     1+ D1D*Q*TEMP + P*DQ1*TEMP*DCSR1 + DAP1 + DB1 + P*Q*PSI13
      D(2)=DVR2*CHI(1)*CHI(3) + VM(1)*CHI(3)*DCHI2 + VM(3)*CHI(1)*DCHI2
     1+ D2D*Q*TEMP + P*DQ2*TEMP*DCSR2 + DAP2
      D(3)=VM(2)*CHI(1)*DCHI3 + VM(1)*CHI(2)*DCHI3 + DVR3*CHI(1)*CHI(2)
     1+ D3D*Q*TEMP + P*DQ3*TEMP*DCSR3 + P*Q*PSI13 + DAP3 + DB3
      DO 222 I=1,3
222   D(I)=D(I)/HART*BOHR
      J=J+1
      GO TO 98
  113 RETURN
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the OH2 ',
     *                'Schinke-Lester potential energy surface')
601   FORMAT(/,2X,T5,'Step size for the numerical derivatives = ',E10.5) 
602   FORMAT(2X,T5,'The Schinke-Lester surface has analytical ',
     *             'derivatives if R1 and R2 are ',
     *       /,2X,T5,'less than 50 bohr, otherwise the derivatives ',
     *               'are computed numerically.')
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
       END                                                                   
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /SATOCM/ DE(3), XRE(3), BET(3)
         COMMON /NDERCM/ STEP
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NDER /1/
         DATA NSURF /0/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters; the energy parameters are in
C   eV and the lengths are in Angstroms.
         DATA DE/4.035D0,4.1662D0,4.035D0/
         DATA XRE/0.9923D0,0.7590D0,0.9923D0/
         DATA BET/2.4113D0,2.0682D0, 2.4113D0/  
C   Initialize the step size used for the numerical derivatives of the
C   energy with respect to the coordinates.
         DATA STEP/5.0D-7/
C
         END
C*****
