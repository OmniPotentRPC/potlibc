      SUBROUTINE PREPOT                             
C
C   System:          H3
C   Functional Form: Double many-body expansion
C   Common Name:     DMBE
C   Interface:       3-1V
C   Number of derivatives: 0
C   References:      A. J. C. Varandas, F. B. Brown, C. A. Mead, 
C                    D. G. Truhlar, and N. C. Blais
C                    J. Chem. Phys. 86, 6258 (1987).
C    B. Liu & P. Siegbahn, J. Chem. Phys., 68, 2457, 1978;
C    D.G. Truhlar & C.J. Horowitz, J. Chem. Phys., 68, 2466, 1978;
C    A.J.C. Varandas, Mol. Phys., 53, 1303, 1984;
C    A.J.C. Varandas, J. Mol. Struct. (Theochem), 21, 401, 1985;
C
C    Note: This is a vectorized version of H3DMBE surface, no derivatives
C 
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy 
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohrs
C      derivatives - hartrees/bohr
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the H
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Coordinates:
C      Internal, Definition: R(1) = R(first H-second H)
C                            R(2) = R(second H-third H)
C                            R(3) = R(third H-first H)
C
C***********************************************************************
C  Calculates double many-body expansion of the H3 potential.
C        POTE1 - energy of surface 1.
C        POTE2 - energy of surface 2.
C        COUP  - nonadiabatic coupling.  
C        Note:  COUP not yet coded.
C  Modified for vectorized L**2 reactive code by M.Mladenovic 11/10/87
     IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXP=100) 
      COMMON/TIKI/ ITIKI, ITIKIP
      DIMENSION VPOTE1(NI), RR(NI,3)
C             
C  Common block for H3 potential parameters, set in BLOCK DATA H3
C                                                          
      DIMENSION CHH(3), CK0(3), CK1(3), DAMPA(3), DAMPB(3), XPAR(15)
C      COMMON /H3DMCM/ ALPH2, ALPHA5, ALPH0, ALPH1, AL0, AL1, AL2, AL3,
C     *   AZ2, BETA1, BETA2, BETA3, BET0, BET1, CD0, CD1, CHH(3), CK0(3),
C     *   CK1(3), DAMPA(3), DAMPB(3), HFD, HFA1, HFA2, HFA3, HFGAM, H2RM,
C     *   H2RMT, H2R0, H2TA, H2TB1, H2TB2, H2TB3, H2TB4, PI, RHOL, SQRT3,
C     *   XPAR(15)
      DIMENSION R(MAXP,3)              
      DIMENSION PER(MAXP), PER2(MAXP), QCOORD(MAXP), RHO(MAXP),  
     *S(MAXP), S2(MAXP), CPHI3(MAXP),GAMMA(MAXP),BETA(MAXP),GAM2(MAXP),
     *BET2(MAXP),T12(MAXP),T1(MAXP),AS(MAXP),T2(MAXP),BR(MAXP)
C                          
C  Common block for 2-body correlation energies and damping factors
C  (computed in H3COR2)
C     
      DIMENSION A1(MAXP), A2(MAXP), A3(MAXP)
      DIMENSION AT1(MAXP), AT2(MAXP), AT3(MAXP)
      DIMENSION B1(MAXP,3), B2(MAXP,3), B3(MAXP,3)
      DIMENSION CE2(MAXP), CE3(MAXP)
      DIMENSION C13(MAXP), C23(MAXP), C33(MAXP)
      DIMENSION C113(MAXP), C213(MAXP), C313(MAXP)
      DIMENSION C123(MAXP), C223(MAXP), C323(MAXP)
      DIMENSION C133(MAXP), C233(MAXP), C333(MAXP)
      DIMENSION COS1(MAXP), COS2(MAXP), COS3(MAXP)
      DIMENSION CORRS(MAXP,3), CORRT(MAXP,3), DAMP(MAXP,3,3) 
      DIMENSION DR1(MAXP), DR2(MAXP), DR3(MAXP), D1(MAXP), D2(MAXP)
      DIMENSION DR12(MAXP), DR23(MAXP), DR31(MAXP)
      DIMENSION DTS1(MAXP,3), DTS2(MAXP,3), DTS3(MAXP,3) 
      DIMENSION DTT1(MAXP,3), DTT2(MAXP,3), DTT3(MAXP,3)
      DIMENSION EB1T(MAXP), EB3T(MAXP)                     
      DIMENSION EB4A(MAXP), EB4B(MAXP), T9(MAXP), T10(MAXP)
      DIMENSION ENT1(MAXP), ENT2(MAXP), ENT3(MAXP)
      DIMENSION EQ(MAXP), EQP1(MAXP), EQP2(MAXP), EQP3(MAXP) 
      DIMENSION ES1(MAXP), ES2(MAXP), ES3(MAXP)        
      DIMENSION ET1(MAXP), ET2(MAXP), ET3(MAXP) 
      DIMENSION EXCH(MAXP), VLEPS(MAXP), VLEPS2(MAXP)
      DIMENSION EXP1(MAXP), EXP2(MAXP), EXP3(MAXP), EXPV3(MAXP)
      DIMENSION F1(MAXP), F11(MAXP), F2(MAXP), F22(MAXP)
      DIMENSION F3(MAXP), F4(MAXP), FF(MAXP)
      DIMENSION G1(MAXP,3), G2(MAXP,3), G3(MAXP,3)                    
      DIMENSION H1(MAXP,3), H2(MAXP,3), H3(MAXP,3)
      DIMENSION HYTAN1(MAXP,3), HYTAN2(MAXP,3), HYTAN3(MAXP,3)
      DIMENSION OMF(MAXP)
      DIMENSION POMS1(MAXP,3), POMS2(MAXP,3), POMS3(MAXP,3) 
      DIMENSION POMT1(MAXP,3), POMT2(MAXP,3), POMT3(MAXP,3)  
      DIMENSION POTE1(MAXP), POTE2(MAXP)
      DIMENSION QSUM(MAXP)
      DIMENSION R12(MAXP), R22(MAXP), R32(MAXP)
      DIMENSION RI(MAXP,3), RRI(MAXP)
      DIMENSION RI6(MAXP,3), RI8(MAXP,3), RI10(MAXP,3)
      DIMENSION RVA11(MAXP), RVA12(MAXP), RVA13(MAXP), RVA14(MAXP)    
      DIMENSION SGNT1(MAXP,3), SGNT2(MAXP,3), SGNT3(MAXP,3)
      DIMENSION T4(MAXP), T5(MAXP), T6(MAXP), T7(MAXP), T8(MAXP)
      DIMENSION TDP1(MAXP,3), TDP2(MAXP,3), TDP3(MAXP,3)
      DIMENSION TF1(MAXP), TF2(MAXP), TF3(MAXP)
      DIMENSION TH1(MAXP,3), TH2(MAXP,3), TH3(MAXP,3)
      DIMENSION TTS1(MAXP,3), TTS2(MAXP,3), TTS3(MAXP,3) 
      DIMENSION TTT1(MAXP,3), TTT2(MAXP,3), TTT3(MAXP,3)
      DIMENSION TW1(MAXP), TW2(MAXP), TW3(MAXP)
      DIMENSION TT(MAXP), RVA(MAXP), EXPVA(MAXP), RVA1(MAXP), VA(MAXP)
      DIMENSION VII(MAXP), VIII(MAXP)
      DIMENSION VC(MAXP), VAN(MAXP) 
      DIMENSION WB(MAXP), WB2(MAXP), EB1(MAXP), EB2(MAXP), EB4(MAXP)
      DIMENSION WE1(MAXP), WE2(MAXP), WE3(MAXP)
      DIMENSION XJ1(MAXP), XJ2(MAXP), XJ3(MAXP)
      DIMENSION XT(MAXP,3), X(MAXP,3) 
      DIMENSION ZT1(MAXP), ZT2(MAXP), ZT3(MAXP)
C     COMMON /POTCM/ Q(3),EE,DQ(3)
      SAVE RHOL, ALPH0, ALPH1, BET0, BET1, ALPHA5, BETA1, BETA2, BETA3,
     *ALPH2, AL0, AL1, AL2, AL3, AZ2, CD0, CD1, CHH, CK0, CK1, HFD, 
     *HFA1, HFA2, HFA3, HFGAM, H2RM, H2R0, H2RMT, H2TA,H2TB1, H2TB2, 
     *H2TB3, H2TB4, PI, SQRT3, XPAR, DAMPA, DAMPB, DASY
      DATA RHOL /2.4848D0/
      DATA ALPH0, ALPH1, BET0, BET1 /25.9528D0, 1.1868D0, 15.7381D0,
     *   0.09729D0/
      DATA ALPHA5, BETA1, BETA2, BETA3 / 8.2433033D-03, 0.53302897D0,
     *   0.39156612D-01, 0.69996945D0/
      DATA ALPH2 /4.735364D-1/
      DATA AL0, AL1, AL2, AL3 / 0.45024472D+01, -0.62467617D+01,
     *   0.40966542D+01, 0.21813012D+01/
      DATA AZ2 /4.453649D-4/
      DATA CD0, CD1 /6.333404D-03, -1.726839D-03/
      DATA CHH /  6.499027D0,      1.243991D+02,  3285.8D0/
      DATA CK0 / -1.372843D-01, -1.638459D-01, -1.973814D-01/
      DATA CK1 /  1.011204D0,      9.988099D-01,  9.399411D-01/
      DATA HFD, HFA1, HFA2, HFA3, HFGAM /0.15796326D0, 2.1977034D0,
     *   1.2932502D0, 0.64375666D0, 2.1835071D0/
      DATA H2RM, H2R0, H2RMT /1.401D0, 6.928203D0, 7.82D0/
      DATA H2TA, H2TB1, H2TB2, H2TB3, H2TB4 /0.448467D0, -0.056687D0,
     *   0.246831D0, -0.018419D0, 0.000598D0/
      DATA PI /3.14159265358979D0/
      DATA SQRT3 /1.73205080756887D0/
      DATA XPAR /-0.9286579D-02,  0.2811592D-03, -0.4665659D-05,
     *            0.2069800D-07,  0.2903613D+02, -0.2934824D+01,
     *            0.7181886D0,     -0.3753218D0,     -0.1114538D+01,
     *            0.2134221D+01, -0.4164343D0,      0.2022584D0,
     *           -0.4662687D-01, -0.4818623D+02,  0.2988468D0/
C
C  Common block for coordinates (computed in H3COOR)
C                   
      DATA DASY /0.174474112D0/
C       DATA DASY /0.0/ 
      NEXP = 4
      DO 10 ID = 1,3
         NEXP = NEXP + 2                   
         DAMPA(ID) = ALPH0/FLOAT(NEXP)**ALPH1
         DAMPB(ID) = BET0*EXP(-BET1*FLOAT(NEXP))
   10 CONTINUE
      RETURN
C
      ENTRY POT(RR,VPOTE1,NI,NSURF)                                          
      ITIKI = ITIKI + 1
      ITIKIP = ITIKIP + NI
      RNULA=1.0E-22
      NIM=NI/MAXP
      IF(NIM*MAXP.EQ.NI) NIM=NIM-1
      NIM=NIM+1
      II=0
      DO 2000 N=1,NIM
         NP=MIN0(NI-II,MAXP)
            DO 15 I=1,3
         DO 20 IN=1,NP 
               R(IN,I)=RR(IN+II,I) 
   20    CONTINUE
   15       CONTINUE
C     
C Set up symmetry coordinates    
C     
C Set up 2-body correlation energies and dispersion damping terms   
C
C Get Leps potentials                            
C
C Get VA term                     
C
C Get VII term 
C
C Get VIII term 
C
C Get 3-body correlation            
C
C 2-body correlation term
C                  
C
C
C      SUBROUTINE H3COOR (R1, R2, R3)
C**********************************************************************
C   Calculates D3H symmetry coordinates 
C**********************************************************************
C                               
      DO 110 I=1,NP
         PER(I)=R(I,1)+R(I,2)+R(I,3)
         R12(I)=R(I,1)*R(I,1)
         R22(I)=R(I,2)*R(I,2)
         R32(I)=R(I,3)*R(I,3)
         RI(I,1) = 1.0D0/R(I,1)
         RI(I,2) = 1.0D0/R(I,2)
         RI(I,3) = 1.0D0/R(I,3)
  110 CONTINUE 
      DO 120 I=1,NP
         PER2(I)=PER(I)*PER(I)
         QCOORD(I) = R12(I) + R22(I) + R32(I) 
         GAMMA(I)=2.0D0*R12(I)-R22(I)-R32(I) 
         BETA(I) = SQRT3*(R22(I)-R32(I)) 
  120 CONTINUE
      DO 130 I=1,NP
         RHO(I) = SQRT(QCOORD(I)/3.0D0)
         GAM2(I) = GAMMA(I) * GAMMA(I)
         BET2(I) = BETA(I) * BETA(I)
  130 CONTINUE 
      DO 140 I=1,NP   
         T12(I) = BET2(I) + GAM2(I)
  140 CONTINUE
      DO 150 I=1,NP
         T1(I) = SQRT(T12(I))
  150 CONTINUE
      DO 160 I=1,NP
         S(I) = T1(I)/QCOORD(I)
  160 CONTINUE     
      DO 170 I=1,NP 
         AS(I)=ABS(S(I))
         S2(I) = S(I) * S(I)
         T2(I) = T1(I)*T12(I)
         BR(I) = GAMMA(I)*(3.0D0*BET2(I) - GAM2(I))
  170 CONTINUE
C      DO 180 I=1,NPC                                                   
C          BR(I) = CVMGT(0.0,BR(I),AS(I).LT.RNULA)
C         T2(I) = CVMGT(1.0,T2(I),AS(I).LT.RNULA)
C  180 CONTINUE                       
C      DO 190 I=1,NP                            
C         CPHI3(I) = BR(I)/T2(I)
C  190 CONTINUE
C                      
      DO 180 I=1,NP                
         IF (S(I) .EQ. 0.0D0) THEN 
C
C   For S=0, CPHI3 should not be used anywhere but set to zero anyway.
C                 
         CPHI3(I) = 0.0D0
      ELSE
         CPHI3(I) = GAMMA(I)*(3.0D0*BET2(I) - GAM2(I))/T2(I)
      END IF                                       
  180 CONTINUE 
C
C      SUBROUTINE H3COR2 (R1,R2,R3) 
C***********************************************************************
C     Calculates 2-body correlation energies for singlet and triplet
C     states of H2, the damping factors for the dispersion terms       
C***********************************************************************
C
C     Loop over three coordinates
C
C     Loop over terms in dispersion expansion                         
C
C     singlet
C     Calculate damping function for the nth dispersion coefficient 
C
C     store damping factors (including CHH coeff and 1/R**NEXP) for
C     later use in computing the 3-body correlation energy.
C
C     triplet
C     Calculate damping function for the nth dispersion coefficient
C
      DENOMS = H2RM + 2.5D0*H2R0
      DENOMT = H2RMT + 2.5D0*H2R0
         DO 410 J=1,3 
      DO 420 I=1,NP
            CORRS(I,J) = 0.0D0
            CORRT(I,J) = 0.0D0    
            X(I,J) = 2.0D0*R(I,J)/DENOMS
            XT(I,J) = 2.0D0*R(I,J)/DENOMT
  420 CONTINUE
  410    CONTINUE
         DO 430 J=1,3
      DO 440 I=1,NP
            TTS1(I,J)=X(I,J)*(DAMPA(1)+X(I,J)*DAMPB(1))
            TTS2(I,J)=X(I,J)*(DAMPA(2)+X(I,J)*DAMPB(2))
            TTS3(I,J)=X(I,J)*(DAMPA(3)+X(I,J)*DAMPB(3)) 
            TTT1(I,J)=XT(I,J)*(DAMPA(1)+XT(I,J)*DAMPB(1))
            TTT2(I,J)=XT(I,J)*(DAMPA(2)+XT(I,J)*DAMPB(2))
            TTT3(I,J)=XT(I,J)*(DAMPA(3)+XT(I,J)*DAMPB(3)) 
  440 CONTINUE
  430    CONTINUE
         DO 450 J=1,3
      DO 460 I=1,NP
            RI6(I,J)=RI(I,J)**6
            RI8(I,J)=RI(I,J)**8
            RI10(I,J)=RI(I,J)**10
            DTS1(I,J)=(1.0D0-EXP(-TTS1(I,J)))**6
            DTS2(I,J)=(1.0D0-EXP(-TTS2(I,J)))**8
            DTS3(I,J)=(1.0D0-EXP(-TTS3(I,J)))**10
            DTT1(I,J)=(1.0D0-EXP(-TTT1(I,J)))**6
            DTT2(I,J)=(1.0D0-EXP(-TTT2(I,J)))**8
            DTT3(I,J)=(1.0D0-EXP(-TTT3(I,J)))**10
  460 CONTINUE 
  450 CONTINUE
         DO 470 J=1,3            
      DO 480 I=1,NP
            POMS1(I,J)=RI6(I,J)*DTS1(I,J) 
            POMS2(I,J)=RI8(I,J)*DTS2(I,J)
            POMS3(I,J)=RI10(I,J)*DTS3(I,J)
            POMT1(I,J)=RI6(I,J)*DTT1(I,J) 
            POMT2(I,J)=RI8(I,J)*DTT2(I,J)
            POMT3(I,J)=RI10(I,J)*DTT3(I,J)
  480 CONTINUE     
  470    CONTINUE
      DO 486 J=1,3 
         DO 485 I=1,NP
            DAMP(I,1,J)=CHH(1)*POMS1(I,J)
            DAMP(I,2,J)=CHH(2)*POMS2(I,J)
            DAMP(I,3,J)=CHH(3)*POMS3(I,J)
  485    CONTINUE
  486 CONTINUE
       DO 490 J=1,3                                                   
      DO 500 I=1,NP
       CORRS(I,J) = - DAMP(I,1,J) - DAMP(I,2,J) - DAMP(I,3,J)
       CORRT(I,J)=-CHH(1)*POMT1(I,J)-CHH(2)*POMT2(I,J)-CHH(3)*POMT3(I,J) 
  500 CONTINUE     
  490  CONTINUE 
C
C      SUBROUTINE H3COR3 (R, CE3) 
C***********************************************************************
C     Calculates 3-body correlation energy 
C***********************************************************************
C
C     Loop over terms in dispersion expansion; dispersion damping factors
C     are computed in H3COR2 and passed through COMMON /H3CRCM/.
C
C     Loop over 3 coordinates
C
C     Calculate G function                       
C
      DO 510 I=1,NP
         CE3(I) = 0.0D0
  510 CONTINUE
      DO 525 ID=1,3
      DO 520 I=1,NP
         G1(I,ID) = 1.0D0 + CK0(ID)*EXP(-CK1(ID)*(R(I,1)-H2RM)) 
         G2(I,ID) = 1.0D0 + CK0(ID)*EXP(-CK1(ID)*(R(I,2)-H2RM)) 
         G3(I,ID) = 1.0D0 + CK0(ID)*EXP(-CK1(ID)*(R(I,3)-H2RM)) 
         TH1(I,ID)=CK1(ID)*R(I,1)
         TH2(I,ID)=CK1(ID)*R(I,2)
         TH3(I,ID)=CK1(ID)*R(I,3)
  520 CONTINUE
  525 CONTINUE
C      DO 535 ID=1,3
C      DO 530 I=1,NP
C         TH1(I,ID)=CVMGT(-TH1(I,ID),TH1(I,ID),TH1(I,ID).LT.RNULA)
C         TH2(I,ID)=CVMGT(-TH2(I,ID),TH2(I,ID),TH2(I,ID).LT.RNULA)
C         TH3(I,ID)=CVMGT(-TH3(I,ID),TH3(I,ID),TH3(I,ID).LT.RNULA)
C         SGNT1(I,ID)=CVMGT(-1.0D0,1.0D0,TH1(I,ID).LT.RNULA)      
C         SGNT2(I,ID)=CVMGT(-1.0D0,1.0D0,TH2(I,ID).LT.RNULA)      
C         SGNT3(I,ID)=CVMGT(-1.0D0,1.0D0,TH3(I,ID).LT.RNULA)            
C  530 CONTINUE
C  535 CONTINUE
       DO 987 ID=1,3
          DO 987 I=1,NP
             IF(TH1(I,ID).LT.RNULA) THEN
             TH1(I,ID)=-TH1(I,ID)
             SGNT1(I,ID)=-1.0D0
             TH2(I,ID)=-TH2(I,ID)
             SGNT2(I,ID)=-1.0D0
             TH3(I,ID)=-TH3(I,ID)
             SGNT3(I,ID)=-1.0D0                                        
            ELSE
             SGNT1(I,ID)=1.0D0
             SGNT2(I,ID)=1.0D0
             SGNT3(I,ID)=1.0D0
             END IF
  987 CONTINUE
C
      DO 545 ID=1,3
      DO 540 I=1,NP
         B1(I,ID)=EXP(-2.0D0*TH1(I,ID))                                        
         B2(I,ID)=EXP(-2.0D0*TH2(I,ID))                                        
         B3(I,ID)=EXP(-2.0D0*TH3(I,ID))
  540 CONTINUE                       
  545 CONTINUE      
      DO 555 ID=1,3
      DO 550 I=1,NP
         HYTAN1(I,ID)=SGNT1(I,ID)*(1.0D0-B1(I,ID))/(1.0D0+B1(I,ID))       
         HYTAN2(I,ID)=SGNT2(I,ID)*(1.0D0-B2(I,ID))/(1.0D0+B2(I,ID))       
         HYTAN3(I,ID)=SGNT3(I,ID)*(1.0D0-B3(I,ID))/(1.0D0+B3(I,ID))
  550 CONTINUE
  555 CONTINUE           
      DO 565 ID=1,3
      DO 560 I=1,NP
         H1(I,ID)=HYTAN1(I,ID)**6               
         H2(I,ID)=HYTAN2(I,ID)**6               
         H3(I,ID)=HYTAN3(I,ID)**6
  560 CONTINUE               
  565 CONTINUE
      DO 575 ID=1,3
      DO 570 I=1,NP
         TDP1(I,ID)=1.0D0-0.5D0*(G2(I,ID)*H3(I,ID)+G3(I,ID)*H2(I,ID))  
         TDP2(I,ID)=1.0D0-0.5D0*(G3(I,ID)*H1(I,ID)+G1(I,ID)*H3(I,ID))  
         TDP3(I,ID)=1.0D0-0.5D0*(G1(I,ID)*H2(I,ID)+G2(I,ID)*H1(I,ID))
  570 CONTINUE                                                     
  575 CONTINUE
      DO 580 I=1,NP
      C113(I)=TDP1(I,1)*DAMP(I,1,1) 
      C213(I)=TDP2(I,1)*DAMP(I,1,2)   
      C313(I)= TDP3(I,1)*DAMP(I,1,3)
      C123(I)=TDP1(I,2)*DAMP(I,2,1) 
      C223(I)= TDP2(I,2)*DAMP(I,2,2) 
      C323(I)= TDP3(I,2)*DAMP(I,2,3)
      C133(I)=TDP1(I,3)*DAMP(I,3,1) 
      C233(I)= TDP2(I,3)*DAMP(I,3,2) 
      C333(I)= TDP3(I,3)*DAMP(I,3,3)
  580 CONTINUE                  
      DO 585 I=1,NP
         C13(I)=C113(I)+C213(I)+C313(I)
         C23(I)=C123(I)+C223(I)+C323(I)
         C33(I)=C133(I)+C233(I)+C333(I) 
  585 CONTINUE
      DO 590 I=1,NP
         CE3(I) = C13(I) + C23(I) + C33(I)                
  590 CONTINUE               
C
C      SUBROUTINE H3LEPS (R1,R2,R3, VLEPS, VLEPS2) 
C***********************************************************************
C     Calculates 3-body extended-Hartree-Fock energy defined by a
C     LEPS-type function.
C***********************************************************************
C   VLEPS  IS LEPS LOWER SURFACE.
C   VLEPS2 IS LEPS UPPER SURFACE.
C
C     Calculate switching term for the Hartree-Fock component of the
C     diatomic triplet state function.  The notation of Thompson et al.
C     (J. Chem. Phys., 82,5597,1985; and references therein) is used
C     throughout.
C                                
C
C     Loop over 3 coordinates
C
C     Calculate extended-Hartree-Fock curve for ground singlet state of
C     H2 [A.J.C. Varandas & J.D. Silva, J. Chem. Soc. Faraday II
C     (submitted)] 
C
C     Compute the HFACE (Hartree-Fock-approximate correlation energy)
C     potential for the lowest triplet state of H2 [A.J.C. Varandas &
C     J. Brandao Mol. Phys.,45,1982,857] without the 2-body correlation
C     energy.
C
C
C     Add in triplet and subtract singlet 2-body correlation terms
C
C     Calculate effective diatomic triplet state curve.
C
C     Triplet energy
C
C     Contribution to coulomb term
C
C     Contribution to exchange term
      DO 610 I=1,NP
         DR1(I) = R(I,1) - H2RM
         DR2(I) = R(I,2) - H2RM 
         DR3(I) = R(I,3) - H2RM 
         D1(I) = -AZ2*QCOORD(I)
         QSUM(I) = 0.0D0
         D2(I) = 1.0D0 + S(I)*S2(I)*CPHI3(I)
      ET1(I)=R(I,1)*(H2TB1+R(I,1)*(H2TB2+R(I,1)*(H2TB3+R(I,1)*H2TB4))) 
      ET2(I)=R(I,2)*(H2TB1+R(I,2)*(H2TB2+R(I,2)*(H2TB3+R(I,2)*H2TB4))) 
      ET3(I)=R(I,3)*(H2TB1+R(I,3)*(H2TB2+R(I,3)*(H2TB3+R(I,3)*H2TB4))) 
         TW1(I) = EXP(-AL3*R(I,1))
         TW2(I) = EXP(-AL3*R(I,2))
         TW3(I) = EXP(-AL3*R(I,3))
  610 CONTINUE    
      DO 620 I=1,NP
         F1(I) = EXP(D1(I)*QCOORD(I))
         F2(I) = EXP(D1(I)*QCOORD(I)*D2(I))
         TF1(I) = EXP(-HFGAM*DR1(I))
         TF2(I) = EXP(-HFGAM*DR2(I))
         TF3(I) = EXP(-HFGAM*DR3(I))
         AT1(I) = H2TA*RI(I,1)*EXP(-ET1(I))
         AT2(I) = H2TA*RI(I,2)*EXP(-ET2(I))
         AT3(I) = H2TA*RI(I,3)*EXP(-ET3(I))
         WE1(I) = TW1(I)*(AL0+R(I,1)*(AL1+R(I,1)*AL2))
         WE2(I) = TW2(I)*(AL0+R(I,2)*(AL1+R(I,2)*AL2))
         WE3(I) = TW3(I)*(AL0+R(I,3)*(AL1+R(I,3)*AL2))
  620 CONTINUE
      DO 630 I=1,NP
         F11(I) = 1.0D0 - F1(I) 
         F22(I) = 1.0D0 - F2(I) 
         F3(I) = 2.0D0*D1(I)*F1(I)
         F4(I) = D1(I)*F2(I)
      ES1(I)=-HFD*(1.0D0+DR1(I)*(HFA1+DR1(I)*(HFA2+DR1(I)*HFA3)))*TF1(I)
      ES2(I)=-HFD*(1.0D0+DR2(I)*(HFA1+DR2(I)*(HFA2+DR2(I)*HFA3)))*TF2(I)
      ES3(I)=-HFD*(1.0D0+DR3(I)*(HFA1+DR3(I)*(HFA2+DR3(I)*HFA3)))*TF3(I)
         ZT1(I) = AT1(I) + CORRT(I,1) - CORRS(I,1)
         ZT2(I) = AT2(I) + CORRT(I,2) - CORRS(I,2)
         ZT3(I) = AT3(I) + CORRT(I,3) - CORRS(I,3)
  630 CONTINUE                                                         
C      DO 640 I=1,NP
C         OMF(I) = CVMGT(F11(I),F22(I),AS(I).LT.RNULA)
C         FF(I) = CVMGT(F1(I),F2(I),AS(I).LT.RNULA)
C  640 CONTINUE
      DO 876 I=1,NP
         IF(AS(I).LT.RNULA) THEN
         OMF(I)=F11(I)
         FF(I)=F1(I)
         ELSE  
         OMF(I)=F22(I)
         FF(I)=F2(I)
      END IF
  876 CONTINUE
      DO 650 I=1,NP 
         ENT1(I) = FF(I)*WE1(I)+OMF(I)*ZT1(I)
         ENT2(I) = FF(I)*WE2(I)+OMF(I)*ZT2(I)
         ENT3(I) = FF(I)*WE3(I)+OMF(I)*ZT3(I)
  650 CONTINUE
      DO 660 I=1,NP
         QSUM(I)=0.5D0*(ES1(I)+ES2(I)+ES3(I)+ENT1(I)+ENT2(I)+ENT3(I)) 
         XJ1(I)=0.5D0*(ES1(I)-ENT1(I))
         XJ2(I)=0.5D0*(ES2(I)-ENT2(I))
         XJ3(I)=0.5D0*(ES3(I)-ENT3(I))
  660 CONTINUE
      DO 670 I=1,NP
         A1(I)=(XJ1(I)-XJ2(I))**2
         A2(I)=(XJ2(I)-XJ3(I))**2
         A3(I)=(XJ3(I)-XJ1(I))**2
  670 CONTINUE
      DO 680 I=1,NP
         EXCH(I) = SQRT(0.5D0*(A1(I)+A2(I)+A3(I)))
  680 CONTINUE
      DO 690 I=1,NP
         VLEPS(I) = QSUM(I)-EXCH(I)
C         VLEPS2(I) = QSUM(I)+EXCH(I)
  690 CONTINUE               
C
C      SUBROUTINE H3VA (R1, R2, R3, VA )      
C***********************************************************************
C     Calculates Va correction energy                               
C***********************************************************************
C          
      DO 710 I=1,NP
         TT(I) = ALPHA5*PER2(I)
         RVA(I) = (R(I,1)-R(I,2))*(R(I,2)-R(I,3))*(R(I,3)-R(I,1))
  710 CONTINUE
      DO 720 I=1,NP
         EXPVA(I) = EXP(-TT(I)*PER(I))
         RVA1(I)=RVA(I)*RVA(I) 
  720 CONTINUE
      DO 730 I=1,NP
         RVA11(I)=RVA1(I)*XPAR(1)
         RVA12(I)=RVA1(I)*RVA1(I)*XPAR(2)
         RVA13(I)=RVA1(I)*RVA1(I)*RVA1(I)*XPAR(3)
         RVA14(I)=RVA1(I)*RVA1(I)*RVA1(I)*RVA1(I)*XPAR(4) 
  730 CONTINUE
      DO 740 I=1,NP
         VAN(I)=RVA11(I)+RVA12(I)+RVA13(I)+RVA14(I)
  740 CONTINUE
      DO 750 I=1,NP
         VA(I) = VAN(I)*EXPVA(I)
  750 CONTINUE                             
C      DO 760 I=1,NP
C         VA(I)=CVMGT(0.0,VAN(I),EXPVA(I).LT.RNULA)
C  760 CONTINUE
      DO 765 I=1,NP
         IF(EXPVA(I).LT.RNULA) THEN
         VAN(I)=0.0
       END IF
  765 CONTINUE               
C      RETURN                                   
C      END
C
C      SUBROUTINE H3VII (R1, R2, R3, E )     
C***********************************************************************
C     Calculates VII correction energy                             
C***********************************************************************
C
C   Compute B1 function.
C             
      DO 810 I=1,NP
         DR12(I) = R(I,1) - R(I,2)
         DR23(I) = R(I,2) - R(I,3)
         DR31(I) = R(I,3) - R(I,1)
         EXP1(I) = EXP(-BETA1*PER(I))
         EXP3(I) = EXP(-BETA3*PER(I))
         EXP2(I) = EXP(-BETA2*PER2(I))
  810 CONTINUE
      DO 820 I=1,NP
         COS1(I) = 0.5D0*RI(I,2)*RI(I,3)*(R12(I)-R22(I)-R32(I))
         COS2(I) = 0.5D0*RI(I,1)*RI(I,3)*(R22(I)-R12(I)-R32(I))
         COS3(I) = 0.5D0*RI(I,1)*RI(I,2)*(R32(I)-R12(I)-R22(I))
         EB1T(I) = (XPAR(5) + XPAR(6)*PER(I))*EXP1(I)
         EB3T(I) = (XPAR(14) + XPAR(15)*PER2(I))*EXP3(I)
         T6(I) = XPAR(10)*EXP1(I)
         T5(I) = XPAR(13)*EXP2(I)
         T4(I) = XPAR(12)*EXP1(I)
         EQ(I) = DR12(I)*DR12(I) + DR23(I)*DR23(I) + DR31(I)*DR31(I)
         EQP1(I) = 2.0D0*(DR12(I) - DR31(I))
         EQP2(I) = 2.0D0*(-DR12(I) + DR23(I))
         EQP3(I) = 2.0D0*(-DR23(I) + DR31(I))
         RRI(I) = RI(I,1) + RI(I,2) + RI(I,3)
  820 CONTINUE
      DO 830 I=1,NP
         WB(I) = 1.0D0+COS1(I)+COS2(I)+COS3(I)
  830 CONTINUE
      DO 840 I=1,NP
         T7(I) = WB(I)*XPAR(11)*EXP2(I)
         WB2(I) = WB(I)*WB(I)
         EB1(I) = WB(I)*(EB1T(I) + EB3T(I))
         EB4B(I) = WB(I)*(T4(I) + T5(I))
  840 CONTINUE 
      DO 850 I=1,NP
         EB4A(I) = WB(I)*(T6(I) + T7(I))
         EB2(I)=WB2(I)*(XPAR(7)+WB(I)*(XPAR(8)+WB(I)*XPAR(9)))*EXP2(I)
  850 CONTINUE
      DO 855 I=1,NP
         EB4(I) = EB4A(I)*RRI(I) + EB4B(I)*EQ(I)
  855 CONTINUE
      DO 860 I=1,NP
         VII(I) = EB1(I) + EB2(I) + EB4(I)
  860 CONTINUE
C
C      SUBROUTINE H3VIII (R1, R2, R3, VIII )       
C**********************************************************************
C     Calculates VIII correction energy                         
C***********************************************************************
C                                
      DO 910 I=1,NP
         T8(I) = ALPH2*(RHO(I)-RHOL)**2
         T9(I)=1.0D0 + S(I)*(S2(I)*CPHI3(I))
         T10(I) = CD0 + CD1*RHO(I)  
  910 CONTINUE
      DO 920 I=1,NP
        EXPV3(I) = DEXP(-1.0*T8(I)) 
  920 CONTINUE 
      DO 930 I=1,NP
         VIII(I) = S2(I)*T9(I)*T10(I)*EXPV3(I)
  930 CONTINUE
C
C UKOLIKO JE S ILI EXP = 0, POSTAVLJAM FUNKCIJU = 0
C 
C      DO 940 I=1,NP
C         VIII(I) = CVMGT(0.0D0,VIII(I),AS(I).LT.RNULA)
C  940 CONTINUE
C      DO 950 I=1,NP
C         VIII(I) = CVMGT(0.0D0,VIII(I),EXPV3(I).LT.RNULA)
C  950 CONTINUE                  
      DO 432 I=1,NP
         IF(AS(I).LT.RNULA) THEN
         VIII(I)=0.0
         END IF
         IF(EXPV3(I).LT.RNULA) THEN
         VIII(I)=0.0
         END IF
 432  CONTINUE               
      DO 1000 I=1,NP
         CE2(I) = CORRS(I,1)+CORRS(I,2)+CORRS(I,3)
 1000 CONTINUE 
      DO 1010 I=1,NP                                                  
         VC(I) = VA(I) + VII(I) + VIII(I) + CE2(I) + CE3(I) + DASY 
 1010 CONTINUE
      DO 1020 I=1,NP
         POTE1(I) = VLEPS(I) + VC(I) 
C         POTE2(I) = VLEPS2(I) + VC(I) 
 1020 CONTINUE 
      DO 1030 I=1,NP
         VPOTE1(I+II) = POTE1(I)
 1030 CONTINUE 
      DO 4000 J=1,NP
         T = 0.0
         IF(R(J,1).LE.T.OR.R(J,2).LE.T.OR.R(J,3).LE.T) THEN
         VPOTE1(J)=1.0D+06
         POTE1(J)=1.0D+06
         END IF
 4000 CONTINUE 
         II=II+NP
 2000 CONTINUE 
      RETURN
      END
C             BLOCK DATA H3
C***********************************************************************
C     Data for the DMBE H3 surface.
C***********************************************************************
C      IMPLICIT REAL*8 (A-H,O-Z)
C      COMMON /H3DMCM/ ALPH2, ALPHA5, ALPH0, ALPH1, AL0, AL1, AL2, AL3,
C     *   AZ2, BETA1, BETA2, BETA3, BET0, BET1, CD0, CD1, CHH(3), CK0(3),
C     *   CK1(3), DAMPA(3), DAMPB(3), HFD, HFA1, HFA2, HFA3, HFGAM, H2RM,
C     *   H2RMT, H2R0, H2TA, H2TB1, H2TB2, H2TB3, H2TB4, PI, RHOL, SQRT3,
C     *   XPAR(15)
C
C      END 

     
