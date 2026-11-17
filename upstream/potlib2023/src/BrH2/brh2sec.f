C
      SUBROUTINE PREPOT
C
C   System:          BrH2
C   Functional form: extended-LEPS plus three-center term
C   Common name:     BrH2SEC
C   Interface:       3-XS
C   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao
C                    J. Phys. Chem. 99, 207 (1995)
C
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
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Br
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Br-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Br-H)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and 
C        derivatives of the ground electronic state energy with respect 
C        to the coordinates.
C      /PT2CM/ NSURF, NDER, NDUM(21)
C        passes the control flags where
C        NDER  = 0 => no derivatives are computed
C              = 1 => derivatives of the energy for the ground electronic 
C                     state with respect to the coordinates are computed
C        NDUM  - not used
C      /PT4CM/ IPTPRT, IDUM(19)
C        IPTPRT passes the FORTRAN unit number used for potential output
c        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the 
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to 
C        Re(AB), the equilibrium bond length for the AB diatomic.  
C        In this potential the AB valley represents H infinitely far from
C        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the terms
C        EASYBC and EASYAC represent the energies in the H2 and the BrH 
C        valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1 
C      IPTPRT           6
C
C*****
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ R(3), PE, DPEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
C
         PARAMETER(CANGAU = 0.52917706D0)
         PARAMETER(CHARKC = 627.5095D0)
         PARAMETER(R2 = 1.41421356D0)
C
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /V3CCOM/ V3CP1, V3CP2, V3CP3, V3CP4, V3CP5, V3CP6,
     *                   V3C, DV3C(3)
         COMMON /LENGCM/ X(3), COUL(3), EXCH(3), RAD
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /LDERCM/ DXDR(3,3), DEQ(3), DERAD(3)
         COMMON /PARMCM/ DE(3), RE(3)
C
         SAVE
C
         WRITE (IPTPRT,100) DE, RE, BETA, 
     1                 ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,
     2                 V3CP1, V3CP2, V3CP3, V3CP4, V3CP5
C                    
C   Convert the constants to atomic units
C
       DO 10 I = 1, 3 
             DE(I)   = DE(I) / CHARKC
             RE(I)   = RE(I) / CANGAU
             BETA(I) = BETA(I) * CANGAU
 10    CONTINUE
C
C    Initialize the energy in the three asymptotic valleys
       EASYAB = DE(1)
       EASYBC = DE(2)
       EASYAC = DE(3)
C
100   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'BrH2 SEC potential energy surface',
     *      //,1X,T5,'Parameters:',
     *       /,2X,T5,'Bond:',T37,'Br-H',T48,'H-H',T58,'H-Br',
     *       /,2X,T5,'Dissociation energies:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Equilibrium bond lengths:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Morse betas:',
     *            T35,F10.5,T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHH: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for ZHBr: ',T45,F10.5,T55,F10.5,
     *       /,2X,T5,'Parameters for V3c: ',5(T45,F10.5,T55,F10.5,
     *       //,1X,'*****'))
C
      RETURN
C
      ENTRY POT
C
C   Check the value of NDER
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
C   Compute bond dependent Sato parameters
C
         CALL SATO
C
C   Compute the LEPS constants
C
         DO 20 I = 1, 3
               ZPO(I)   = 1.0D0+Z(I)
               OP3Z(I)  = 1.0D0+3.0D0*Z(I)
               TOP3Z(I) = 2.0D0*OP3Z(I)
               ZP3(I)   = Z(I)+3.0D0
               TZP3(I)  = 2.0D0*ZP3(I)
               DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
 20      CONTINUE
C
         PE = 0.D0
C
         DO 30 I = 1, 3
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
               COUL(I) = DO4Z(I) * (ZP3(I)*X(I) - TOP3Z(I)) * X(I)
               EXCH(I) = DO4Z(I) * (OP3Z(I)*X(I) - TZP3(I)) * X(I)
               PE      = PE + COUL(I)
 30      CONTINUE
C
         RAD = SQRT((EXCH(1) - EXCH(2)) ** 2 +
     1              (EXCH(2) - EXCH(3)) ** 2 +
     2              (EXCH(3) - EXCH(1)) ** 2 )
C
         PE  = PE - (RAD / R2) + DE(2)
C
C   Compute the three-center term
         CALL V3CFCN
C
C   Sum all the terms that make up the energy
         PE = PE + V3C
C
C   Compute the derivatives if NDER = 1
C
         IF (NDER .NE. 1) GO TO 700
C
C   DCOUL computes the derivatives of the coulombic terms w.r.t. R(I)
C   DEXCH computes the derivatives of the exchange terms w.r.t. R(I) and
C         also computes the derivative of the term RAD w.r.t. R(i)
C
C   DCOUL must be called before DEXCH becaues this subprogram computes
C   the DXDR terms needed by both subprograms.
C
C   Compute DEQ = dE/dQ * [dQ/dX * dX/dR + dQ/dZ * dZ/dR]
         CALL DCOUL
C   Compute DERAD = dRAD/dJ * [dJ/dX * dX/dR + dJ/dZ * dZ/dR]
         CALL DEXCH
C
C   Put all the partial derivatives together:
C      DEDR  = DEQ + DERAD
C      DEQ   = dE/dQ * DQ/DR
C      DERAD = dE/dRAD * DRAD/DR
C
         DO 40 I = 1,3
               DPEDR(I) = DEQ(I) - DERAD(I)/R2
40       CONTINUE
C
C   Add in the derivatives of V3C w.r.t. R(I)
C
         DO 50 I = 1, 3
               DPEDR(I) = DPEDR(I) + DV3C(I)
50       CONTINUE
C
700   CONTINUE
C
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C
         RETURN
         END
C****
C
         SUBROUTINE SATO
C
C   This subroutine determines the Sato parameters as a function of the
C   bond lengths.
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
         COMMON /PT1CM/ R(3), PE, DPEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         SAVE
C
      DO 10 I = 1, 3
      DO 10 J = 1, 3
10          DZ(I,J) = 0.0D0
C
      R1P2 = R(1) * R(1)
      R2P2 = R(2) * R(2)
      R3P2 = R(3) * R(3)
      R2I  = 1.0D0/R(2)
      ZCHI = 0.0D0
C
C   CAPR IS THE Br TO H2 DISTANCE
      CAPR2 = 0.5D0 * (R1P2 + R3P2 - 0.5D0 * R2P2)
      IF (CAPR2 .LE. 0.D0) THEN
          CAPR2 = 0.0D0
          CAPR  = 0.0D0
      ELSE
          CAPR = SQRT(CAPR2)
C   Compute Chi, the angle between CAPR and RHH
          COSCHI = 0.5D0*(R3P2 - R1P2)/(R(2)*CAPR)
          SCHIP2 = 1.0D0 - COSCHI*COSCHI
          ZCHI = CAPR2*(ZHBRP1*SCHIP2 + ZHBRP2*SCHIP2*SCHIP2)
      ENDIF
C
      Z(1) = 0.1675D0 + ZCHI
      Z(3) = Z(1)
      Z(2) = ZHHP1*R2P2/(ZHHP2+R2P2)
C
C   Compute the derivative of the Sato parameters with respect to 
C   bond lengths.
C   The HH Sato parameter depends on R(2) only
C
      IF (NDER .NE. 1) GO TO 700
      TMP   = R(2)/(ZHHP2 + R2P2) 
      DZ(2,2)  = 2.0D0*Z(2)*(R2I - TMP)
C
C   Compute the derivative of ZHBr with respect R(i) if CHI is
C   is between 0 and 90 degrees.
      IF (CAPR .EQ. 0.D0) GO TO 100
      CAPRI = 1.0D0/CAPR
C   Compute the derivative of R(Br-HH) with respect to R(i)
      DCAPR1 =  0.5D0*R(1)*CAPRI
      DCAPR2 = -0.25D0*R(2)*CAPRI
      DCAPR3 =  0.5D0*R(3)*CAPRI
C   Compute the derivative of CHI with respect to R(i)
C   These derivatives a actually sin(CHI)*dCHI/dR(i)
      TMP1   = COSCHI*CAPRI
      TMP2   = 1.0D0/(R(2)*CAPR)
      DCHI1  = R(1)*TMP2 + TMP1*DCAPR1
      DCHI2  = TMP1*(CAPR + R(2)*DCAPR2)*R2I
      DCHI3  = TMP1*DCAPR3 - R(3)*TMP2
C   Compute the derivative of ZHBr with respect to CAPR and CHI
       DZCHIR = 2.0D0*ZCHI*CAPRI
       DZCHIC = 2.0D0*CAPR2*COSCHI*(ZHBRP1 + 2.0D0*ZHBRP2*SCHIP2)
C   Compute the derivative of ZHBr with respect to R(i)
       DZ(1,1) = DZCHIR*DCAPR1 + DZCHIC*DCHI1
       DZ(1,2) = DZCHIR*DCAPR2 + DZCHIC*DCHI2
       DZ(1,3) = DZCHIR*DCAPR3 + DZCHIC*DCHI3
       DZ(3,1) = DZ(1,1)
       DZ(3,2) = DZ(1,2)
       DZ(3,3) = DZ(1,3)
C
100    CONTINUE
700    CONTINUE
C
       RETURN
       END
C*****
C
         SUBROUTINE V3CFCN
C
C   This subroutine computes the energy and the derivatives of the 
C   three-center term for the H-Br-H barrier as a function of the geometry.
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
         COMMON /PT1CM/ R(3), PE, DPEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /V3CCOM/ V3CP1, V3CP2, V3CP3, V3CP4, V3CP5, V3CP6,
     *                   V3C, DV3C(3)
         SAVE
C
         A1  =  V3CP1
         A2  =  V3CP2
         A3  =  V3CP3
         A4  =  V3CP4
         A5  =  V3CP5
         A6  =  V3CP6
C
         R1P3   = R(1) + R(3)
         R1M3   = R(1) - R(3)
         R1P3M2 = R1P3 - R(2)
         R1M3P2 = R1M3*R1M3
         R132P2 = R1P3M2*R1P3M2
C
C   Compute the bond angle from the three internal coordinates
C
         DIST  = R(1)*R(1) + R(3)*R(3) - R(2)*R(2)
         COSTH = DIST/(2.0D0*R(1)*R(3))
         SINTP2 = 1.0D0 - COSTH*COSTH
         SINTP6 = SINTP2*SINTP2*SINTP2
         SINTP8 = SINTP6*SINTP2
C
         FCN1  = EXP(-A2*R1P3)
         FCN2  = EXP(-A3*R1M3P2)
         FCN3  = EXP(-A4*R132P2)
         FCN4  = 1.D0 + A5*SINTP2 + A6*SINTP8
C
         V3C = A1*FCN1*FCN2*FCN3*FCN4
C
C   Compute the derivatives of V3C with respect to the internal coordinates
C   if NDER = 1
C
         IF (NDER .NE. 1) GO TO 700
C
C   Compute the partial derivatives of the FCN terms
C
      DFCN1 = -A2
      DFCN2 = -2.0D0*A3*R1M3
      DFCN3 = -2.0D0*A4*R1P3M2
      DFCN4 = 2.0D0*A5*COSTH + 
     *        8.0D0*A6*COSTH*SINTP6
C
C   Compute the partial derivative of THETA w.r.t. R(i)
C   DTHDR(i) = SIN(THETA)*DTHDR(I)
C
         DTHDR1  = COSTH/R(1) - 1.0D0/R(3)
         DTHDR2  = R(2)/(R(1)*R(3))
         DTHDR3  = COSTH/R(3) - 1.0D0/R(1)
         V3CDTH  = A1*FCN1*FCN2*FCN3*DFCN4
C
         DV3C(1) = V3C * (DFCN1 + DFCN2 + DFCN3) + V3CDTH*DTHDR1
         DV3C(2) = -V3C * DFCN3 + V3CDTH*DTHDR2
         DV3C(3) = V3C * (DFCN1 - DFCN2 + DFCN3) + V3CDTH*DTHDR3
C
700      CONTINUE
C
         RETURN
         END
C*****
C
      SUBROUTINE DCOUL
C
C    This subprogram computes the contribution to the derivative of
C    the energy with respect to R(i) from the coulombic terms.
C    That is:
C            DE/DR(j) = DE/DQ(i) * DQ(i)/DR(j)
C
C    In this subprogram the following convention is assumed for the 
C    numbering of the internal coordinates:
C    For Br + AB -> ABr + B 
C                   R(1) = R(Br-A)
C                   R(2) = R(A-B)
C                   R(3) = R(Br-B)
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /LENGCM/ X(3), COUL(3), EXCH(3), RAD
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /LDERCM/ DXDR(3,3), DEQ(3), DERAD(3)
C
         DIMENSION DQDX(3), DQDZ(3)
C
         SAVE
C
         DO 10 I = 1, 3
               DEQ(I)  = 0.D0
               DQDX(I) = 0.D0
               DQDZ(I) = 0.D0
         DO 10 J = 1, 3
               DXDR(I,J) = 0.D0
10       CONTINUE
C
C   First compute the partial derivatives dQ/dX, dQ/dZ, and 
C   the derivative of X w.r.t. R(i), DX/DR. 
C
         DO 20 I = 1, 3
            DQDX(I)   = DO4Z(I)*(TZP3(I)*X(I) - TOP3Z(I))
            DQDZ(I)   = DO4Z(I)*X(I)*(X(I) - 6.0D0) - COUL(I)/ZPO(I)
            DXDR(I,I) = -BETA(I)*X(I)
20       CONTINUE
C
C   Compute the contribution to the derivate of E w.r.t. R(i) from
C   the coulombic terms, i.e., dE/dR(i) = dQ/dX * DX/DR(i) + dQ/dZ *DZ/DR(i)
C
         DO 30 I = 1, 3
         DO 30 J = 1, 3
               DEQ(I) = DEQ(I) + DQDX(J)*DXDR(J,I) + 
     *                           DQDZ(J)*DZ(J,I)
30       CONTINUE
C
         RETURN
         END
C****
C
      SUBROUTINE DEXCH
C
C    This subprogram computes the contribution to the derivative of
C    the energy with respect to R(i) from the exchange terms.
C    That is:
C            DE/DR(j) = DE/DJ(i) * DJ(i)/DR(j)
C    where
C            E(LEPS)  = Q1 + Q2 + Q3 - RAD/SQRT(2)
C
C    In this subprogram the following convention is assumed for the 
C    numbering of the internal coordinates:
C    For Br + AB -> ABr + B 
C                   R(1) = R(Br-A)
C                   R(2) = R(A-B)
C                   R(3) = R(Br-B)
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /LENGCM/ X(3), COUL(3), EXCH(3), RAD
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /LDERCM/ DXDR(3,3), DEQ(3), DERAD(3)
C
         DIMENSION DJDX(3), DJDZ(3), DJDR(3,3), DRADDJ(3)
C
         SAVE
C
         DO 10 I = 1, 3
               DERAD(I) = 0.D0
               DJDX(I)  = 0.D0
               DJDZ(I)  = 0.D0
10       CONTINUE
C
         SUMJ = 0.0D0
C
C   First compute the partial derivatives dJ/dX, dJ/dZ
C
         DO 20 I = 1, 3
            SUMJ      = SUMJ + EXCH(I)
            DJDX(I)   = DO4Z(I)*(TOP3Z(I)*X(I) - TZP3(I))
            DJDZ(I)   = DO4Z(I)*X(I)*(3.0D0*X(I) - 2.0D0) - 
     *                  EXCH(I)/ZPO(I)
20       CONTINUE
C
C   Compute the partial derivatives of the exchange term, J, w.r.t. R(i).
C   The outer loop, I, is over the bond lengths, the inner loop is over
C   each X and Z term.
C
         DO 30 I = 1, 3
         DO 30 J = 1, 3
               DJDR(J,I) = DJDX(J)*DXDR(J,I) + DJDZ(J)*DZ(J,I)
30       CONTINUE
C
C   Compute the partial derivatives of the term RAD w.r.t. J(I).
C
         DO 40 I = 1, 3
               DRADDJ(I) = (3.0D0*EXCH(I) - SUMJ)/RAD 
40       CONTINUE
C
         DO 60 I = 1, 3
         DO 70 J = 1, 3
               DERAD(I) = DERAD(I) + DRADDJ(J)*DJDR(J,I)
70       CONTINUE
60       CONTINUE
C
         RETURN
         END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /V3CCOM/ V3CP1, V3CP2, V3CP3, V3CP4, V3CP5, V3CP6,
     *                   V3C, DV3C(3)
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /PARMCM/ DE(3), RE(3)
C 
         SAVE
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER, NSURF /6, 1, 0/
         DATA NDUM /21*0/ 
C
C   Initialize the potential parameters
         DATA DE/90.4473D0, 109.559D0, 90.4473D0/
         DATA RE/1.41443D0, 0.74144D0, 1.41443D0/
         DATA BETA/1.81088D0, 1.94198D0, 1.81088D0/
         DATA ZHHP1, ZHHP2 /2.91230499D-1, 4.99160974D0/
         DATA ZHBRP1, ZHBRP2 /3.84262775D-3, 2.87642885D-3/
         DATA V3CP1, V3CP2 /1.67457884D0, 7.09518875D-1/
         DATA V3CP3, V3CP4 /1.87582679D-1, 2.43517074D-1/
         DATA V3CP5, V3CP6 /2.04145206D-1, 2.57833022D0/
C
         END
C*****
