      SUBROUTINE PREPOT
C
C   System:          FH2
C   Functional Form: eLEPS plus three-center term
C   Common Name:     No. 5A
C   Interface:       3-2V
C   Number of derivatives: 1
C   References:      R. Steckler, D.G. Truhlar, and B.C. Garrett, 
C                    J. Chem. Phys. vol. 82, pp. 5499-5505 (1985).
C
C    Note: This is a vectorized version of the surface.
C          This version makes AB and AC Sato parameters a function of PHI.
C          PHI is angle of BC-to-A vector with BC axis.
C          BC is a homonuclear diatom.  The BC Sato parameter is made a
C          function of the R(BC) distance and the bend angle PHI.
C          A three center energy term is also added to correct the HFH
C          barrier height.
C
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
C   Surfaces: 
C      ground electronic state
C
C   Zer0 of energy is F infinitely removed from H2 at R=RE.
C
c	Modified by GCL on 90-06-20 to run with the GNVP code
c	Note: This routine still has the analytic derivatives.
c
      IMPLICIT double precision (A-H,O-Z)
      LOGICAL LCOL, LDERIV, LZERO
      CHARACTER*80 DATAF
      DIMENSION DM(3), RE(3), BETA(3)
      DIMENSION Z(3), Q(3), XJ(3), DQ(3,3), DJ(3,3), DZ(3,3)
      DIMENSION R(3)
      dimension rr(nt,3), ee(nt)
      COMMON /PT1CM/ RSV(3),V,DV(3)
      DATA TDEG /57.29577951D0/
      DATA DM,RE,BETA/141.196,109.449,141.196,.9170,.7419,.9170,2.2187,
     *   1.9420,2.2187/
      DATA A1,A2,A,B,C,D,F,G/0.0395D0,0.201D0,1.26D0,1.60D0,0.0D0,0.0D0,
     *  0.0D0,0.2D0/
      DATA BET,BETP,ZFH0,ZHH0,ZHHST/-0.101168D0,-0.46183D0,0.170D0,
     *    0.120D0,0.23107D0/
      DATA DEG1,DEG2,DEG3,DEG4,DEG3P/10.0D0,20.0D0,60.0D0,90.0D0,41.0D0/
      DATA C1,C2,C3,C4,C5,C6/1.15D-4,-6.0D-6,0.016D0,0.0005D0,0.00069D0,
     *    3.308D-2/
      DATA C3C,ALF,ALFP,ALFT,P3C,Q3C/0.270063D0,0.18D0,0.078D0,2.14D0,
     *   0.95D0,0.48D0/
C   EXLARG IS LN(LARGEST FLOATING POINT NUMBER ALLOWED)
      DATA EXLARG/80.D0/
      save lderiv,dm,re,beta,ZHHDIF,OMP3C,OMQ3C,HA2,G2
      save /pt1cm/
C
C   SET POTENTIAL ENERGY SURFACE PARAMETERS
C        LDERIV = TRUE IF DERIVATIVES ARE TO BE COMPUTED
      LDERIV = .FALSE.
C
      WRITE (6,600) DM, RE, BETA,
     *   A, ZFH0, A1, DEG1, A2, DEG2, B, DEG3, C, DEG4, D, DEG3P,
     *   F, C1, G, C2, BET, C3, BETP, C4, C5, C6
      WRITE (6,601) C3C, ALF, ALFP, ALFT, P3C, Q3C
      DO 10 I = 1,3
         DM(I) = DM(I) / 627.5095D0
         RE(I) = RE(I) / 0.52917706D0
         BETA(I) = BETA(I) * 0.52917706D0
   10 CONTINUE
C   USEFUL CONSTANTS
      ZHHDIF = ZHHST - ZHH0
      OMP3C = 1.D0 - P3C
      OMQ3C = 1.D0 - Q3C
      HA2 = 0.5D0 * A2
      G2 = G * G           
C
      RETURN
C                            
      ENTRY POT(rr, ee, nt, nsurf)
c                            
      do 21 ms = 1, nt
         do 22 msz = 1, 3
            r(msz) = rr(ms,msz)
   22  continue
C   TEST R'S FOR PHYSICAL VALUES
c      DO 20 I = 1,3
c         R(I) = RSV(I)
c   20 CONTINUE
c      DO 30 I1 = 1,3
c         I2 = MOD(I1,3) + 1
c         I3 = MOD(I2,3) + 1
c         T = R(I2) + R(I3)
c         IF (R(I1) .GT. T) R(I1) = T
c   30 CONTINUE
      R12 = R(1) * R(1)
      R22 = R(2) * R(2)
      R32 = R(3) * R(3)
C   CAPR IS THE F TO H2 DISTANCE
      CAPR2 = 0.5D0 * (R12 + R32 - 0.5D0 * R22)
      IF (CAPR2 .LT. 0.D0) CAPR2 = 0.0D0
      CAPR = SQRT(CAPR2)
C   PHI IS THE ANGLE BETWEEN CAPR AND RHH
C       PHI IS RESTRICTED TO BE BETWEEN O AND 90 DEGREES
      T = 0.5D0 * (R32 - R12)
      IF (ABS(T) .GT. 1.D-7) T = T / (R(2) * CAPR)
      T2 = 1.0D0
      SGNPHI = SIGN(T2,T)
      COSPHI = SGNPHI * T
      COSPHI = MIN(COSPHI,T2)
      PHI = ACOS(COSPHI)
      SINP2 = 1.D0 - COSPHI * COSPHI
      IF(SINP2 .LT. 0.D0) SINP2 = 0.D0
      SINPHI = SQRT(SINP2)
C   LCOL IS SET TRUE IF COLLINEAR GEOMETRY
      LCOL = ABS(SINPHI) .LT. 1.D-6 .OR. ABS(CAPR) .LT. 1.D-6
C   COMPUTE COLLINEAR SATO PARAMETERS
      Z(1) = ZFH0
      T = A * (R(2) - B)
      IF (.NOT.(T .LT. EXLARG)) GO TO 35
         T = EXP(T)
         T1 = 1.D0 / T
         HSECH = 1.D0/(T + T1)
         ZHH = A1 + HA2 * ((T - T1) * HSECH + 1.D0)
         GO TO 36
   35 CONTINUE
         HSECH = 0.D0
         ZHH = A1 + A2
   36 CONTINUE
      T = R(2) - F
      EXX = C * EXP(D * T * T)
      ZHH = ZHH - EXX
      Z(2) = ZHH
      IF (.NOT.LDERIV) GO TO 50
         DO 40 I = 1,3
            DO 40 J = 1,3
               DZ(I,J) = 0.D0
   40    CONTINUE
C  ONLY THE HH SATO PARAMETER HAS A NONZERO DERIVATIVE - W.R.T. R2
         DZ(2,2) = 2.D0 * A * A2 * HSECH * HSECH
     *      - 2.D0 * D * T * EXX
   50 CONTINUE
C   COMPUTE NONCOLLINEAR CONTRIBUTIONS TO SATO PARAMETERS
      IF (LCOL) GO TO 160
C   COMPUTE Z(PHI) FOR FH SATO PARAMETER
C   DZPHI IS THE DERIVATIVE OF Z(PHI) W.R.T. PHI
         DEG = TDEG * PHI
         IF (.NOT.(DEG.LT.0. .OR. DEG.GT.90.D0)) GO TO 60
            WRITE(6,6000) DEG
            STOP
   60    CONTINUE
         IF (.NOT.(DEG .LE. DEG1)) GO TO 70
            ZPHI = 0.D0
            DZPHI = 0.D0
            GO TO 140
   70    CONTINUE
            IF (.NOT.(DEG .LE. DEG2)) GO TO 80
               DEGDIF = DEG - DEG1
               DD2 = DEGDIF * DEGDIF
               ZPHI = DD2 * (C1 + C2 * DEGDIF)
               IF (LDERIV) DZPHI = TDEG *
     *            DEGDIF * (2.D0 * C1 + 3.D0 * C2 * DEGDIF)
               GO TO 130
   80       CONTINUE
               IF (.NOT.(DEG .LE. DEG3)) GO TO 90
                  DEGDIF = DEG - DEG3P
                  ZPHI = C3 + C4 * DEGDIF
                  IF (LDERIV) DZPHI = C4 * TDEG
                  GO TO 120
   90          CONTINUE
                  IF (.NOT.(DEG .LE. DEG4)) GO TO 100
                     ZPHI = C5 + C6 * SINPHI * SINPHI
                     IF (LDERIV) DZPHI = 2.D0 * C6 * SINPHI * COSPHI
                     GO TO 110
  100             CONTINUE
                     WRITE (6,6001) DEG
                     STOP
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
C   COMPUTE SCALE FACTOR TO DAMP OUT CORRECTION AT LARGE DISTANCE
         SCALE1 = CAPR2 / (G2 + CAPR2)
         Z(1) = Z(1) + SCALE1 * ZPHI
C   COMPUTE SCALE FACTOR FOR HH SATO PARAMETER
         SCALE2 = (ZHH - ZHH0) / ZHHDIF
C   COMPUTE PHI DEPENDENT PART OF HH SATO PARAMETER
         CPHI = (BET + BETP * SINP2) * SINP2
C   COMPUTE HH SATO PARAMTER
         Z(2) = Z(2) + SCALE1 * SCALE2 * CPHI
         IF (.NOT.LDERIV) GO TO 150
C   COMPUTE DERIVATIVES
C   FIRST, DERIVATIVES OF PHI W.R.T. R1, R2, R3
            T = 1.D0 / CAPR
            T1 = 0.5D0 * COSPHI * T
            T2 = SGNPHI / R(2)
            DPHI1 = R(1) * (T1 + T2)
            DPHI3 = R(3) * (T1 - T2)
            T1 = 1.D0 / SINPHI
            T2 = T1 * T
            DPHI1 = DPHI1 * T2
            DPHI3 = DPHI3 * T2
            DPHI2 = R(2) * COSPHI * T1 * (1.D0 / R22 - 0.25D0 / CAPR2)
            T1 = (1.D0 - SCALE1) / CAPR2
            T2 = T1 * ZPHI
C   NEXT, DERIVATIVES OF FH SATO PARAMETERS
            DZ(1,1) = SCALE1 * (T2 * R(1) + DZPHI * DPHI1)
            DZ(1,2) = SCALE1 * (-0.5D0 * T2 * R(2) + DZPHI * DPHI2)
            DZ(1,3) = SCALE1 * (T2 * R(3) + DZPHI * DPHI3)
            T1 = T1 * CPHI
C   DERIVATIVE OF THE PHI DEPENDENT PART OF HH SATO W.R.T. PHI
            DCPHI = 2.D0 * (BET + 2.D0 * BETP * SINP2) * SINPHI * COSPHI
            T2 = SCALE1 * SCALE2
C   LAST, DERIVATIVES OF HH SATO PARAMETER
            DZ(2,1) = T2 * (T1 * R(1) + DCPHI * DPHI1)
            DZ(2,2) = DZ(2,2) * (1.D0 + SCALE1 * CPHI / ZHHDIF)
     *         + T2 * (-0.5D0 * T1 * R(2) + DCPHI * DPHI2)
            DZ(2,3) = T2 * (T1 * R(3) + DCPHI * DPHI3)
  150    CONTINUE
  160 CONTINUE
C   COMPUTE POTENTIAL FROM LEPS FORM
      V = 0.D0
      Z(3) = Z(1)
      IF (.NOT.LDERIV) GO TO 170
         DO 165 I = 1,3
            DZ(3,I) = DZ(1,I)
            DV(I) = 0.D0
            Q(I) = 0.D0
            XJ(I) = 0.D0
            DO 162 J = 1,3
               DQ(I,J) = 0.D0
               DJ(I,J) = 0.D0
  162       CONTINUE
  165    CONTINUE
  170 CONTINUE
      LZERO = .TRUE.
      DO 200 I = 1,3
         EX = EXP(BETA(I) * (RE(I) - R(I)))
         LZERO = LZERO .AND. EX.EQ.0.0D0
         IF (EX.EQ.0.0D0) GO TO 195
            ZP1 = Z(I) + 1.D0
            ZP3 = ZP1 + 2.D0
            OP3Z = 1.D0 + 3.D0 * Z(I)
            T = 0.25D0 * DM(I) * EX / ZP1
            T1 = ZP3 * EX
            T2 = OP3Z * EX
C   Q AND XJ ARE THE COULUMB AND EXCHANGE PARTS, INDEX I RUNS OVER THE
C   COORDINATES OF WHICH THEY ARE EXPLICIT FUNCTIONS
            Q(I) = T * (T1 - 2.D0  * OP3Z)
            XJ(I) = T * (T2 - 2.D0  * ZP3)
C   PUT SUM OF Q'S IN V
            V = V + Q(I)
            IF (.NOT.LDERIV) GO TO 190
C   DERIVATIVE OF Q AND ZJ W.R.T. R1, R2, R3
C   ELEMENT I,J IS THE DERIVATIVE OF Q(I) W.R.T. RJ
               T3 = 2.D0 * T * (EX + 2.D0) / ZP1
               DO 180 J = 1,3
                  DQ(I,J) = -T3 * DZ(I,J)
                  DJ(I,J) = -DQ(I,J)
  180          CONTINUE
               T = 2.D0*T*BETA(I)
               DQ(I,I) = DQ(I,I) - T * (T1 - OP3Z)
               DJ(I,I) = DJ(I,I) - T * (T2 - ZP3)
  190       CONTINUE
  195    CONTINUE
  200 CONTINUE
      IF (LZERO) GO TO 270
         XJX = 0.D0
C   COMPUTE THE TOTAL EXCHANGE TERM
         DO 230 I1 = 1,3
            I2 = MOD(I1,3) + 1
            T = XJ(I1) - XJ(I2)
            XJX = XJX + T * T
            IF (.NOT.LDERIV) GO TO 220
C   PUT THE DERIVATIVES OF THE TOTAL EXCHANGE TERM IN DV
               DO 210 J = 1,3
                  DV(J) = DV(J) + T * (DJ(I1,J) - DJ(I2,J))
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
         XJX = SQRT(0.5D0 * XJX)
C   COMPUTE THE LEPS PART OF THE POTENTIAL ENERGY
         V = V - XJX
         IF (.NOT.LDERIV) GO TO 260
C   COMPUTE DERIVATIVES OF THE LEPS PART
            IF(ABS(XJX) .GT. 1.D-14) T = 0.5D0 / XJX
            DO 250 J = 1,3
C   THE DERIVATIVE OF THE EXCHANGE PART MUST BE DIVIDED BY 2*XJX
               DV(J) = -T * DV(J)
               DO 240  I = 1,3
C   THEN ADD IN THE DERIVATIVE OF Q
                  DV(J) = DV(J) + DQ(I,J)
  240          CONTINUE
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
C   DM(2) ADDED TO PUT ZERO AT REQ FOR REACTANTS
      V = V + DM(2)
C   COMPUTE THE 3-CENTER TERM
      RSUM = R(1) + R(3)
      RDIF = R(1) - R(3)
      RDIF2 = RDIF * RDIF
      T1 = -ALFP * RSUM
      EX1 = C3C * EXP(T1 * RSUM)
      T2 = -ALF * RDIF
      EX2 = OMP3C * EXP(T2 * RDIF)
      T3 = -ALFT * RDIF2
      EX3 = P3C * EXP(T3 * RDIF2)
      R1I = 1.D0 / R(1)
      R3I = 1.D0 / R(3)
      T4 = R1I*R3I
C   COSTH IS THE ANGLE BETWEEN R1 AND R3
      COSTH = 0.5D0 * (R12 + R32 - R22) * T4
      T5 = OMQ3C * COSTH
      T = Q3C + T5 * COSTH
C   E3C CONTAINS THE THREE CENTER CORRECTION TO V
      E3C = (EX2 + EX3) * T * EX1
      V = V + E3C
      IF (.NOT.LDERIV) GO TO 280
C   COMPUTE DERIVATIVE OF THE 3C TERM
C   FIRST, DERIVATIVE OF COSTH
         DCOS1 = R3I - COSTH*R1I
         DCOS2 = - R(2)*T4
         DCOS3 = R1I - COSTH*R3I
         T4 = EX2 + EX3
         T5 = T4 * T5
C   NOW, DERIVATIVES OF E3C
         DE3C1 = 2.D0*((T2 * EX2 + 2.D0 * T3 * RDIF * EX3 + T1 * T4) 
     *      * T + T5 * DCOS1) * EX1
         DE3C2 = 2.D0 * T5 * DCOS2 * EX1
         DE3C3 = 2.D0*((-T2 * EX2 - 2.D0 * T3 * RDIF * EX3 + T1 * T4)
     *      * T + T5 * DCOS3) * EX1
         DV(1) = DV(1) + DE3C1
         DV(2) = DV(2) + DE3C2
         DV(3) = DV(3) + DE3C3
  280 CONTINUE
         ee(ms) = v
  21   continue
      RETURN
  500 FORMAT(1X,L1)
  501 FORMAT (3F20.10)
  502 FORMAT (4F20.10)
  600 FORMAT (/' LEPS POTENTIAL ENERGY SURFACE PARAMETERS'
     * /' BOND',20X,2HAB,8X,2HBC,8X,2HAC /' DISS. ENERGIES',5X,
     * 3F10.5,'  KCAL/MOL' /' EQUILIBRIUM',8X,3F10.5,'  ANGSTROMS'
     * /' MORSE BETA',9X,3F10.5,'  ANSTROMS''-1'
     * /' HH SATO PARAMETER FIT', T41, ' FH SATO PARAMETER FIT'
     * /' A=',T10,1PE13.5,  T41,' ZFH0=',T50,E13.5
     * /' A1=',T10,E13.5,   T41,' DEG1=',T50,1PE13.5
     * /' A2=',T10,E13.5,   T41,' DEG2=',T50,E13.5
     * /' B=',T10,E13.5,    T41,' DEG3=',T50,E13.5
     * /' C=',T10,E13.5,    T41,' DEG4=',T50,E13.5
     * /' D=',T10,E13.5,    T41,' DEG3P=',T50,E13.5
     * /' F=',T10,E13.5,    T41,' C1=',T50,E13.5
     * /' G=',T10,E13.5,    T41,' C2=',T50,E13.5
     * /' BET=',T10,E13.5,  T41,' C3=',T50,E13.5
     * /' BETP=',T10,E13.5, T41,' C4=',T50,E13.5
     * / T41,' C5=',T50,E13.5
     * /T41,' C6=',T50,E13.5)
  601 FORMAT(' PARAMETERS FOR THREE-CENTER TERM'/
     * ' C3C=',T10,E13.5/' ALF=',T10,E13.5/' ALFP=',T10,E13.5/
     * ' ALFT=',T10,E13.5/' P3C=',T10,E13.5/' Q3C=',T10,E13.5/)
  900 FORMAT(1X,A80)
 6000 FORMAT(' IN FH2N4 POTENTIAL, PHI=',F10.5,' DEGREES, BUT MUST BE',
     *   ' GREATER THAN ZERO.')
 6001 FORMAT(' IN FH2N4 POTENTIAL, PHI=',F10.5,' DEGREES, BUT MUST BE',
     *   ' LESS THAN 90 DEGREES.')
 6100 FORMAT(' PROBLEM WITH POTENTIAL DATA FILE, CHECK FILE NAME')
 6101 FORMAT(' POTENTIAL DATA FROM FILE',T40,A80)
      END

