C
      SUBROUTINE PREPOT
C
C   System:          H3
C   Common name:     PK2
C   Reference:       R. N. Porter and M. Karplus
C                    J. Chem. Phys. 40, 1105 (1964).  
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PPARMS.
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
         COMMON /PT1CM/ R(3),E,DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /PARMCM/ D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM 
C
C   Conversion constants
C   CEV (eV) converts from Hartree to kcal/mol
          PARAMETER (CEV = 27.21161D0)
          PARAMETER (HARTRE = 627.5095D0)
C
      DIMENSION EX(3),S(3),SS(3),E1(3),E3(3),SUME(3),DIFE(3),Q(3),FJ(3),
     1 DE1(3),DE3(3),DQ(3),DJ(3),DS(3)
      DIMENSION RR(3),DELJ(3)
C   CONVERT TO ATOMIC UNITS
      D1 = D1/CEV
      D3 = D3/CEV
      EPS = EPS/CEV
      DEL = DEL/CEV
C   Echo the potential paramters
         WRITE (IPRT, 600) D1, D3, EPS, DEL, RE, ALPH, BETA, FKAP, GAM
C
C    Set the values of the classical energy in the three asymptotic valleys
             EASYAB = D1
             EASYBC = EASYAB
             EASYAC = EASYAB
C
      RETURN
C
      ENTRY POT
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
      DO 12 I = 1,3
      GAMR = GAM*R(I)
      EXPG = EXP(-GAMR)
      ZETA = 1.0D0 + FKAP*EXPG
      RHO = R(I)*ZETA
      RR(I) = 1.0D0/R(I)
      EX(I) = DEL*EXP(-R(I)-R(I))
      EXPR = EXP(-RHO)
      S(I) = (1.0D0+RHO+(RHO*RHO)/3.0D0)*EXPR
      SS(I) = S(I)*S(I)
      REL = RE-R(I)
      EXA = EXP(ALPH*REL)
      EXB = EXP(BETA*REL)
      E1(I) = D1*(1.0D0-EXA)**2-D1
      E3(I) = D3*(1.0D0+EXB)**2-D3
      SUME(I) = E1(I)+E3(I)
      DIFE(I) = E1(I)-E3(I)
      Q(I) = 0.5D0*(SUME(I)+SS(I)*DIFE(I))
      FJ(I) = 0.5D0*(DIFE(I)+SS(I)*SUME(I))
      DE1(I) = 2.0D0*D1*ALPH*EXA*(1.0D0-EXA)
      DE3(I) = 2.0D0*D3*BETA*EXB*(-1.0D0-EXB)
      DS(I) = -RHO*(1.0D0+RHO)*EXPR*(ZETA+GAMR*(1.0D0-ZETA))/3.0D0
      DJ(I) = S(I)*SUME(I)*DS(I)+0.5D0*((1.0D0+SS(I))*DE1(I)-(1.0D0-SS(I        
     V))*
     1   DE3(I))
      DELJ(I) = (1.0D0+RR(I))*EX(I)
   12 DQ(I) = S(I)*DIFE(I)*DS(I)+0.5D0*((1.0D0+SS(I))*DE1(I)+(1.0D0-SS(I        
     V))*
     1   DE3(I))
      DJ(1) = DJ(1) + 2.0D0*S(1)*DS(1)*(DELJ(2)+DELJ(3))
      DJ(2) = DJ(2) + 2.0D0*S(2)*DS(2)*(DELJ(1)+DELJ(3))
      DJ(3) = DJ(3) + 2.0D0*S(3)*DS(3)*(DELJ(1)+DELJ(2))
      FJ(1) = FJ(1) + SS(1)*(DELJ(2)+DELJ(3))
      FJ(2) = FJ(2) + SS(2)*(DELJ(1)+DELJ(3))
      FJ(3) = FJ(3) + SS(3)*(DELJ(1)+DELJ(2))
      S12 = SS(1)-SS(2)
      S23 = SS(2)-SS(3)
      S13 = S12+S23
      FJ12 = FJ(1)-FJ(2)
      FJ23 = FJ(2)-FJ(3)
      FJ13 = FJ12+FJ23
      DJ12 = -SS(2)*EX(1)*(2.0D0+2.0D0*RR(1)+RR(1)*RR(1))
      DJ13 = -SS(3)*EX(1)*(2.0D0+2.0D0*RR(1)+RR(1)*RR(1))
      DJ21 = -SS(1)*EX(2)*(2.0D0+2.0D0*RR(2)+RR(2)*RR(2))
      DJ23 = -SS(3)*EX(2)*(2.0D0+2.0D0*RR(2)+RR(2)*RR(2))
      DJ31 = -SS(1)*EX(3)*(2.0D0+2.0D0*RR(3)+RR(3)*RR(3))
      DJ32 = -SS(2)*EX(3)*(2.0D0+2.0D0*RR(3)+RR(3)*RR(3))
      DJ112 = DJ(1)-DJ12
      DJ113 = DJ(1)-DJ13
      DJ212 = DJ21-DJ(2)
      DJ223 = DJ(2)-DJ23
      DJ313 = DJ31-DJ(3)
      DJ323 = DJ32-DJ(3)
      DJ123 = DJ12-DJ13
      DJ213 = DJ21-DJ23
      DJ312 = DJ31-DJ32
      S123 = S(1)*S(2)*S(3)
      QT = Q(1)+Q(2)+Q(3)-EPS*S123
      CS123 = 1.0D0-S123
      A = CS123*CS123-0.5D0*(S12*S12+S23*S23+S13*S13)
      B = -QT*CS123+0.5D0*(FJ12*S12+FJ23*S23+FJ13*S13)
      C = QT*QT-0.5D0*(FJ12*FJ12+FJ23*FJ23+FJ13*FJ13)
      A1 = -(DS(1)+DS(1))*(CS123*S(2)*S(3)+S(1)*(S12+S13))
      A2 = -(DS(2)+DS(2))*(CS123*S(1)*S(3)+S(2)*(S23-S12))
      A3 = -(DS(3)+DS(3))*(CS123*S(1)*S(2)-S(3)*(S13+S23))
      F = CS123*EPS+QT
      B1 = -CS123*DQ(1)+DS(1)*(S(2)*S(3)*F+S(1)*(FJ12+FJ13))
     1       +0.5D0*(S12*DJ112+S23*DJ123+S13*DJ113)
      B2 = -CS123*DQ(2)+DS(2)*(S(1)*S(3)*F+S(2)*(FJ23-FJ12))
     1       +0.5D0*(S12*DJ212+S23*DJ223+S13*DJ213)
      B3 = -CS123*DQ(3)+DS(3)*(S(1)*S(2)*F-S(3)*(FJ13+FJ23))
     1       +0.5D0*(S12*DJ312+S23*DJ323+S13*DJ313)
      F = QT+QT
      C1 = F*(DQ(1)-EPS*S(2)*S(3)*DS(1))-(FJ12*DJ112+FJ23*DJ123
     1     +FJ13*DJ113)
      C2 = F*(DQ(2)-EPS*S(1)*S(3)*DS(2))-(FJ12*DJ212+FJ23*DJ223
     1     +FJ13*DJ213)
      C3 = F*(DQ(3)-EPS*S(2)*S(1)*DS(3))-(FJ12*DJ312+FJ23*DJ323
     1     +FJ13*DJ313)
      BAC = B*B-A*C
      RBAC = SQRT(BAC)
      E = -B - RBAC                                                     02JUL83
      IF(E.NE.0.D0) E = E/A                                             02JUL83
      T = 2.0D0*RBAC                                                    02JUL83
      DEDR(1) = E*E*A1+(E+E)*B1+C1                                      02JUL83
      IF(DEDR(1).NE.0.D0) DEDR(1) = DEDR(1)/T                           02JUL83
      DEDR(2) = E*E*A2+(E+E)*B2+C2                                      02JUL83
      IF(DEDR(2).NE.0.D0) DEDR(2) = DEDR(2)/T                           02JUL83
      DEDR(3) = E*E*A3+(E+E)*B3+C3                                      02JUL83
      IF(DEDR(3).NE.0.D0) DEDR(3) = DEDR(3)/T
      E = E + D1
C
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',
     *                'potential energy surface PK2',
     *       //,2X,T5,'Potential energy surface parameters:',
     *        /,2X,T15,'D1',  T25,F15.8,T50,'D3',  T60,F15.8,
     *        /,2X,T15,'EPS', T25,F15.8,T50,'DEL', T60,F15.8,
     *        /,2X,T15,'RE',  T25,F15.8,T50,'ALPH',T60,F15.8,
     *        /,2X,T15,'BETA',T25,F15.8,T50,'FKAP',T60,F15.8,
     *        /,2X,T15,'GAM', T25,F15.8)
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not supported in this ',
     *              'version of the potential.')
C
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PARMCM/ D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM 
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/
C
C   Initialze the parameters for the potential
      DATA D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM  / 4.7466D0,1.9668D0,
     1       1.40083D0,1.04435D0,1.000122D0,-17.5D0,28.2D0,.60D0,0.65D0/
         END
C*****
