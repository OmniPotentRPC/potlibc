C
      SUBROUTINE PREPOT
C
C   System:          IH2
C   Functional form: 
C   Common name:     
C   Reference:       L. Raff, L. Stivers, R. N. Porter, D. L. Thompson, L. B. Sims, 
C                    J. Chem. Phys., Vol. 52, 3449 (1970).
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PPARMS
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
         COMMON /PARMCM/ D1(3),D3(3),ALPH(3),BETA(3),RE(3),
     *                   A(3),C(3),SIG(3),RC(3)
         DIMENSION E1(3),E3(3),FJ(3),DE1(3),DE3(3),Q(3),DJ(3),DQ(3)
C   Echo the potential name and the potential parameters to the 

         WRITE (IPRT, 600) D1, D3, ALPH, BETA, RE 
         WRITE (IPRT, 610) A, C, SIG, RC
C
C   Convert the potential parameters to atomic units.
      DO 5 I = 1,3
      D1(I) = D1(I)/27.21161D0
      D3(I) = D3(I)/27.21161D0
    5 C(I) = C(I)/27.21161D0
C
C   Initialize the asymptotic energy values
         EASYAB = D1(1)
         EASYBC = D1(2)
         EASYAC = D1(3)
C
      RETURN
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
      E=QT-EXCH
      E = E + D1(2)
      IF (NDER .EQ. 1) THEN
          DO 3 I=1,3
3              DEDR(I)=DQ(I)-0.5D0*DJ(I)*(3.0D0*FJ(I)-SUMJ)/EXCH
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
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PARMCM/ D1(3),D3(3),ALPH(3),BETA(3),RE(3),
     *                   A(3),C(3),SIG(3),RC(3)
C
C   Initialize the flags for the potential
         DATA NDER /1/
         DATA NSURF /0/
         DATA NDUM /21*0/
         DATA IPRT /6/
C   The energy values are in eV.
         DATA D1 / 3.194D0, 4.7466D0, 3.194D0/
         DATA D3 / 1.44399D0, 1.9668D0, 1.44399D0/
         DATA ALPH / 0.9468D0, 1.04435D0, 0.9468D0/
         DATA BETA / 0.794D0, 1.0001D0, 0.794D0/
C   The lengths are in bohr.
         DATA RE / 3.032D0, 1.402D0, 3.032D0/
         DATA A / -3.60481D0, 1.0D0, -3.60481D0/ 
         DATA C / 73599.62D0, 25.5530D0, 73599.62D0/ 
         DATA SIG / 2.47075D0, 1.67564D0, 2.47075D0/
         DATA RC / 4.25D0, 1.6D0, 4.25D0/
C
         END
C*****
