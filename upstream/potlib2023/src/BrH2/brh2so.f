C
      subroutine prepot                                                 4s19m94
C
C   System:          BrH2
C   Functional form: extended-LEPS plus three-center term
C   Common name:     BrH2SEC
C   Interface:       3-1V
C   Number of electronic surfaces: 2
C   Number of derivatives: 0
C   Reference:       G. C. Lynch, D. G. Truhlar, F. B. Brown, and J.-g. Zhao
C                    J. Phys. Chem. 99, 207 (1995)
C
C   Note:            This is a vectorized version of the BRH2SEC surface
C                    without derivatives
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy 
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground and first excited electronic states, diabatic representation
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Br
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Br-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Br-H)
C
C
      implicit double precision (a-h,o-z)
      dimension e(*),r(3,*)                                             4s19m94
      call prepot11                                                     4s19m94
      call prepot12                                                     4s19m94
      call prepot22                                                     4s19m94
      return                                                            4s19m94
      entry pot(r,e,nt,nsurf)                                           4s19m94
      if(nsurf.eq.1)then                                                4s19m94
       call pot11(r,e,nt)                                               4s19m94
      elseif (nsurf.eq.2) then                                          4s19m94
       call pot12(r,e,nt)                                               4s19m94
      elseif (nsurf.eq.3) then                                          4s19m94
       call pot22(r,e,nt)                                               4s19m94
      else                                                              4s19m94
       stop                                                             4s19m94
      endif                                                             4s19m94 
      return                                                            4s19m94
      end                                                               4s19m94
      SUBROUTINE PREPOT11
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER (IUNIT = 6)
C
         PARAMETER(CANGAU = 0.52917706D0)
         PARAMETER(CHARKC = 627.5095D0)
C
         COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
         COMMON /V3CCOM/ V3CP1, V3CP2, V3CP3, V3CP4, V3CP5, V3CP6,
     *                   V3C, DV3C(3)
         COMMON /LENGCM/X(3), COUL(3), EXCH(3), RAD
         COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
         COMMON /LDERCM/ DXDR(3,3), DEQ(3), DERAD(3)
         COMMON /DECOM/ DE(3)
C
         DIMENSION D(3), RE(3), BETAA(3), B(3),r(3),rr(n,3),e(*)
C
         SAVE
C
         DATA D/90.4473D0, 109.559D0, 90.4473D0/
         xlamb=.4568d0/27.21161d0
         A1  =  V3CP1
         A2  =  V3CP2
         A3  =  V3CP3
         A4  =  V3CP4
         A5  =  V3CP5
         A6  =  V3CP6
         DATA RE/1.41443D0, 0.74144D0, 1.41443D0/
         DATA BETAA/1.81088D0, 1.94198D0, 1.81088D0/
         DATA R2/1.41421356D0/
C
         ZHHP1  =  2.91230499D-1
         ZHHP2  =  4.99160974D0
         ZHBRP1 =  3.84262775D-3
         ZHBRP2 =  2.87642885D-3
         V3CP1  =  1.67457884D0
         V3CP2  =  7.09518875D-1
         V3CP3  =  1.87582679D-1
         V3CP4  =  2.43517074D-1
         V3CP5  =  2.04145206D-1
         V3CP6  =  2.57833022D0
C
         WRITE (IUNIT,100) D, RE, BETAA, 
     1                 ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,
     2                 V3CP1, V3CP2, V3CP3, V3CP4, V3CP5
C                    
C   Convert the constants to atomic units
C
         DO 10 I = 1, 3 
               DE(I)   = D(I) / CHARKC
               RE(I)   = RE(I) / CANGAU
               BETA(I) = BETAA(I) * CANGAU
 10      CONTINUE
C
C
100   FORMAT(/,2X,T5,'PREPOT has been called for the BrH2 potential',
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
C
         ENTRY POT11(rr,e,n)
C
      do 7777 k=1,n
       r(1)=rr(k,1)
        r(2)=rr(k,2)
        r(3)=rr(k,3)

C   Compute bond dependent Sato parameters
      R1P2 = R(1) * R(1)
      R2P2 = R(2) * R(2)
      R3P2 = R(3) * R(3)
      R2I  = 1.0D0/R(2)
      ZCHI = 0.0D0
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
               DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
               B(I) = BETA(I)*DO4Z(I)*2.0D0
 20      CONTINUE
         PE = 0.D0
         DO 30 I = 1, 3
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
               COUL(I) = DO4Z(I) * (ZP3(I)*X(I) - TOP3Z(I)) * X(I)
               EXCH(I) = DO4Z(I) * (OP3Z(I)*X(I) - TZP3(I)) * X(I)
               PE      = PE + COUL(I)
 30      CONTINUE
         RAD = SQRT((EXCH(1) - EXCH(2)) ** 2 +
     &              (EXCH(2) - EXCH(3)) ** 2 +
     &              (EXCH(3) - EXCH(1)) ** 2 )
         PE  = PE - (RAD / R2) + DE(2)
C   Compute the three-center term
         R1P3   = R(1) + R(3)
         A1  =  V3CP1
         A2  =  V3CP2
         A3  =  V3CP3
         A4  =  V3CP4
         A5  =  V3CP5
         A6  =  V3CP6
         R1M3   = R(1) - R(3)
         R1P3M2 = R1P3 - R(2)
         R1M3P2 = R1M3*R1M3
         R132P2 = R1P3M2*R1P3M2
C   Compute the bond angle from the three internal coordinates
         DIST  = R(1)*R(1) + R(3)*R(3) - R(2)*R(2)
         COSTH = DIST/(2.0D0*R(1)*R(3))
         SINTP2 = 1.0D0 - COSTH*COSTH
         SINTP6 = SINTP2*SINTP2*SINTP2
         SINTP8 = SINTP6*SINTP2
         FCN4  = 1.D0 + A5*SINTP2 + A6*SINTP8
         V3C = A1*FCN4*exp(-(a2*r1p3+a3*r1m3p2+a4*r132p2))
C   Sum all the terms that make up the energy
         PE = PE + V3C
         e(k)=pe+xlamb/3.d0
7777  continue
c
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
         RETURN
         END
C
      SUBROUTINE PREPOT12
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension rr(n,3),e(*)
      save
      xlamb=0.4568d0/27.21161d0
      RETURN
C
      ENTRY POT12(rr,e,n)
      do 7777 k=1,n
       e(k)=-dsqrt(2.d0)*xlamb/3.d0
7777  continue
      RETURN
      END
C
      SUBROUTINE PREPOT22
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IUNIT = 6)
C
      PARAMETER(CANGAU = 0.52917706D0)
      PARAMETER(CHARKC = 627.5095D0)
C
      COMMON /SATOCM/ ZHHP1, ZHHP2, ZHBRP1, ZHBRP2, 
     *                   Z(3), DZ(3,3)
      COMMON /LENGCM/X(3), COUL(3), EXCH(3), RAD
      COMMON /LSATCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),
     *                   DO4Z(3), BETA(3)
      COMMON /DECOM/ DE(3)
C
      DIMENSION D(3), RE(3), BETAA(3), B(3),rr(n,3),e(*),r(3)
C
      SAVE
C
      DATA D/90.4473D0, 109.559D0, 90.4473D0/
      DATA RE/1.41443D0, 0.74144D0, 1.41443D0/
      DATA BETAA/1.81088D0, 1.94198D0, 1.81088D0/
      DATA R2/1.41421356D0/
C
      data ccon /0.5d0/
      xlamb=0.4568d0/27.21161d0
      ZHHP1  =  2.91230499D-1
      ZHHP2  =  4.99160974D0
      ZHBRP1 =  3.84262775D-3
      ZHBRP2 =  2.87642885D-3

      WRITE (IUNIT,100) D, RE, BETAA, 
     1                 ZHHP1, ZHHP2, ZHBRP1, ZHBRP2,
     2                 V3CP1, V3CP2, V3CP3, V3CP4, V3CP5
C                    
C   Convert the constants to atomic units
C
      DO 10 I = 1, 3 
         DE(I)   = D(I) / CHARKC
         RE(I)   = RE(I) / CANGAU
         BETA(I) = BETAA(I) * CANGAU
 10   CONTINUE
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
C
      ENTRY POT22(rr,e,n)
      do 7777 k=1,n
       r(1)=rr(k,1)
       r(2)=rr(k,2)
       r(3)=rr(k,3)
       y2=exp(-beta(2)*(r(2)-re(2)))
       y1=exp(-beta(1)*(r(1)-re(1)))             
       y3=exp(-beta(3)*(r(3)-re(3)))             
       vmorse=de(2)*(y2**2-2.d0*y2)
       vanti1=de(1)*(y1**2+2.d0*y1)
       vanti3=de(3)*(y3**2+2.d0*y3)
       pe=de(2)+vmorse+(vanti1+vanti3)*ccon+2.d0*xlamb/3.d0
       e(k)=pe
7777  continue
c
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
      RETURN
      END

