c
      subroutine prepot
c System:          F2H
c Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
c Common name:     JOT-II
c Reference:       N. Jonathon, S. Okuda, and D. Timlin, Mol. Phys.
c                  24, 1143 (1972).
c
c PREPOT must be called once before any calls to POT are made.  It should
c never be called more than once.
c The potential parameters are included in DATA statements.
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

      implicit double precision (a-h,o-z)

      common /pt1cm/ r(3), energy, dedr(3)
      common /pt2cm/ nsurf, nder, ndum(21)
      common /pt5cm/ easyab, easybc, easyac

      dimension de(3),re(3),b(3),s(3),hde(3),hpf(3),pfs(3),pft(3)
      save de,re,b,s,cevau,canau,hde,hpf,pfs,pft
c
c Constants used by the code.  de is the dissociation energy (eV),
c re is the equilibrium bond length (bohr), b is the Morse beta
c parameter (bohr**-1), and s is the Sato parameter (unitless) for
c the appropriate diatom (1 = I-H, 2 = H-I', 3 = I-I').  The other
c constants may be found in the literature reference.  cevau is the
c conversion factor between eV and hartrees.
c
      data de     / 6.12d0,    1.63d0,  6.12d0 /
      data re     / 0.9171d0,  1.418d0, 0.9171d0 /
      data b      / 2.219d0,   2.920d0, 2.219d0 /
      data s      / 0.d0,     -0.35d0,  0.d0 /
      data cevau  / 27.2113961d0 /
      data canau  / 0.529177249d0 /
c
c Write header potential parameters to FORTRAN unit 6.
c
      write(6,*)
      write(6,*) 'PREPOT has been called for the HF2 JOT-II potential'
      write(6,*) 'energy surface.  Scalar version (optimized) including'
      write(6,*) 'analytic first derivatives in internal coordinates.'
      write(6,*) 'Last modified 7 February 1997 (TCA).'
      write(6,1000) de(1),de(2),de(3),re(1),re(2),re(3),b(1),b(2),b(3),
     &              s(1),s(2),s(3),cevau
c
c Calculate first derivatives (default).  If no first derivatives are
c desired, nder should be set equal to zero after prepot is called.
c
      nder=1
c
c Convert to atomic units and calculate some useful constants.
c
      do 10 i=1,3
        de(i)=de(i)/cevau
        re(i)=re(i)/canau
        b(i)=b(i)*canau
        pf=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
        hde(i)=0.5d0*de(i)
        hpf(i)=0.5d0*pf
        pfs(i)=-de(i)*b(i)
        pft(i)=-pf*b(i)
10    continue

      easyab=de(1)
      easybc=de(2)
      easyac=de(3)

 1000 format(/,10x,'H-I  ',7x,'I-I'' ',7x,'H-I'' ',
     &       /,1x,'De',5x,3(f8.4,4x),'eV',
     &       /,1x,'Re',5x,3(f8.4,4x),'angstrom',
     &       /,1x,'B ',5x,3(f8.4,4x),'angstrom**-1',
     &       /,1x,'S ',5x,3(f8.4,4x),
     &      //,1x,'1.00 hartree = ',f11.7,' eV',
     &       /,1x,'1.00 angstrom = ',f11.9,' bohr',/)

      return

      entry pot
c
c compute coulomb and exchange terms
c
        y1=dexp(b(1)*(re(1)-r(1)))
        hs1=hde(1)*(y1-2.d0)*y1
        ht1=hpf(1)*(y1+2.d0)*y1
        wq1=hs1+ht1
        wj1=hs1-ht1
        y2=dexp(b(2)*(re(2)-r(2)))
        hs2=hde(2)*(y2-2.d0)*y2
        ht2=hpf(2)*(y2+2.d0)*y2
        wq2=hs2+ht2
        wj2=hs2-ht2
        y3=dexp(b(3)*(re(3)-r(3)))
        hs3=hde(3)*(y3-2.d0)*y3
        ht3=hpf(3)*(y3+2.d0)*y3
        wq3=hs3+ht3
        wj3=hs3-ht3
c
c compute potential energy
c
        rad=dsqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
        energy=de(2)+wq1+wq2+wq3-rad
c
c  compute first derivatives if necessary
c
        if (nder .eq. 1) then
          radi=1.d0/rad
          hradi=0.5d0*radi
          hds1dr1=pfs(1)*(y1-1.d0)*y1
          hds2dr2=pfs(2)*(y2-1.d0)*y2
          hds3dr3=pfs(3)*(y3-1.d0)*y3
          hdt1dr1=pft(1)*(y1+1.d0)*y1
          hdt2dr2=pft(2)*(y2+1.d0)*y2
          hdt3dr3=pft(3)*(y3+1.d0)*y3
          dq1dr1=hds1dr1+hdt1dr1
          dq2dr2=hds2dr2+hdt2dr2
          dq3dr3=hds3dr3+hdt3dr3
          dj1dr1=hds1dr1-hdt1dr1
          dj2dr2=hds2dr2-hdt2dr2
          dj3dr3=hds3dr3-hdt3dr3
          dedr(1)=dq1dr1-hradi*dj1dr1*(2.d0*wj1-wj2-wj3)
          dedr(2)=dq2dr2-hradi*dj2dr2*(2.d0*wj2-wj1-wj3)
          dedr(3)=dq3dr3-hradi*dj3dr3*(2.d0*wj3-wj1-wj2)
        end if

      return

      end
