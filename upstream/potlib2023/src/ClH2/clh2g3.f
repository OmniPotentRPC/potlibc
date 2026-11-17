c
      subroutine prepot
c Last modified on 3 August 1995 by TCA.
c
c System:          ClH2
c Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)
c                  plus E3C (a three-center term)
c                  and modifications to improve the saddle point
c                  bend potential; based on GQQ surface
c Common name:     G3
c Reference:       Allison, T. C.; Lynch, G. C.; Truhlar, D. G.;
c                  Gordon, M. S. J. Phys. Chem., to be submitted.
c
c PREPOT must be called once before any calls to POT are made.  It should
c never be called more than once.
c The potential parameters are included in DATA statements.
c
c The coordinates are defined as follows:
c           R(1) = R(H-Cl)
c           R(2) = R(H-H')
c           R(3) = R(H'-Cl)
c
c The zero of energy for this potential occurs when the Cl atom is
c "infinitely" far from the the H2 diatom, and the H2 diatom is at
c its equilibrium seperation.
c
c
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
      implicit double precision (a-h,o-z)
      common /pt1cm/ r(3), energy, dedr(3)
      common /pt2cm/ nsurf, nder, ndum(21)
      common /pt5cm/ easyab, easybc, easyac
c
      dimension de(3),re(3),b(3),s(3),pf(3),pfs(3),pfd(3)
      save de,re,b,s,rj,ap,at,bp,q2,q4,ckcau,cangau,c1,c2
c
c Constants used by the code.  de is the dissociation energy (kcal/mol),
c re is the equilibrium bond length (angstroms), b is the Morse beta
c parameter (angstroms**-1), and s is the sato parameter (unitless) for
c the appropriate diatom (1 = H-Cl, 2 = H=H', 3 = H'-Cl).  The other
c constants may be found in the literature reference.  ckcau is the
c conversion factor between kcal/mol and hartrees and cangau is the
c conversion factor between angstroms and bohr.
c
      data de     / 106.447d0, 109.458d0, 106.447d0 /
      data re     / 1.2732d0, 0.74127d0, 1.2732d0 /
      data b      / 1.8674d0, 1.9413d0, 1.8674d0 /
      data s      / 0.1835d0, 0.167d0, 0.1835d0 /
      data rj     / 0.0758016022d0 /
      data ap     / 0.0008627355d0 /
      data at     / 0.2981969860d0 /
      data bp     / 0.1439106543d0 /
      data q2     / 0.6940323070d0 /
      data q4     / 1.6963092005d0 /
      data c1     / 3.1053877618397071175758280546d+0 /
      data c2     / -1.8942608491155350536449828878d+0 /
      data ckcau  / 627.5095d0 /
      data cangau / 0.52917706d0 /
c
c Write header potential parameters to FORTRAN unit 6.
c
      write(6,*)
      write(6,*) 'PREPOT has been called for the ClH2 potential energy'
      write(6,*) 'surface G3.  Scalar version (optimized) including'
      write(6,*) 'analytic first derivatives in internal coordinates.'
      write(6,*) 'Last modified 3 August 1995 (TCA).'
      write(6,1000) de(1),de(2),de(3),re(1),re(2),re(3),b(1),b(2),b(3),
     &              s(1),s(2),s(3),rj,ap,at,bp,q2,q4,c1,c2,ckcau,cangau
c
c Calculate first derivatives (default).  If no first derivatives are
c desired, nder should be set equal to zero after prepot is called.
c
      nder=1
c
c Convert to atomic units and calculate some useful constants.
c
      do 10 i=1,3
        de(i)=de(i)/ckcau
        re(i)=re(i)/cangau
        b(i)=b(i)*cangau
        pf(i)=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))
        pfs(i)=de(i)+pf(i)
        pfd(i)=de(i)-pf(i)
10    continue
c
      easyab=de(1)
      easybc=de(2)
      easyac=de(3)
c
1000  format(/,10x,'H-Cl ',7x,'H-H'' ',7x,'H''-Cl',
     &       /,1x,'De',5x,3(f7.3,5x),'kcal/mol',
     &       /,1x,'Re',5x,f8.4,4x,f9.5,3x,f8.4,4x,'kcal/mol',
     &       /,1x,'B ',5x,3(f8.4,4x),'angstrom**-1',
     &       /,1x,'S ',5x,3(f8.4,4x),
     &      //,1x,'rj = ',f11.7,3x,'hartrees',
     &       /,1x,'ap = ',f13.9,1x,'angstrom**-4',
     &       /,1x,'at = ',f10.6,4x,'angstrom**-2',
     &       /,1x,'bp = ',f10.6,4x,'angstrom**-2',
     &       /,1x,'q2 = ',f10.6,
     &       /,1x,'q4 = ',f9.5,
     &       /,1x,'c1 = ',f9.5,
     &       /,1x,'c2 = ',f9.5,
     &      //,1x,'1 hartree = ',f8.4,' kcal/mol',
     &       /,1x,'1 angstrom = ',f10.8,' bohr',/)
c
      return
c
      entry pot
c
      r1i=1.d0/r(1)
      r2i=1.d0/r(2)
      r3i=1.d0/r(3)
      r1s=r(1)*r(1)
      r2s=r(2)*r(2)
      r3s=r(3)*r(3)
c
c Compute three-center term.
c
      f1=dexp(-ap*(r(1)+r(3))**4)
      f2=dexp(-at*(r(1)-r(3))**2)
      f3=dexp(-bp*(r(1)+r(3)-r(2))**2)
      ctheta=(r1s+r3s-r2s)*(0.5d0*r1i*r3i)
      if (ctheta .lt. -1.d0) then
        ctheta=-1.d0
      else if (ctheta .gt. 1.d0) then
        ctheta=1.d0
      end if
      stheta=dsin(dacos(ctheta))
      st2=stheta*stheta
      f4=1.d0+q2*st2+q4*st2*st2
c
c Compute modification to bending potential.
c
      x=(r(1)-r(3))*r2i
      x2=x*x
      g5=(0.5d0*(1.d0+ctheta))**5
      g=g5*(0.5d0*(1.d0+ctheta))
      fmod=1.d0+g*(c1*(1.d0-x2)+c2*(1.d0-x2*x2))
c
c Compute coulomb and exchange terms.
c
      y1=dexp(b(1)*(re(1)-r(1)))
      wq1=(0.5d0*(de(1)+pf(1)*fmod)*y1-pfd(1))*y1
      wj1=(0.5d0*(de(1)-pf(1)*fmod)*y1-pfs(1))*y1
      y2=dexp(b(2)*(re(2)-r(2)))
      wq2=(0.5d0*pfs(2)*y2-pfd(2))*y2
      wj2=(0.5d0*pfd(2)*y2-pfs(2))*y2
      y3=dexp(b(3)*(re(3)-r(3)))
      wq3=(0.5d0*(de(3)+pf(3)*fmod)*y3-pfd(3))*y3
      wj3=(0.5d0*(de(3)-pf(3)*fmod)*y3-pfs(3))*y3
c
c Compute potential energy.
c
      rad=dsqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))
      energy=de(2)+wq1+wq2+wq3-rad+f4*(rj*f1*f2*f3)
c
      if (nder .eq. 1) then
c
c Compute contributions to the derivatives.
c
        dfrr=(1.d0-x2)*(c1+c2*(1.d0+x2))*3.d0*g5
        dfrs=2.d0*x*g*(c1+2.d0*c2*x2)*r2i
        dfrt1=r3i-ctheta*r1i
        dfrt2=-r(2)*(r1i*r3i)
        dfrt3=r1i-ctheta*r3i
        dfr1=dfrr*dfrt1-dfrs
        dfr2=dfrr*dfrt2+dfrs*x
        dfr3=dfrr*dfrt3+dfrs
c
        y1s=y1*y1
        y3s=y3*y3
c
        dj1r1=y1*(de(1)*b(1)*(1.d0-y1)
     &        -pf(1)*(-b(1)*(1.d0+fmod*y1)+0.5d0*dfr1*y1))
        dj1r2=-0.5d0*pf(1)*dfr2*y1s
        dj1r3=-0.5d0*pf(1)*dfr3*y1s
        dj2r2=b(2)*y2*(pfs(2)-pfd(2)*y2)
        dj3r1=-0.5d0*pf(3)*dfr1*y3s
        dj3r2=-0.5d0*pf(3)*dfr2*y3s
        dj3r3=y3*(de(3)*b(3)*(1.d0-y3)
     &        -pf(3)*(-b(3)*(1.d0+fmod*y3)+0.5d0*dfr3*y3))
c
        df1r1=-2.d0*at*f1*(r(1)-r(3))
        df2r1=-4.d0*ap*f2*(r(1)+r(3))**3
        df3r1=-2.d0*bp*f3*(r(1)+r(3)-r(2))
        df4rr=2.d0*ctheta*(q2+2.d0*q4*st2)
c
        f123=f1*f2*f3
        f124=f1*f2*f4
        f134=f1*f3*f4
        f234=f2*f3*f4
        wj12d=wj1-wj2
        wj23d=wj2-wj3
        wj31d=wj3-wj1
        hradi=0.5d0/rad
c
c Assemble derivatives.
c
        dedr(1)=y1*(de(1)*b(1)*(1.d0-y1)+pf(1)*(-b(1)*(1.d0+fmod*y1)
     &   +0.5d0*dfr1*y1))+0.5d0*pf(3)*dfr1*y3s
     &   -hradi*(wj12d*dj1r1-wj23d*dj3r1+wj31d*(dj3r1-dj1r1))
     &   +rj*(f124*df3r1+f134*df2r1+f234*df1r1-f123*df4rr*dfrt1)
        dedr(2)=0.5d0*pf(1)*dfr2*y1s+b(2)*y2*(pfd(2)-pfs(2)*y2)
     &   +0.5d0*pf(3)*dfr2*y3s
     &   -hradi*(wj12d*(dj1r2-dj2r2)+wj23d*(dj2r2-dj3r2)
     &   +wj31d*(dj3r2-dj1r2))-rj*(f123*df4rr*dfrt2+f124*df3r1)
        dedr(3)=0.5d0*pf(1)*dfr3*y1s+y3*(de(3)*b(3)*(1.d0-y3)
     &   +pf(3)*(-b(3)*(1.d0+fmod*y3)+0.5d0*dfr3*y3))
     &   -hradi*(wj12d*dj1r3-wj23d*dj3r3+wj31d*(dj3r3-dj1r3))
     &   +rj*(f124*df3r1+f134*df2r1-f234*df1r1-f123*df4rr*dfrt3)
      end if
c
      return
c
      end
