C
      subroutine prepot
C
C   System:          F2H2
C   Common name:     SQSBDE
C   Number of derivatives: 0
C   Number of bodies: 4
C   Number of electronic surfaces: 1
C   Interface:       4-1S
C   Data file:       ga.in.sqsbde
C   References:      M. Quack, M. Suhm Mol. Phys. 69 (1990) P791
C                    Kofranek et al Chem. Phys. 121 (1988) 137
C                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507
C
C   Protocol:
C      PREPOT - initializes the potential's variables 
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy 
C
C   Units: 
C      energies    - cm-1
C      coordinates - distances: bohr; angles: degrees
C
C   Surfaces: 
C      ground electronic state
C
C
C M. A. SUHM  SEPTEMBER 1990
C
C INPUT : SIX INTERNAL COORDINATES 
C --R1,R2,R IN ATOMIC UNITS (BOHR)
C --TH1,TH2,TAU IN DEGREES
C DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791
C OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN CM-1
C 
C THE PARAMETERS OF THE POTENTIAL ARE READ IN FROM UNIT 1
C (FILENAME GA.IN)
C
C THE 6-D POTENTIAL IS BASED ON 1070 AB INITIO POINTS FROM
C KOFRANEK ET AL (CHEM PHYS 121(1988) 137) UP TO 42000 CM-1
C AND A PERTURBATION TREATMENT BY RIJKS AND WORMER
C (JCP 90(1989) 6507) AS WELL AS MONOMER We,WeXe FROM 
C HUBER/HERZBERG. 29 FREE, 15 MONOMER OR AB INITIO
C BASED AND 18 CONSTRAINED PARAMETERS ARE USED.
C SQSBDE IS AN EMPIRICALLY REFINED VERSION WHICH HAS THE
C CORRECT EXPERIMENTAL HYDROGEN BOND WAVENUMBER OF
C 1062CM-1 (MILLER, ACC. CHEM. RES. 23 (1990) P10) AND
C B ROTATIONAL CONSTANT. FOR THIS PURPOSE, THE AB INITIO
C RAB COORDINATES WERE SCALED BY 1/1.035.
C THE RMS OF THE WEIGHTED FIT (W=1 FOR R1,2=/=1.7374 OR
C R1,2=1.7374 AND E<2220CM-1, W=0.1 FOR R1,2=1.7374 AND 
C 2220<E<3210CM-1, W=0.001 FOR R1,2=1.7374 AND E>3210CM-1)
C IS 32.8CM-1. 
C
C E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS
C THE ABSOLUTE MINIMUM OF THE POTENTIAL
C
C
      subroutine prepot
C
C M. A. SUHM  SEPTEMBER 1989
C
C preparatory calculations for gssqs-routine
C
      implicit double precision (a-h,o-z)
      dimension fac(0:2*4)
      dimension fct(401),efct(401)
      parameter(ma=135)
      common/ang/ g(0:4,0:4,0:2*4,0:4),
     *            pf(0:4,0:4)
      data pi/3.1415926535898/
C
C prefactors for polynomials
C
      open(unit=4,file='ga.in.sqsbde',status='old')
      do 200 i=1,ma
       read(4,*) a(i),da(i)
  200 continue
C
      fac(0)=1.
      do 101 i=1,2*4
      fac(i)=fac(i-1)*i
  101 continue
      do 102 l=0,4
      do 103 m=0,l
      pf(l,m)=(-1)**m*sqrt((2*l+1)*fac(l-m)/(4*pi*fac(l+m)))
  103 continue
  102 continue
C 
C 3j-symbols for potential energy surface
C (abc/d-d0)=g(a,b,c,d)
C
C
      fct(1)=1.0d0
      efct(1)=0.0d0
      do 104 i=1,400
      x=i
      j=i+1
      efct(j)=efct(i)
      fct(j)=fct(i)*x
    1  if(fct(j).ge.10.0d0) then
        fct(j)=0.1d0*fct(j)
        efct(j)=efct(j)+1.0d0
        goto 1
       endif
  104 continue
C
C
      m2=0
      m3=0
      do 108 j1=0,4*2,2
      do 107 j2=0,4*2,2
      do 106 j3=iabs(j1-j2),iabs(j1+j2),4
      do 105 m1=0,min(j1,j2),2
      m2=-m1
      i1=max0(0,(j2-j3-m1)/2,(j1-j3+m2)/2)
      i2=min0((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)
      ia=i1
      ib=(j1+j2-j3)/2-i1
      ic=(j1-m1)/2-i1
      id=(j2+m2)/2-i1
      ie=(j3-j2+m1)/2+i1
      if=(j3-j1-m2)/2+i1
      sum=(-1)**i1/(fct(ia+1)*fct(ib+1)*fct(ic+1)*fct(id+1)*
     *    fct(ie+1)*fct(if+1))
      esum=-efct(ia+1)-efct(ib+1)-efct(ic+1)-efct(id+1)-
     *    efct(ie+1)-efct(if+1)
      nmax=i2-i1
      f=1d0
      if(nmax.gt.0) then
      do 109 i=nmax,1,-1
      x=(ib-i+1)*(ic-i+1)*(id-i+1)
      x=x/((ia+i)*(ie+i)*(if+i))
      f=1d0-x*f
  109 continue
      endif
      sum=sum*f
      f=fct((j1+j2-j3)/2+1)*fct((j1-j2+j3)/2+1)*
     * fct((j2-j1+j3)/2+1)*fct((j1+m1)/2+1)*fct((j1-m1)/2+1)*
     * fct((j2+m2)/2+1)*fct((j2-m2)/2+1)*fct((j3+m3)/2+1)*
     * fct((j3-m3)/2+1)/fct((j1+j2+j3)/2+2)
      e=efct((j1+j2-j3)/2+1)+efct((j1-j2+j3)/2+1)+
     * efct((j2-j1+j3)/2+1)+efct((j1+m1)/2+1)+efct((j1-m1)/2+1)+
     * efct((j2+m2)/2+1)+efct((j2-m2)/2+1)+efct((j3+m3)/2+1)+
     * efct((j3-m3)/2+1)-efct((j1+j2+j3)/2+2)
      cg3j=sum*dsqrt(f)*10d0**(esum+0.5d0*e)*(-1)**((j2-j1-m3)/2)
      g(j1/2,j2/2,j3/2,m1/2)=cg3j
  105 continue
  106 continue
  107 continue
  108 continue
      return
      end
C
C M. A. SUHM  SEPTEMBER 1989
C
      subroutine gsrsqs(r,ac,ar)
      parameter (ndmax=1100)
      implicit double precision (a-h,o-z)
      dimension ar(0:4,0:4,0:2*4)
      dimension ac(3,0:4,0:4,0:2*4)
      dab1=r-3.d0
      dab2=r-4.d0
      dab3=r-4.9d0
      do 301 j1=0,4,1
      do 302 j2=j1,4,1
      do 303 j3=(j2-j1),(j1+j2),2
      ar(j1,j2,j3)=ac(1,j1,j2,j3)*exp(-0.20d0*dab1**2)+
     *            ac(2,j1,j2,j3)*exp(-0.35d0*dab2**2)+
     *            ac(3,j1,j2,j3)*exp(-0.80d0*dab3**2)
      ar(j2,j1,j3)=(-1)**(j1+j2)*ar(j1,j2,j3)
  303 continue
  302 continue
  301 continue
      return
      end
C
C M. A. SUHM  SEPTEMBER 1989
C
      subroutine pot(r1,r2,r,th1,th2,tau,e)
      parameter(ma=135)
      implicit double precision(a-h,o-z)
      dimension x(4),y(4),z(4)
      dimension p1(0:4,0:4),p2(0:4,0:4),cs(4)
      dimension ar(0:4,0:4,0:2*4)
      dimension ac(3,0:4,0:4,0:2*4)
      dimension vp(0:4,0:4,0:2*4)
      common/ang/ g(0:4,0:4,0:2*4,0:4),
     *            pf(0:4,0:4)
      common/par/ a(ma)
      data amh/1837.152756d0/amf/34631.9678d0/
      phicon=0.01745329252
      pi=4.d0*atan(1.d0)
      pif=1./(2.*sqrt(pi))
      fac=amh/(amh+amf)
      x(3)=0.d0-r1*cos(th1)*fac
      y(3)=0.d0-r1*sin(th1)*fac
      z(3)=0.d0
      x(1)=r1*cos(th1)*(1.-fac)
      y(1)=r1*sin(th1)*(1.-fac)
      z(1)=0.d0
      x(4)=r-r2*cos(th2)*fac
      y(4)=0.d0-r2*sin(th2)*cos(tau)*fac
      z(4)=0.d0-r2*sin(th2)*sin(tau)*fac
      x(2)=r+r2*cos(th2)*(1.-fac)
      y(2)=0.d0+r2*sin(th2)*cos(tau)*(1.-fac)
      z(2)=0.d0+r2*sin(th2)*sin(tau)*(1.-fac)
      r12=sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)
      r23=sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)
      r14=sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)
      r34=sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)
      do 101 j=1,4
      cs(j)=cos(j*tau)
  101 continue
      s=sin(th1)
      c=cos(th1)
      s2=s*s
      s3=s2*s
      s4=s3*s
      s5=s4*s
      c2=c*c
      c3=c2*c
      c4=c3*c
      c5=c4*c
      p1(0,0)=pf(0,0)
      p1(1,0)=pf(1,0)*c
      p1(2,0)=pf(2,0)*0.5*(3*c2-1.)
      p1(3,0)=pf(3,0)*0.5*(5*c3-3*c)
      p1(4,0)=pf(4,0)*0.125*(35*c4-30*c2+3.)
      p1(1,1)=pf(1,1)*s
      p1(2,1)=pf(2,1)*3*s*c
      p1(3,1)=pf(3,1)*1.5*s*(5*c2-1)
      p1(4,1)=pf(4,1)*2.5*s*(7*c3-3*c)
      p1(2,2)=pf(2,2)*3*s2
      p1(3,2)=pf(3,2)*15*s2*c
      p1(4,2)=pf(4,2)*7.5*s2*(7*c2-1)
      p1(3,3)=pf(3,3)*15*s3
      p1(4,3)=pf(4,3)*105*s3*c
      p1(4,4)=pf(4,4)*105*s4
      s=sin(th2)
      c=cos(th2)
      s2=s*s
      s3=s2*s
      s4=s3*s
      s5=s4*s
      c2=c*c
      c3=c2*c
      c4=c3*c
      c5=c4*c
      p2(0,0)=pf(0,0)
      p2(1,0)=pf(1,0)*c
      p2(2,0)=pf(2,0)*0.5*(3*c2-1.)
      p2(3,0)=pf(3,0)*0.5*(5*c3-3*c)
      p2(4,0)=pf(4,0)*0.125*(35*c4-30*c2+3.)
      p2(1,1)=pf(1,1)*s
      p2(2,1)=pf(2,1)*3*s*c
      p2(3,1)=pf(3,1)*1.5*s*(5*c2-1)
      p2(4,1)=pf(4,1)*2.5*s*(7*c3-3*c)
      p2(2,2)=pf(2,2)*3*s2
      p2(3,2)=pf(3,2)*15*s2*c
      p2(4,2)=pf(4,2)*7.5*s2*(7*c2-1)
      p2(3,3)=pf(3,3)*15*s3
      p2(4,3)=pf(4,3)*105*s3*c
      p2(4,4)=pf(4,4)*105*s4
      icu=31
      do 102 j1=0,4,1
      do 103 j2=j1,4,1
      do 104j3=(j2-j1),(j1+j2),2
      ac(1,j1,j2,j3)=a(icu)
      ac(2,j1,j2,j3)=a(icu+35)
      ac(3,j1,j2,j3)=a(icu+70)
      icu=icu+1
  104 continue
  103 continue
  102 continue
      call gsrsqs(r,ac,ar)
      do 105 j1=0,4,1
      do 106 j2=0,4,1
      do 107 j3=abs(j2-j1),(j2+j1),2
      vp(j1,j2,j3)=g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)
      do 108 m=1,min(j1,j2),1
      vp(j1,j2,j3)=vp(j1,j2,j3)+(-1)**m*2*g(j1,j2,j3,m)*
     *            cs(m)*p1(j1,m)*p2(j2,m)
  108 continue
      vp(j1,j2,j3)=vp(j1,j2,j3)*(-1)**(j1-j2)*(2*j3+1)*pif
  107 continue
  106 continue
  105 continue
C
      e=0.0d0
      do 109 j1=0,4,1
      do 110 j2=0,4,1
      do 111 j3=abs(j2-j1),(j2+j1),2
      e=e+vp(j1,j2,j3)*ar(j1,j2,j3)
  111 continue
  110 continue
  109 continue
C
c     e=e+a(1)
      if(r.lt.7.2)then
      e=e+vp(0,0,0)*exp(-((7.2-r)/r)**2)*a(2)/r**6
      else
      e=e+vp(0,0,0)*a(2)/r**6
      endif
      if(r.lt.8.4)then
      e=e+vp(0,0,0)*exp(-((8.4-r)/r)**2)*a(3)/r**8
      else
      e=e+vp(0,0,0)*a(3)/r**8
      endif
      if(r.lt.9.6)then
      e=e+vp(0,0,0)*exp(-((9.6-r)/r)**2)*a(4)/r**10
      else
      e=e+vp(0,0,0)*a(4)/r**10
      endif
      if(r.lt.9.5) then
      e=e+vp(1,1,2)*a(5)/r**3*exp(-((9.5-r)/r)**2)
      else
      e=e+vp(1,1,2)*a(5)/r**3
      endif
      if(r.lt.10.5) then
      e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/r**4*exp(-((10.5-r)/r)**2)
      else
      e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/r**4
      endif
      if(r.lt.11.5) then
      e=e+vp(2,2,4)*a(7)/r**5*exp(-((11.5-r)/r)**2)
      e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/r**5*exp(-((11.5-r)/r)**2)
      else
      e=e+vp(2,2,4)*a(7)/r**5
      e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/r**5
      endif
      e=e+a(9)*exp(-0.5*r12)
      e=e+a(10)*exp(-1.5*r12)
      e=e+a(13)*exp(-1.5*r34)
      e=e+a(17)*(exp(-10*r14)+exp(-10*r23))
      if(r.lt.4.5)then
      e=e+a(18)/r**6*((3*(cos(th1))**2+1)+
     *               (3*(cos(th2))**2+1))*
     *                exp(-((4.5-r)/r)**2)
      else
      e=e+a(18)/r**6*((3*(cos(th1))**2+1)+
     *               (3*(cos(th2))**2+1))
      endif
C
C monomer terms
C
      e=e+47635d0*(1.-a(19)*exp(-0.1*(r23-3.)**2))*
     *         (1.-exp(-1.1953*
     *         (1.d0-a(26)*(1.d0-tanh(a(27)*(r14-3.)))
     *                                                      )*
     *         (r1-1.7374*(1.d0
     *         -a(20)*
     *         (1.-tanh(a(21)*(r14-3.)))
     *                                   ))))**2
      e=e+47635d0*(1.-a(19)*exp(-0.1*(r14-3.)**2))*
     *         (1.-exp(-1.1953*
     *         (1.d0-a(26)*(1.d0-tanh(a(27)*(r23-3.)))
     *                                                      )*
     *         (r2-1.7374*(1.d0
     *         -a(20)*
     *         (1.-tanh(a(21)*(r23-3.)))
     *                                   ))))**2
c
c     next line corrects for asymptotic diatomic potential which is
c     added in separately.
c
      e=e-(ptntl(r1-1.7374d0)+ptntl(r2-1.7374d0))*219474.63067d0
      e=e+a(28)*(exp(-4*(r14-r1)**2)+exp(-4*(r23-r2)**2))
      e=e+a(22)*exp(-4*r14)*exp(-4*r23)
      e=e+a(23)
      e=e+a(24)*(exp(-r14)+exp(-r23))
     *        *(exp(-15*r1)+exp(-15*r2))
      e=e/219474.63067d0
      return
      end
      FUNCTION PTNTL(X)
C
C     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING
C     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.
C     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND
C     IN CM-1 IN THE ICP REFERENCE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
      PARAMETER (AM=1.1953D0,DM=0.21704,ONE=1.D0)
C
      PTNTL=DM*(ONE-EXP(-AM*X))**2
C
      RETURN
      END      
