      subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Embedded atom
c   Common name:                   CleR
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  For a recent
c           discussion of this potential, see:
c           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c           and D. G. Truhlar, "Analytic Potentials Energy Functions
c           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  References: F. Cleri and V. Rosato, Phys. Rev. B 48, 22 (1993). 
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system (output)
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system (input)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximal number of atoms (input)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
c
c  calculate the potential energy
c
      v=0d0

      a=0.1221/autoev
      si=1.316/autoev
      p=8.612
      q=2.516
      a0=4.050/autoang

      r0=a0/dsqrt(2d0)

      u=0d0
          do 2 i=1,natom-1
          do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      rel=rr/r0
      ur=0d0
      ur=a*dexp(-p*(rel-1d0))
    3 u=u+ur
    2 continue

      ud=0d0
          do 5 i=1,natom
      udw=0d0
          do 4 j=1,natom
      urr=0d0
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      rel=rr/r0
          if(j.eq.i) go to 41
      urr=(si**2)*dexp(-2d0*q*(rel-1d0))
      udw=udw+urr  
   41 continue
    4 continue
      udn=dsqrt(udw)
    5 ud=ud-udn
      v=v+u+ud
      return
      end
