      subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Embedded atom
c   Common name:                   PapKEP
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
c  References: N. I. Papanicolaou, G. C. Kallinteris, G. A. Evangelakis, 
c   D. A. Papaconstantopoulos, Computational Materials Science, 
c   17, 224, (2000).
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

c
      implicit real*8(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      dimension x(maxatom),y(maxatom),z(maxatom)

      r0 = 2.8119d0 / autoang
      a  = 0.0334d0 / autoev
      p  = 14.6147d0
      ep =  0.7981d0 / autoev
      q  = 1.1120d0

      v = 0.d0
      do 2 i=1,natom
      t1 = 0.d0
      t2 = 0.d0
      do 3 j=1,natom
      if (i.eq.j) go to 3
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

      dr = rr/r0 - 1.d0
      t1 = t1 + dexp(-1.d0*p*dr)
      t2 = t2 + dexp(-2.d0*q*dr)
    3 continue
      v = v + a*t1 - ep*dsqrt(t2)
    2 continue

      return
      end
