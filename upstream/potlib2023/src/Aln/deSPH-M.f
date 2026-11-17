c
      subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Morse
c   Common name:                   deSPH/M
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
c  References: P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase
c              J. Phys. Chem., 99, 8147 (1995).
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
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      parameter(autokpm=627.5095d0)
      dimension x(maxatom),y(maxatom),z(maxatom)

      re = 2.834d0  / autoang
      de = 12.933d0 / autokpm
      a0 = 0.5282d0 * autoang
      a2 = 4.390d-3 * autoang**3
      a3 = 5.109d0  * autoang**4

      v=0.d0
c     calculate energy of two-body interaction
      do 2 i=1,natom-1
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dr = rr-re
      aa = a0 + a2*dr**2 + a3*dr**3
      xmorse = de*(1.d0-dexp(-aa*dr))**2-de
      v=v+xmorse
    3 continue
    2 continue
  
      return
      end
