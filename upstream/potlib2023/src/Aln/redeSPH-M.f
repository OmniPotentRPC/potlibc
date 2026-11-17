      subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Morse
c   Common name:                   redeSPH/M
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Many-body aluminum potential energy function.  The parameters 
c           for this PEF have been parameterized against a data set of                             
c           aluminum cluster energies and bulk data.
c
c  References: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,                          
c              and D. G. Truhlar, "Analytic Potentials Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c              P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase
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

c  calculate the potential energy
      re = 3.28236443575964820  /autoang
      de = 0.928187591597459671 /autoev
      b0 = 0.855886663409868076 *autoang
      b2 = 0.755251587689301340 *autoang**3
      b3 = 1.47093307278944785  *autoang**4

      v=0.d0
c  1. calculate energy of two-body interactions Morse
      do 2 i=1,natom-1
      do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      dr = rr-re
      beta = b0 + b2*dr**2 + b3*dr**3
      xmorse = de*(1.d0-dexp(-beta*dr))**2-de
      v=v+xmorse
    3 continue
    2 continue
  
      return
      end
