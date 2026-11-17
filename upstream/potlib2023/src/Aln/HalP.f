             subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Lennard-Jones
c   Common name:                   HalP
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c   Notes:  Many-body aluminum potential energy function.  For a recent
c           discussion of this potential, see:
c           A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c           and D. G. Truhlar, "Analytic Potentials Energy Functions 
c           for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
c
c  References: T. Halicioglu and G.M. Pound, Phys. Stat. Sol. (b)
c		  30, 619 (1975). 
c
c  Units:
c	energies    - hartrees
c	coordinates - bohrs
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


      sig=2.620d0/autoang
      de=4.d0*0.3922d0/autoev
      en=6.000d0
      em=12.000d0
c      re=sig*(em/en)**(1.d0/(em-en))

      v1=0.d0
      do 1 i=1,natom
      do 1 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      vem=(sig/rr)**em
      ven=(sig/rr)**en
      v1=v1+vem-ven
    1 continue
      v=de*v1
      return
      end
