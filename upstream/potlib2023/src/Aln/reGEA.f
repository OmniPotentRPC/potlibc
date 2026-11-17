      subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Embedded atom
c   Common name:                   reGEA
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

c
c  calculate the potential energy

c new parameters
        a= 1.82725940400586229
        b= 0.701319003419638487
        cc= 0.260371275036638983
        v0= 90.0439667806546140

c Original
c      a=1.589647708
c      b=0.866262749
c      cc=0.4865568924
c      v0=44.8337922

      v1=0d0
          do 4 i=1,natom
      ud=0d0
          do 5 j=1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
          if(j.eq.i) go to 41
      udd=0d0
      udd=dexp(-b*rr)
      ud=ud+udd
   41 continue
    5 continue
      v1=v1+dsqrt(ud)
    4 continue
      u=0d0
          do 2 i=1,natom-1
          do 3 j=i+1,natom
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)                                                       
      uu=0d0
      uu=dexp(-a*rr)
      u=u+uu
    3 continue
    2 continue
      v=u*v0-v1*cc
      return
      end
c
