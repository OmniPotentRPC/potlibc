             subroutine pot(x,y,z,v,natom,maxatom)

c   System:                        Aln
c   Functional form:               Embedded atom
c   Common name:                   Gol
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
c  References: H. Gollisch, Surface Science 166 (1986), 87.
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

c  calculate and tabulate the distances between all atoms

      aa = 350.415d0  / autoev
      bb = 10.3812d0  / autoev
      a  = 2.7338d0 *  autoang
      b  = 1.3681 * autoang
      en = 0.6d0

      v=0.d0
      do i=1,natom
      u=0.d0
      ud=0.d0    
      do j=1,natom
      if(i.ne.j) then
       dx=x(i)-x(j)
       dy=y(i)-y(j)
       dz=z(i)-z(j)
       rr=dsqrt(dx*dx+dy*dy+dz*dz)
       u=u+dexp(-rr*a)
       ud=ud+dexp(-rr*b)
      endif
      enddo
      v=v+aa*u-bb*ud**en
      enddo
      return
      end
