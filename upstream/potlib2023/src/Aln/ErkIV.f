      subroutine pot(x,y,z,v,natom,maxatom)
c   System:                        Aln
c   Functional form:               Heuristic
c   Common name:                   ErkIV
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
c  References: S. Erkoc, Phys. Stat. Sol. (b) 161 (1990), 211. 
c	       S. Erkoc, Empirical Potential Energy Functions 
c	          Used in the Simulations of Materials Properties,  
c		  Annual Reviews of Computational Physics IX, 
c                 D. Stauffer, ed., World Scientific, Singapore 2001, 
c		  pp. 1-103.
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
      dimension x(natom),y(natom),z(natom),dist(natom,natom)
      dimension u(natom,natom)
      PARAMETER(autoang = 0.5291277249d0, autoeV = 27.2113961d0)

      en = 2.073d0
      r0 = 2.47d0/autoang
      aa = 6.20d0/autoev
      alpha = dlog(2.d0)
      c2 = 1.d0        ! A = c2*aa
      c3 = -1.553508d0 ! B = c3/c2

      vpair  = 0.0d0
      vmany  = 0.0d0
      r02   = r0**2
      en2  = en*2d0

c  pair potential, energy is denoted as vpair
      DO i=1,natom-1
       DO j=i+1,natom
        dx=x(i)-x(j)
        dy=y(i)-y(j)
        dz=z(i)-z(j)
        rr=sqrt(dx*dx+dy*dy+dz*dz)
        dist(i,j)=rr
        dist(j,i)=rr
        r1  = rr/r0
        r2  = r0/rr
        temp1=dexp(-2.d0*alpha*r1**2)
        temp2=dexp(-1.d0*alpha*r1**2)
        u(i,j) = aa*((r2**en2)*temp1 - (r2**en)*temp2)
        u(j,i) = u(i,j) 
        vpair = u(i,j) + vpair
       END DO
      END DO

      IF (natom .GE. 3) THEN
      DO i=1, natom-2
       DO j=i+1,natom-1
        DO k=j+1,natom
        fijk = dexp(-(dist(i,k)**2 + dist(j,k)**2)/r02)
        fikj = dexp(-(dist(i,j)**2 + dist(k,j)**2)/r02)
        fjki = dexp(-(dist(j,i)**2 + dist(k,i)**2)/r02)
        Wijk = (u(i,j)*fijk+u(i,k)*fikj+u(j,k)*fjki)
        vmany = Wijk + vmany
        END DO
       END DO
      END DO
      END IF

      v = c2*vpair + c3*vmany


      END
