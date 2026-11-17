C   System:          H3plus
C   Functional form: DIM
C   Common name:     H3+ KBNN
C   Number of derivatives: 0
C   Number of bodies: 3
C   Number of electronic surfaces: 3
C   Interface: 3-1V
C
C   References: V. G. Ushakov, K. Nobusada, and V. I. Osherov, 
c               Phys. Chem. Chem. Phys., Vol. 3, 63 (2001);
c               H. Kamisaka, W. Bian, K. Nobusada, and H. Nakamura,
c               J. Chem. Phys. Vol. 116, 654 (2002).
C
c   Notes: This is a 3x3 diabatic representation of the h3+ system.
c          This potential routine was coded from the functional forms and 
c          parameters as described in the publications listed above 
c          (see the References section), but it may not be exactly
c          equivalent to the potential that was used in the dynamics 
c          calculations of the Kamisaka paper (private communication, 2002).
c
c     Coded May 28, 2002 by A. G. Anderson, A. W. Jasper, 
c                        and D. G. Truhlar
c
c     definitions for parameters:
c     r -> array with dimensions R(nt,3)
c          R(i,1) is distance between atoms 1 and 2 (A-B)
c          R(i,2) is distance between atoms 2 and 3 (B-C)
c          R(i,3) is distance between atoms 1 and 3 (A-C)
c          where 1 <= i <= nt
c          and distances are in atomic units
c     De -> array of the derivatives of the diabatic potential
c     nt -> number of geometries to be computed. in this program, nt is 
c          always = 1
c     e ->  array dimension nt containing the diabatic potential energy(ies)
c     nsurf -> which surface we want the energy(ies) for. 
c              nsurf = 1: U11
c              nsurf = 2: U12
c              nsurf = 3: U22
c              nsurf = 4: U13
c              nsurf = 5: U23
c              nsurf = 6: U33
c
c     This potential and uses the functional forms of 
c
c     V. G. Ushakov, K. Nobusada, and V. I. Osherov, 
c     Phys. Chem. Chem. Phys., Vol. 3, 63 (2001)
c
c     and the updated parameters in
c
c     H. Kamisaka, W. Bian, K. Nobusada, and H. Nakamura,
c     J. Chem. Phys. Vol. 116, 654 (2002).

*******************************************
      subroutine prepot

      implicit double precision(a-h,o-z)
      dimension r(nt,3),De(3,nt),e(nt)

c     Table 2: fitting parameters of eqns 10 and 11
      double precision c1(3),c2(3),c3(3),c4(3),c5(3),c6(3),c7(3)
      double precision c8(3),c9(3),c10(3),c11(3),c12(3),c13(3),c14(3)
      double precision p(3,3)
      double precision c(3,3)
      double precision osg(3),tsg(3),tsu(3)
      double precision h(3),g(3),u(3)
      double precision potential(3,3)

c     Table 2: fitting params for eqns 10 and 11
      data c1  /-2.03221d0, -0.51486d0, -3.14355d0/
      data c2  / 1.46569d0,  0.08741d0,  1.39001d0/
      data c3  / 1.66412d0, -0.12057d0, -0.11632d0/
      data c4  / 0.27272d0, -0.02033d0, -0.28688d0/
      data c5  / 0.06060d0,  0.08057d0, -0.42105d0/
      data c6  /-2.25049d0, -0.37742d0,  0.44070d0/
      data c7  / 0.01945d0,  0.05521d0,  0.46790d0/
      data c8  /-0.05753d0,  0.01096d0,  0.05605d0/
      data c9  /-0.01288d0, -0.00318d0, -0.05596d0/
      data c10 / 0.57271d0,  0.00257d0,  0.10729d0/
      data c11 /-0.07831d0,  0.00774d0,  0.16146d0/
      data c12 / 0.50812d0,  0.00916d0, -0.09256d0/
      data c13 /-0.07156d0, -0.02141d0, -0.15563d0/
      data c14 / 0.98170d0,  0.74157d0,  0.79817d0/

c     Table 1: diatomic parameters for eqns 7-9
      data d_h/-0.17445d0/, d_g/0.10272d0/, d_u/0.62871d0/
      data b_h/ 2.09607d0/, b_g/0.71935d0/, b_u/0.87080d0/
      data r_h/ 1.40104d0/, r_g/2.00313d0/, r_u/0.34405d0/

      data a/1.13804d0/, b/0.52961d0/

      save

      write(6,*)'Prepot called for h3+KBNN surface'
      write(6,*)'Version finalized May 28, 2002'

      return

*******************************************

      entry pot(r,e,de,nt,nsurf)

      do ii=1,nt

      r1 = r(ii,1)
      r2 = r(ii,2)
      r3 = r(ii,3)

c     equation 11 (referring to the Ushakov paper)
      do L=1,3
       p(L,1) = c1(L)+c2(L)*r1+c4(L)*r1*r1+c5(L)*(r2+r3)**2+c6(L)*r2*r3+
     &          c7(L)*r1*(r2+r3)+c8(L)*r1**3+c9(L)*(r2+r3)**3+ 
     &          c10(L)*r2*r3*(r2+r3)+c11(L)*r1*(r2+r3)**2+ 
     &          c12(L)*r1*r2*r3+c13(L)*r1*r1*(r2+r3)+c3(L)*(r2+r3)

       p(L,2) = c1(L)+c2(L)*r2+c4(L)*r2*r2+c5(L)*(r3+r1)**2+c6(L)*r3*r1+
     &          c7(L)*r2*(r3+r1)+c8(L)*r2**3+c9(L)*(r3+r1)**3+ 
     &          c10(L)*r3*r1*(r3+r1)+c11(L)*r2*(r3+r1)**2+ 
     &          c12(L)*r2*r3*r1+c13(L)*r2*r2*(r3+r1)+c3(L)*(r3+r1)

       p(L,3) = c1(L)+c2(L)*r3+c4(L)*r3*r3+c5(L)*(r1+r2)**2+c6(L)*r1*r2+
     &          c7(L)*r3*(r1+r2)+c8(L)*r3**3+c9(L)*(r1+r2)**3+ 
     &          c10(L)*r1*r2*(r1+r2)+c11(L)*r3*(r1+r2)**2+ 
     &          c12(L)*r3*r1*r2+c13(L)*r3*r3*(r1+r2)+c3(L)*(r1+r2)
      enddo

c     equation 10
      do L=1,3
       c(L,1) = p(L,1)*dexp(-1.0d0*c14(L)*c14(L)*(r2+r3))
       c(L,2) = p(L,2)*dexp(-1.0d0*c14(L)*c14(L)*(r3+r1))
       c(L,3) = p(L,3)*dexp(-1.0d0*c14(L)*c14(L)*(r1+r2))
      enddo

c     equations 7,8,9
      do k=1,3
       osg(k) = dexp(-1.0d0*b_h*(r(ii,k)-r_h))*d_h*
     &          (1+b_h*(r(ii,k)-r_h)+a*(r(ii,k)-r_h)**2+
     &          b*(r(ii,k)-r_h)**3)
       tsg(k) = d_g*( dexp(-2.0d0*b_g*(r(ii,k)-r_g)) - 
     &          2.0d0*dexp(-1.0*b_g*(r(ii,k)-r_g)) + 
     &          0.0015d0*dexp(-8.0d0*b_g*(r(ii,k)-r_g)) )
       tsu(k) = d_u*( dexp(-2.0d0*b_u*(r(ii,k)-r_u)) + 
     &          2.0d0*dexp(-1.0d0*b_u*(r(ii,k)-r_u)) )
      enddo

c     eqns 4,5,6
      do i=1,3
       h(i) = osg(i) + c(1,i)
       g(i) = tsg(i) + c(2,i)
       u(i) = tsu(i) + c(3,i)
      enddo

c     this 'potential' 3X3 is the potential matrix
      potential(1,1) = h(1) + 0.5d0*( g(2) + u(2) + g(3) + u(3) )
      potential(2,2) = h(2) + 0.5d0*( g(1) + u(1) + g(3) + u(3) )
      potential(3,3) = h(3) + 0.5d0*( g(1) + u(1) + g(2) + u(2) )
      potential(1,2) = 0.5d0*( g(3) - u(3) )
      potential(1,3) = 0.5d0*( g(2) - u(2) )
      potential(2,1) = 0.5d0*( g(3) - u(3) )
      potential(2,3) = 0.5d0*( g(1) - u(1) )
      potential(3,1) = 0.5d0*( g(2) - u(2) )
      potential(3,2) = 0.5d0*( g(1) - u(1) )

      if (nsurf .eq. 1) e(ii) = potential(1,1) - d_h
      if (nsurf .eq. 2) e(ii) = potential(1,2)
      if (nsurf .eq. 3) e(ii) = potential(2,2) - d_h
      if (nsurf .eq. 4) e(ii) = potential(1,3)
      if (nsurf .eq. 5) e(ii) = potential(2,3)
      if (nsurf .eq. 6) e(ii) = potential(3,3) - d_h
      if (nsurf .gt. 6 .or. nsurf .lt. 1) then
               e(ii) = 0.d0
               write(6,*)'WARN:  NSURF = ',nsurf
      endif

c derivatives (zero for now)
       do i=1,3
        de(i,ii) = 0.d0
       enddo

c end nt loop
      enddo

      return
      end
*******************************************
