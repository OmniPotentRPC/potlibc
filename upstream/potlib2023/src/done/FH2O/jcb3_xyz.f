
      subroutine xyzTOjcb3(ct,rt)
      implicit none  
                     !   D  F
!           /|\         /
!     A  /\  |         /
!     H    \ |        /
!           \|       /        
!  ----------\------------------/ C  H
!             \ 
!              \  
!               \
!             O  B
!     
      real*8 ct(3,4),xAB(3),xABC(3),rt(6) !ct in HOH + F= HO + HF order
      real*8 xvec(3,6),crop1(3),crop2(3)
      real*8 mass(4)

      mass(1)=1.00783d0
      mass(2)=15.99491d0
      mass(3)=1.00783d0
      mass(4)=18.99840d0
      
      xAB(:)=(mass(1)*ct(:,1)+mass(2)*ct(:,2))/(mass(1)+mass(2))
      xABC(:)=(mass(1)*ct(:,1)+mass(2)*ct(:,2)+mass(3)*ct(:,3))/
     $ (mass(1)+mass(2)+mass(3))

      xvec(:,1)=ct(:,1)-ct(:,2)
      xvec(:,2)=ct(:,3)-xAB(:)
      xvec(:,3)=ct(:,4)-xABC(:)
      
      rt(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1))) 
      rt(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2))) 
      rt(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3))) 

      rt(4)=dot_product(xvec(:,1),xvec(:,2))/rt(1)/rt(2)
      rt(5)=dot_product(xvec(:,2),xvec(:,3))/rt(2)/rt(3)

      rt(4)=dacos(max(-1.0d0,min(1.0d0,rt(4))))
      rt(5)=dacos(max(-1.0d0,min(1.0d0,rt(5))))
      
      call crossp(xvec(:,1),xvec(:,2),crop1) 
      call crossp(xvec(:,2),xvec(:,3),crop2) 
      
      rt(6)=dot_product(crop1,crop2)/dsqrt(dot_product(crop1,crop1))/
     $dsqrt(dot_product(crop2,crop2))
      
      rt(6)=dacos(max(-1.0d0,min(1.0d0,rt(6))))

      rt(6)=dacos(-1.d0)-rt(6) 

      return
      end subroutine xyzTOjcb3


!******************************************************************************
c     subroutine crossp(v1,v2,v3)
c     real(kind=8)::v1(1:3),v2(1:3),v3(1:3)
c     ! cross product of two 3-dimensional vectors
c     v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
c     v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
c     v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
c     return
c     end subroutine crossp

!******************************************************************************
      subroutine jcb3tocart(rint,cart)
      implicit none
      real*8 rint(6),cart(3,4)
      real*8 r1,r2,r3,th1,th2,phia,sth1,sth2,sth3,cth1,cth2,cth3
      real*8 wa,wb,wc,wd,abm,cdm,ama,amb,amc,amd
!-->-->rint(1-6)=r1,r2,r3,th1,th2,phi
!-->-->r1 is vector A->B, r2 is a vector pointing from the COM of AB to
!C
!-->-->r3 is a vector pointing from the COM of ABC to D
!-->-->th1 is the angle for r1 and r2
!-->-->th2 is the angle for r2 and r3
!-->-->phia is the diheral angle for r1*r2 plane and r2*r3 plane.
!-->-->cart(3,1-4)=cartesian coordinates of four atoms A,B,C,D.
      real*8 temp1,temp2,temp3

      ama=1.00783d0
      amb=15.99491d0
      amc=1.00783d0
      amd=18.99840d0

      r1=rint(1)
      r2=rint(2)
      r3=rint(3)
      th1=rint(4)
      th2=rint(5)
      phia=rint(6)
!-->-->ABC+D Jacobi coordinates 
!-->-->            Schematic diagram
c          X     
c         /|\                 - D  F
c          |                -   /
c          |              -    /
c     H A  |            -r3   /
c        \ |          -      /
c         \th1      -th2    /
c----------e---------------C-------------->Z 
c        r1|\      r2    /   H
c          | \         /
c          |  \      /
c                B O             
!c          |   
c               
!-->-->origin is the center of mass of AB.
!-->-->r2 lies on Z positive axis, AB is in the XZ positive plane
!-->-->Xa is always positive and Xb is always negative.
!-->-->Xc is always positive
      abm=r1/(ama+amb)
      wa=amb*abm
      wb=ama*abm
      sth1=dsin(th1)
      cth1=dcos(th1)
      sth2=dsin(th2)
      cth2=dcos(th2)
      cth3=dcos(phia)
      sth3=dsin(phia)
      cart(1,1)=wa*sth1
      cart(2,1)=0d0
      cart(3,1)=wa*cth1
      cart(1,2)=-wb*sth1
      cart(2,2)=0d0
      cart(3,2)=-wb*cth1
      cart(1,3)=0d0
      cart(2,3)=0d0
      cart(3,3)=r2
      cart(1,4)=r3*sth2*cth3
      cart(2,4)=r3*sth2*sth3
      cart(3,4)=r3*cth2+amc*r2/(ama+amb+amc)
      end subroutine jcb3tocart 

!******************************************************************************

      subroutine cart_diJacobi(c1,r)
! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
!  c in angstrom and ABCD, r in angstrom and radians
      implicit none
      real*8 r(1:6),c(1:3,1:4),Mass(1:4),c1(1:3,1:4)
      real*8 crop1(1:3),crop2(1:3)
      real*8 vec(1:3,1:9),dotpp(1:9),bohr
      real*8 dotp,costh1,costh2,AMass,BMass,CMass,DMass,costh3

c     c(:,1)=c1(:,4)
c     c(:,2)=c1(:,1)
c     c(:,3)=c1(:,3)
c     c(:,4)=c1(:,2)
      c=c1  ! c in HOHF order angstrom
      bohr=1.0d0
      AMass=2.0141d0
      BMass=15.99491d0
      CMass=AMASS
      DMass=18.9984d0
      vec(:,1)=c(:,1)-c(:,2) !  vector B->A OH
      r(1)    =dsqrt( sum(vec(:,1)**2) )/bohr
      vec(:,2)=c(:,2)+vec(:,1)*AMass/(AMass+BMass) !vector BX (X COM of BA)
      vec(:,3)=c(:,3)-c(:,4) ! vector D->C OC
      r(2)    =dsqrt( sum(vec(:,3)**2) )/bohr
      vec(:,4)=c(:,4)+vec(:,3)*CMass/(CMass+DMass) !vector Y(Y COM of CD)
      vec(:,5)=vec(:,4)-vec(:,2) !vecor R X->Y
      r(3)    =dsqrt( sum(vec(:,5)**2) )/bohr
!-
      vec(:,6)=c(:,3)-c(:,2) ! vector  O->C
      vec(:,7)=c(:,2)-c(:,3) ! vector  C->O
      vec(:,8)=c(:,4)-c(:,3) ! vector  C->O'
!     dompp(5)=dsqrt(dot_product(vec(:,5),vec(:,5)))
!     dotpp(6)=dsqrt(dot_product(vec(:,6),vec(:,6)))
!     dotpp(7)=dsqrt(dot_product(vec(:,7),vec(:,7)))

      dotpp(1)=dsqrt(dot_product(vec(:,1),vec(:,1)))
      dotpp(3)=dsqrt(dot_product(vec(:,3),vec(:,3)))
      dotpp(5)=dsqrt(dot_product(vec(:,5),vec(:,5)))

      if (dotpp(1)*dotpp(5).lt. 1.0d-16)then
        costh1=0.0d0
      else
        costh1=dot_product(vec(:,1),vec(:,5))/(dotpp(1)*dotpp(5))
      endif
      if (dotpp(3)*dotpp(5).lt. 1.0d-16)then
        costh2=0.0d0
      else
        costh2=dot_product(vec(:,3),vec(:,5))/(dotpp(3)*dotpp(5))
      endif
!  note the round off errors
      r(4) = dacos( max(-1.0d0,min(1.0d0,costh1))) ! in radians
      r(5) = dacos( max(-1.0d0,min(1.0d0,costh2))) ! in radians

      if (dabs(dcos(r(4))*dcos(r(5))).lt.1.0d-16)then
c       write(*,*) '   ******'
c       write(*,'(3f12.6)')c
c      write(*,*)'some 3-atom collinear, dihedral angle cant be defined'
c       write(*,*) '   ******'
!       stop
!       return
      else
! compute the Jacobi dihedral angle, based on the norm of the plane
       call crossp(vec(:,1), vec(:,5),crop1) ! crossproduct for the norm
       call crossp(vec(:,3), vec(:,5),crop2) !
       if(dsqrt( sum(crop1(:)**2))*dsqrt(sum(crop2(:)**2)).lt.1.0d-10)
     &  then
        r(6) = 0.0d0
       else
       costh3=dot_product(crop1(:),crop2(:))/(dsqrt(sum(crop1(:)**2))
     $ *dsqrt(sum(crop2(:)**2)))
!    write(99999,'(f20.10)') costh3
       r(6) = dacos( max(-1.0d0,min(1.0d0,costh3)))
       endif
      endif
      return
      end subroutine cart_diJacobi 
!******************************************************************************

      subroutine diJacobi_cart(r,c)
! convert the diatom-diatom Jacobi coordinates to XYZ coordinates
! r in angstrom and radians,  c in bohr and ABCD
      implicit none
      real*8 r(1:6),c(1:3,1:4),Mass(1:4)
      real*8 cm(3,2),vec(1:3,1:6), PI
      real*8 XA,XB,YC,YD,Cxsin,Dxsin,dotp
      real*8 RX1,RX2,RX3,RX4,RX5,RX6
      real*8 AMass,BMass,CMass,DMass

      PI=dacos(-1.0d0)
      
      AMass=1.00783d0
      BMass=15.99491d0
      CMass=AMASS
      DMass=18.99834d0

      Mass(1)=AMass
      Mass(2)=BMass
      Mass(3)=CMass
      Mass(4)=DMass
!     diatomic-diatomic AB+CD Jacobi is
!     r1, r2, R, theta1, theta2, phi
!     r1, distance AB
!     r2, distance CD
!     R,  distance between the centers of mass AB(X) and CD(Y)
!     theta1, angle from vector B->A to vector X->Y
!     theta2, angle from vector D->C to vector X->Y
!     phi, dihedral angle between the plane AXY and plane XYD    
!
!   A/- r1                D\
!     \                     \ r2
!      \ theta1    R         \
!      X\-------------------->\ Y--->  
!        \                     \ th2
!       B \                    -/ C

!   put ABY on xy plane, X origin point, X->Y (R)  on the x axis
!   and also put A in the positive y direction
!             ^y axis
!        A    |
!        /\-  |
!        r1\  |
!           \ |
!            \|     R
!    ---------X--------->Y-->x axis
!              \
!               \
!                B
      c=0.0d0
      XA=Bmass/(Amass+Bmass)*r(1)
      XB=Amass/(Amass+Bmass)*r(1)
!    A
      c(1,1)= XA*dcos(r(4))
      c(2,1)= XA*dsin(r(4))
      c(3,1)= 0.0d0
!    B
      c(1,2)= XB*dcos(PI-r(4))
      c(2,2)=-XB*dsin(PI-r(4))
      c(3,2)=0.d0

!     put ABY on xy plane, X origin point
!               ^y
!          A    |
!          /\-  |            
!          r1\  |               \D
!             \ |                \ |
!              \|     R           \|
!      ---------X----------------->Y-->x
!                \                 |\
!                 \                |-\/C
!                  B
!      C in +z and -z phase is identical, so put it in +z
      YC=Dmass/(Cmass+Dmass)*r(2)
      YD=Cmass/(Cmass+Dmass)*r(2)
      Cxsin = YC*dsin(r(5))
      Dxsin =-YD*dsin(r(5))
!     C
      c(1,3)=r(3)+YC*dcos(r(5))
      c(2,3)= Cxsin*dcos(r(6))
      c(3,3)= Cxsin*dsin(r(6))
!     D
      c(1,4)=r(3)-YD*dcos(r(5))
      c(2,4)= Dxsin*dcos(r(6))
      c(3,4)= Dxsin*dsin(r(6))

      return
      end subroutine diJacobi_cart 

!******************************************************************************
!******************************************************************************
      subroutine crossp(v1,v2,v3)
      real*8 v1(1:3),v2(1:3),v3(1:3)
      ! cross product of two 3-dimensional vectors
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end subroutine crossp
      
