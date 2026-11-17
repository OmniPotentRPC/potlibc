ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      subroutine fh2o_soc_nn(xt,vg,socg,vgsoc)
      use nnparam
      implicit none
      
      real*8 xt(3,4),xt2(3,4),xn(3,4),xvec(3,6),rint(6),rt(6)
      real*8 tr1basis(nrbasissoc),tr2basis(nrbasissoc) 
      real*8 r1function,r2function,rcom,the1,the2,phi
      real*8 p2function,tp2basis(npbasissoc)
      integer n,k,l,m,j,nrcom
      real*8 vg,vgsoc,socg,socs,vpesa,vpesb,vpesc
      real*8 socrg,socrs,socpg,socps
      real*8 r1switch,r2switch,p2switch,switch

      call fh2oNN(xt,vg,vpesa,vpesb,vpesc) !xt in HHFO and angstrom
      vg=vg*23.0605d0 
      vg=min(vg,100.0d0)
      if(vg.ge.100.0d0)then
        socg=0.0d0
        vgsoc=vg+socg/349.75d0
        return
      endif
      
      xt2=xt

      call orderxnt(xt2,xn) !xn 1st H closer to F

      xvec(:,1)=xn(:,2)-xn(:,1) !  H2->H1
      xvec(:,2)=xn(:,3)-xn(:,1) !  F->H1
      xvec(:,3)=xn(:,4)-xn(:,1) !  O->H1
      xvec(:,4)=xn(:,3)-xn(:,2) !  H2->F
      xvec(:,5)=xn(:,4)-xn(:,2) !  H2->O
      xvec(:,6)=xn(:,4)-xn(:,3) !  F-O 
     
      rt(1)=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
      rt(2)=dsqrt(dot_product(xvec(:,2),xvec(:,2)))
      rt(3)=dsqrt(dot_product(xvec(:,3),xvec(:,3))) !the breaking bond
      rt(4)=dsqrt(dot_product(xvec(:,4),xvec(:,4)))
      rt(5)=dsqrt(dot_product(xvec(:,5),xvec(:,5)))
      rt(6)=dsqrt(dot_product(xvec(:,6),xvec(:,6)))

cccccccc  reactant 
       xt2(:,1)=xt(:,1)
       xt2(:,2)=xt(:,4)
       xt2(:,3)=xt(:,2)
       xt2(:,4)=xt(:,3)

       call xyzTOjcb3(xt2,rint) !-->xt2 in HOHF order,angstrom, rt in angstrom and radians
       rcom=rint(3)
       the2=rint(5)
       phi=rint(6) 

       call soc_react_ground(rcom,the2,phi,socrg)

cccccccc  product       
       xt2(:,1)=xn(:,2)
       xt2(:,2)=xn(:,4)
       xt2(:,3)=xn(:,1)
       xt2(:,4)=xn(:,3)

       call cart_diJacobi(xt2,rint) !xt2 in HOHF order,angstrom, rt in angstrom and radians
       rcom=rint(3)
       the1=rint(4)
       the2=rint(5)
       phi=rint(6) 

       call soc_product_ground(rcom,the1,the2,phi,socpg)

       switch=(1.d0-dtanh(20.0d0*(rt(3)-1.6d0)))/2.0d0

       socg=socrg*switch+socpg*(1.0d0-switch)
       
      vgsoc=vg+socg/349.75d0

      vgsoc=min(100.0d0,vgsoc)
      
      return
      end subroutine fh2o_soc_nn
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine orderxnt(x1,x2)
      implicit none
      real*8 x1(3,4),x2(3,4),xvec(3,2),r1,r2,xt(3,4)
      xt=x1
      xvec(:,1)=xt(:,1)-xt(:,3)
      xvec(:,2)=xt(:,2)-xt(:,3)
      r1=dsqrt(dot_product(xvec(:,1),xvec(:,1)))
      r2=dsqrt(dot_product(xvec(:,2),xvec(:,2)))

      if(r1.gt.r2)then
        x2(:,1)=xt(:,2)  !-->1st H closer to F
        x2(:,2)=xt(:,1)
        x2(:,3)=xt(:,3)
        x2(:,4)=xt(:,4)
      else
        x2=xt
      endif

      return
      end subroutine orderxnt

      real*8 function r1function(r,n)
      implicit none
      real*8, parameter:: r0=1.8d0,c=0.45d0
      real*8 r
      integer n

      if(n.eq.0)then
        r1function=1.0d0
        return
      else
        r1function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function r1function
     
     
      real*8 function r2function(r,n)
      implicit none
      real*8,parameter:: r0=2.0d0,c=0.42d0
      real*8 r
      integer n

      if(n.eq.0)then
        r2function=1.0d0
        return
      else
        r2function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function r2function
     
      real*8 function p2function(r,n)
      implicit none
      real*8,parameter:: r0=2.3d0,c=0.40d0
      real*8 r
      integer n

      if(n.eq.0)then
        p2function=1.0d0
        return
      else
        p2function=dtanh(dble(n)*c*(r-r0))
        return
      endif

      end function p2function
     

*********************SOC**************************
      subroutine soc_react_ground(rcom,t2,phi,socrg)
      use nnparam
      implicit none
      real*8 rcom,t2,phi,socrg,r1function
      real*8 tr1basis(nrbasissoc)
      integer j,n,k,l,nrcom

      tr1basis=0.0d0
      j=0
      do n=0,6
       do k=0,2
        do l=0,6
         
         j=j+1
         nrcom=n
         if(l.le.3)then
          tr1basis(j)=r1function(rcom,nrcom)*dcos(dble(k)*phi)*
     $                dcos(dble(l)*t2)
         else
          tr1basis(j)=r1function(rcom,nrcom)*dcos(dble(k)*phi)*
     $                 dsin(dble(l-3)*t2)
         endif
        enddo
       enddo
      enddo
      if(j.ne.nrbasissoc)then
       write(*,*)'error for nrbasissoc'
       stop
      endif
      socrg=dot_product(tr1basis(:),r1soccoef(:))

      return
      end subroutine soc_react_ground

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine soc_product_ground(rcom,t1,t2,phi,socpg)    
      use nnparam
      implicit none
      real*8 rcom,t1,t2,phi,socpg,socps,p2function
      real*8 tp2basis(npbasissoc)
      integer j,n,k,l,m,nrcom

      tp2basis=0.0d0
      j=0
      do n=0,5
       do k=0,2
        do l=0,3
         do m=0,3 
          j=j+1
          nrcom=n
          tp2basis(j)=p2function(rcom,nrcom)*dcos(dble(k)*phi)*
     $ dcos(dble(l)*t1)*dcos(dble(m)*t2)
         enddo
        enddo
       enddo
      enddo
      if(j.ne.npbasissoc)then
        write(*,*)'error for nrbasissoc'
        stop
      endif
       
      socps=dot_product(tp2basis(:),p2soccoef(:))
      socpg=-1.0d0*socps

      return
      end subroutine soc_product_ground
