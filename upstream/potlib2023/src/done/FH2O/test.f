
      program main
      use nnparam
      implicit none
      integer i,j,k,l,m,n,n1
      character*100 asymb
      real*8 v1,v2,v3,vg,socg,vgsoc,soc1

      real*8 c1(3,4),c2(3,4),jcb(6),ct(3,4)

      call pes_init
      jcb=(/1.02d0,4.9d0,3.0d0,40.0d0,0.0d0,0.0d0/)
      jcb(4:6)=jcb(4:6)*radian

      open(21,file='out-stationary-casci',status='old')
      open(23,file='stationary.pes')
      write(23,'(a5,5a16)')'No.','Eab','SOC ab','EFit 1','EFit2',
     % 'SOCfit'
      do i=1,6 
       read(21,*)n1
       read(21,*)v1,soc1
       read(21,*) asymb, c1(:,1)
       read(21,*) asymb, c1(:,2)
       read(21,*) asymb, c1(:,3)
       read(21,*) asymb, c1(:,4) ! c1 in HHFO 

       if(i.eq.5)then
          c2(:,1)=c1(:,2)
          c2(:,2)=c1(:,4)
          c2(:,3)=c1(:,1)
          c2(:,4)=c1(:,3) ! c2 in HOH...F order
          call xyzTOjcb3(c2,jcb)
          write(*,'(a2,3f15.6)')'H',c2(:,1)
          write(*,'(a2,3f15.6)')'O',c2(:,2)
          write(*,'(a2,3f15.6)')'H',c2(:,3)
          write(*,'(a2,3f15.6)')'F',c2(:,4)
          write(*,'(6f10.5)')jcb(1:3),jcb(4:6)/dacos(-1.d0)*180.0d0
!-->   0.96900   1.04766   2.03327 106.21027  33.12256  64.17181
          call jcb3tocart(jcb,ct)
          call xyzTOjcb3(c2,jcb)
          write(*,'(6f10.5)')jcb(1:3),jcb(4:6)/dacos(-1.d0)*180.0d0
       endif
       
       call fh2o_soc_nn(c1,v2,socg,vgsoc)
       v1=(v1-vpesmin)*627.5095d0

       c2(:,1)=c1(:,2)
       c2(:,2)=c1(:,1)
       c2(:,3)=c1(:,3)
       c2(:,4)=c1(:,4) ! c2 in HH..FO order
       call fh2o_soc_nn(c2,v3,socg,vgsoc)

       write(23,'(i5,6f16.6)')i,v1,soc1,v2,v3,socg
      enddo

      end program main

