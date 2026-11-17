      program main
      implicit none
      real*8 c(3,5),vpes,vab
      integer ndata,i,itmp,ns
      real*8 c2(3,5),rccl,rmses,maxe
      character*8 asymb
      
      call pes_init()

      ndata=2
      open(11,file='test-h3o2.xyz',status='old')
      do i=1,ndata
       read(11,*)itmp
       read(11,*)vab
       vab=(vab+152.040123400851d0)*27.2114d0
       read(11,*)asymb,c(:,1)
       read(11,*)asymb,c(:,2)
       read(11,*)asymb,c(:,3)
       read(11,*)asymb,c(:,4)
       read(11,*)asymb,c(:,5)

       call h3o2pipNN(c,vpes)
       write(6,'(3f20.10)')vab,vpes,vpes-vab	   

      enddo

      end  
