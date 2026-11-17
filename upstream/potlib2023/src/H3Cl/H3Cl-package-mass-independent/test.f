        include 'pes-H3Cl.f'
        implicit real*8(a-h,o-z)
        dimension r(6),x(3,4)
        character fn*20,fn1*30
        call pes_init

        open(1,file='test.in',status='old')
        open(2,file='test.out',status='unknown')
        read(1,'(A30)')fn1 
        read(1,*)fn,x(:,1)
        read(1,*)fn,x(:,2)
        read(1,*)fn,x(:,3)
        read(1,*)fn,x(:,4)
        call pesh3cl(x,vtot)
        write(2,'(A,f15.8)')fn1,vtot
        close(1)
        close(2)
        end

 
