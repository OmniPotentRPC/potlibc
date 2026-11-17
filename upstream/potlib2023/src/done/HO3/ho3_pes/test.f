        include "ho3_pes.f"

        program main
        implicit none
        integer :: n 
        character*2 :: atom(4)
        character*9 :: geom
        double precision :: q(3,4),vpes,vpes_kcal,v_r

        call pes_init

        n=0 
        open(51,file='test.xyz',status='old')
        open(61,file='test.dat',status='replace')

        write(61,'(9H=========,20H ===================,20H =============
     $======)')
        write(61,'(A9,2A20)')'Geom','Vpes/eV','Vpes_rel/kcal*mol-1'
        write(61,'(9H=========,20H ===================,20H =============
     $======)')

        do while(.true.)
           read(51,*,end=100)
           read(51,*)geom
           read(51,*)atom(1),q(:,1) !--> H
           read(51,*)atom(2),q(:,2) !--> O
           read(51,*)atom(3),q(:,3) !--> O
           read(51,*)atom(4),q(:,4) !--> O
           call ho3_pes_interface(q,vpes)
           n=n+1 
           if (n.eq.1) then
              v_r=vpes
           endif
           if (n.ne.1) then
              vpes_kcal=(vpes-v_r)*23.0605d0 ! eV to kcal/mol
              write(61,'(A9,f20.8,f20.2)')geom,vpes,vpes_kcal
           endif
        enddo
100     continue
        write(61,'(9H=========,20H ===================,20H =============
     $======)')
        close(51)
        close(61)

        end program
