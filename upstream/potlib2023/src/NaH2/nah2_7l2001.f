c   System:             NaH2                                                    
c   Common name:        7l                                                      
c   Number of electronic surfaces: 2                                            
c   Number of derivatives: 0                                                    
C   Interface:       potlib2001
C                                                                               
c   Note:	This version contains long range forces,                              
c		no derivatives                                                               
C                                                                               
c This is NaH2 potential matrix 7l, finalized on 5/17/99                        
c This version contains long range forces.                                      
C                                                                               
c  The potential is called by first issuing a call to PREPOT which initializes  
c  all surfaces.  The potential is called by calling POT with the parameters    
c  r, e, nt, nsurf.  r(nt,3) and ei(nt) are double precision arrays of          
c  dimension nt which is an integer representing the number of geometries       
c  input.  The integer nsurf selects the desired surface.  nsurf=1 is the       
c  lower diabatic surface, nsurf=2 is the coupling surface, nsurf=3 is          
c  the upper diabatic surface, nsurf = 4 selects the lower adiabatic            
c  surface, and nsurf = 5 selects the upper adiabatic surface                   
                                                                                
c                                                                               
c ********************************************************************          
c ********************************************************************          
C                                                                               
      subroutine prepot                                                         
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE22CM/ ret10,bet10,bet10b,f10,                               
     +                    det190,bet190,ret190,                                 
     +                    t1alp,t1rho,det1c,bet1c,ret1c,                        
     +                    det20,bet20,ret20,                                    
     +                    det290,bet290,ret290,                                 
     +                    t2alp,t2rho,                                          
     +                    des290,s2alp22,s2rho22,deh222,s20alp,s20rho,          
     +                    c2a,c2b,c2c,una2p,                                    
     +                    rnaie22,rnapol1,rnaquad,rqq,rqh,rb,                   
     +                    balp22,brho22,rd22,                                   
     +                    h2hc1,h2hc2,h2hc3,h2hc4,h2hc5,                        
     +                    h2qc1,h2qc2,h2qc3,h2qc4,h2qc5,                        
     +                    denah22,renah22,bf22,b022,gamnah22                    
C                                                                               
C      common /com_para/ g_myrank,g_nprocs                                      
C      integer g_myrank, g_nprocs                                               
C                                                                               
C      if (g_myrank .eq. 0) then                                                
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
        WRITE(NFLAG(18),*) 'PREPOT called for NaH2 potential matrix 7l.'        
        WRITE(NFLAG(18),*) 'Final version.  Last modified 17 May 1999.'         
        WRITE(NFLAG(18),*) 'Long-range forces are included.'                    
        WRITE(NFLAG(18),*) 'NFLAG(4)=0=>SELECT 11,12,22 DIA. PES VIA',          
     x                     ' NASURF'                                            
        WRITE(NFLAG(18),*) '        =1=>SELECT ADIA. PES (MINUS ROOT)'          
        WRITE(NFLAG(18),*) '            PUT IN GROUND STATE PES',               
     X                                ' (=> NASURF(1,1)=1)'                     
        WRITE(NFLAG(18),*) '        =2=>SELECT ADIA. PES (PLUS ROOT)'           
        WRITE(NFLAG(18),*) '            PUT IN EXCITED STATE PES',              
     X                                ' (=> NASURF(2,2)=1)'                     
        WRITE(NFLAG(18),123) NFLAG(4),NASURF(1,1),NASURF(1,2),                  
     x                                NASURF(2,2)                               
123     FORMAT(' INPUT NFLAG(4),NASURF11,12,22 VALUES:',4I5)                    
        WRITE(NFLAG(18),*) ' '                                                  
C                                                                               
C      endif                                                                    
C                                                                               
        call prepot11                                                           
        call prepot12                                                           
        call prepot22                                                           
C                                                                               
C                                                                               
      EZERO(1)=DEH2                                                             
      EZERO(2)=DEH222                                                           
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
C                                                                               
      INDEXES(1) = 11                                                           
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
c                                                                               
c ********************************************************************          
c ********************************************************************          
C                                                                               
      SUBROUTINE POT                                                            
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
      IF(nflag(4).eq.0) then                                                    
         if (nasurf(1,1) .gt. 0) call pot11                                     
         if (nasurf(2,2) .gt. 0) call pot22                                     
         if (nasurf(1,2)+nasurf(2,1) .gt. 0) call pot12                         
C                                                                               
      ELSEIF (nflag(4) .eq. 1) then                                             
         call pot11                                                             
         call pot12                                                             
         call pot22                                                             
            ENGYGS = 0.5d0*(ENGYGS+ENGYES(1)-                                   
     +                dsqrt((ENGYGS-ENGYES(1))**2+                              
     +                4.d0*ENGYIJ(1)**2))                                       
C                                                                               
      ELSEIF (nflag(4) .eq. 2) then                                             
         call pot11                                                             
         call pot12                                                             
         call pot22                                                             
            ENGYES(1) = 0.5d0*(ENGYGS+ENGYES(1)+                                
     +                  dsqrt((ENGYGS-ENGYES(1))**2+                            
     +                  4.d0*ENGYIJ(1)**2))                                     
      end if                                                                    
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
      end                                                                       
C                                                                               
c ********************************************************************          
c ********************************************************************          
c  U11 surface - lower diabatic surface                                         
c                                                                               
      subroutine prepot11                                                       
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
C misc constant                                                                 
C                                                                               
      rac2=1.d0/dsqrt(2.d0)                                                     
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE POT11                                                          
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
C misc constants                                                                
C                                                                               
      rac2=1.d0/dsqrt(2.d0)                                                     
C                                                                               
c======================================================================         
c useful quantities calculated first, r1^2, R2^2, R3^2, atom to diatom          
c center of mass distances, the H-Na-H bond angle and the Jacobian              
c angle chi.                                                                    
C=======================================================================        
C                                                                               
      r1s=r(1)*r(1)                                                             
      r2s=r(2)*r(2)                                                             
      r3s=r(3)*r(3)                                                             
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))                                 
      rnags=rnag*rnag                                                           
      rnagc=rnags*rnag                                                          
C                                                                               
      rmid = 0.5d0*(r(1)+r(3))                                                  
      rave = 0.5d0*(r(1)-r(3))                                                  
C                                                                               
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up             
C                                                                               
      cos1 = (r1s+r2s-r3s)/(2.0d0*r(1)*r(2))                                    
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(3)*r(2))                                    
      cosg = 0.5d0*(cos1-cos2)                                                  
C                                                                               
      frac=0.958d0                                                              
C                                                                               
C H to NaH center of masses                                                     
C                                                                               
      rnag1= dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))                  
      rnag3= dsqrt(dabs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))                  
C                                                                               
      rnag1c=rnag1*rnag1*rnag1                                                  
      rnag3c=rnag3*rnag3*rnag3                                                  
C                                                                               
C the shortest of any two distances . . .                                       
C                                                                               
      rsm1=r(3)+(r(2)-r(3))*                                                    
     +     (0.5d0+0.5d0*dtanh(0.5d0*(r(3)-r(2))))                               
C                                                                               
      rsm2=r(3)+(r(1)-r(3))*                                                    
     +     (0.5d0+0.5d0*dtanh(0.5d0*(r(3)-r(1))))                               
C                                                                               
       rsm3=r(1)+(r(2)-r(1))*                                                   
     +      (0.5d0+0.5d0*dtanh(0.5d0*(r(1)-r(2))))                              
C                                                                               
      rexp1=dexp(-r(1))                                                         
      rexp2=dexp(-r(2))                                                         
      rexp3=dexp(-r(3))                                                         
C                                                                               
C checks on angles . . .                                                        
C                                                                               
      if (rnag.eq.0.0d0) then                                                   
          cstsq=1.0d0                                                           
      else                                                                      
          cstsq=(0.5d0*(r3s-r1s)/(r(2)*rnag))**2                                
      endif                                                                     
C                                                                               
      cosa=(r1s+r3s-r2s)/(2*r(1)*r(3))                                          
C                                                                               
      coscos = cstsq*(0.5d0+0.5d0*cosa)*(1.0d0-dexp(-10.0d0*rnags))             
C                                                                               
      if (rnag1 .eq. 0.0d0) then                                                
          cstsq1 = 1.0d0                                                        
      else                                                                      
          cstsq1= 0.25d0*((r3s-r2s-r1s+2.0d0*r3s*frac)/                         
     +            (r(1)*rnag1))**2                                              
      endif                                                                     
C                                                                               
      if (rnag3 .eq. 0.0d0) then                                                
          cstsq3 = 1.0d0                                                        
      else                                                                      
          cstsq3= 0.25d0*((r1s-r2s-r3s+2.0d0*r1s*frac)/                         
     +            (r(3)*rnag3))**2                                              
      endif                                                                     
C                                                                               
C============================================                                   
C This section contains the diatomic curves                                     
C============================================                                   
C                                                                               
C The first NaH singlet                                                         
C                                                                               
      xes1 = dexp(-bes1*(r(1)-res1))                                            
      s1 = des1*onethd*xes1*(xes1+2.0d0)                                        
C                                                                               
C The H2 singlet                                                                
C                                                                               
C=========================================================                      
C new singlet, designed to minic lower adiabatic behavior                       
C=========================================================                      
C                                                                               
      beh2 = bf11*(b011-(r(2)/gamh2a11)+(r(2)/gamh2b11)**2)/                    
     +       (bf11-(r(2)/gamh2a11)+(r(2)/gamh2b11)**2)                          
C                                                                               
      sswt = 0.5d0*(1.0d0-dtanh((rmid-s2rho)/s2alp))*cosg**2                    
      deh2m = deh2 + deh2*(des20 - 1.0d0)*sswt                                  
C                                                                               
      xeh2 = dexp(-beh2*(r(2)-reh211))                                          
      s2 = deh2m*xeh2*(xeh2-2.0d0)                                              
C                                                                               
C The second NaH singlet                                                        
C                                                                               
      xes3 = dexp(-bes1*(r(3)-res1))                                            
      s3 = des1*onethd*xes3*(xes3+2.0d0)                                        
C                                                                               
C================================================================               
C================================================================               
C                                                                               
C  H2 ionization energy:                                                        
C                                                                               
      h2ie=(h2i1+h2i2*r(2))*dexp(-h2i3*r(2))+                                   
     +     r2s*(h2i4+h2i5*r(2))*dexp(-h2i6*r(2))-s2a+hie                        
C                                                                               
C  H2 polarizability:                                                           
C                                                                               
C   al1 is the parallel polarizability for h2                                   
C   al2 is the perpendicaular polarizability for h2                             
C   alt is the anisotropic? polarizability                                      
C   alb is the mean polarizability                                              
C                                                                               
      al1=dexp(h2a1*(h2a2-r(2)))*(h2a3+h2a4*r(2)+h2a5*r2s+                      
     +    h2a6*r(2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s*r(2)+                         
     +    h2a9*r2s*r2s*r2s)+2.0d0*hpol                                          
C                                                                               
      al2=dexp(h2b1*(h2b2-r(2)))*(h2b3+h2b4*r(2)+h2b5*r2s+                      
     +    h2b6*r(2)*r2s+h2b7*r2s*r2s+h2b8*r2s*r2s*r(2)+                         
     +    h2b9*r2s*r2s*r2s)+2.0d0*hpol                                          
C                                                                               
      alt=twothd*(al1-al2)                                                      
      alb=onethd*al1+twothd*al2                                                 
C                                                                               
      alhh=alb+alt*(1.5d0*cstsq-0.5d0)                                          
C                                                                               
C  The NaH + H dispersion and NaH + H induction terms must go on the lower diaba
C  so that the adiabat is fitted correctly in the region of the coupling        
C                                                                               
C NaH properties                                                                
C                                                                               
      rnahie1=(hnai1+hnai2*r(1))/                                               
     +        (hnai3+hnai4*r(1)+hnai5*r1s+hnai6*r(1)*r1s)+rnaie                 
C                                                                               
      rnahie3=(hnai1+hnai2*r(3))/                                               
     +        (hnai3+hnai4*r(3)+hnai5*r3s+hnai6*r(3)*r3s)+rnaie                 
C                                                                               
      rnahmu1=(hnau1*r(1)+hnau2*r1s)/                                           
     +        (hnau3+hnau4*r(1)+hnau5*r1s+hnau6*r(1)*r1s)                       
C                                                                               
      rnahmu3=(hnau1*r(3)+hnau2*r3s)/                                           
     +        (hnau3+hnau4*r(3)+hnau5*r3s+hnau6*r(3)*r3s)                       
C                                                                               
      rnahpar1=rnapol+hpol+dexp(-hnaa1*r(1))*                                   
     +         (hnaa2+hnaa3*r(1)+hnaa4*r1s+hnaa5*r1s*r(1))                      
C                                                                               
      rnahpar3=rnapol+hpol+dexp(-hnaa1*r(3))*                                   
     +         (hnaa2+hnaa3*r(3)+hnaa4*r3s+hnaa5*r3s*r(3))                      
C                                                                               
      rnahperp1=rnapol+hpol+dexp(-hnab1*r(1))*                                  
     +          (hnab2+hnab3*r(1)+hnab4*r1s+hnab5*r1s*r(1))                     
C                                                                               
      rnahperp3=rnapol+hpol+dexp(-hnab1*r(3))*                                  
     +          (hnab2+hnab3*r(3)+hnab4*r3s+hnab5*r3s*r(3))                     
C                                                                               
      rnahpol1=rnahperp1+(rnahpar1-rnahperp1)*cstsq1                            
      rnahpol1b=onethd*rnahpar1+twothd*rnahperp1                                
C                                                                               
      rnahpol3=rnahperp3+(rnahpar3-rnahperp3)*cstsq3                            
      rnahpol3b=onethd*rnahpar3+twothd*rnahperp3                                
C                                                                               
C induction forces                                                              
C dipole induced dipole                                                         
C                                                                               
      div1=1.0d0/(rnag1c*rnag1c*rnag1c+rd)                                      
      rmuimu1=-rnahmu1**2*hpol*(1.5d0*cstsq1+0.5d0)*rnag1c*div1                 
C                                                                               
      div3=1.0d0/(rnag3c*rnag3c*rnag3c+rd)                                      
      rmuimu3=-rnahmu3**2*hpol*(1.5d0*cstsq3+0.5d0)*rnag3c*div3                 
C                                                                               
C dispersion forces                                                             
C                                                                               
      c690=-1.5d0*rnaie*h2ie*rnapol*al2/                                        
     +     ((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb90))                              
C                                                                               
      c690=c690*f90*rnagc                                                       
C                                                                               
      c60=-1.5d0*rnaie*h2ie*rnapol*al1/                                         
     +    ((rnaie+h2ie)*(rnagc*rnagc*rnagc+rb0))                                
C                                                                               
      c60=c60*f0*rnagc                                                          
C                                                                               
      c6=c690+(c60-c690)*cstsq*(1.0d0-dexp(-39.0625d0*rnag*rnagc))              
C                                                                               
      c61=-1.5d0*rnahie1*hie*rnag1c/(rnahie1+hie)*rnahpol1*hpol*div1            
C                                                                               
      dexp1=dexp(-39.0625d0*rnag1c*rnag1)                                       
C                                                                               
      c61=c61*(1.0d0+(rnahpol1b/rnahpol1-1.0d0)*dexp1)                          
C                                                                               
      c63=-1.5d0*rnahie3*hie*rnag3c/(rnahie3+hie)*rnahpol3*hpol*div3            
C                                                                               
      dexp3=dexp(-39.0625d0*rnag3c*rnag3)                                       
C                                                                               
      c63=c63*(1.0d0+(rnahpol3b/rnahpol3-1.0d0)*dexp3)                          
C                                                                               
      cind1=rmuimu1*(1.0d0+(1.0d0/(1.5d0*cstsq1+0.5d0)-1.0d0)*dexp1)            
C                                                                               
      cind3=rmuimu3*(1.0d0+(1.0d0/(1.5d0*cstsq3+0.5d0)-1.0d0)*dexp3)            
C                                                                               
      c6=c6*cevau*0.5d0*(1.0d0-dtanh(balp*(r(2)-brho)))                         
C                                                                               
      c6swt1=0.5d0*(1.0d0-dtanh(dalp*(r(1)-drho)))                              
      c6swt3=0.5d0*(1.0d0-dtanh(dalp*(r(3)-drho)))                              
C                                                                               
      c61=c61*c6swt1*cevau                                                      
      c63=c63*c6swt3*cevau                                                      
      cind1=cind1*c6swt1                                                        
      cind3=cind3*c6swt3                                                        
C                                                                               
      elec=c61+c63+cind1+cind3+c6                                               
C                                                                               
                                                                                
C===============================================                                
C here is where it all gets put together . . .                                  
C===============================================                                
C                                                                               
      E = s1 + s2 + s3 + EZERO(1)                                               
C                                                                               
      elec = elec/cevau                                                         
C                                                                               
      E = E+elec                                                                
C                                                                               
C Here we will try a post-LEPS correction term                                  
C                                                                               
      xect = dexp(-bect*(rmid-rect))                                            
      vmid = dect*xect*(xect-2.0d0)*                                            
     +       dexp(-bect2*rave**2)*(0.5d0*(1.0d0-cosa))**2                       
C                                                                               
      ENGYGS=(E+vmid)*cevau                                                     
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C ********************************************************************          
C ********************************************************************          
C                                                                               
C U12 coupling surface . . .                                                    
C                                                                               
      subroutine prepot12                                                       
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE12CM/ ap,alp,rd0,bet,r20                                    
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT12                                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      logical leven                                                             
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE12CM/ ap,alp,rd0,bet,r20                                    
C                                                                               
C set leven true for an even U12 surface, and false for an odd                  
C surface                                                                       
                                                                                
      leven=.true.                                                              
C                                                                               
      theta = 15.6d0*dacos(-1.0d0)/180.0d0                                      
C                                                                               
      r1s=r(1)**2                                                               
      r2s=r(2)**2                                                               
      r3s=r(3)**2                                                               
C                                                                               
C Here we do the function evaluations . . .                                     
C                                                                               
      rnagsq=0.5d0*r1s+0.5d0*r3s-0.25d0*r2s                                     
      if (rnagsq .gt. 0.0d0) then                                               
          rnag=dsqrt(rnagsq)                                                    
      else                                                                      
          rnag=0.0d0                                                            
      endif                                                                     
C                                                                               
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up             
C                                                                               
      cos1 = (r1s+r2s-r3s)/(2.0d0*r(1)*r(2))                                    
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(3)*r(2))                                    
      cosg = 0.5d0*(cos1-cos2)                                                  
C                                                                               
      gc  = dcos(theta)*rnag + dsin(theta)*r(2)                                 
      gc2 = dsin(theta)*rnag - dcos(theta)*r(2)                                 
      E = ap*cosg*dexp(-(1.0d0/alp)**4*(gc - rd0)**4)*                          
     +    dexp(-(1.0d0/bet)**2*(gc2 - r20)**2)                                  
C                                                                               
      ENGYIJ(1) = E*cevau                                                       
C                                                                               
C leven equal to true for an even surface, false for an odd one                 
C                                                                               
      if (leven) then                                                           
          ENGYIJ(1)=dsqrt( ENGYIJ(1)**4/                                        
     +                     (ENGYIJ(1)**2+1.0d-8) )                              
      endif                                                                     
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C ********************************************************************          
C ********************************************************************          
C                                                                               
C The upper diabatic surface  U22                                               
C                                                                               
      subroutine prepot22                                                       
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE12CM/ ap,alp,rd0,bet,r20                                    
C                                                                               
         COMMON /PRE22CM/ ret10,bet10,bet10b,f10,                               
     +                    det190,bet190,ret190,                                 
     +                    t1alp,t1rho,det1c,bet1c,ret1c,                        
     +                    det20,bet20,ret20,                                    
     +                    det290,bet290,ret290,                                 
     +                    t2alp,t2rho,                                          
     +                    des290,s2alp22,s2rho22,deh222,s20alp,s20rho,          
     +                    c2a,c2b,c2c,una2p,                                    
     +                    rnaie22,rnapol1,rnaquad,rqq,rqh,rb,                   
     +                    balp22,brho22,rd22,                                   
     +                    h2hc1,h2hc2,h2hc3,h2hc4,h2hc5,                        
     +                    h2qc1,h2qc2,h2qc3,h2qc4,h2qc5,                        
     +                    denah22,renah22,bf22,b022,gamnah22                    
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT22                                                          
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE12CM/ ap,alp,rd0,bet,r20                                    
C                                                                               
         COMMON /PRE22CM/ ret10,bet10,bet10b,f10,                               
     +                    det190,bet190,ret190,                                 
     +                    t1alp,t1rho,det1c,bet1c,ret1c,                        
     +                    det20,bet20,ret20,                                    
     +                    det290,bet290,ret290,                                 
     +                    t2alp,t2rho,                                          
     +                    des290,s2alp22,s2rho22,deh222,s20alp,s20rho,          
     +                    c2a,c2b,c2c,una2p,                                    
     +                    rnaie22,rnapol1,rnaquad,rqq,rqh,rb,                   
     +                    balp22,brho22,rd22,                                   
     +                    h2hc1,h2hc2,h2hc3,h2hc4,h2hc5,                        
     +                    h2qc1,h2qc2,h2qc3,h2qc4,h2qc5,                        
     +                    denah22,renah22,bf22,b022,gamnah22                    
C                                                                               
C                                                                               
C================================================================               
C This section calcluates useful quantities - the A to BC center                
C of mass distance, the angles between A-BC and BC, the shorest                 
C of any of the three diatomic distances, and the bond angle                    
C H-Na-H                                                                        
C================================================================               
C                                                                               
      rmid = 0.5d0*(r(1)+r(3))                                                  
      rave = 0.5d0*(r(1)-r(3))                                                  
C                                                                               
      r1s=r(1)*r(1)                                                             
      r2s=r(2)*r(2)                                                             
      r3s=r(3)*r(3)                                                             
C                                                                               
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))                                 
      rnags=rnag*rnag                                                           
      rnagc=rnags*rnag                                                          
C                                                                               
      if (rnag .eq. 0.0d0) then                                                 
          cstsq=0.0d0                                                           
      else                                                                      
          cstsq=(0.5d0*(r3s-r1s)/(r(2)*rnag))**2                                
      endif                                                                     
C                                                                               
      cosa=(r1s+r3s-r2s)/(2*r(1)*r(3))                                          
C                                                                               
      coscos = cstsq*(0.5d0+0.5d0*cosa)*(1.0d0-dexp(-10.0d0*rnags))             
C                                                                               
C cosg = 0.5*(cos1-cos2) looks like cos(chi), but it never blows up             
C                                                                               
      cos1 = (r1s+r2s-r3s)/(2.0d0*r(1)*r(2))                                    
      cos2 = (r3s+r2s-r1s)/(2.0d0*r(3)*r(2))                                    
      cosg = 0.5d0*(cos1-cos2)                                                  
C                                                                               
      frac=0.958d0                                                              
      rnag1= dsqrt(dabs(r2s+frac*(r3s-r2s-r1s)+r1s*frac*frac))                  
      rnag3= dsqrt(dabs(r2s+frac*(r1s-r2s-r3s)+r3s*frac*frac))                  
C                                                                               
      rnag1c=rnag1*rnag1*rnag1                                                  
      rnag3c=rnag3*rnag3*rnag3                                                  
C                                                                               
      if (rnag1 .eq. 0.0d0) then                                                
          cstsq1 = 1.0d0                                                        
      else                                                                      
         cstsq1= 0.25d0*((r3s-r2s-r1s+2.0d0*r3s*frac)/                          
     +           (r(1)*rnag1))**2                                               
      endif                                                                     
C                                                                               
      if (rnag3 .eq. 0.0d0) then                                                
          cstsq3 = 1.0d0                                                        
      else                                                                      
          cstsq3= 0.25d0*((r1s-r2s-r3s+2.0d0*r1s*frac)/                         
     +            (r(3)*rnag3))**2                                              
      endif                                                                     
C                                                                               
C=========================================================                      
C The first NaH singlet                                                         
C=========================================================                      
C                                                                               
C new singlet, designed to minic lower adiabatic behavior                       
C                                                                               
      beta = bf22*(b022+(r(1)/gamnah22)**8)/                                    
     +       (bf22+(r(1)/gamnah22)**8)                                          
      xenah = dexp(-beta*(r(1)-renah22))                                        
C                                                                               
      s1 = denah22*xenah*(xenah-2.0d0)                                          
C                                                                               
C The first NaH triplet                                                         
C                                                                               
      g10 = r(1)-ret10                                                          
      t10 = f10*dexp(-bet10*g10)+dexp(-bet10b*g10)                              
      t10 = t10/(1.0d0+f10)                                                     
C                                                                               
      xet190 = dexp(-bet190*(r(1)-ret190))                                      
      t190 = det190*onethd*xet190*(xet190 + 2.0d0)                              
C                                                                               
      xet1c = dexp(-bet1c*(r(1)-ret1c))                                         
      t1c = det1c*onethd*xet1c*(xet1c + 2.0d0)                                  
C                                                                               
      t1m = t190 +(t1c-t190)*                                                   
     +      0.5d0*(1.0d0+dtanh(t1alp*(r(2)-t1rho)))                             
C                                                                               
      t1 = t1m + (t10-t1m)*cosg**2                                              
C                                                                               
C The H2 triplet                                                                
C                                                                               
      xet20 = dexp(-bet20*(r(2)-ret20))                                         
      t20 = det20*onethd*xet20*(xet20 + 2.0d0)                                  
C                                                                               
      xet290 = dexp(-bet290*(r(2)-ret290))                                      
      t290 = det290*onethd*xet290*(xet290 + 2.0d0)                              
C                                                                               
      t2 = t20 + (t290-t20)*(1.0d0-cosg**2)*                                    
     +     0.5d0*(1.0d0-dtanh(t2alp*(rmid-t2rho)))                              
C                                                                               
C The H2 singlet                                                                
C new singlet, designed to minic lower adiabatic behavior                       
C The parameters are the same as those used in prepot11/pot11                   
C                                                                               
      beh2 = bf11*(b011-(r(2)/gamh2a11)+(r(2)/gamh2b11)**2)/                    
     +       (bf11-(r(2)/gamh2a11)+(r(2)/gamh2b11)**2)                          
C                                                                               
      sswt = 0.5d0*(1.0d0-dtanh(s2alp22*(rmid-s2rho22)))                        
      deh2m = deh222 + deh222*(des290 - 1.0d0)*sswt                             
C                                                                               
      xeh2 = dexp(-beh2*(r(2)-reh2))                                            
      s2a = deh2m*xeh2*(xeh2-2.0d0) + una2p                                     
C                                                                               
      s2=0.5d0*(s2a+t2-dtanh((s2a-t2)/0.1d0)*(s2a-t2))                          
C                                                                               
      xeh2 = dexp(-1.34d0*beh2*(r(2)-reh2))                                     
      s20b = (deh222-una2p)*xeh2*(xeh2-2.0d0)                                   
C                                                                               
      s2 = s2 + (s20b - s2)*cosg**2*                                            
     +     0.5d0*(1.0d0-dtanh(s20alp*(rmid-s20rho)))                            
C                                                                               
C The second NaH singlet                                                        
C new singlet, designed to minic lower adiabatic behavior                       
C                                                                               
      beta = bf22*(b022+(r(3)/gamnah22)**8)/                                    
     +       (bf22+(r(3)/gamnah22)**8)                                          
      xenah = dexp(-beta*(r(3)-renah22))                                        
C                                                                               
      s3 = denah22*xenah*(xenah-2.0d0)                                          
C                                                                               
C The second NaH triplet                                                        
C                                                                               
      xet3 = dexp(-bet1*(r(3)-ret1))                                            
      t30 = det1*onethd*xet3*(xet3 + 2.0d0)                                     
C                                                                               
      g30 = r(3)-ret10                                                          
      t30 = f10*dexp(-bet10*g30)+dexp(-bet10b*g30)                              
      t30 = t30/(1.0d0+f10)                                                     
C                                                                               
      xet390 = dexp(-bet190*(r(3)-ret190))                                      
      t390 = det190*onethd*xet390*(xet390 + 2.0d0)                              
C                                                                               
      xet3c = dexp(-bet1c*(r(3)-ret1c))                                         
      t3c = det1c*onethd*xet3c*(xet3c + 2.0d0)                                  
C                                                                               
      t3m = t390 + (t3c - t390)*                                                
     +      0.5d0*(1.0d0+dtanh(t1alp*(r(2)-t1rho)))                             
C                                                                               
      t3 = t3m + (t30 - t3m)*cosg**2                                            
C                                                                               
C===============================================================                
C this long section contains data pertinent to long range forces                
C===============================================================                
C                                                                               
C  H2 ionization energy:                                                        
C                                                                               
      h2ie=(h2i1+h2i2*r(2))*dexp(-h2i3*r(2))+                                   
     +     r2s*(h2i4+h2i5*r(2))*dexp(-h2i6*r(2))-s2a+hie                        
C                                                                               
C  H2 polarizability:                                                           
C                                                                               
C   al1 is the parallel polarizability for h2                                   
C   al2 is the perpendicaular polarizability for h2                             
C   alt is the anisotropic polarizability                                       
C   alb is the mean polarizability                                              
C                                                                               
      al1=dexp(h2a1*(h2a2-r(2)))*(h2a3+h2a4*r(2)+h2a5*r2s+                      
     +    h2a6*r(2)*r2s+h2a7*r2s*r2s+h2a8*r2s*r2s*r(2)+                         
     +    h2a9*r2s*r2s*r2s)+2.0d0*hpol                                          
C                                                                               
      al2=dexp(h2b1*(h2b2-r(2)))*(h2b3+h2b4*r(2)+h2b5*r2s+                      
     +    h2b6*r(2)**3+h2b7*r2s*r2s+h2b8*r2s*r2s*r(2)+                          
     +    h2b9*r2s*r2s*r2s)+2.0d0*hpol                                          
C                                                                               
      alt=twothd*(al1-al2)                                                      
      alb=onethd*al1+twothd*al2                                                 
      alhh=alb+alt*(1.5d0*cstsq-0.5d0)                                          
C                                                                               
C Here follow properties necesary for the long range forces in                  
C the Na + H2 channel:                                                          
C                                                                               
C H2 properties:                                                                
C                                                                               
      h2hex=dexp(-h2hc1*r(2))*                                                  
     +      (h2hc2*r(2)+h2hc3*r2s+h2hc4*r(2)*r2s+h2hc5*r2s*r2s)                 
C                                                                               
      h2quad=dexp(-h2qc1*r(2))*                                                 
     +       (h2qc2*r(2)+h2qc3*r2s+h2qc4*r(2)*r2s+h2qc5*r2s*r2s)                
C                                                                               
C NaH properties                                                                
C                                                                               
      rnahie1=(hnai1+hnai2*r(1))/                                               
     +        (hnai3+hnai4*r(1)+hnai5*r1s+hnai6*r(1)*r1s)+                      
     +        rnaie22+una2p                                                     
C                                                                               
      rnahie3=(hnai1+hnai2*r(3))/                                               
     +        (hnai3+hnai4*r(3)+hnai5*r3s+hnai6*r(3)*r3s)+                      
     +        rnaie22+una2p                                                     
C                                                                               
      rnahmu1=(hnau1*r(1)+hnau2*r1s)/                                           
     +        (hnau3+hnau4*r(1)+hnau5*r1s+hnau6*r(1)*r1s)                       
C                                                                               
      rnahmu3=(hnau1*r(3)+hnau2*r3s)/                                           
     +        (hnau3+hnau4*r(3)+hnau5*r3s+hnau6*r(3)*r3s)                       
C                                                                               
      rnahpar1=159.267d0+hpol+dexp(-hnaa1*r(1))*                                
     +         (hnaa2+hnaa3*r(1)+hnaa4*r1s+hnaa5*r1s*r(1))                      
C                                                                               
      rnahpar3=159.267d0+hpol+dexp(-hnaa1*r(3))*                                
     +         (hnaa2+hnaa3*r(3)+hnaa4*r3s+hnaa5*r3s*r(3))                      
C                                                                               
      rnahperp1=159.267d0+hpol+dexp(-hnab1*r(1))*                               
     +          (hnab2+hnab3*r(1)+hnab4*r1s+hnab5*r1s*r(1))                     
C                                                                               
      rnahperp3=159.267d0+hpol+dexp(-hnab1*r(3))*                               
     +          (hnab2+hnab3*r(3)+hnab4*r3s+hnab5*r3s*r(3))                     
C                                                                               
      rnahpol1=rnahperp1+(rnahpar1-rnahperp1)*cstsq1                            
      rnahpol1b=onethd*rnahpar1+twothd*rnahperp1                                
C                                                                               
      rnahpol3=rnahperp3+(rnahpar3-rnahperp3)*cstsq3                            
      rnahpol3b=onethd*rnahpar3+twothd*rnahperp3                                
C                                                                               
C electrostatic forces                                                          
C                                                                               
C quadrupole quadrupole                                                         
C                                                                               
      qqint=3.0d0-7.0d0*cstsq                                                   
      qq=0.75d0*rnaquad*h2quad/(rnagc*rnag*rnag+rqq)*qqint                      
      dexp2=dexp(-39.0625d0*rnagc*rnag)                                         
      qq=qq*(1.0d0+(-4.0d0*pi/(3.0d0*qqint)-1.0d0)*dexp2)                       
C                                                                               
C quadrupole hexadecapole                                                       
C                                                                               
      qhint=-45.0d0-147.0d0*cstsq*cstsq+162.0d0*cstsq                           
      qh=0.3125d0*rnaquad*h2hex/(rnagc*rnagc*rnag+rqh)*(qhint)                  
      qh=qh*(1.0d0+(-40.8d0*pi/qhint-1.0d0)*dexp2)                              
C                                                                               
C induction forces                                                              
C dipole induced dipole                                                         
C                                                                               
      div1=1.0d0/(rnag1c*rnag1c*rnag1c+rd22)                                    
      rmuimu1=-rnahmu1**2*hpol*(1.5d0*cstsq1+0.5d0)*rnag1c*div1                 
C                                                                               
      div3=1.0d0/(rnag3c*rnag3c*rnag3c+rd22)                                    
      rmuimu3=-rnahmu3**2*hpol*(1.5d0*cstsq3+0.5d0)*rnag3c*div3                 
C                                                                               
C dispersion forces                                                             
C                                                                               
      c6=-1.5d0*rnaie22*h2ie/(rnaie22+h2ie)*rnapol1*alhh/                       
     +   (rnagc*rnagc+rb)                                                       
      c6=c6*(1.0d0+(alb/alhh-1.0d0)*dexp2)                                      
C                                                                               
      c61=-1.5d0*rnahie1*hie/(rnahie1+hie)*rnahpol1*hpol*rnag1c*div1            
      dexp1=dexp(-39.0625d0*rnag1c*rnag1)                                       
      c61=c61*(1.0d0+(rnahpol1b/rnahpol1-1.0d0)*dexp1)                          
C                                                                               
      c63=-1.5d0*rnahie3*hie/(rnahie3+hie)*rnahpol3*hpol*rnag3c*div3            
      dexp3=dexp(-39.0625d0*rnag3c*rnag3)                                       
      c63=c63*(1.0d0+(rnahpol3b/rnahpol3-1.0d0)*dexp3)                          
C                                                                               
      cind1=rmuimu1*(1.0d0+(1.0d0/(1.5d0*cstsq1+0.5d0)-1.0d0)*dexp1)            
C                                                                               
      cind3=rmuimu3*(1.0d0+(1.0d0/(1.5d0*cstsq3+0.5d0)-1.0d0)*dexp3)            
C                                                                               
      elec=(qq+qh+c6*cevau)*0.5d0*                                              
     +     (1.0d0-dtanh(balp22*(r(2)-brho22)))                                  
C                                                                               
      c6swt1=0.5d0*(1.0d0-dtanh(dalp*(r(1)-drho)))                              
      c6swt3=0.5d0*(1.0d0-dtanh(dalp*(r(3)-drho)))                              
C                                                                               
      c61=c61*c6swt1*cevau                                                      
      c63=c63*c6swt3*cevau                                                      
      cind1=cind1*c6swt1                                                        
      cind3=cind3*c6swt3                                                        
C                                                                               
      elec=elec+c61+c63+cind1+cind3                                             
C                                                                               
C===============================================================                
C                                                                               
      coul1=0.5d0*(s1+t1)                                                       
      coul2=0.5d0*(s2+t2)                                                       
      coul3=0.5d0*(s3+t3)                                                       
C                                                                               
      exch1=0.5d0*(s1-t1)                                                       
      exch2=0.5d0*(s2-t2)                                                       
      exch3=0.5d0*(s3-t3)                                                       
C                                                                               
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2                      
      cplg2=c2a*dexp(-c2b*w-c2c*(r(1)+r(2)+r(3)))                               
      E=coul1+coul2+coul3+ EZERO(2)  -                                          
     +                dsqrt(w+cplg2*cplg2)/                                     
     +                dsqrt(2.0d0)                                              
C                                                                               
      elec = elec/cevau                                                         
C                                                                               
      ENGYES(1)=(E+elec)*cevau                                                  
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C***************************************************************************    
C**************************************************************************     
C                                                                               
      BLOCK DATA PTPACM                                                         
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
         PARAMETER (onethd=1.0d0/3.0d0)                                         
         PARAMETER (twothd=2.0d0/3.0d0)                                         
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
         COMMON /PRE11CM/ des1,bes1,res1,                                       
     +                    reh2,bf,b0,gamh2a,gamh2b,deh2,                        
     +                    des20,s2alp,s2rho,                                    
     +                    dect,bect,rect,bect2,                                 
     +                    h2a1,h2a2,h2a3,h2a4,h2a5,h2a6,h2a7,h2a8,h2a9,         
     +                    h2b1,h2b2,h2b3,h2b4,h2b5,h2b6,h2b7,h2b8,h2b9,         
     +                    h2i1,h2i2,h2i3,h2i4,h2i5,h2i6,                        
     +                    rnaie,rnapol,hpol,hie,                                
     +                    hnai1,hnai2,hnai3,hnai4,hnai5,hnai6,                  
     +                    hnau1,hnau2,hnau3,hnau4,hnau5,hnau6,                  
     +                    hnaa1,hnaa2,hnaa3,hnaa4,hnaa5,                        
     +                    hnab1,hnab2,hnab3,hnab4,hnab5,                        
     +                    balp,brho,rb0,rb90,f90,f0,rd,dalp,drho,               
     +                    reh211,bf11,b011,gamh2a11,gamh2b11                    
C                                                                               
         COMMON /PRE12CM/ ap,alp,rd0,bet,r20                                    
C                                                                               
         COMMON /PRE22CM/ ret10,bet10,bet10b,f10,                               
     +                    det190,bet190,ret190,                                 
     +                    t1alp,t1rho,det1c,bet1c,ret1c,                        
     +                    det20,bet20,ret20,                                    
     +                    det290,bet290,ret290,                                 
     +                    t2alp,t2rho,                                          
     +                    des290,s2alp22,s2rho22,deh222,s20alp,s20rho,          
     +                    c2a,c2b,c2c,una2p,                                    
     +                    rnaie22,rnapol1,rnaquad,rqq,rqh,rb,                   
     +                    balp22,brho22,rd22,                                   
     +                    h2hc1,h2hc2,h2hc3,h2hc4,h2hc5,                        
     +                    h2qc1,h2qc2,h2qc3,h2qc4,h2qc5,                        
     +                    denah22,renah22,bf22,b022,gamnah22                    
C                                                                               
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                               
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,1,0/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
C NaH singlet                                                                   
C                                                                               
      DATA des1   / 0.55118d0/                                                  
      DATA bes1   / 1.05118d0/                                                  
      DATA res1   / 3.38583d0/                                                  
C                                                                               
C H2 singlet                                                                    
C                                                                               
      DATA reh2   / 1.401d0/                                                    
      DATA bf     / 1.627d0/                                                    
      DATA b0     / 1.194d0/                                                    
      DATA gamh2a / 2.323d0/                                                    
      DATA gamh2b / 3.126d0/                                                    
      DATA deh2   / 4.74806d0/                                                  
C                                                                               
      DATA reh211    /  1.401d0/                                                
      DATA bf11      /  1.627d0/                                                
      DATA b011      /  1.194d0/                                                
      DATA gamh2a11  /  2.323d0/                                                
      DATA gamh2b11  /  3.126d0/                                                
C                                                                               
      DATA renah22   /  3.566044d0/                                             
      DATA bf22      /  0.864d0/                                                
      DATA b022      /  0.594d0/                                                
      DATA gamnah22  /  7.376d0/                                                
      DATA denah22   /  1.97134878d0/                                           
C                                                                               
      DATA des20  / 0.996d0/                                                    
      DATA s2alp  / 3.0d0/                                                      
      DATA s2rho  / 7.0d0/                                                      
C                                                                               
C H-Na-H term                                                                   
C                                                                               
      DATA dect   / 2.24409d0/                                                  
      DATA bect   / 0.57874d0/                                                  
      DATA rect   / 3.46457d0/                                                  
      DATA bect2  / 0.87143d0/                                                  
C                                                                               
C data blocks for the long range forces                                         
C                                                                               
      DATA h2a1   /   2.38029900d0/                                             
      DATA h2a2   /   1.47400334d0/                                             
      DATA h2a3   /  -0.226331771d0/                                            
      DATA h2a4   /  -2.17339562d0/                                             
      DATA h2a5   /   9.54597998d0/                                             
      DATA h2a6   / -23.3366749d0/                                              
      DATA h2a7   /  24.7115748d0/                                              
      DATA h2a8   / -12.8786102d0/                                              
      DATA h2a9   /   2.75711758d0/                                             
C                                                                               
      DATA h2b1   /   1.35164094d0/                                             
      DATA h2b2   /  -0.254404564d0/                                            
      DATA h2b3   / -10.6627055d0/                                              
      DATA h2b4   / -14.0112406d0/                                              
      DATA h2b5   /  -7.80660398d0/                                             
      DATA h2b6   /   5.47697753d0/                                             
      DATA h2b7   /  -6.95875664d0/                                             
      DATA h2b8   /   3.64334771d0/                                             
      DATA h2b9   /  -0.492505596d0/                                            
C                                                                               
      DATA h2i1   / 160.484d0/                                                  
      DATA h2i2   / -20.397d0/                                                  
      DATA h2i3   /   3.913d0/                                                  
      DATA h2i4   /   0.777d0/                                                  
      DATA h2i5   /  -7.583/                                                    
      DATA h2i6   /   1.501/                                                    
C                                                                               
      DATA rnaie  /   5.139d0/                                                  
      DATA rnapol / 159.267d0/                                                  
      DATA hpol   /   4.5d0/                                                    
      DATA hie    /  13.598d0/                                                  
C                                                                               
      DATA hnai1  /-193.341d0/                                                  
      DATA hnai2  /  16.968d0/                                                  
      DATA hnai3  / -77.120d0/                                                  
      DATA hnai4  /  -1.024d0/                                                  
      DATA hnai5  /   6.009d0/                                                  
      DATA hnai6  /  -1.241d0/                                                  
C                                                                               
      DATA hnau1  / 119.880d0/                                                  
      DATA hnau2  /  -0.620d0/                                                  
      DATA hnau3  /  31.335d0/                                                  
      DATA hnau4  /  70.414d0/                                                  
      DATA hnau5  / -13.016d0/                                                  
      DATA hnau6  /   0.8681d0/                                                 
C                                                                               
      DATA hnaa1  /   0.67d0/                                                   
      DATA hnaa2  / -92.229d0/                                                  
      DATA hnaa3  /  51.01019d0/                                                
      DATA hnaa4  /-130.37826d0/                                                
      DATA hnaa5  /   6.40285d0/                                                
C                                                                               
      DATA hnab1  /   0.34333d0/                                                
      DATA hnab2  / -92.229d0/                                                  
      DATA hnab3  / -48.54867d0/                                                
      DATA hnab4  / -17.30523d0/                                                
      DATA hnab5  /   0.99821d0/                                                
C                                                                               
      DATA balp   /   2.07234d0/                                                
      DATA brho   /   2.0d0/                                                    
      DATA rb0    /   5.1998697814229d7/                                        
      DATA rb90   /   5.0d8/                                                    
      DATA f90    /   1.0d0/                                                    
      DATA f0     /   1.0d0/                                                    
      DATA rd     /   9.1538750508008d8/                                        
      DATA dalp   /   1.43307d0/                                                
      DATA drho   /   5.71654d0/                                                
C                                                                               
      DATA ap     /   0.835d0/                                                  
      DATA alp    /   4.0d0/                                                    
      DATA rd0    /   4.0d0/                                                    
      DATA bet    /   0.5d0/                                                    
      DATA r20    /  -1.3d0/                                                    
C                                                                               
C NaH triplets                                                                  
C                                                                               
      DATA ret10  /   3.22835d0/                                                
      DATA bet10  /   0.57874d0/                                                
      DATA bet10b /   2.44094d0/                                                
      DATA f10    /   9.44882d0/                                                
      DATA det190 /   4.32441d0/                                                
      DATA bet190 /   2.43701d0/                                                
      DATA ret190 /   1.95480d0/                                                
      DATA t1alp  /   0.30778d0/                                                
      DATA t1rho  /   3.21850d0/                                                
      DATA det1c  /   4.92126d0/                                                
      DATA bet1c  /   0.65748d0/                                                
      DATA ret1c  /   2.71654d0/                                                
C                                                                               
C H2 triplets                                                                   
C                                                                               
      DATA det20  /   3.29921d0/                                                
      DATA bet20  /   2.64567d0/                                                
      DATA ret20  /   2.44882d0/                                                
      DATA det290 /   3.97638d0/                                                
      DATA bet290 /   3.73937d0/                                                
      DATA ret290 /   0.96024d0/                                                
      DATA t2alp  /   1.0d0/                                                    
      DATA t2rho  /   9.0d0/                                                    
C                                                                               
C H2 singlets                                                                   
C                                                                               
      DATA des290   /   0.91433d0/                                              
      DATA s2alp22  /   0.66778d0/                                              
      DATA s2rho22  /   5.30472d0/                                              
      DATA deh222   /   4.74806d0/                                              
      DATA s20alp   /   0.5d0/                                                  
      DATA s20rho   /   6.0d0/                                                  
C                                                                               
C misc constants                                                                
C                                                                               
      DATA c2a    /   1.38661d0/                                                
      DATA c2b    /   0.27362d0/                                                
      DATA c2c    /   0.03764d0/                                                
      DATA una2p  /   2.1037d0/                                                 
C                                                                               
C data blocks for the long range forces                                         
C                                                                               
      DATA rnaie22  /   3.0353/                                                 
      DATA rnapol1  / 290.4d0/                                                  
      DATA rnaquad  / -16.8992d0/                                               
      DATA rqq      /   9.766021100d4/                                          
      DATA rqh      /   2.355466062927d6/                                       
      DATA rb       /   9.3119383371714d5/                                      
      DATA balp22   /   1.93701d0/                                              
      DATA brho22   /   1.4d0/                                                  
      DATA rd22     /   9.1351724748364d8/                                      
C                                                                               
      DATA h2hc1    /   0.677d0/                                                
      DATA h2hc2    /   0.378d0/                                                
      DATA h2hc3    /   0.281d0/                                                
      DATA h2hc4    /  -1.022d0/                                                
      DATA h2hc5    /   0.665d0/                                                
C                                                                               
      DATA h2qc1    /   1.532d0/                                                
      DATA h2qc2    /  -0.04635d0/                                              
      DATA h2qc3    /   1.826d0/                                                
      DATA h2qc4    /  -2.508d0/                                                
      DATA h2qc5    /   1.856d0/                                                
C                                                                               
      END                                                                       
                                                                                
                                                                                
