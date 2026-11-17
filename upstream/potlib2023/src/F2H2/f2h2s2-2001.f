C                                                                               
      subroutine prepot                                                         
C                                                                               
C   System:          F2H2                                                       
C   Common name:     S2                                                         
C   Number of derivatives: 0                                                    
C   Number of bodies: 4                                                         
C   Number of electronic surfaces: 1                                            
C                                                                               
C   Interface:       potlib2001
C                                                                               
C   Note:            This surface is based on the SQSBDE surface                
C                    by M.A. Suhm.  It was refit to give better                 
C                    tunneling splitting for the ground state.                  
C                                                                               
C   Protocol:                                                                   
C                                                                               
C      PREPOT - initializes the potential's variables                           
C               must be called once before any calls to POT                     
C      POT    - driver for the evaluation of the energy                         
C                                                                               
C   Units:                                                                      
C                                                                               
C      energies    - cm-1                                                       
C      coordinates - distances: bohr; angles: degrees                           
C                                                                               
C   Surfaces:                                                                   
C                                                                               
C      ground electronic state                                                  
C                                                                               
C M. A. SUHM  SEPTEMBER 1990                                                    
C                                                                               
C INPUT : SIX INTERNAL COORDINATES                                              
C --R1,R2,R IN ATOMIC UNITS (BOHR)                                              
C --TH1,TH2,TAU IN RADIANS                                                      
C DEFINITION AS IN QUACK+SUHM MOL.PHYS. 69 (1990) P791                          
C OUTPUT : HF DIMER (SQSBDE) POTENTIAL ENERGY IN HARTREE                        
C                                                                               
C THE PARAMETERS OF THE POTENTIAL ARE READ IN FROM UNIT 1                       
C (FILENAME GA.IN)                                                              
C                                                                               
C THE 6-D POTENTIAL IS BASED ON 1070 AB INITIO POINTS FROM                      
C KOFRANEK ET AL (CHEM PHYS 121(1988) 137) UP TO 42000 CM-1                     
C AND A PERTURBATION TREATMENT BY RIJKS AND WORMER                              
C (JCP 90(1989) 6507) AS WELL AS MONOMER We,WeXe FROM                           
C HUBER/HERZBERG. 29 FREE, 15 MONOMER OR AB INITIO                              
C BASED AND 18 CONSTRAINED PARAMETERS ARE USED.                                 
C SQSBDE IS AN EMPIRICALLY REFINED VERSION WHICH HAS THE                        
C CORRECT EXPERIMENTAL HYDROGEN BOND WAVENUMBER OF                              
C 1062CM-1 (MILLER, ACC. CHEM. RES. 23 (1990) P10) AND                          
C B ROTATIONAL CONSTANT. FOR THIS PURPOSE, THE AB INITIO                        
C RAB COORDINATES WERE SCALED BY 1/1.035.                                       
C THE RMS OF THE WEIGHTED FIT (W=1 FOR R1,2=/=1.7374 OR                         
C R1,2=1.7374 AND E<2220CM-1, W=0.1 FOR R1,2=1.7374 AND                         
C 2220<E<3210CM-1, W=0.001 FOR R1,2=1.7374 AND E>3210CM-1)                      
C IS 32.8CM-1.                                                                  
C                                                                               
C E(1.7445,1.7404,5.1440,9.0007,64.1361,180)=-1.307CM-1 IS                      
C THE ABSOLUTE MINIMUM OF THE POTENTIAL                                         
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
C preparatory calculations for gssqs-routine                                    
C                                                                               
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
      parameter (ma=135)                                                        
C                                                                               
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                           
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)                         
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
                                                                                
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
C   References:      W.C. Necoechea, D.G. Truhlar Chem. Phys. Lett. 248,        
C                    (1995) 182                                                 
C                    M. Quack, M. Suhm Mol. Phys. 69 (1990) P791                
C                    Kofranek et al Chem. Phys. 121 (1988) 137                  
C                    Rijks and Wormer J. Chem. Phys. 90(1989) 6507              
C                                                                               
       REF(1)='W.C. Necoechea, D.G. Truhlar,'                                   
       REF(2)='Chem. Phys. Lett. 248, 182(1995)'                                
       REF(3)='M. Quack, M. Suhm, Mol. Phys. 69, P791(1990)'                    
       REF(4)='Kofranek et al., Chem. Phys. 121, 137(1988)'                     
       REF(5)='Rijks and Wormer, J. Chem. Phys. 90, 6507(1989)'                 
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 9                                                            
      INDEXES(3) = 1                                                            
      INDEXES(4) = 9                                                            
C                                                                               
      IRCTNT=3                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
C prefactors for polynomials                                                    
C                                                                               
      fac(0)=1.                                                                 
      do 101 i=1,2*4                                                            
         fac(i)=fac(i-1)*i                                                      
  101 continue                                                                  
      do 102 l=0,4                                                              
         do 103 m=0,l                                                           
            pf(l,m)=(-1)**m*sqrt((2*l+1)*fac(l-m)/(4*pi*fac(l+m)))              
  103    continue                                                               
  102 continue                                                                  
C                                                                               
C 3j-symbols for potential energy surface                                       
C (abc/d-d0)=g(a,b,c,d)                                                         
C                                                                               
C                                                                               
      fct(1)=1.0d0                                                              
      efct(1)=0.0d0                                                             
      do 104 i=1,400                                                            
         x=i                                                                    
         j=i+1                                                                  
         efct(j)=efct(i)                                                        
         fct(j)=fct(i)*x                                                        
    1    if(fct(j).ge.10.0d0) then                                              
           fct(j)=0.1d0*fct(j)                                                  
           efct(j)=efct(j)+1.0d0                                                
           goto 1                                                               
         endif                                                                  
  104 continue                                                                  
C                                                                               
C                                                                               
      m2=0                                                                      
      m3=0                                                                      
      do 108 j1=0,4*2,2                                                         
         do 107 j2=0,4*2,2                                                      
            do 106 j3=iabs(j1-j2),iabs(j1+j2),4                                 
               do 105 m1=0,min(j1,j2),2                                         
                  m2=-m1                                                        
                  i1=max0(0,(j2-j3-m1)/2,(j1-j3+m2)/2)                          
                  i2=min0((j1+j2-j3)/2,(j1-m1)/2,(j2+m2)/2)                     
                  ia=i1                                                         
                  ib=(j1+j2-j3)/2-i1                                            
                  ic=(j1-m1)/2-i1                                               
                  id=(j2+m2)/2-i1                                               
                  ie=(j3-j2+m1)/2+i1                                            
                  if=(j3-j1-m2)/2+i1                                            
                  sum=(-1)**i1/(fct(ia+1)*fct(ib+1)*fct(ic+1)*fct(id+1)*        
     *                fct(ie+1)*fct(if+1))                                      
                  esum=-efct(ia+1)-efct(ib+1)-efct(ic+1)-efct(id+1)-            
     *                 efct(ie+1)-efct(if+1)                                    
                  nmax=i2-i1                                                    
                  f=1.0d0                                                       
                  if(nmax.gt.0) then                                            
                     do 109 i=nmax,1,-1                                         
                        x=(ib-i+1)*(ic-i+1)*(id-i+1)                            
                        x=x/((ia+i)*(ie+i)*(if+i))                              
                        f=1.0d0-x*f                                             
  109                continue                                                   
                  endif                                                         
                  sum=sum*f                                                     
                  f=fct((j1+j2-j3)/2+1)*fct((j1-j2+j3)/2+1)*                    
     *              fct((j2-j1+j3)/2+1)*fct((j1+m1)/2+1)*                       
     *              fct((j1-m1)/2+1)*                                           
     *              fct((j2+m2)/2+1)*fct((j2-m2)/2+1)*fct((j3+m3)/2+1)*         
     *              fct((j3-m3)/2+1)/fct((j1+j2+j3)/2+2)                        
                  e=efct((j1+j2-j3)/2+1)+efct((j1-j2+j3)/2+1)+                  
     *              efct((j2-j1+j3)/2+1)+efct((j1+m1)/2+1)+                     
     *              efct((j1-m1)/2+1)+                                          
     *              efct((j2+m2)/2+1)+efct((j2-m2)/2+1)+                        
     *              efct((j3+m3)/2+1)+                                          
     *              efct((j3-m3)/2+1)-efct((j1+j2+j3)/2+2)                      
                  cg3j=sum*dsqrt(f)*10.0d0**(esum+0.5d0*e)*                     
     *                 (-1)**((j2-j1-m3)/2)                                     
                  g(j1/2,j2/2,j3/2,m1/2)=cg3j                                   
  105          continue                                                         
  106       continue                                                            
  107    continue                                                               
  108 continue                                                                  
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
      subroutine gsrsqs                                                         
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      parameter (ndmax=1100)                                                    
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
C                                                                               
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                           
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)                         
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
C      dab1=rcm1-3.0d0                                                          
C      dab2=rcm1-4.0d0                                                          
C      dab3=rcm1-4.9d0                                                          
C                                                                               
      dab1=R(1)-3.0d0                                                           
      dab2=R(1)-4.0d0                                                           
      dab3=R(1)-4.09d0                                                          
      do 301 j1=0,4,1                                                           
         do 302 j2=j1,4,1                                                       
            do 303 j3=(j2-j1),(j1+j2),2                                         
               ar(j1,j2,j3)=ac(1,j1,j2,j3)*exp(-0.20d0*dab1**2)+                
     *         ac(2,j1,j2,j3)*exp(-0.35d0*dab2**2)+                             
     *         ac(3,j1,j2,j3)*exp(-0.80d0*dab3**2)                              
               ar(j2,j1,j3)=(-1)**(j1+j2)*ar(j1,j2,j3)                          
  303       continue                                                            
  302    continue                                                               
  301 continue                                                                  
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
C M. A. SUHM  SEPTEMBER 1989                                                    
C                                                                               
      subroutine pot                                                            
C                                                                               
      implicit double precision(a-h,o-z)                                        
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
      parameter (ma=135)                                                        
C                                                                               
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                           
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)                         
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
      dimension x(4),y(4),z(4)                                                  
      dimension p1(0:4,0:4),p2(0:4,0:4),cs(4)                                   
      dimension vp(0:4,0:4,0:2*4)                                               
C                                                                               
      data amh/1837.152756d0/amf/34631.9678d0/                                  
C                                                                               
C     PUT COORDINATES IN PROPER ARRAYS                                          
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C*********************************                                              
C                                *                                              
C   Assign entries from R-array  *                                              
C   to variables used internally *                                              
C   by the POT subroutine        *                                              
C                                *                                              
C*********************************                                              
C                                                                               
      rcm=R(1)                                                                  
      r1=R(2)                                                                   
      r2=R(3)                                                                   
      th1=R(4)                                                                  
      th2=R(5)                                                                  
      tau=R(6)                                                                  
C         rcm=2.5d0                                                             
C         r1=0.9d0                                                              
C         r2=0.9d0                                                              
C         th1=1.04719755120d0                                                   
C         th2=1.04719755120d0                                                   
C         tau=0.78539816340d0                                                   
      phicon=0.01745329252d0                                                    
      pif=1.0d0/(2.0d0*sqrt(pi))                                                
      fac2=amh/(amh+amf)                                                        
      x(3)=0.0d0-r1*cos(th1)*fac2                                               
      y(3)=0.0d0-r1*sin(th1)*fac2                                               
      z(3)=0.0d0                                                                
      x(1)=r1*cos(th1)*(1.0d0-fac2)                                             
      y(1)=r1*sin(th1)*(1.0d0-fac2)                                             
      z(1)=0.0d0                                                                
      x(4)=rcm-r2*cos(th2)*fac2                                                 
      y(4)=0.0d0-r2*sin(th2)*cos(tau)*fac2                                      
      z(4)=0.0d0-r2*sin(th2)*sin(tau)*fac2                                      
      x(2)=rcm+r2*cos(th2)*(1.0d0-fac2)                                         
      y(2)=0.0d0+r2*sin(th2)*cos(tau)*(1.0d0-fac2)                              
      z(2)=0.0d0+r2*sin(th2)*sin(tau)*(1.0d0-fac2)                              
      r12=sqrt((x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2)                    
      r23=sqrt((x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2)                    
      r14=sqrt((x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2)                    
      r34=sqrt((x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2)                    
      do 101 j=1,4                                                              
      cs(j)=cos(j*tau)                                                          
  101 continue                                                                  
      s=sin(th1)                                                                
      c=cos(th1)                                                                
      s2=s*s                                                                    
      s3=s2*s                                                                   
      s4=s3*s                                                                   
      s5=s4*s                                                                   
      c2=c*c                                                                    
      c3=c2*c                                                                   
      c4=c3*c                                                                   
      c5=c4*c                                                                   
      p1(0,0)=pf(0,0)                                                           
      p1(1,0)=pf(1,0)*c                                                         
      p1(2,0)=pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)                                    
      p1(3,0)=pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)                                  
      p1(4,0)=pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)                       
      p1(1,1)=pf(1,1)*s                                                         
      p1(2,1)=pf(2,1)*3.0d0*s*c                                                 
      p1(3,1)=pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)                                  
      p1(4,1)=pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)                                
      p1(2,2)=pf(2,2)*3.0d0*s2                                                  
      p1(3,2)=pf(3,2)*15.0d0*s2*c                                               
      p1(4,2)=pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)                                 
      p1(3,3)=pf(3,3)*15.0d0*s3                                                 
      p1(4,3)=pf(4,3)*105.0d0*s3*c                                              
      p1(4,4)=pf(4,4)*105.0d0*s4                                                
      s=sin(th2)                                                                
      c=cos(th2)                                                                
      s2=s*s                                                                    
      s3=s2*s                                                                   
      s4=s3*s                                                                   
      s5=s4*s                                                                   
      c2=c*c                                                                    
      c3=c2*c                                                                   
      c4=c3*c                                                                   
      c5=c4*c                                                                   
      p2(0,0)=pf(0,0)                                                           
      p2(1,0)=pf(1,0)*c                                                         
      p2(2,0)=pf(2,0)*0.5d0*(3.0d0*c2-1.0d0)                                    
      p2(3,0)=pf(3,0)*0.5d0*(5.0d0*c3-3.0d0*c)                                  
      p2(4,0)=pf(4,0)*0.125d0*(35.0d0*c4-30.0d0*c2+3.0d0)                       
      p2(1,1)=pf(1,1)*s                                                         
      p2(2,1)=pf(2,1)*3.0d0*s*c                                                 
      p2(3,1)=pf(3,1)*1.5d0*s*(5.0d0*c2-1.0d0)                                  
      p2(4,1)=pf(4,1)*2.5d0*s*(7.0d0*c3-3.0d0*c)                                
      p2(2,2)=pf(2,2)*3.0d0*s2                                                  
      p2(3,2)=pf(3,2)*15.0d0*s2*c                                               
      p2(4,2)=pf(4,2)*7.5d0*s2*(7.0d0*c2-1.0d0)                                 
      p2(3,3)=pf(3,3)*15.0d0*s3                                                 
      p2(4,3)=pf(4,3)*105.0d0*s3*c                                              
      p2(4,4)=pf(4,4)*105.0d0*s4                                                
      icu=31                                                                    
      do 102 j1=0,4,1                                                           
         do 103 j2=j1,4,1                                                       
            do 104 j3=(j2-j1),(j1+j2),2                                         
               ac(1,j1,j2,j3)=a(icu)                                            
               ac(2,j1,j2,j3)=a(icu+35)                                         
               ac(3,j1,j2,j3)=a(icu+70)                                         
               icu=icu+1                                                        
  104       continue                                                            
  103    continue                                                               
  102 continue                                                                  
C                                                                               
      call gsrsqs                                                               
C                                                                               
      do 105 j1=0,4,1                                                           
         do 106 j2=0,4,1                                                        
            do 107 j3=abs(j2-j1),(j2+j1),2                                      
               vp(j1,j2,j3)=g(j1,j2,j3,0)*p1(j1,0)*p2(j2,0)                     
               do 108 m=1,min(j1,j2),1                                          
                  vp(j1,j2,j3)=vp(j1,j2,j3)+(-1)**m*2*g(j1,j2,j3,m)*            
     *            cs(m)*p1(j1,m)*p2(j2,m)                                       
  108          continue                                                         
               vp(j1,j2,j3)=vp(j1,j2,j3)*(-1)**(j1-j2)*(2*j3+1)*pif             
  107       continue                                                            
  106    continue                                                               
  105 continue                                                                  
C                                                                               
      e=0.0d0                                                                   
C                                                                               
      do 109 j1=0,4,1                                                           
         do 110 j2=0,4,1                                                        
            do 111 j3=abs(j2-j1),(j2+j1),2                                      
               e=e+vp(j1,j2,j3)*ar(j1,j2,j3)                                    
  111       continue                                                            
  110    continue                                                               
  109 continue                                                                  
C                                                                               
c     e=e+a(1)                                                                  
C                                                                               
      if(rcm.lt.7.2d0)then                                                      
         e=e+vp(0,0,0)*exp(-((7.2d0-rcm)/rcm)**2)*a(2)/rcm**6                   
      else                                                                      
         e=e+vp(0,0,0)*a(2)/rcm**6                                              
      endif                                                                     
      if(rcm.lt.8.4d0)then                                                      
         e=e+vp(0,0,0)*exp(-((8.4d0-rcm)/rcm)**2)*a(3)/rcm**8                   
      else                                                                      
         e=e+vp(0,0,0)*a(3)/rcm**8                                              
      endif                                                                     
      if(rcm.lt.9.6d0)then                                                      
         e=e+vp(0,0,0)*exp(-((9.6d0-rcm)/rcm)**2)*a(4)/rcm**10                  
      else                                                                      
         e=e+vp(0,0,0)*a(4)/rcm**10                                             
      endif                                                                     
      if(rcm.lt.9.5d0) then                                                     
         e=e+vp(1,1,2)*a(5)/rcm**3*exp(-((9.5d0-rcm)/rcm)**2)                   
      else                                                                      
         e=e+vp(1,1,2)*a(5)/rcm**3                                              
      endif                                                                     
      if(rcm.lt.10.5d0) then                                                    
         e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4*                                 
     *     exp(-((10.5d0-rcm)/rcm)**2)                                          
      else                                                                      
         e=e+(vp(2,1,3)-vp(1,2,3))*a(6)/rcm**4                                  
      endif                                                                     
      if(rcm.lt.11.5d0) then                                                    
         e=e+vp(2,2,4)*a(7)/rcm**5*exp(-((11.5d0-rcm)/rcm)**2)                  
         e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5*                                 
     *     exp(-((11.5d0-rcm)/rcm)**2)                                          
      else                                                                      
         e=e+vp(2,2,4)*a(7)/rcm**5                                              
         e=e+(vp(3,1,4)+vp(1,3,4))*a(8)/rcm**5                                  
      endif                                                                     
      e=e+a(9)*exp(-0.5d0*r12)                                                  
      e=e+a(10)*exp(-1.5d0*r12)                                                 
      e=e+a(13)*exp(-1.5d0*r34)                                                 
      e=e+a(17)*(exp(-10.0d0*r14)+exp(-10.0d0*r23))                             
      if(rcm.lt.4.5d0)then                                                      
         e=e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+                         
     *     (3.0d0*(cos(th2))**2+1.0d0))*                                        
     *     exp(-((4.5d0-rcm)/rcm)**2)                                           
      else                                                                      
         e=e+a(18)/rcm**6*((3.0d0*(cos(th1))**2+1.0d0)+                         
     *     (3.0d0*(cos(th2))**2+1.0d0))                                         
      endif                                                                     
C                                                                               
C monomer terms                                                                 
C                                                                               
      e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r23-3.0d0)**2))*                   
     *         (1.0d0-exp(-1.1953d0*                                            
     *         (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r14-3.0d0)))                     
     *                                                      )*                  
     *         (r1-1.7374d0*(1.0d0                                              
     *         -a(20)*                                                          
     *         (1.0d0-tanh(a(21)*(r14-3.0d0)))                                  
     *                                        ))))**2                           
      e=e+47635.0d0*(1.0d0-a(19)*exp(-0.1d0*(r14-3.0d0)**2))*                   
     *         (1.0d0-exp(-1.1953d0*                                            
     *         (1.0d0-a(26)*(1.0d0-tanh(a(27)*(r23-3.0d0)))                     
     *                                                      )*                  
     *         (r2-1.7374d0*(1.0d0                                              
     *         -a(20)*                                                          
     *         (1.0d0-tanh(a(21)*(r23-3.0d0)))                                  
     *                                        ))))**2                           
c                                                                               
c     next line corrects for asymptotic diatomic potential which is             
c     added in separately.                                                      
c                                                                               
      e=e-(ptntl(r1-1.7374d0)+ptntl(r2-1.7374d0))*219474.63067d0                
      e=e+a(28)*(exp(-4.0d0*(r14-r1)**2)+exp(-4.0d0*(r23-r2)**2))               
      e=e+a(22)*exp(-4.0d0*r14)*exp(-4.0d0*r23)                                 
      e=e+a(23)                                                                 
      e=e+a(24)*(exp(-r14)+exp(-r23))*                                          
     *  (exp(-15.0d0*r1)+exp(-15.0d0*r2))                                       
C                                                                               
C   Units conversion:  Conversion to atomic units                               
C                                                                               
      e=e/219474.63067d0                                                        
C                                                                               
      ENGYGS=e                                                                  
C                                                                               
C                                                                               
      CALL EUNITZERO                                                            
C                                                                               
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
      end                                                                       
C                                                                               
      FUNCTION PTNTL(X)                                                         
C                                                                               
C     HF DIATOMIC POTENTIAL CURVE FOR SQSBDE SURFACE, CALCULATED USING          
C     PARAMETERS GIVEN IN JCP, VOL. 95(5), P. 31.                               
C     NOTICE THAT THE PARAMETER DM IS GIVEN IN ATOMIC UNITS HERE AND            
C     IN CM-1 IN THE JCP REFERENCE.                                             
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
      parameter (ma=135)                                                        
      PARAMETER (AM=1.1953D0,DM=0.21704D0,ONE=1.0D0)                            
C                                                                               
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                           
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)                         
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
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
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
                                                                                
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
C      SAVE                                                                     
C                                                                               
      PTNTL=DM*(ONE-EXP(-AM*X))**2                                              
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      BLOCK DATA PTPACM                                                         
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
      parameter (ma=135)                                                        
C                                                                               
      COMMON/PT3CM/ EZERO(ISURF+1)                                              
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON/POTCM/  A(MA),DA(MA),fac(0:2*4),                                   
     +               fct(401),efct(401),                                        
     +               ar(0:4,0:4,0:2*4),                                         
     +               ac(3,0:4,0:4,0:2*4)                                        
      common/ang/    g(0:4,0:4,0:2*4,0:4),pf(0:4,0:4)                           
C                                                                               
       DATA NASURF /1,35*0/                                                     
       DATA NDER /0/                                                            
       DATA NFLAG /1,1,15*0,6,0,0/                                              
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/4,0,0/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /4/                                                           
C                                                                               
      DATA A(1),DA(1)  / 1558.0D0,  0.0D0/                                      
      DATA A(2),DA(2)  / -202900000.0D0,  0.0D0/                                
      DATA A(3),DA(3)  / -3703000000.0D0,  0.0D0/                               
      DATA A(4),DA(4)  / -11350000000.0D0,  0.0D0/                              
      DATA A(5),DA(5)  / -1736693.094999999D0,  0.0D0/                          
      DATA A(6),DA(6)  / -4416487.722000003D0,  0.0D0/                          
      DATA A(7),DA(7)  / 16012389.50999999D0,  0.0D0/                           
      DATA A(8),DA(8)  / -6389444.854999989D0,  0.0D0/                          
      DATA A(9),DA(9)  / -866.86737555133D0,  117.5897899275601D0/              
      DATA A(10),DA(10)  / 129789.3634701688D0,  1937.563168749279D0/           
      DATA A(11),DA(11)  / 0.0D0,  0.0D0/                                       
      DATA A(12),DA(12)  / 0.0D0,  0.0D0/                                       
      DATA A(13),DA(13)  / 2622184.756658763D0,  29726.63014604442D0/           
      DATA A(14),DA(14)  / 0.0D0,  0.0D0/                                       
      DATA A(15),DA(15)  / 0.0D0,  0.0D0/                                       
      DATA A(16),DA(16)  / 0.0D0,  0.0D0/                                       
      DATA A(17),DA(17)  / 184254701843.3105D0,  17553092934.44629D0/           
      DATA A(18),DA(18)  / -320000.0D0,  0.0D0/                                 
      DATA A(19),DA(19)  / 3.326438943860243D-2,  1.353645657747024D-3/         
      DATA A(20),DA(20)  / -4.576887206887037D-3,  1.865017706495818D-4/        
      DATA A(21),DA(21)  / 0.4097306660396534D0,  2.777020547964315D-2/         
      DATA A(22),DA(22)  / 32623551943283.75D0,  3630087056915.047D0/           
      DATA A(23),DA(23)  / 0.0D0,  0.0D0/                                       
      DATA A(24),DA(24)  / 71754497894202.0D0,  2525041881368.016D0/            
      DATA A(25),DA(25)  / 0.0D0,  0.0D0/                                       
      DATA A(26),DA(26)  / 5.032387872479127D-2,  1.860431838521345D-3/         
      DATA A(27),DA(27)  / 1.46832057915281D0,  9.71324803167759D-2/            
      DATA A(28),DA(28)  / 3384.631559286528D0,  85.55678229458908D0/           
      DATA A(29),DA(29)  / 0.0D0,  0.0D0/                                       
      DATA A(30),DA(30)  / 0.0D0,  0.0D0/                                       
      DATA A(31),DA(31)  / -160267.101272529D0,  4055.338776717705D0/           
      DATA A(32),DA(32)  / 0.0D0,  0.0D0/                                       
      DATA A(33),DA(33)  / 0.0D0,  0.0D0/                                       
      DATA A(34),DA(34)  / 0.0D0,  0.0D0/                                       
      DATA A(35),DA(35)  / 0.0D0,  0.0D0/                                       
      DATA A(36),DA(36)  / -16268.89465981041D0,  541.1396773442175D0/          
      DATA A(37),DA(37)  / 0.0D0,  0.0D0/                                       
      DATA A(38),DA(38)  / 0.0D0,  0.0D0/                                       
      DATA A(39),DA(39)  / 2802.576101972474D0,  56.80382831682482D0/           
      DATA A(40),DA(40)  / 0.0D0,  0.0D0/                                       
      DATA A(41),DA(41)  / 0.0D0,  0.0D0/                                       
      DATA A(42),DA(42)  / 0.0D0,  0.0D0/                                       
      DATA A(43),DA(43)  / 0.0D0,  0.0D0/                                       
      DATA A(44),DA(44)  / 0.0D0,  0.0D0/                                       
      DATA A(45),DA(45)  / 0.0D0,  0.0D0/                                       
      DATA A(46),DA(46)  / 2860.063943395988D0,  44.36333848084973D0/           
      DATA A(47),DA(47)  / 0.0D0,  0.0D0/                                       
      DATA A(48),DA(48)  / 0.0D0,  0.0D0/                                       
      DATA A(49),DA(49)  / 0.0D0,  0.0D0/                                       
      DATA A(50),DA(50)  / 0.0D0,  0.0D0/                                       
      DATA A(51),DA(51)  / 0.0D0,  0.0D0/                                       
      DATA A(52),DA(52)  / 509.1020221188392D0,  18.73758800462167D0/           
      DATA A(53),DA(53)  / 0.0D0,  0.0D0/                                       
      DATA A(54),DA(54)  / 0.0D0,  0.0D0/                                       
      DATA A(55),DA(55)  / 0.0D0,  0.0D0/                                       
      DATA A(56),DA(56)  / 0.0D0,  0.0D0/                                       
      DATA A(57),DA(57)  / 0.0D0,  0.0D0/                                       
      DATA A(58),DA(58)  / 0.0D0,  0.0D0/                                       
      DATA A(59),DA(59)  / 0.0D0,  0.0D0/                                       
      DATA A(60),DA(60)  / 0.0D0,  0.0D0/                                       
      DATA A(61),DA(61)  / 0.0D0,  0.0D0/                                       
      DATA A(62),DA(62)  / 0.0D0,  0.0D0/                                       
      DATA A(63),DA(63)  / 0.0D0,  0.0D0/                                       
      DATA A(64),DA(64)  / 0.0D0,  0.0D0/                                       
      DATA A(65),DA(65)  / 0.0D0,  0.0D0/                                       
      DATA A(66),DA(66)  / 0.0D0,  0.0D0/                                       
      DATA A(67),DA(67)  / 0.0D0,  0.0D0/                                       
      DATA A(68),DA(68)  / -3467.234170178781D0,  52.91056956195189D0/          
      DATA A(69),DA(69)  / 0.0D0,  0.0D0/                                       
      DATA A(70),DA(70)  / 0.0D0,  0.0D0/                                       
      DATA A(71),DA(71)  / 0.0D0,  0.0D0/                                       
      DATA A(72),DA(72)  / 0.0D0,  0.0D0/                                       
      DATA A(73),DA(73)  / 2304.632858038298D0,  44.02327215328341D0/           
      DATA A(74),DA(74)  / 0.0D0,  0.0D0/                                       
      DATA A(75),DA(75)  / 0.0D0,  0.0D0/                                       
      DATA A(76),DA(76)  / 0.0D0,  0.0D0/                                       
      DATA A(77),DA(77)  / 0.0D0,  0.0D0/                                       
      DATA A(78),DA(78)  / 0.0D0,  0.0D0/                                       
      DATA A(79),DA(79)  / 0.0D0,  0.0D0/                                       
      DATA A(80),DA(80)  / 0.0D0,  0.0D0/                                       
      DATA A(81),DA(81)  / 0.0D0,  0.0D0/                                       
      DATA A(82),DA(82)  / 0.0D0,  0.0D0/                                       
      DATA A(83),DA(83)  / 0.0D0,  0.0D0/                                       
      DATA A(84),DA(84)  / -644.8150830210325D0,  18.1151965381963D0/           
      DATA A(85),DA(85)  / 0.0D0,  0.0D0/                                       
      DATA A(86),DA(86)  / 0.0D0,  0.0D0/                                       
      DATA A(87),DA(87)  / 0.0D0,  0.0D0/                                       
      DATA A(88),DA(88)  / 0.0D0,  0.0D0/                                       
      DATA A(89),DA(89)  / 0.0D0,  0.0D0/                                       
      DATA A(90),DA(90)  / 0.0D0,  0.0D0/                                       
      DATA A(91),DA(91)  / 0.0D0,  0.0D0/                                       
      DATA A(92),DA(92)  / 0.0D0,  0.0D0/                                       
      DATA A(93),DA(93)  / 0.0D0,  0.0D0/                                       
      DATA A(94),DA(94)  / 0.0D0,  0.0D0/                                       
      DATA A(95),DA(95)  / 0.0D0,  0.0D0/                                       
      DATA A(96),DA(96)  / 0.0D0,  0.0D0/                                       
      DATA A(97),DA(97)  / 0.0D0,  0.0D0/                                       
      DATA A(98),DA(98)  / 0.0D0,  0.0D0/                                       
      DATA A(99),DA(99)  / 0.0D0,  0.0D0/                                       
      DATA A(100),DA(100)  / 0.0D0,  0.0D0/                                     
      DATA A(101),DA(101)  / 21489.15561258642D0,  989.768118754484D0/          
      DATA A(102),DA(102)  / 1895.153684423094D0,  111.7819362295832D0/         
      DATA A(103),DA(103)  / 0.0D0,  0.0D0/                                     
      DATA A(104),DA(104)  / 768.3954452456346D0,  18.43253058647917D0/         
      DATA A(105),DA(105)  / -395.2301759253078D0,  12.84348126946628D0/        
      DATA A(106),DA(106)  / 2041.09421622329D0,  246.1094667134039D0/          
      DATA A(107),DA(107)  / -2422.952146150754D0,  70.83017811025411D0/        
      DATA A(108),DA(108)  / 0.0D0,  0.0D0/                                     
      DATA A(109),DA(109)  / 0.0D0,  0.0D0/                                     
      DATA A(110),DA(110)  / -297.3729718893483D0,  19.35802742071735D0/        
      DATA A(111),DA(111)  / 0.0D0,  0.0D0/                                     
      DATA A(112),DA(112)  / 0.0D0,  0.0D0/                                     
      DATA A(113),DA(113)  / 39.46815694764132D0,  11.34523118392423D0/         
      DATA A(114),DA(114)  / 0.0D0,  0.0D0/                                     
      DATA A(115),DA(115)  / 720.5552230122485D0,  33.69084291038598D0/         
      DATA A(116),DA(116)  / 0.0D0,  0.0D0/                                     
      DATA A(117),DA(117)  / 0.0D0,  0.0D0/                                     
      DATA A(118),DA(118)  / 0.0D0,  0.0D0/                                     
      DATA A(119),DA(119)  / 0.0D0,  0.0D0/                                     
      DATA A(120),DA(120)  / 0.0D0,  0.0D0/                                     
      DATA A(121),DA(121)  / 0.0D0,  0.0D0/                                     
      DATA A(122),DA(122)  / 0.0D0,  0.0D0/                                     
      DATA A(123),DA(123)  / 0.0D0,  0.0D0/                                     
      DATA A(124),DA(124)  / 0.0D0,  0.0D0/                                     
      DATA A(125),DA(125)  / 0.0D0,  0.0D0/                                     
      DATA A(126),DA(126)  / 0.0D0,  0.0D0/                                     
      DATA A(127),DA(127)  / 0.0D0,  0.0D0/                                     
      DATA A(128),DA(128)  / 0.0D0,  0.0D0/                                     
      DATA A(129),DA(129)  / 0.0D0,  0.0D0/                                     
      DATA A(130),DA(130)  / 0.0D0,  0.0D0/                                     
      DATA A(131),DA(131)  / 0.0D0,  0.0D0/                                     
      DATA A(132),DA(132)  / 0.0D0,  0.0D0/                                     
      DATA A(133),DA(133)  / 0.0D0,  0.0D0/                                     
      DATA A(134),DA(134)  / 0.0D0,  0.0D0/                                     
      DATA A(135),DA(135)  / 0.0D0,  0.0D0/                                     
C                                                                               
       END                                                                      
