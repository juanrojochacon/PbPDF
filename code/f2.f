      IMPLICIT NONE

      double precision  X,q,f2,f2c,f2b,fl,flc,flb,q2
      double precision mc,mt,mb,asmur,mur,fr2
      integer AAA,ZZZ,ipn,i,iordmstw,pset,isetdssz,j,isetmstw
      integer dim_f2
      parameter (dim_f2=18) !number of data of the experiment
      double precision Q2_Au_139(dim_f2),X_Au_139(dim_f2),
     1                 Au_139(dim_f2),ES_Au_139(dim_f2),
     1                 EE_Au_139(dim_f2),y(dim_f2)
      double precision f2p(0:40,dim_f2),f2n(0:40,dim_f2),
     1                 f2d(0:40,dim_f2)
      double precision f2dssz(51,dim_f2),f2pbdssz(51,dim_f2)
      double precision f2eps09(31,dim_f2),f2pbeps09(31,dim_f2)
      double precision f2finaldssz(51,dim_f2),f2finaleps09(31,dim_f2)
c      double precision ersf2n(dim_f2),erf2n(dim_f2) 
c      double precision ersf2p(dim_f2),erf2p(dim_f2)
c      double precision ersf2d(dim_f2),erf2d(dim_f2)
      double precision ersf2dssz(dim_f2),ersf2pbdssz(dim_f2)
      double precision ersf2eps09(dim_f2),ersf2pbeps09(dim_f2) 
      double precision erf2dssz(dim_f2),erf2pbdssz(dim_f2)
      double precision erf2eps09(dim_f2),erf2pbeps09(dim_f2)
      double precision factoreps09(31,dim_f2),factordssz(51,dim_f2)
      double precision ersfactordssz(dim_f2),erfactordssz(dim_f2)
      double precision ersfactoreps09(dim_f2),erfactoreps09(dim_f2)
      double precision ersf2finaldssz(dim_f2),erf2finaldssz(dim_f2)
      double precision ersf2finaleps09(dim_f2),erf2finaleps09(dim_f2)
      double precision stateps09(dim_f2),statdssz(dim_f2) 
      double precision systeps09(dim_f2),systdssz(dim_f2)   
      COMMON/A/AAA,ZZZ 
      COMMON/iordCommon/iordmstw
      common/pseteps09/pset 
      common/psetdssz/isetdssz 
      common/psetmstw/isetmstw

      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mc,mb,ASMUR,alphaSMZ,alphaSorder,alphaSnfmax



      iordmstw=1 ! no other value available for nPDFs

     	    FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2
	    MUR = 1.D0                ! input mu_r in GeV
	    ASMUR = 0.49128D0         ! input value of alpha_s at mu_r
	    MC = 1.4D0                ! charm quark mass
	    MB = 4.75D0               ! bottom quark mass
	    MT = 175d10               ! top quark mass
	    
	    

        CALL WATE96
	    CALL INITALPHAS(IORDmstw, FR2, mur, ASMUR, MC, MB, MT)

      
c      ipn=1 ! 1 = proton, 2 = neutron, 3 = nucleus DSSZ, 4 = nucleus EPS09 
      
      isetdssz=0
      pset=1
c     reading data      
      open(20, file='nmched.data',status='unknown')
      read(20,*)
      read(20,*)
      do i=1,dim_f2,1
       read(20,*) X_Au_139(i),Q2_Au_139(i),y(i),Au_139(i),EE_Au_139(i)
     1 ,ES_Au_139(i)
      enddo
      close(20)   
      open (21,FILE='nmched_f2_2.data',STATUS='unknown')
      open (22,FILE='nmched_f2_factor.data',STATUS='unknown')
      
       
c     F2 proton (mstw2008)     
      ipn=1
      do j=0,40
       isetmstw=j
       do i=1,dim_f2,1     
        x=X_Au_139(i) 
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2p(j,i)=f2
       enddo
      enddo  

c      do i=1,dim_f2,1
c       ersf2p(i)=0.d0
c       do j=1,20
c        ersf2p(i)=ersf2p(i)+((f2p(2*j,i)-
c     1  f2p(2*j-1,i))**2.d0)
c       enddo
c       erf2p(i)=dsqrt(ersf2p(i))/2.d0
c      enddo   

   
c     F2 neutron (mstw2008)     
      ipn=2 
      do j=0,40
       isetmstw=j
       do i=1,dim_f2,1     
        x=X_Au_139(i)
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2n(j,i)=f2
       enddo 
      enddo   
  
c      do i=1,dim_f2,1
c       ersf2n(i)=0.d0
c       do j=1,20
c        ersf2n(i)=ersf2n(i)+((f2n(2*j,i)-
c     1  f2n(2*j-1,i))**2.d0)
c       enddo
c       erf2n(i)=dsqrt(ersf2n(i))/2.d0
c      enddo      

c     F2 deuterium
      do j=0,40
       do i=1,dim_f2,1 
        f2d(j,i)=(f2n(j,i)+f2p(j,i))/2.d0
       enddo  
      enddo 

c      do i=1,dim_f2,1 
c       erf2d(i)=(erf2n(i)+erf2p(i))/2.d0
c      enddo  

      isetmstw=0   !central value (when calling dssz and eps, we are calling ratios so we need to multiply by the pdf central value)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     F2 nucleus A (the one we want to convert)

      AAA = 4  !1 for ipn=1,2. Change if IPN=3 or IPN=4
      ZZZ = 2 !1 for ipn=1, 0 for ipn=2. change if IPN=3 or IPN=4 
      
c     F2 Au DSSZ
      ipn=3
           
      do j=1,51
       if (j.eq.1) then
        ISETdssz=0
       elseif (j.le.26) then
        ISETdssz=j-1
       else
        ISETdssz=-j+26
       endif
       call DSSZINI(ISETdssz)
       do i=1,dim_f2,1     
        x=X_Au_139(i) 
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2dssz(j,i)=f2
       enddo
      enddo

      do i=1,dim_f2,1
       ersf2dssz(i)=0.d0
       do j=1,25
        ersf2dssz(i)=ersf2dssz(i)+((f2dssz(j+1,i)
     1   -f2dssz(j+26,i))**2.d0)
       enddo
       erf2dssz(i)=dsqrt(ersf2dssz(i))/2.d0
      enddo
 
c     F2 Au eps09    
      ipn=4
      
      do j=1,31
       pset=j
       do i=1,dim_f2,1     
        x=X_Au_139(i)
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2eps09(j,i)=f2
       enddo  
      enddo

      do i=1,dim_f2,1
       ersf2eps09(i)=0.d0
       do j=1,15
        ersf2eps09(i)=ersf2eps09(i)+((f2eps09(2*j,i)-
     1  f2eps09(2*j+1,i))**2.d0)
       enddo
       erf2eps09(i)=dsqrt(ersf2eps09(i))/2.d0
      enddo   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     F2(Pb)

      AAA = 208 !  1 for ipn=1,2, change if IPN=3 or IPN=4
      ZZZ = 82  !  1 for ipn=1, 0 for ipn=2, change if IPN=3 or IPN=4
         
c     F2 Pb DSSZ
      ipn=3
           
      do j=1,51
       if (j.eq.1) then
        ISETdssz=0
       elseif (j.le.26) then
        ISETdssz=j-1
       else
        ISETdssz=-j+26
       endif
       call DSSZINI(ISETdssz)     
       do i=1,dim_f2,1     
        x=X_Au_139(i)
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2pbdssz(j,i)=f2
       enddo
      enddo
      
      do i=1,dim_f2,1
       ersf2pbdssz(i)=0.d0
       do j=1,25
        ersf2pbdssz(i)=ersf2pbdssz(i)+((f2pbdssz(j+1,i)-
     1   f2pbdssz(j+26,i))**2.d0)
       enddo
       erf2pbdssz(i)=dsqrt(ersf2pbdssz(i))/2.d0
      enddo
 
c     F2 Pb eps09    
      ipn=4
      
      do j=1,31
      pset=j
c     pset=1 central value
       do i=1,dim_f2,1     
        x=X_Au_139(i)
        q2=Q2_Au_139(i)
        q=dsqrt(q2)
        call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
        f2pbeps09(j,i)=f2
       enddo  
      enddo
      do i=1,dim_f2,1
       ersf2pbeps09(i)=0.d0
       do j=1,15
        ersf2pbeps09(i)=ersf2pbeps09(i)+((f2pbeps09(2*j,i)-
     1  f2pbeps09(2*j+1,i))**2.d0)
       enddo
       erf2pbeps09(i)=dsqrt(ersf2eps09(i))/2.d0
      enddo 


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Final result
      
      do j=1,31
       do i=1,dim_f2,1
        factoreps09(j,i)= f2d(0,i)*f2pbeps09(j,i)/f2eps09(j,i)         
        f2finaleps09(j,i)= Au_139(i)*factoreps09(j,i)
       enddo
      enddo

      do i=1,dim_f2,1
       ersfactoreps09(i)=0.d0
       do j=1,15
        ersfactoreps09(i)=ersfactoreps09(i)+((factoreps09(2*j,i)-
     1  factoreps09(2*j+1,i))**2.d0)
       enddo
       erfactoreps09(i)=dsqrt(ersfactoreps09(i))/2.d0
       erf2finaleps09(i)=erfactoreps09(i)*f2finaleps09(1,i)
       stateps09(i)=factoreps09(1,i)*EE_Au_139(i)
       systeps09(i)=factoreps09(1,i)*ES_Au_139(i)
      enddo   
      
      do j=1,51
       do i=1,dim_f2,1  
        factordssz(j,i)=f2d(0,i)*f2pbdssz(j,i)/f2dssz(j,i)     
        f2finaldssz(j,i)= Au_139(i)*factordssz(j,i)   
       enddo
      enddo

      do i=1,dim_f2,1
       ersfactordssz(i)=0.d0
       do j=1,25
        ersfactordssz(i)=ersfactordssz(i)+((factordssz(2*j,i)-
     1  factordssz(2*j+1,i))**2.d0)
       enddo
       erfactordssz(i)=dsqrt(ersfactordssz(i))/2.d0
       erf2finaldssz(i)=erfactordssz(i)*f2finaldssz(1,i)
       statdssz(i)=factordssz(1,i)*EE_Au_139(i)
       systdssz(i)=factordssz(1,i)*ES_Au_139(i)
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
      
       write(21,*) '  X     ','Q^2  ',' y  ','  F2(Pb)EPS09 ',
     1  'S(F2(Pb)EPS09) ',
     1'statEPS09  ',' systEPS09  ', '  F2(Pb)DSSZ ','  S(F2(Pb)DSSZ) ',
     1 ' statDSSZ  ','  systDSSZ  '
       write(21,*) '-------------------------------------------------',
     1 '-------------------------------------------------------------',
     1 '---------'
       do i=1,dim_f2,1 
       if(q2_au_139(i).GE.1.69d0) then  !Starting point of EPS09 evolution
         if(i.EQ.1) then 
          write (21,1002) X_au_139(i),Q2_au_139(i),y(i),
     1    f2finaleps09(1,i),
     1    erf2finaleps09(i),stateps09(i),systeps09(i),f2finaldssz(1,i),
     1    erf2finaldssz(i),statdssz(i),systdssz(i)

         else
          write (21,1001) X_au_139(i),Q2_au_139(i),y(i),
     1    f2finaleps09(1,i),
     1    erf2finaleps09(i),stateps09(i),systeps09(i),f2finaldssz(1,i),
     1    erf2finaldssz(i),statdssz(i),systdssz(i)
       end if
      end if
      enddo
1002  format(f6.4,f6.1,f6.2,8e13.6)
1001  format(f6.4,f6.1,f6.2,8e13.6)
      close(21)


      write(22,*) '  X     ','Q^2  ', 'F2pbeps09', 'S(F2pbeps09)',
     1 ' F2Aeps09', ' S(F2Apbeps09)', ' factoreps09', ' erfactoreps09',
     1  ' F2pbdssz', ' S(F2pbdssz)',' F2Adssz', ' S(F2Adssz)',  
     1  ' factordssz',' erfactordssz'
      do i=1,dim_f2,1 
       if(q2_au_139(i).GE.1.69d0) then
       write (22,1005) X_au_139(i),Q2_au_139(i),f2pbeps09(1,i), 
     1 erf2pbeps09(i),f2eps09(1,i),erf2eps09(i),factoreps09(1,i),
     1 erfactoreps09(i),f2pbdssz(1,i),erf2pbdssz(i),f2dssz(1,i),
     1 erf2dssz(i),factordssz(1,i),erfactordssz(i)
       endif
      enddo
1005  format(f6.4,f6.1,12e13.6)
      close(22)
      end program
