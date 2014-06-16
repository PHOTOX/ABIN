      subroutine shake(x,y,z,vx,vy,vz,iq,iv) 
      use mod_array_size
      use mod_general
      use mod_system, ONLY: am,NShake,IShake1,IShake2,dshake,shake_tol
      !implicit real*8(a-h,o-z)
      implicit none
      real*8 :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      integer,intent(in) :: iq,iv
      real*8  :: mi,mj,mij,agama
      integer :: iw,i,j,ixshake,iiter,maxcycle,itest
      real*8  :: dijiter2,xij,yij,zij,rij2
      real*8  :: xdotij,ydotij,zdotij,dot
      real*8  :: xiiter,yiiter,ziiter,xjiter,yjiter,zjiter
      !TODO implicit none
      
      maxcycle=1000
      do iw=1,nwalk

      if(iw.eq.2.and.istage.eq.2) exit   !shake only centroid variable i.e. first normal mode

! PS&DHchange: SHAKE algorithm implemented, iteratively solved until convergence
! criterion agama < shake_tol are met for all bonds.
! We stop program if we do not converge after 1000 steps.
! Velocities recalculated in a primitive way by method of finite differences.

       if(iq.eq.1)then

       Do IIter=1,maxcycle

!-------variable which controls convergence,0=converged, 1=not converged
        itest=0

        Do ixshake=1,NShake
        i=IShake1(ixshake)
        j=IShake2(ixshake)
        mi=am(i)
        mj=am(j)
        mij=((1/mi))+(1/(mj))
        xiiter=x(i,iw)
        yiiter=y(i,iw)
        ziiter=z(i,iw)
        xjiter=x(j,iw)
        yjiter=y(j,iw)
        zjiter=z(j,iw)

        dijiter2=(xiiter-xjiter)**2+(yiiter-yjiter)**2 + (ziiter-zjiter)**2
        agama=dijiter2-dshake(ixshake)
        agama=agama/(4*dijiter2*mij)

        if(abs(agama).gt.shake_tol)then
         itest=1
! Rapaport, p.275, here with masses(small modification)
         x(i,iw)=xiiter-agama*(xiiter-xjiter)/mi
         y(i,iw)=yiiter-agama*(yiiter-yjiter)/mi
         z(i,iw)=ziiter-agama*(ziiter-zjiter)/mi
         x(j,iw)=xjiter+agama*(xiiter-xjiter)/mj
         y(j,iw)=yjiter+agama*(yiiter-yjiter)/mj
         z(j,iw)=zjiter+agama*(ziiter-zjiter)/mj
        endif

! ixshake loop
        enddo

        if(itest.eq.0)then
         exit
        endif
! IITER loop
        enddo
        
        if (iiter.ge.maxcycle)then
          write(*,*)'Shake not converged after',maxcycle,'iterations. Exiting...'      
          stop
        endif
        
        endif


! velocity update
        if(iv.eq.1)then
       Do IIter=1,maxcycle

!-------variable which controls convergence,0=converged, 1=not converged
        itest=0
       
        Do ixshake=1,NShake
        i=IShake1(ixshake)
        j=IShake2(ixshake)
        mi=am(i)
        mj=am(j)
        mij=((1/mi))+(1/(mj))
        
        xij=x(j,iw)-x(i,iw)
        yij=y(j,iw)-y(i,iw)
        zij=z(j,iw)-z(i,iw)
        xdotij=vx(j,iw)-vx(i,iw)
        ydotij=vy(j,iw)-vy(i,iw)
        zdotij=vz(j,iw)-vz(i,iw)

        dot=xij*xdotij+yij*ydotij+zij*zdotij
        rij2=xij**2+yij**2+zij**2
        agama=-dot/(2*rij2*mij)

        if(abs(agama).gt.shake_tol)then
           itest=1
           vx(i,iw)=vx(i,iw)-agama*xij/mi
           vy(i,iw)=vy(i,iw)-agama*yij/mi
           vz(i,iw)=vz(i,iw)-agama*zij/mi
           vx(j,iw)=vx(j,iw)+agama*xij/mj
           vy(j,iw)=vy(j,iw)+agama*yij/mj
           vz(j,iw)=vz(j,iw)+agama*zij/mj
        endif

! ixshake loop
        enddo

        if(itest.eq.0)then
           exit
        endif
! IITER loop
        enddo
        if (iiter.ge.maxcycle)then
           write(*,*)'Velocity shake not converged after',maxcycle,'iterations. Exiting...'      
           stop
        endif
        endif


      enddo
   

      return  
      end

