!----- Velocity verlet ALGORITHM                              Daniel Hollas,2012
!at this moment does not contain shake routines, which are only in respashake
!function

! GLE and NHC thermostats available at the moment

      subroutine verletstep(x,y,z,px,py,pz,amt,amg,dt,eclas,fxc,fyc,fzc)
      use mod_array_size
      use mod_general
      use mod_nhc, ONLY:inose,imasst
      use mod_gle
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8 fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
      real*8 px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
      real*8  :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
      real*8  :: dt,eclas,equant
      real*8  :: ekin_p
      integer :: iat,iw
      integer :: iq,iv


!---THERMOSTATS------------------
      if(inose.eq.1)then
       if (imasst.eq.1)then
        call shiftNHC_yosh_mass(px,py,pz,amt,dt/2.0d0)
       else
        call shiftNHC_yosh(px,py,pz,amt,dt/2.0d0)
       endif
      endif

      if (inose.eq.3)  call wn_step(px,py,pz,amt)
      if (inose.eq.2)then
       langham=langham+ekin_p(px,py,pz)
       call gle_step(px,py,pz,amt)
       langham=langham-ekin_p(px,py,pz)
      endif

       call shiftP (px,py,pz,fxc,fyc,fzc,amt,dt/2.0d0)

       if(conatom.gt.0)then
        do iw=1,nwalk
         do iat=1,conatom
          px(iat,iw)=0.0
          py(iat,iw)=0.0
          pz(iat,iw)=0.0
         enddo
        enddo
       endif

      call shiftX(x,y,z,px,py,pz,amt,dt)

      call force_clas(fxc,fyc,fzc,x,y,z,eclas)

      call shiftP (px,py,pz,fxc,fyc,fzc,amt,dt/2.0d0)

     !---THERMOSTATS------------------

      if (inose.eq.3)  call wn_step(px,py,pz,amt)
      if (inose.eq.2)then
       langham=langham+ekin_p(px,py,pz)
       call gle_step(px,py,pz,amt)
       langham=langham-ekin_p(px,py,pz)
      endif

      if(inose.eq.1)then
       if (imasst.eq.1)then
        call shiftNHC_yosh_mass(px,py,pz,amt,dt/2.0d0)
       else
        call shiftNHC_yosh(px,py,pz,amt,dt/2.0d0)
       endif
      endif


      end


