! initial version                         P. Slavicek a kol., Mar 25, 2009
!
!--------------------------------------------------------------------------
      SUBROUTINE temperature(px,py,pz,amt,dt,eclas)
      use mod_array_size
      use mod_general
      use mod_estimators, ONLY:est_temp_cumul,entot_cumul
      use mod_system, ONLY: nshake
      use mod_nhc, ONLY: inose,nhcham,calc_nhcham
      use mod_gle, ONLY:langham
      implicit none
      integer :: iw,iat
      real*8,intent(in)  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
      real*8,intent(in)  :: amt(npartmax,nwalkmax)
      real*8,intent(in)  :: dt,eclas
      real*8  :: est_temp,temp1,ekin_mom
      real*8  :: it2
      
      it2=it/ncalc
      ekin_mom=0.0d0
      do iw=1,nwalk
       do iat=1,natom
        temp1=px(iat,iw)**2+py(iat,iw)**2+pz(iat,iw)**2
        temp1=0.5d0*temp1/amt(iat,iw)
        ekin_mom=ekin_mom+temp1
       enddo
      enddo

      est_temp=2*ekin_mom/(dime*natom*nwalk-nshake-f-dime*conatom)
      est_temp_cumul=est_temp_cumul+est_temp


      if(modulo(it,nwrite).eq.0.and.inose.eq.1)then
       call calc_nhcham()
       write(2,'(F15.2,2F10.2,E20.10)')it*dt*autofs,est_temp*autok,est_temp_cumul*autok/it2,nhcham+ekin_mom+eclas
      endif

      if(modulo(it,nwrite).eq.0.and.inose.ne.1.and.ipimd.ne.2)then
       if(inose.eq.2)then
        write(2,'(F15.2,2F10.2,E20.10)')it*dt*autofs,est_temp*autok,est_temp_cumul*autok/it2,langham+ekin_mom+eclas
       else
        write(2,'(F15.2,2F10.2)')it*dt*autofs,est_temp*autok,est_temp_cumul*autok/it2
       endif
      endif

      if(ipimd.eq.0)then
       entot_cumul=entot_cumul+eclas+ekin_mom
       if(modulo(it,nwrite).eq.0)then
        write(1,'(F15.2,4E20.10)')it*dt*autofs,eclas,ekin_mom,eclas+ekin_mom,entot_cumul/it2
       endif
      endif

      if(modulo(it,nwrite).eq.0.and.ipimd.eq.2)then
       write(1,'(F15.2,3E20.10)')it*dt*autofs,eclas,ekin_mom,eclas+ekin_mom
      endif
      
      RETURN
      END


      real*8 function ekin_v (vx,vy,vz)
      use mod_array_size
      use mod_general
      use mod_system, ONLY:am
      implicit none
      real*8,intent(in)  :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8  :: temp1,ekin_mom 
      integer :: iw,iat

      ekin_mom=0.0d0

      do iw=1,nwalk
       do iat=1,natom
        temp1=vx(iat,iw)**2+vy(iat,iw)**2+vz(iat,iw)**2
        temp1=0.5d0*temp1*am(iat)
        ekin_mom=ekin_mom+temp1
       enddo
      enddo

      ekin_v=ekin_mom

      RETURN
      END
!
      real*8 function ekin_p (px,py,pz)
      use mod_array_size
      use mod_general
      use mod_system, ONLY:am
      implicit none
      real*8,intent(in)  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
      real*8  :: temp1,ekin_mom 
      integer :: iw,iat

      ekin_mom=0.0d0

      do iw=1,nwalk
       do iat=1,natom
        temp1=px(iat,iw)**2+py(iat,iw)**2+pz(iat,iw)**2
        temp1=0.5d0*temp1/am(iat)
        ekin_mom=ekin_mom+temp1
       enddo
      enddo

      ekin_p=ekin_mom

      RETURN
      END
