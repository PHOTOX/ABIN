! initial version                         P. Slavicek a kol., Mar 25, 2009
!
!--------------------------------------------------------------------------
module mod_kinetic
   use mod_const, ONLY: DP
   implicit none
   private
   public :: temperature, ekin_v, ekin_p, est_temp_cumul, entot_cumul
   real(DP) :: est_temp_cumul=0.0d0,entot_cumul=0.0d0
   save
   contains
   subroutine temperature(px,py,pz,amt,dt,eclas)
      use mod_const, ONLY:AUtoFS,AUtoK
      use mod_general
      use mod_system, ONLY: dime, f, conatom
      use mod_nhc, ONLY: inose,nhcham,calc_nhcham
      use mod_gle, ONLY:langham
      use mod_shake, only: nshake
      implicit none
      integer :: iw,iat
      real(DP),intent(in)  :: px(:,:),py(:,:),pz(:,:)
      real(DP),intent(in)  :: amt(:,:)
      real(DP),intent(in)  :: dt,eclas
      real(DP)  :: est_temp,temp1,ekin_mom
      real(DP)  :: it2
      
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
      end subroutine temperature


      real(DP) function ekin_v (vx,vy,vz)
      use mod_general, ONLY: nwalk, natom
      use mod_system, ONLY:am
      implicit none
      real(DP),intent(in)  :: vx(:,:),vy(:,:),vz(:,:)
      real(DP)  :: temp1,ekin_mom 
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
      END function ekin_v
!
      real(DP) function ekin_p (px,py,pz)
      use mod_general, ONLY: nwalk, natom
      use mod_system, ONLY:am
      implicit none
      real(DP),intent(in)  :: px(:,:),py(:,:),pz(:,:)
      real(DP)  :: temp1,ekin_mom 
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
      END function ekin_p

end module mod_kinetic
