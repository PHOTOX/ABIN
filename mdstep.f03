
module mdstep
   use mod_system, ONLY: conatom, constrainP
   use mod_kinetic, ONLY: ekin_p
   implicit none
   private
   public :: verletstep,respastep,respashake,force_clas,force_quantum

   contains

   SUBROUTINE shiftX (rx,ry,rz,px,py,pz,mass,dt)
   real*8,intent(inout) :: rx(:,:),ry(:,:),rz(:,:)
   real*8,intent(in)   :: px(:,:),py(:,:),pz(:,:)
   real*8,intent(in)    :: mass(:,:)
   real*8,intent(in)    :: dt

      RX = RX + dt * PX/MASS
      RY = RY + dt * PY/MASS
      RZ = RZ + dt * PZ/MASS

   RETURN
   END

   SUBROUTINE shiftP (px,py,pz,fx,fy,fz,dt)
   real*8,intent(inout)  :: px(:,:),py(:,:),pz(:,:)
   real*8,intent(in)     :: fx(:,:),fy(:,:),fz(:,:)
   real*8,intent(in)     :: dt

      PX = PX + dt*FX
      PY = PY + dt*FY
      PZ = PZ + dt*FZ

   RETURN
   END

!----- Velocity verlet ALGORITHM                              Daniel Hollas,2012
!  At this moment it does not contain shake routines,
!  which are only in respashake function
!  GLE and NHC thermostats available at the moment
   subroutine verletstep(x,y,z,px,py,pz,amt,dt,eclas,fxc,fyc,fzc)
   use mod_nhc, ONLY:inose, imasst, shiftNHC_yosh, shiftNHC_yosh_mass
   use mod_gle, ONLY:langham, gle_step
!   use mod_interfaces, ONLY:force_clas
   real*8,intent(inout) :: x(:,:),y(:,:),z(:,:)
   real*8,intent(inout) :: fxc(:,:),fyc(:,:),fzc(:,:)
   real*8,intent(inout) :: px(:,:),py(:,:),pz(:,:)
   real*8,intent(in)    :: amt(:,:)
   real*8,intent(in)    :: dt
   real*8,intent(inout) :: eclas


   !---THERMOSTATS------------------
   if(inose.eq.1)then
   
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass(px,py,pz,amt,dt/2)
      else
         call shiftNHC_yosh(px,py,pz,amt,dt/2)
      endif
   
   endif
   
   !TODO zakazat inose=3 v init pro stabilni verzi
   if (inose.eq.3)  call wn_step(px,py,pz,amt)
   
   if (inose.eq.2)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   endif
   !--END OF THERMOSTATS
   
   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   if(conatom.gt.0) call constrainP (px,py,pz)
   
   call shiftX(x,y,z,px,py,pz,amt,dt)
   
   call force_clas(fxc,fyc,fzc,x,y,z,eclas)

   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   if(conatom.gt.0) call constrainP (px,py,pz)
   !---THERMOSTATS------------------
   
   if (inose.eq.3)  call wn_step(px,py,pz,amt)
   
   if (inose.eq.2)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   endif
   
   
   if(inose.eq.1)then
   
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass(px,py,pz,amt,dt/2)
      else
         call shiftNHC_yosh(px,py,pz,amt,dt/2)
      endif
   
   endif
   
   
   end
   
   !----- RESPA ALGORITHM                                              10.12.2012
   ! General algorithm of the propagation:
   ! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  
   ! further info is before subroutine respa_shake
   subroutine respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
              fxc,fyc,fzc,fxq,fyq,fzq)
   use mod_general, ONLY: istage
   use mod_nhc, ONLY:inose,imasst,shiftNHC_yosh,shiftNHC_yosh_mass
   use mod_gle, ONLY:langham, gle_step
!   use mod_interfaces, ONLY:force_quantum,force_clas
   real*8,intent(inout)  :: x(:,:),y(:,:),z(:,:)
   real*8,intent(inout)  :: fxc(:,:),fyc(:,:),fzc(:,:)
   real*8,intent(inout)  :: fxq(:,:),fyq(:,:),fzq(:,:)
   real*8,intent(inout)  :: px(:,:),py(:,:),pz(:,:)
   real*8,intent(in)     :: amg(:,:),amt(:,:)
   real*8,intent(in)     :: dt
   real*8,intent(inout)  :: eclas,equant
   integer               :: iabin
   
   if (inose.eq.2)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   end if
   
   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2*nabin))
      else
         call shiftNHC_yosh(px,py,pz,amt,dt/(2*nabin))
      end if
   end if
   
   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   if(conatom.gt.0) call constrainP (px,py,pz)
   
   
   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,-dt/(2*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,-dt/(2*nabin))
      end if
   end if
   
   
   do iabin=1,nabin
   
      if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif
   end if

   if (inose.eq.2.and.iabin.ne.1)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   end if

! if(istage.eq.2)then
!  call XtoU(transpx,transpy,transpz,px,py,pz)
! endif
!----quantum forces propagated using normal modes 
!---- GLE must use physical masses, i.e. cartesian coordinates

   call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))
   if(conatom.gt.0) call constrainP (px,py,pz)
 
   call shiftX(x,y,z,px,py,pz,amt,dt/nabin)

   call force_quantum(fxq,fyq,fzq,x,y,z,amg,equant)

   call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))
   if(conatom.gt.0) call constrainP (px,py,pz)

   if (inose.eq.2.and.iabin.ne.nabin)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   end if

   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif
   endif
   
   ! iabin loop
   enddo
   
   
   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,-dt/(2*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,-dt/(2*nabin))
      endif
   endif
   
   call force_clas(fxc,fyc,fzc,x,y,z,eclas)
   
   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   if(conatom.gt.0) call constrainP (px,py,pz)
   
   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif
   endif
   
   if (inose.eq.2)then

      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   
   endif
   
   
   end
   
   
!----- RESPA ALGORITHM WITH RATTLE                                            10.12.2012
! General algorithm of the propagation:
! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  

! Implementation of shake+RESPA according to www.chim.unifi.it/orac/MAN/node9.html

! some discussion concerning the Shake+nose-hoover may be found in: 
! p.199,M.Tuckermann, Statistical mechanics.. (2010)
! basically it says that NHC should not copromise constraints. 
! However, this doesn't apply to massive thermostatting that we use for PIMD!!!
! This is also obvious from the behaviour of pinyMD, which doesn't allow massive
! thermostatting with constraints.

! Therefore, we use different thermostating scheme here (subroutine shiftNHC_yosh).
! User must define "molecules", whose atoms share constraints. 
! Atoms which are not part of any constraint can be molecules themselves.
! We then append each "molecule" with its own thermostat.
! Definition of molecules are controlled by variables nmolt,natmolt nshakemol.
! Note that molecules must be in sequential order.
! For global thermostattting (i.e. dynamics is not disturbed that much), specify the system as one molecule. But don't do this with PIMD!!

! I.e for system of Chloride in 5 TIP3P waters, there will be 6
! molecules and 6 NHC chains. 


subroutine respashake(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
                 fxc,fyc,fzc,fxq,fyq,fzq)
      use mod_general, ONLY: istage
      use mod_nhc, ONLY:inose,shiftNHC_yosh,shiftNHC_yosh_mass
      use mod_system, ONLY:nshake
      use mod_fftw3
!      use mod_interfaces, ONLY:force_quantum,force_clas,shake,utox,xtou,qtox,xtoq
      real*8,intent(inout) :: x(:,:),y(:,:),z(:,:)
      real*8,intent(inout) :: px(:,:),py(:,:),pz(:,:)
      real*8,intent(in)    :: amg(:,:),amt(:,:)
      real*8,intent(in)    :: dt
      real*8,intent(inout) :: eclas,equant
      real*8,intent(inout) :: fxc(:,:),fyc(:,:),fzc(:,:)
      real*8,intent(inout) :: fxq(:,:),fyq(:,:),fzq(:,:)
      real*8 vx(:,:),vy(:,:),vz(:,:)
      real*8 transx(:,:),transy(:,:),transz(:,:)
      real*8 transxv(:,:),transyv(:,:),transzv(:,:)
      integer :: iabin,iq,iv


      if(inose.eq.1)then
        call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif

      call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
      if(conatom.gt.0) call constrainP (px,py,pz)

! RATTLE HERE!
      if(nshake.ge.1)then
      iq=0 ;iv=1
      !TODO: do this conversion inside the shake routine
      vx=px/amt
      vy=py/amt
      vz=pz/amt
      if(istage.eq.1)then
       call QtoX(x,y,z,transx,transy,transz)
       call QtoX(vx,vy,vz,transxv,transyv,transzv)
       call shake(transx,transy,transz,transxv,transyv,transzv,iq,iv) 
       call XtoQ(transxv,transyv,transzv,vx,vy,vz)
      else
       call shake(x,y,z,vx,vy,vz,iq,iv) 
      endif
      px=vx*amt
      py=vy*amt
      pz=vz*amt
      endif
!------END OF RATTLE     


      if(inose.eq.1)then
       call shiftNHC_yosh (px,py,pz,amt,-dt/(2*nabin))
      endif

!------RESPA LOOP
      do iabin=1,nabin
     
       if(inose.eq.1)then
        call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
       endif

       call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))
       if(conatom.gt.0) call constrainP (px,py,pz)

!--------RATTLE HERE!
   if(nshake.ge.1)then
   iq=0 ;iv=1
   vx=px/amt
   vy=py/amt
   vz=pz/amt
   if(istage.eq.1)then
    call QtoX(x,y,z,transx,transy,transz)
    call QtoX(vx,vy,vz,transxv,transyv,transzv)
    call shake(transx,transy,transz,transxv,transyv,transzv,iq,iv) 
    call XtoQ(transxv,transyv,transzv,vx,vy,vz)
!    upravujeme pouze rychlosti, tak asi nepotrebujeme transformovat zpatky
!    souradnice      
    call XtoQ(transx,transy,transz,x,y,z)
   else
    call shake(x,y,z,vx,vy,vz,iq,iv) 
   endif
    px=vx*amt
    py=vy*amt
    pz=vz*amt
   endif
!------END OF RATTLE     

!CONSTRAINING ATOMS
    if(conatom.gt.0) call constrainP (px,py,pz)

    call shiftX(x,y,z,px,py,pz,amt,dt/nabin)

!------SHAKE , iq=1, iv=0
   if(NShake.ge.1)then 
    iq=1
    iv=0
    if(istage.eq.1)then
     call QtoX(x,y,z,transx,transy,transz)
     call shake(transx,transy,transz,transxv,transyv,transzv,iq,iv) 
     call XtoQ(transx,transy,transz,x,y,z)
    else
     call shake(x,y,z,vx,vy,vz,iq,iv) 
    endif
   endif
!------END OF SHAKE     

    call force_quantum(fxq,fyq,fzq,x,y,z,amg,equant)

    call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))
    if(conatom.gt.0) call constrainP (px,py,pz)

   if(inose.eq.1)then
     call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
   endif

   enddo
!-----END OF RESPA LOOP

   if(inose.eq.1)then
    call shiftNHC_yosh (px,py,pz,amt,-dt/(2*nabin))
   endif

   call force_clas(fxc,fyc,fzc,x,y,z,eclas)

   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   if(conatom.gt.0) call constrainP (px,py,pz)

   if(inose.eq.1)then
     call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
   endif


   end

end module mdstep
