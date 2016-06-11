
module mod_mdstep
   use mod_const,    only: DP
   use mod_kinetic,  only: ekin_p
   use mod_transform
   implicit none
   private
   public :: verletstep, respastep, respashake

   contains

   subroutine shiftX (rx, ry, rz, px, py, pz, mass, dt)
   real(DP),intent(inout) :: rx(:,:), ry(:,:), rz(:,:)
   real(DP),intent(in)    :: px(:,:), py(:,:), pz(:,:)
   real(DP),intent(in)    :: mass(:,:)
   real(DP),intent(in)    :: dt

   RX = RX + dt * PX / MASS
   RY = RY + dt * PY / MASS
   RZ = RZ + dt * PZ / MASS

   end subroutine shiftX


   subroutine shiftP (px, py, pz, fx, fy, fz, dt)
   use mod_general, only: nwalk, natom
   use mod_system,  only: conatom, constrainP
   real(DP),intent(inout)  :: px(:,:), py(:,:), pz(:,:)
   real(DP),intent(in)     :: fx(:,:), fy(:,:), fz(:,:)
   real(DP),intent(in)     :: dt

   PX = PX + dt * FX
   PY = PY + dt * FY
   PZ = PZ + dt * FZ

   if(conatom.gt.0) call constrainP (px, py, pz)

   end subroutine shiftP


   subroutine thermostat(px,py,pz,amt,dt)
   use mod_nhc
   use mod_gle
   real(DP),intent(inout)  :: px(:,:),py(:,:),pz(:,:)
   real(DP),intent(in)     :: amt(:,:)
   real(DP),intent(in)     :: dt

   if(inose.eq.1)then
   
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass(px,py,pz,amt,dt)
      else
         call shiftNHC_yosh(px,py,pz,amt,dt)
      endif
   
   endif
   
   if (inose.eq.3)  call wn_step(px,py,pz,amt)
   
   if (inose.eq.2)  call gle_step(px,py,pz,amt)

   end subroutine thermostat
           

!--VELOCITY VERLET INTEGRATOR WITH THERMOSTATS
!  At this moment it does not contain shake routines,
!  which are only in respashake function.
!  Works with NHC, GLE and PIGLET thermostats.
!  Contains propagation of normal modes according to:
!  eq 23 from J. Chem. Phys. 133, 124104 2010
   subroutine verletstep(x,y,z,px,py,pz,amt,dt,eclas,fxc,fyc,fzc)
   use mod_general, ONLY: pot, ipimd, inormalmodes
   use mod_nhc, ONLY:inose
   use mod_interfaces, only: force_clas
   real(DP),intent(inout) :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) :: fxc(:,:),fyc(:,:),fzc(:,:)
   real(DP),intent(inout) :: px(:,:),py(:,:),pz(:,:)
   real(DP),intent(in)    :: amt(:,:)
   real(DP),intent(in)    :: dt
   real(DP),intent(inout) :: eclas


   if(inose.gt.0) call thermostat(px, py, pz, amt, dt/2)
   
   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   
   if(inormalmodes.eq.1)then
      ! Warning, initial hack, passing amt here
      call propagate_nm(x, y, z, px, py, pz, amt, dt)
   else
      call shiftX(x, y, z, px, py, pz, amt, dt)
   end if
   
   call force_clas(fxc,fyc,fzc,x,y,z,eclas)

   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)

   if(inose.gt.0) call thermostat(px, py, pz, amt, dt/2)
   
   end subroutine verletstep

   
   !----- RESPA ALGORITHM                                              10.12.2012
   ! General algorithm of the propagation:
   ! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  
   ! further info is before subroutine respa_shake
   subroutine respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
              fxc,fyc,fzc,fxq,fyq,fzq)
   use mod_general,  only: nabin
   use mod_nhc,      only: inose,imasst,shiftNHC_yosh,shiftNHC_yosh_mass
   use mod_gle,      only: gle_step
   use mod_interfaces, only: force_clas, force_quantum
   real(DP),intent(inout)  :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)  :: fxc(:,:),fyc(:,:),fzc(:,:)
   real(DP),intent(inout)  :: fxq(:,:),fyq(:,:),fzq(:,:)
   real(DP),intent(inout)  :: px(:,:),py(:,:),pz(:,:)
   real(DP),intent(in)     :: amg(:,:),amt(:,:)
   real(DP),intent(in)     :: dt
   real(DP),intent(inout)  :: eclas,equant
   integer                 :: iabin
   
   call thermostat(px, py, pz, amt, dt/(2*nabin))
   
   call shiftP (px, py, pz, fxc, fyc, fzc, dt/2)
   
   do iabin = 1, nabin
   
      if (inose.gt.0.and.iabin.ne.1) call thermostat(px, py, pz, amt, dt/(2*nabin))

      call shiftP(px, py, pz, fxq, fyq, fzq, dt/(2*nabin))
  
      call shiftX(x, y, z, px, py, pz, amt, dt/nabin)
 
      call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)
 
      call shiftP(px, py, pz, fxq, fyq, fzq, dt/(2*nabin))
 
      if (inose.gt.0.and.iabin.ne.nabin) call thermostat(px, py, pz, amt, dt/(2*nabin))
   
   enddo
   
   call force_clas(fxc,fyc,fzc,x,y,z,eclas)
   
   call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)
   
   call thermostat(px, py, pz, amt, dt/(2*nabin))
   
   end subroutine respastep
   
   
!----- RESPA ALGORITHM WITH RATTLE                                            10.12.2012
! General algorithm of the propagation:
! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  

! Implementation of shake+RESPA according to www.chim.unifi.it/orac/MAN/node9.html

! some discussion concerning the Shake+nose-hoover may be found in: 
! p.199,M.Tuckermann, Statistical mechanics.. (2010)
! basically it says that NHC should not copromise constraints. 
! However, this doesn't apply to massive thermostatting that we use for PIMD!!!

! Therefore, we use different thermostating scheme here (subroutine shiftNHC_yosh).
! User must define "molecules", whose atoms share constraints. 
! Atoms which are not part of any constraint can be molecules themselves.
! We then append each "molecule" with its own thermostat.
! Definition of molecules are controlled by variables nmolt,natmolt nshakemol.
! Note that molecules must be in sequential order.
! For global thermostattting (i.e. dynamics is not disturbed that much), specify the system as one molecule.
! But don't do this with PIMD!!

! I.e for system of Chloride in 5 TIP3P waters, there will be 6
! molecules and 6 NHC chains. 


subroutine respashake(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
                 fxc,fyc,fzc,fxq,fyq,fzq)
      use mod_general, ONLY: nabin
      use mod_nhc, ONLY:inose,shiftNHC_yosh,shiftNHC_yosh_mass
      use mod_shake, only: shake, nshake
      use mod_interfaces, only: force_clas, force_quantum
      real(DP),intent(inout) :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(inout) :: px(:,:),py(:,:),pz(:,:)
      real(DP),intent(in)    :: amg(:,:),amt(:,:)
      real(DP),intent(in)    :: dt
      real(DP),intent(inout) :: eclas,equant
      real(DP),intent(inout) :: fxc(:,:),fyc(:,:),fzc(:,:)
      real(DP),intent(inout) :: fxq(:,:),fyq(:,:),fzq(:,:)
      integer  :: iabin


      if(inose.eq.1)then
        call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif

      call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)

! RATTLE HERE!
      if(nshake.ge.1) call shake(x,y,z,px,py,pz,amt,0,1) 


      if(inose.eq.1)then
       call shiftNHC_yosh (px,py,pz,amt,-dt/(2*nabin))
      endif

!------RESPA LOOP
      do iabin=1,nabin
     
       if(inose.eq.1)then
        call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
       endif

       call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))

!--------RATTLE HERE!
   if(nshake.ge.1) call shake(x,y,z,px,py,pz,amt,0,1) 

    call shiftX(x,y,z,px,py,pz,amt,dt/nabin)

!------SHAKE , iq=1, iv=0
   if(NShake.gt.0) call shake(x,y,z,px,py,pz,amt,1,0) 

    call force_quantum(fxq,fyq,fzq,x,y,z,amg,equant)

    call shiftP (px,py,pz,fxq,fyq,fzq,dt/(2*nabin))

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

   if(inose.eq.1)then
     call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
   endif


   end subroutine respashake

end module mod_mdstep
