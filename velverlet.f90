!----- Velocity verlet ALGORITHM                              Daniel Hollas,2012
!at this moment does not contain shake routines, which are only in respashake
!function

! GLE and NHC thermostats available at the moment

subroutine verletstep(x,y,z,px,py,pz,amt,dt,eclas,fxc,fyc,fzc)
use mod_array_size
use mod_general
use mod_system, ONLY: constrainP
use mod_nhc, ONLY:inose,imasst,shiftNHC_yosh,shiftNHC_yosh_mass
use mod_gle
use mod_interfaces, ONLY:shiftP,shiftX,force_clas,ekin_p
implicit none
real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8,intent(inout) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
real*8,intent(inout) :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
real*8,intent(in)    :: amt(npartmax,nwalkmax)
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

