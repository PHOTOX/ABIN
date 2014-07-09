!----- RESPA ALGORITHM                                              10.12.2012
! General algorithm of the propagation:
! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  

subroutine respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
           fxc,fyc,fzc,fxq,fyq,fzq)
use mod_array_size
use mod_general
use mod_nhc, ONLY:inose,imasst,shiftNHC_yosh,shiftNHC_yosh_mass
use mod_gle
use mod_interfaces, ONLY:shiftP,shiftX,force_quantum,force_clas
use mod_kinetic, ONLY:ekin_p
use mod_system, ONLY:constrainP, conatom
implicit none
real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8,intent(inout)  :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
real*8,intent(inout)  :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
real*8,intent(inout)  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
real*8,intent(in)     :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
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


