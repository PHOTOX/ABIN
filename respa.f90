!----- RESPA ALGORITHM                                              10.12.2012
! General algorithm of the propagation:
! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)  

subroutine respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
           fxc,fyc,fzc,fxq,fyq,fzq)
use mod_array_size
use mod_general
use mod_nhc, ONLY:inose,imasst
use mod_gle
implicit none
real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8,intent(inout)  :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
real*8,intent(inout)  :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
real*8,intent(inout)  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
real*8,intent(in)     :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
real*8                :: dt,eclas,equant,ekin_p
integer               :: iabin,iat,iw

if (inose.eq.2)then
   langham=langham+ekin_p(px,py,pz)
   call gle_step(px,py,pz,amt)
   langham=langham-ekin_p(px,py,pz)
end if

if(inose.eq.1)then
   if (imasst.eq.1)then
      call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2.0*nabin))
   else
      call shiftNHC_yosh(px,py,pz,amt,dt/(2.0*nabin))
   end if
end if

call shiftP (px,py,pz,fxc,fyc,fzc,amt,dt/2.0d0)


if(inose.eq.1)then
   if (imasst.eq.1)then
      call shiftNHC_yosh_mass (px,py,pz,amt,-dt/(2.0*nabin))
   else
      call shiftNHC_yosh (px,py,pz,amt,-dt/(2.0*nabin))
   end if
end if


do iabin=1,nabin

   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2.0*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,dt/(2.0*nabin))
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

   call shiftP (px,py,pz,fxq,fyq,fzq,amt,dt/(2.0d0*nabin))


   if(conatom.gt.0)then
      do iw=1,nwalk
         do iat=1,conatom
            px(iat,iw)=0.0
            py(iat,iw)=0.0
            pz(iat,iw)=0.0
         enddo
      enddo
   endif
 
   call shiftX(x,y,z,px,py,pz,amt,dt/nabin)

   call force_quantum(fxq,fyq,fzq,x,y,z,amg,equant)

   call shiftP (px,py,pz,fxq,fyq,fzq,amt,dt/(2.0*nabin))

   if (inose.eq.2.and.iabin.ne.nabin)then
      langham=langham+ekin_p(px,py,pz)
      call gle_step(px,py,pz,amt)
      langham=langham-ekin_p(px,py,pz)
   end if

   if(inose.eq.1)then
      if (imasst.eq.1)then
         call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2.0*nabin))
      else
         call shiftNHC_yosh (px,py,pz,amt,dt/(2.0*nabin))
      endif
   endif

! iabin loop
enddo


if(inose.eq.1)then
   if (imasst.eq.1)then
      call shiftNHC_yosh_mass (px,py,pz,amt,-dt/(2.0*nabin))
   else
      call shiftNHC_yosh (px,py,pz,amt,-dt/(2.0*nabin))
   endif
endif

call force_clas(fxc,fyc,fzc,x,y,z,eclas)

call shiftP (px,py,pz,fxc,fyc,fzc,amt,dt/2.0d0)

if(inose.eq.1)then
   if (imasst.eq.1)then
      call shiftNHC_yosh_mass (px,py,pz,amt,dt/(2.0*nabin))
   else
      call shiftNHC_yosh (px,py,pz,amt,dt/(2.0*nabin))
   endif
endif

if (inose.eq.2)then

   langham=langham+ekin_p(px,py,pz)
   call gle_step(px,py,pz,amt)
   langham=langham-ekin_p(px,py,pz)

endif


end


