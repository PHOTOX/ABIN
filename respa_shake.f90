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
      use mod_array_size
      use mod_general
      use mod_nhc, ONLY:inose,shiftNHC_yosh,shiftNHC_yosh_mass
      use mod_system, ONLY:nshake
      use mod_fftw3
      use mod_interfaces, ONLY:shiftP,shiftX,force_quantum,force_clas,shake,utox,xtou,qtox,xtoq
      implicit none
      real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(inout) :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
      real*8,intent(in)    :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
      real*8,intent(in)    :: dt
      real*8,intent(inout) :: eclas,equant
      real*8,intent(inout) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
      real*8,intent(inout) :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
      real*8 transxv(npartmax,nwalkmax),transyv(npartmax,nwalkmax),transzv(npartmax,nwalkmax)
      integer :: iabin,iq,iv,iat,iw


      if(inose.eq.1)then
        call shiftNHC_yosh (px,py,pz,amt,dt/(2*nabin))
      endif

       call shiftP (px,py,pz,fxc,fyc,fzc,dt/2)

! RATTLE HERE!
      if(nshake.ge.1)then
      iq=0
      iv=1
      do iat=1,natom
       do iw=1,nwalk
        vx(iat,iw)=px(iat,iw)/amt(iat,iw)
        vy(iat,iw)=py(iat,iw)/amt(iat,iw)
        vz(iat,iw)=pz(iat,iw)/amt(iat,iw)
       enddo
      enddo
      if(istage.eq.1)then
       call QtoX(x,y,z,transx,transy,transz)
       call QtoX(vx,vy,vz,transxv,transyv,transzv)
       call shake(transx,transy,transz,transxv,transyv,transzv,iq,iv) 
       call XtoQ(transxv,transyv,transzv,vx,vy,vz)
      else
       call shake(x,y,z,vx,vy,vz,iq,iv) 
      endif
      do iat=1,natom
       do iw=1,nwalk
       px(iat,iw)=vx(iat,iw)*amt(iat,iw)
       py(iat,iw)=vy(iat,iw)*amt(iat,iw)
       pz(iat,iw)=vz(iat,iw)*amt(iat,iw)
       enddo
      enddo
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

!--------RATTLE HERE!
      if(nshake.ge.1)then
      iq=0
      iv=1
      do iat=1,natom
       do iw=1,nwalk
        vx(iat,iw)=px(iat,iw)/amt(iat,iw)
        vy(iat,iw)=py(iat,iw)/amt(iat,iw)
        vz(iat,iw)=pz(iat,iw)/amt(iat,iw)
       enddo
      enddo
      if(istage.eq.1)then
       call QtoX(x,y,z,transx,transy,transz)
       call QtoX(vx,vy,vz,transxv,transyv,transzv)
       call shake(transx,transy,transz,transxv,transyv,transzv,iq,iv) 
       call XtoQ(transxv,transyv,transzv,vx,vy,vz)
!      upravujeme pouze rychlosti, tak asi nepotrebujeme transformovat zpatky
!      souradnice      
!      call XtoQ(transx,transy,transz,x,y,z)
      else
       call shake(x,y,z,vx,vy,vz,iq,iv) 
      endif
      do iat=1,natom
       do iw=1,nwalk
       px(iat,iw)=vx(iat,iw)*amt(iat,iw)
       py(iat,iw)=vy(iat,iw)*amt(iat,iw)
       pz(iat,iw)=vz(iat,iw)*amt(iat,iw)
       enddo
      enddo
      endif
!------END OF RATTLE     

!CONSTRAINING ATOMS
       if(conatom.gt.0)then
        do iw=1,nwalk
         do iat=1,conatom
          px(iat,iw)=0.0d0
          py(iat,iw)=0.0d0
          pz(iat,iw)=0.0d0
         enddo
        enddo
       endif
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


      end


