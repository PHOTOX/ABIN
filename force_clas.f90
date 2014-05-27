      subroutine force_clas(fx,fy,fz,x,y,z,energy)
      use mod_array_size
      use mod_general
      use mod_qmmm, ONLY: qmmmtype
      use mod_nab, ONLY: ipbc,wrap,nsnb
      use mod_sbc
      use mod_bag
      use mod_nhc, ONLY: inose
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8 transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
      real*8 fxab(npartmax,nwalkmax),fyab(npartmax,nwalkmax),fzab(npartmax,nwalkmax)
      real*8 sumux(npartmax),sumuy(npartmax),sumuz(npartmax)
      integer :: iat,iw
      real *8 :: energy,eclas

! Initialization
      do iw=1,nwalk
       do iat=1,natom
       fx(iat,iw)=0.0d0
       fy(iat,iw)=0.0d0
       fz(iat,iw)=0.0d0
       fxab(iat,iw)=0.0d0
       fyab(iat,iw)=0.0d0
       fzab(iat,iw)=0.0d0
       enddo
      enddo

      eclas=0.0d0


! Back stage transformation,cartesian coordinates are kept in trans
! matrices(even if staging is OFF!)
      if(istage.eq.1)then
      call QtoX(x,y,z,transx,transy,transz)
      endif
      if(istage.eq.2)then
      call UtoX(x,y,z,transx,transy,transz)
      endif
      if(istage.eq.0)then
        do iat=1,natom
         do iw=1,nwalk
          transx(iat,iw)=x(iat,iw)
          transy(iat,iw)=y(iat,iw)
          transz(iat,iw)=z(iat,iw)
         enddo
        enddo 
      endif
! end of the back stage transformation

! wraping molecules back to the box 
      if (pot.eq.'nab'.and.ipbc.eq.1.and.modulo(it,nsnb).eq.0) call wrap(transx,transy,transz)

!--Here we decide which forces we want. By default we call external program by
! force_abin routine
     SELECT CASE (pot)
        case ("mm")
          call force_LJCoul(transx,transy,transz,fxab,fyab,fzab,eclas)
        case ("harm")
          call force_harmon(transx,transy,transz,fxab,fyab,fzab,eclas)
        case ("2dho")
          call force_2dho(transx,transy,transz,fxab,fyab,fzab,eclas)
        case ("morse")
          call force_morse(transx,transy,transz,fxab,fyab,fzab,eclas)
        case ("guillot")
          call force_guillot(transx,transy,transz,fxab,fyab,fzab,eclas)
        case ("nab")
          call force_nab(transx,transy,transz,fxab,fyab,fzab,eclas)
        case DEFAULT
          call force_abin(transx,transy,transz,fxab,fyab,fzab,eclas)
          eclas=eclas/nwalk
     END SELECT

      if (isbc.eq.1) call force_sbc(transx,transy,transz,fxab,fyab,fzab)
      if (ibag.eq.1) call force_bag(transx,transy,transz,fxab,fyab,fzab)

!---------QMMM SECTION-----------------
     if(iqmmm.eq.1)then
      if(qmmmtype.eq.'nab') call force_nab(transx,transy,transz,fxab,fyab,fzab,eclas)
      if(qmmmtype.eq.'abin') call force_LJCoul(transx,transy,transz,fxab,fyab,fzab,eclas)
     endif

!--------------------------------------

! Classical force without the stage transformation
!  
      if(istage.eq.0)then
       if(inose.eq.2)then  ! pro kvantovy termostat je jiny hamiltonian
        do iw=1,nwalk 
         do iat=1,natom
          fx(iat,iw)=fx(iat,iw)+fxab(iat,iw)
          fy(iat,iw)=fy(iat,iw)+fyab(iat,iw)
          fz(iat,iw)=fz(iat,iw)+fzab(iat,iw)
!          eclas=eclas*nwalk  eclas stejne neovlivnuje dynamiku.....
         enddo 
        enddo
       else
        do iw=1,nwalk 
         do iat=1,natom
          fx(iat,iw)=fx(iat,iw)+fxab(iat,iw)/nwalk
          fy(iat,iw)=fy(iat,iw)+fyab(iat,iw)/nwalk 
          fz(iat,iw)=fz(iat,iw)+fzab(iat,iw)/nwalk 
         enddo 
        enddo
       endif
      endif

!----TRANSFORMING FORCES FROM CARTESIAN TO STAGING or NORMAL MODE COORDS--!
      if(istage.eq.1)then
!----forces are divided by nwalk inside the FXtoFQ routine!---------------!
       call FXtoFQ(fxab,fyab,fzab,fx,fy,fz)
      endif
      if(istage.eq.2)then
       if (idebug.eq.1) call printf(fxab,fyab,fzab)
       call XtoU(fxab,fyab,fzab,fx,fy,fz)
!       if (idebug.eq.1) call printf(fx,fy,fz)
      endif
    
      if(conatom.gt.0)then
       do iw=1,nwalk
        do iat=1,conatom
         fx(iat,iw)=0.0
         fy(iat,iw)=0.0
         fz(iat,iw)=0.0
        enddo
       enddo
      endif
      energy=eclas

!-----for PBC we do wrapping of molecules back to the box    
      if(pot.eq.'nab'.and.ipbc.eq.1)then

! Stage transformation,
      if(istage.eq.1)then
       call XtoQ(transx,transy,transz,x,y,z)
      endif
      if(istage.eq.2)then
       call XtoU(transx,transy,transz,x,y,z)
      endif
! end of the stage transformation
      if(istage.eq.0)then
        do iat=1,natom
         do iw=1,nwalk
          x(iat,iw)=transx(iat,iw)
          y(iat,iw)=transy(iat,iw)
          z(iat,iw)=transz(iat,iw)
         enddo
        enddo 
      endif

      endif

      end
                                        
