
module mod_forces
   use mod_const, only: DP
   implicit none
   private 
   public :: force_clas, force_quantum
   contains

   subroutine force_clas(fx,fy,fz,x,y,z,energy)
      use mod_general
      use mod_qmmm,     only: force_LJCoul
      use mod_nab,      only: ipbc,wrap,nsnb,force_nab
      use mod_sbc,      only: force_sbc, isbc !,ibag
      use mod_system,   only: conatom
      use mod_nhc,      only: inose
      use mod_harmon,   only: force_harmon,force_2dho,force_morse
      use mod_guillot,  only: force_guillot
      use mod_utils,    only: printf
      use mod_transform
      use mod_interfaces, only: force_abin
      real(DP),intent(inout) ::  x(:,:),y(:,:),z(:,:)
      real(DP),intent(inout) ::  fx(:,:),fy(:,:),fz(:,:)
      real(DP),intent(out)   ::  energy
      real(DP)  :: transx(size(x,1),size(x,2))
      real(DP)  :: transy(size(x,1),size(x,2))
      real(DP)  :: transz(size(x,1),size(x,2))
      real(DP)  :: fxab(size(x,1),size(x,2))
      real(DP)  :: fyab(size(x,1),size(x,2))
      real(DP)  :: fzab(size(x,1),size(x,2))
      integer   :: iat,iw
      real (DP) :: eclas

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
!      if (ibag.eq.1) call force_bag(transx,transy,transz,fxab,fyab,fzab)

!---------QMMM SECTION-----------------
!- ONIOM method (iqmmm=1) is called in force_abin
     if(iqmmm.eq.2) call force_nab(transx,transy,transz,fxab,fyab,fzab,eclas)
     if(iqmmm.eq.3) call force_LJCoul(transx,transy,transz,fxab,fyab,fzab,eclas)

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
         fx(iat,iw)=0.0d0
         fy(iat,iw)=0.0d0
         fz(iat,iw)=0.0d0
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

   end subroutine force_clas

   subroutine force_quantum(fx,fy,fz,x,y,z,amg,energy)
      use mod_array_size
      use mod_general
      use mod_nhc, ONLY: temp,inose
      implicit none
      real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(in) :: amg(:,:)
      real(DP),intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
      real(DP),intent(out) :: energy
      real(DP) :: ak(size(x,1),size(x,2))
      real(DP) :: equant
      integer  :: iat,iw,i,j,kplus,kminus

      do iat=1,natom
       do iw=1,nwalk
        fx(iat,iw)=0.0d0
        fy(iat,iw)=0.0d0
        fz(iat,iw)=0.0d0
       enddo
      enddo

! TODO: ak parametry se nemuseji pocitat pokazde znova
! Setting the quantum force constants
! ak is defined is m*P/(beta^2*hbar^2)
        do iw=1,nwalk
         do i=1,natom
          ak(i,iw)=nwalk*amg(i,iw)*TEMP**2
         enddo
        enddo
! for PI+GLE we have different hamiltonian
       if(inose.eq.2)then
        do iw=1,nwalk
         do i=1,natom
          ak(i,iw)=nwalk*ak(i,iw)
         enddo
        enddo
       endif

      equant=0.0d0


! If the staging transformation is not used
      if(istage.eq.0)then
       do j=1,natom 
        do i=1,nwalk 
         kplus=i+1
         kminus=i-1
         if(i.eq.1)then
          kminus=nwalk
         endif
         if(i.eq.nwalk)then
          kplus=1
         endif
         fx(j,i)=(x(j,i)-x(j,kplus))
         fx(j,i)=fx(j,i)+(x(j,i)-x(j,kminus))
         fx(j,i)=-fx(j,i)*ak(j,i)
         fy(j,i)=(y(j,i)-y(j,kplus))
         fy(j,i)=fy(j,i)+(y(j,i)-y(j,kminus))
         fy(j,i)=-fy(j,i)*ak(j,i)
         fz(j,i)=(z(j,i)-z(j,kplus))
         fz(j,i)=fz(j,i)+(z(j,i)-z(j,kminus))
         fz(j,i)=-fz(j,i)*ak(j,i)
         equant=equant+0.5*ak(j,i)*(x(j,i)-x(j,kplus))**2
         equant=equant+0.5*ak(j,i)*(y(j,i)-y(j,kplus))**2
         equant=equant+0.5*ak(j,i)*(z(j,i)-z(j,kplus))**2
        enddo
       enddo
      endif

! If the staging transformation is used
       if(istage.eq.1.or.istage.eq.2)then
        do j=1,natom 
         do i=1,nwalk 
          fx(j,i)=-ak(j,i)*x(j,i)   
          fy(j,i)=-ak(j,i)*y(j,i)       
          fz(j,i)=-ak(j,i)*z(j,i) 
          equant=equant+0.5d0*ak(j,i)*x(j,i)**2
          equant=equant+0.5d0*ak(j,i)*y(j,i)**2
          equant=equant+0.5d0*ak(j,i)*z(j,i)**2
         enddo
        enddo
       endif 

      energy=equant
      
   end subroutine force_quantum
                                        
end module mod_forces
