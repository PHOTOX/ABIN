!---                    by Daniel Hollas,15.3.2012
!--contains spherical potential SBC and elastig bag(TO DO)
MODULE mod_sbc
   use mod_const, only: DP, AUtoM, ME, PI
   implicit none
   private 
   public    :: isbc, rb_sbc, kb_sbc, rho
   public    :: force_sbc, sbc_init
   integer   :: isbc=0  
   real(DP)  :: rb_sbc=-1, kb_sbc=0.02d0, rho=-1
   real(DP),parameter  :: fact=autom**3/me
   real(DP)            :: mass_total=0.0d0
   integer, parameter  :: ibag=0  !elastic bag
   save
   CONTAINS

   SUBROUTINE sbc_init(x,y,z)
      use mod_const,    only: ANG
      use mod_general,  only: natom !,nwalk
      use mod_system,   only: am, names
      real(DP),intent(in) :: x(:,:), y(:,:), z(:,:)
      real(DP) :: r, rmax, xcm, ycm, zcm
      integer  :: iw, iat

      mass_total=0.0d0

      do iat=1,natom
         mass_total=mass_total+am(iat)
      enddo

      rmax=0.0d0

      rb_sbc=rb_sbc*ANG

!---calculation of Center-Of-Mass
      iw=1
      call calc_com(x, y, z, iw ,xcm, ycm, zcm)

      write(*,*)'Center of mass is at [A]:',xcm/ang,ycm/ang,zcm/ang

!-----calculation of size-of-cluster
      do iat=1,natom
       r=(x(iat,iw)-xcm)**2+(y(iat,iw)-ycm)**2+(z(iat,iw)-zcm)**2
       r=sqrt(r)
       if(r.gt.rmax.and.names(iat).ne.'H')then
        rmax=r
       endif
      enddo

!---SETTING SIZE OF THE CLUSTER      
      write(*,*)'Cluster radius is r= ',rmax/ang,'  Angstrom'
      if (rb_sbc.le.0.and.rho.le.0)then
         write(*,*)'Cluster radius for spherical boundary conditions not specified.'
         if(rho.lt.0)then
            write(*,*)'Will be set equal the initial radius of the system.'
            rb_sbc=rmax
         end if
      endif
!---DETERMINATION of cluster size from given density
      if (rho.gt.0)then
         write(*,*)'Calculating cluster radius from given densty.'
         rho=rho*fact !conversion from g/L to atomic units
         rb_sbc=mass_total/rho*3/4/PI
         rb_sbc=rb_sbc**(1/3.)
      endif

      if (rmax.gt.rb_sbc)then
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
         write(*,*)'WARNING: Cluster radius is bigger than rb_sbc.'
         write(*,*)'Is this really what you want?'
         write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
      endif

      write(*,*)'rb_sbc[A]=',rb_sbc/ang

      END SUBROUTINE
!--------------------------------
      SUBROUTINE force_sbc(x,y,z,fx,fy,fz)
      use mod_general,only: natom, nwalk, it, nwrite
      use mod_system, ONLY: names 
      real(DP),intent(in)    :: x(:,:), y(:,:), z(:,:)
      real(DP),intent(inout) :: fx(:,:), fy(:,:), fz(:,:)
      real(DP)  :: r, frb, rmax, xcm, ycm, zcm
      integer   :: iat, iw

      rmax=0.0d0

      iw=1
!---calculation of Center-Of-Mass
      call calc_com(x, y, z, iw ,xcm, ycm, zcm)

!      write(*,*)'COM: ',xcm,ycm,zcm

!---kb_sbc and rb_sbc must be specified in input.in
      do iw=1,nwalk
       do iat=1,natom
        r=(x(iat,iw)-xcm)**2+(y(iat,iw)-ycm)**2+(z(iat,iw)-zcm)**2
        r=sqrt(r)
        if(r.gt.rb_sbc.and.names(iat).ne.'H')then
         frb=-kb_sbc*(r-rb_sbc)
         fx(iat,iw)=fx(iat,iw)+frb*(x(iat,iw)-xcm)/(r)
         fy(iat,iw)=fy(iat,iw)+frb*(y(iat,iw)-ycm)/(r)
         fz(iat,iw)=fz(iat,iw)+frb*(z(iat,iw)-zcm)/(r)
        endif

        if(r.gt.rmax.and.names(iat).ne.'H')then
         rmax=r
        endif

       enddo
      enddo

!      if(idebug.eq.1)then
       if(modulo(it,nwrite).eq.0) write(11,*)it,rmax,mass_total/(4.0d0/3.0d0*pi*rmax**3)/fact
!      endif

      return
   end subroutine

   subroutine calc_com(x, y, z, iw ,xcm, ycm, zcm)
   use mod_general, only: natom
   use mod_system, ONLY: am
   real(DP),intent(out) :: xcm, ycm, zcm
   real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
   integer, intent(in)  :: iw
   integer  :: iat

   xcm=0.0d0
   ycm=0.0d0
   zcm=0.0d0

   do iat=1,natom
      xcm=xcm+x(iat,iw)*am(iat)
      ycm=ycm+y(iat,iw)*am(iat)
      zcm=zcm+z(iat,iw)*am(iat)
   enddo
   xcm=xcm/mass_total
   ycm=ycm/mass_total
   zcm=zcm/mass_total

   end subroutine calc_com

END MODULE mod_sbc


!      MODULE mod_bag
!      real(DP),allocatable :: neigh1(:),neigh2(:),neigh3(:)
!      real(DP) :: some_parameters
!      save
 !     CONTAINS
 !
 !     SUBROUTINE bag_init()
 !     !allocate(neigh1,neig)
 !     END SUBROUTINE
!
!      SUBROUTINE force_bag(x,y,z,fx,fy,fz)
!      use mod_general
!      use mod_system, ONLY:am,names 
!      implicit none
!      real(DP) x(:,:),y(:,:),z(:,:)
!      real(DP) fx(:,:),fy(:,:),fz(:,:)
!      real(DP)  :: r,frb
!      integer :: iat,iw
!
!
!
!
!
!      END SUBROUTINE
!      END MODULE

