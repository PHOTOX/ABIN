!---                    by Daniel Hollas,15.3.2012
!--- contains spherical potential SBC and elastig bag(TO DO)
      MODULE mod_sbc
      use mod_array_size
      real*8  :: xcm=0.0d0,ycm=0.0d0,zcm=0.0d0,mass_total=0.0d0
      real*8  :: rmax=0.0d0
      real*8  :: rb_sbc=-1,kb_sbc=0.02d0,rho=-1
      real*8  :: fact=autom**3/me
      save
      CONTAINS
      SUBROUTINE sbc_init(x,y,z)
      use mod_general, ONLY:natom !,nwalk
      use mod_system, ONLY:am,names
      implicit none
      real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 r
      integer :: iw,iat

      rb_sbc=rb_sbc*ang

!---calculation of Center-Of-Mass
      iw=1
      do iat=1,natom
       xcm=xcm+x(iat,iw)*am(iat)
       ycm=ycm+y(iat,iw)*am(iat)
       zcm=zcm+z(iat,iw)*am(iat)
       mass_total=mass_total+am(iat)
      enddo
      xcm=xcm/mass_total
      ycm=ycm/mass_total
      zcm=zcm/mass_total
      write(*,*)'Center of mass is at [A]:',xcm/ang,ycm/ang,zcm/ang

!-----calculation of size-of-cluster
      iw=1
      do iat=1,natom
       r=(x(iat,iw)-xcm)**2+(y(iat,iw)-ycm)**2+(z(iat,iw)-zcm)**2
       r=dsqrt(r)
       if(r.gt.rmax.and.names(iat).ne.'H')then
        rmax=r
       endif
      enddo

!---SETTING SIZE OF THE CLUSTER      
      write(*,*)'Cluster radius is r= ',rmax/ang,'  Angstrom'
      if (rb_sbc.le.0.and.rho.le.0)then
       write(*,*)'Cluster radius for spherical boundary conditions not specified.'
       write(*,*)'Will be set equal the initial radius of the system.'
       rb_sbc=rmax
      endif
!---DETERMINATION of cluster size from given density
      if (rho.gt.0)then
       write(*,*)'Calculating cluster radius from given densty.'
       rho=rho*fact !conversion from g/L to atomic units
       rb_sbc=mass_total/rho*3/4/pi
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
      use mod_general
      use mod_system, ONLY:am,names 
      implicit none
      real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(inout) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8  :: r,frb
      integer :: iat,iw

      rmax=0.0d0

!---calculation of Center-Of-Mass
      iw=1
      do iat=1,natom
       xcm=xcm+x(iat,iw)*am(iat)
       ycm=ycm+y(iat,iw)*am(iat)
       zcm=zcm+z(iat,iw)*am(iat)
      enddo
      xcm=xcm/mass_total
      ycm=ycm/mass_total
      zcm=zcm/mass_total

!---kb_sbc and rb_sbc must be specified in input.in
      do iw=1,nwalk
       do iat=1,natom
        r=(x(iat,iw)-xcm)**2+(y(iat,iw)-ycm)**2+(z(iat,iw)-zcm)**2
        r=dsqrt(r)
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

      END MODULE mod_sbc


!      MODULE mod_bag
!      use mod_array_size
!      real*8,allocatable :: neigh1(:),neigh2(:),neigh3(:)
!      real*8 :: some_parameters
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
!      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
!      real*8 fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
!      real*8  :: r,frb
!      integer :: iat,iw
!
!
!
!
!
!      END SUBROUTINE
!      END MODULE

