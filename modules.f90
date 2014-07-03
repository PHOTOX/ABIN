!--- File with modules                                created by Daniel Hollas,9.2.2012

!!!!!! We are now using modules to initialize some variables and also
!!!!!! for passing global variables to different subroutines.
!------------------------------------------------------------------------------------



!-----Contains various array limits. Modify here if you need larger arrays (for big systems)
!-----Also contains some physical constants
!----- If you set these too large, compiler will complain...maybe some
!additional flags are needed
      module mod_array_size
      implicit none
      integer,parameter :: NPARTMAX=3000,NWALKMAX=300,MAXCHAIN=10,NSHAKEMAX=1
      integer,parameter :: NBINMAX=3000,NDISTMAX=30
      integer,parameter :: NSTMAX=15,NTRAJMAX=1
      real*8, parameter :: AMU=1823.d0,ANG=1.8897d0,AUTOFS=0.02419d0,PI=3.14159265d0
      real*8, parameter :: AUTOK=3.1577464d5,me=9.109382e-31 !electron mass
      real*8, parameter :: AUTOM=5.2917720859e-11 ! atomic length
      real*8, parameter :: AMBTOAU=0.84d0/15.3067320d0 ! atomic length
      real*8, parameter :: AUTOKCAL=627.51d0,AUTOKK=3.1577322d5,AUTOEV=27.2114 
      save
      end module  
!------------------------------------------------------
!------------------------------------------------------

!--General simulation parameters      
module mod_general
   implicit none
   integer :: it=0,ipimd=0,istage=0,nwalk=1,ihess=0,md=1  
   character(len=10) :: pot='default'
   integer :: imini=0,nstep=1,nabin=50
   integer :: nwrite=1,nwritex=1,ncalc=1,nrest=1,nwritev=0
   integer :: icv=0,irest=0,anal_ext=0,idebug=0
   integer :: irandom=156873,natom=0,pid,conatom=0
   integer :: isbc=0,ibag=0  !spherical boundary conditions,elastic bag
   integer :: iqmmm=0
   real*8  :: gamm=20., gammthr=1e-10 !minthr=1e-15
   integer :: parrespa=0 
   integer :: dime=3,f=3 !dimenze systemu a pocet zakonu zachovani  
   save
end module
      
!------------------------------------------------------
!------------------------------------------------------

!-----Some information about simulated system, especially for distributions and shake
      module mod_system
      use mod_array_size
      implicit none
      real*8  :: am(npartmax)
      character(len=2) :: names(npartmax)
      integer :: inames(npartmax),imass_init=1
!---distributions (distance,angle,dihedral)
      integer :: ndist=0,nbin=1000,dist1(ndistmax),dist2(ndistmax)
      integer :: nang=0,ang1(ndistmax),ang2(ndistmax),ang3(ndistmax)
      integer :: ndih=0,dih1(ndistmax),dih2(ndistmax),dih3(ndistmax),dih4(ndistmax),nbin_ang=180
      real*8  :: dist(nbinmax,ndistmax),dist_ang(nbinmax,ndistmax),dist_dih(nbinmax,ndistmax) 
      real*8  :: xmin=0.5d0,xmax=5.0d0
      real*8  :: shiftdih=360.0d0 ! 0 for (-180,180), 360 for (0,360)
!SHAKE stuff
      real*8  :: dshake(nshakemax),shake_tol=0.001
      integer :: nshake=0,IShake1(nshakemax),IShake2(nshakemax),nmol=1
!      integer :: molfirst(npartmax)  !not in use anymore
      integer :: natmol(npartmax),nshakemol(npartmax)
      save
      CONTAINS
      subroutine mass_init()
      use mod_general, ONLY: natom
      implicit none
      integer :: i
      do i=1,natom
       if(names(i).eq.'H')then
        am(i)=1.00d0
       else if(names(i).eq.'O')then
        am(i)=16.00d0
       else if (names(i).eq.'S')then
        am(i)=32.07d0
       else if (names(i).eq.'SE')then
        am(i)=78.971d0
       else if (names(i).eq.'TE')then
        am(i)=127.60d0
       else if (names(i).eq.'N')then
        am(i)=14.00d0
       else if (names(i).eq.'P')then
        am(i)=30.97d0
       else if (names(i).eq.'AS')then
        am(i)=74.922d0
       else if (names(i).eq.'SB')then
        am(i)=121.76d0
       else if (names(i).eq.'BI')then
        am(i)=208.980d0
       else if (names(i).eq.'F')then
        am(i)=19.0d0
       else if (names(i).eq.'CL')then
        am(i)=35.45d0
       else if (names(i).eq.'BR')then
        am(i)=79.904d0
       else if (names(i).eq.'I')then
        am(i)=126.90447d0
       else if (names(i).eq.'LI')then
        am(i)=6.94d0
       else if (names(i).eq.'NA')then
        am(i)=22.99d0
       else if (names(i).eq.'K')then
        am(i)=39.1d0
       else if (names(i).eq.'BE')then
        am(i)=9.01d0
       else if (names(i).eq.'MG')then
        am(i)=24.305d0
       else if (names(i).eq.'CA')then
        am(i)=40.078d0
       else if (names(i).eq.'B')then
        am(i)=10.81d0
       else if (names(i).eq.'AL')then
        am(i)=26.982d0
       else if (names(i).eq.'C')then
        am(i)=12.00d0
       else if (names(i).eq.'SI')then
        am(i)=28.085d0
       else if (names(i).eq.'GE')then
        am(i)=72.630d0
       else if (names(i).eq.'SN')then
        am(i)=118.710d0
       else if (names(i).eq.'PB')then
        am(i)=207.2d0
       else if (names(i).eq.'HE')then
        am(i)=4.0d0
       else if (names(i).eq.'NE')then
        am(i)=20.18d0
       else if (names(i).eq.'AR')then
        am(i)=39.95d0
       else if (names(i).eq.'KR')then
        am(i)=83.798d0
       else if (names(i).eq.'XE')then
        am(i)=131.29d0
       else if (names(i).eq.'FE')then
        am(i)=55.845d0
       else if (names(i).eq.'TI')then
        am(i)=47.867d0
       else if (names(i).eq.'V')then
        am(i)=47.867d0
       else if (names(i).eq.'CR')then
        am(i)=51.9961d0
       else if (names(i).eq.'MN')then
        am(i)=54.938d0
       else if (names(i).eq.'CO')then
        am(i)=58.933d0
       else if (names(i).eq.'NI')then
        am(i)=58.693d0
       else if (names(i).eq.'CU')then
        am(i)=63.546d0
       else if (names(i).eq.'ZN')then
        am(i)=65.38d0
       else if (names(i).eq.'AG')then
        am(i)=107.8682d0
       else if (names(i).eq.'AU')then
        am(i)=196.967d0
       else if (names(i).eq.'PT')then
        am(i)=195.084d0
       else if (names(i).eq.'CD')then
        am(i)=112.414d0
       else if (names(i).eq.'HG')then
        am(i)=200.592d0
       else if (names(i).eq.'U')then
        am(i)=238.028d0
       else 
        write(*,*)'Unknown atom.',names(i),'Exiting..'
        stop
       endif
      enddo
      end subroutine

      subroutine dist_init()
         implicit none
         dist=0.0d0
         dist_ang=0.0d0
         dist_dih=0.0d0
      end subroutine

      subroutine shake_init(x,y,z)
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 xi,yi,zi,xj,yj,zj
      integer :: ixshake,i,j
        Do ixshake=1,NShake
        i=IShake1(ixshake)
        j=IShake2(ixshake)
        xi=x(i,1)
        yi=y(i,1)
        zi=z(i,1)
        xj=x(j,1)
        yj=y(j,1)
        zj=z(j,1)
        dshake(ixshake)=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
        enddo
      end subroutine
!---potentially useful for guillot and other empirical force fields, because string
!---comparison is very cpu demanding!!
      subroutine inames_guillot()
      use mod_general, ONLY:natom
      implicit none
      integer :: i
      do i=1,natom
       if(names(i).eq.'H')then
        inames(i)=1
       else if(names(i).eq.'O')then
        inames(i)=0
       else if (names(i).eq.'CL')then
        inames(i)=2
       endif
      enddo
     end subroutine

      end module 
!------------------------------------------------------
!------------------------------------------------------
