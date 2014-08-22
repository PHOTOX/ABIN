      module mod_array_size
      implicit none
      integer,parameter :: npartmax=2000,nwalkmax=1
      integer,parameter :: nbinmax=3100,ndistmax=30
      real*8, parameter :: amu=1823.d0,ang=1.89d0,pi=3.14159265
      integer, parameter :: ndist=5,nang=0,ndih=0
      character*2 :: shit
      save
      end module  


program analyze_movie
      use mod_array_size
      implicit none
      real*8  :: x(npartmax,1),y(npartmax,1),z(npartmax,1)
      real*8  :: r,anorm,dx,dbin
      integer :: idist,iw,ipom,ian
      real*8  :: dihmin,dihmax,delta,get_dihedral
      real*8  :: angmin,angmax,alfa,get_angle
      integer :: nbin,natfirst(ndist),natsecond(ndist)
      integer :: angles(3,10),dih(4,10),nbin_ang
      real*8  :: dist(nbinmax,ndist),dist_ang(nbinmax,ndist),dist_dih(nbinmax,ndist) 
      real*8  :: xmin(ndist),xmax(ndist)
      integer :: ios,iat,natom=6,igeom,ngeoms=5000
      real*8  :: shiftdih=360.d0,time
      character*20 :: chgeom

      natfirst=(/ 1,3,4,5,6 /)
      natsecond=(/ 2,1,1,2,2 /)
!     angles = reshape((/ 3,1,4, 5, 2, 6,1,4,2 /), shape(angles))
!     dih = reshape((/ 3,1,2, 5/), shape(dih))
      xmin=(/0.5,0.5,0.5,0.5,0.5/)
      xmax=(/5.5,5.5,5.5,5.5,5.5/)

      nbin=500
      nbin_ang=90
      dihmin=0.0
      dihmax=180.0
      angmin=0.0
      angmax=180.0

      do igeom=1,ngeoms

      if(igeom.lt.10)then
      write(chgeom,'(A,I1)')'FMSTRAJS/Traj.',igeom
      elseif (igeom.lt.100)then
      write(chgeom,'(A,I2)')'FMSTRAJS/Traj.',igeom
      elseif (igeom.lt.1000)then
      write(chgeom,'(A,I3)')'FMSTRAJS/Traj.',igeom
      elseif (igeom.lt.10000)then
      write(chgeom,'(A,I4)')'FMSTRAJS/Traj.',igeom
      endif
!      write(*,*)chgeom
      open(100,file=chgeom)
      read(100,*,IOSTAT=ios)time,x(1,1),y(1,1),z(1,1),x(2,1),y(2,1),z(2,1),x(3,1),y(3,1),z(3,1) &
                               ,x(4,1),y(4,1),z(4,1),x(5,1),y(5,1),z(5,1),x(6,1),y(6,1),z(6,1)
      close(100)

      if(ios.ne.0)then
              write(*,*)'Traj doesn.t exist',igeom
              cycle
      endif
!      write(*,*)natom
!      write(*,*)
      do iat=1,natom
       x(iat,1)=x(iat,1)/ang
       y(iat,1)=y(iat,1)/ang
       z(iat,1)=z(iat,1)/ang
!       write(*,*)'O', x(iat,1),y(iat,1),z(iat,1)
      enddo


      iw=1

      do idist=1,ndist

       r=(x(natfirst(idist),iw)-x(natsecond(idist),iw))**2
       r=r+(y(natfirst(idist),iw)-y(natsecond(idist),iw))**2
       r=r+(z(natfirst(idist),iw)-z(natsecond(idist),iw))**2
       r=dsqrt(r)

       dbin=(xmax(idist)-xmin(idist))/nbin

       ipom=((r)-xmin(idist))/dbin
       if(ipom.gt.nbin.or.ipom.lt.0)then
        write(*,*)'problems with distribution function'
        stop
       endif

       dist(ipom,idist)=dist(ipom,idist)+1.0
      enddo

      dbin=(angmax-angmin)/nbin_ang

      do idist=1,nang

       alfa=get_angle(x,y,z,angles(1,idist),angles(2,idist),angles(3,idist),iw)
       ipom=(alfa-angmin)/dbin

       if(ipom.gt.nbin_ang.or.ipom.lt.0)then
        write(*,*)'problems with angle distribution function'
        write(*,*)'For angle between atoms:',angles(1,idist),angles(2,idist),angles(3,idist)
        write(*,*)'Value of ipom=',ipom
        stop
       endif

       dist_ang(ipom,idist)=dist_ang(ipom,idist)+1.0
      enddo

      dbin=(dihmax-dihmin)/nbin_ang

      do idist=1,ndih
       delta=get_dihedral(x,y,z,dih(1,idist),dih(2,idist),dih(3,idist),dih(4,idist),iw)
       ipom=(delta-dihmin)/dbin

       if(ipom.gt.nbin_ang.or.ipom.lt.0)then
        write(*,*)'problems with dihedral angle distribution function'
        write(*,*)'For dihedral between atoms:',dih(1,idist),dih(2,idist),dih(3,idist),dih(4,idist)
        write(*,*)'Value of ipom=',ipom
        stop
       endif
       dist_dih(ipom,idist)=dist_dih(ipom,idist)+1.0
      enddo

      write(*,*)igeom

      enddo
      close(100)
      !--------PRINTING---------------------

       open(10,file='ang_wigner.dat')
       do idist=1,nang
        anorm=0.0

        do ian=1,nbin_ang
         anorm=anorm+dist_ang(ian,idist)
        enddo
       
        do ian=1,nbin_ang
         write(10,*)ian*dbin+angmin,dist_ang(ian,idist)/(anorm*dbin)
        enddo

        write(10,*)
       enddo
       close(10)

        open(128,file='dist_wigner.dat')
        do idist=1,ndist
        anorm=0.0
        dx=(xmax(idist)-xmin(idist))/nbin

        do ian=1,nbin
         anorm=anorm+dist(ian,idist)
        enddo
        
        do ian=1,nbin
         write(128,*)ian*dx+xmin(idist),dist(ian,idist)/(anorm*dx)
        enddo
        write(128,*)

        enddo

        close(128)

       open(10,file='dih_wigner.dat')
       do idist=1,ndih
        anorm=0.0

        do ian=1,nbin_ang
         anorm=anorm+dist_dih(ian,idist)
        enddo
       
        do ian=1,nbin_ang
         write(10,*)ian*dbin+dihmin,dist_dih(ian,idist)/(anorm*dbin)
        enddo

        write(10,*)
       enddo
       close(10)


      end



      real*8 function get_angle(x,y,z,at1,at2,at3,iw)
      use mod_array_size
      implicit none
      real*8 x(npartmax,1),y(npartmax,1),z(npartmax,1)
      integer :: at1,at2,at3,iw
      real*8  :: vec1x,vec1y,vec1z
      real*8  :: vec2x,vec2y,vec2z

      vec1x=x(at1,iw)-x(at2,iw)
      vec1y=y(at1,iw)-y(at2,iw)
      vec1z=z(at1,iw)-z(at2,iw)
      vec2x=x(at3,iw)-x(at2,iw)
      vec2y=y(at3,iw)-y(at2,iw)
      vec2z=z(at3,iw)-z(at2,iw)
      get_angle=180/pi*acos((vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/ &
      (sqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2)))

      return 
      end


      real*8 function get_dihedral(x,y,z,at1,at2,at3,at4,iw)
      use mod_array_size
      implicit none
      real*8 x(npartmax,1),y(npartmax,1),z(npartmax,1)
      integer :: at1,at2,at3,at4,iw
      real*8  :: vec1x,vec1y,vec1z
      real*8  :: vec2x,vec2y,vec2z
      real*8  :: vec3x,vec3y,vec3z,sign
      real*8  :: norm1x,norm1y,norm1z,norm2x,norm2y,norm2z
      real*8  :: shiftdih=360.d0
      

      vec1x=x(at1,iw)-x(at2,iw)
      vec1y=y(at1,iw)-y(at2,iw)
      vec1z=z(at1,iw)-z(at2,iw)
      vec2x=x(at3,iw)-x(at2,iw)
      vec2y=y(at3,iw)-y(at2,iw)
      vec2z=z(at3,iw)-z(at2,iw)
      vec3x=x(at4,iw)-x(at3,iw)
      vec3y=y(at4,iw)-y(at3,iw)
      vec3z=z(at4,iw)-z(at3,iw)

      norm1x=vec1y*vec2z-vec1z*vec2y
      norm1y=vec1z*vec2x-vec1x*vec2z
      norm1z=vec1x*vec2y-vec1y*vec2x
      norm2x=vec3y*vec2z-vec3z*vec2y
      norm2y=vec3z*vec2x-vec3x*vec2z
      norm2z=vec3x*vec2y-vec3y*vec2x

      sign = norm1x*vec3x+norm1y*vec3y+norm1z*vec3z
      get_dihedral = 180/pi*acos((norm1x*norm2x+norm1y*norm2y+norm1z*norm2z)/ &
      (sqrt(norm1x**2+norm1y**2+norm1z**2)*sqrt(norm2x**2+norm2y**2+norm2z**2)))

      !if (sign.lt.0) get_dihedral = shiftdih-get_dihedral

      return
      end
      
      
      
      
