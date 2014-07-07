!---This is a template file for user-defined analysis function.
!---This function will be called if anal_ext=1 in input(section general).
!---Should the user need something more then coordinates(velocities,forces),
! he/she must also modify  analysis.f90 and possibly also  abin.f90
      module mod_analyze_ext
      use mod_array_size
      real*8  :: xmin=-20.0d0,xmax=20.0d0
      integer,parameter :: nbin=500
      real*8  dist(nbin,nwalkmax)
      save
      contains
      subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: names,am
      implicit none
      real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(in) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8,intent(in) :: amt(npartmax,nwalkmax)
      real*8  :: anorm,dx,dbin
      integer :: iw,ipom,ian,iat
      real*8  :: p(npartmax,nwalkmax)

      iat=1

       dbin=(xmax-xmin)/nbin
       do iw=1,nwalk


       p(iat,iw)=(vx(iat,iw)**2+vy(iat,iw)**2+vz(iat,iw)**2)
       p(iat,iw)=am(iat)*sqrt(p(iat,iw))

       if(dime.eq.1.and.pot.eq.'2dho') p(1,iw)=vx(1,iw)*am(iat)


       ipom=ceiling( ((p(iat,iw))-xmin)/dbin )
       if(ipom.gt.nbin.or.ipom.le.0)then
        write(*,*)'problems with p distribution function'
        write(*,*)'p=',p(iat,iw)
        stop
       endif

       dist(ipom,iw)=dist(ipom,iw)+1.0
       enddo

      if(modulo(it,nwrite).eq.0)then
!--print output every nwrite steps              
        open(128,file='distp.dat')
        do iw=1,nwalk
        anorm=0.0d0
        dx=(xmax-xmin)/nbin

        do ian=1,nbin
         anorm=anorm+dist(ian,iw)
        enddo
        
        do ian=1,nbin
         write(128,*)ian*dx+xmin,dist(ian,iw)/(anorm*dx)
        enddo
        write(128,*)

        enddo

        close(128)
      endif
   
     end subroutine analyze_ext

     end module
                                        


