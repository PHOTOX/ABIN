!-This is a template file for user-defined analysis function.
!-This function will be called if anal_ext=1 in input(section general).
!-Should the user need something more then coordinates(velocities,forces),
! he/she must also modify  analysis.f90 and possibly also  abin.f90
module mod_analyze_ext
   use mod_const, only: DP
   use mod_array_size, only: nbinmax
   implicit none
   private
   public :: analyze_ext
   real(DP)  :: xmin=-20.0d0,xmax=20.0d0
   integer,parameter :: nbin=500
   real(DP), allocatable :: dist(:,:)
   save
   contains

   subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
      use mod_general,  only: it, natom, nwalk, pot, nwrite
      use mod_system,   only: names, am, dime
      use mod_utils,    only: abinerror
      real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(in) :: amt(:,:)
      real(DP)  :: p(size(vx,1),size(vx,2))
      real(DP)  :: anorm,dx,dbin
      integer :: iw,ipom,ian,iat

      if (.not.allocated( dist )) allocate( dist(nbinmax, nwalk) )

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
        call abinerror('analyze_ext')
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
                                        


