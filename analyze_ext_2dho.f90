!---This is a template file for user-defined analysis function.
!---This function will be called if anal_ext=1 in input(section general).
!---Should the user need something more then coordinates(velocities,forces),
! he/she must also modify  analysis.f90 and possibly abin.f90
module mod_analyze_ext
   use mod_const, only: DP
   use mod_array_size, only: nbinmax
   implicit none
   private 
   public :: analyze_ext
!----force constants for displaced oscillator      
   real(DP)  :: k1p=0.334d0,k2p=0.00334d0
!----displacement
   real(DP)  :: dd1=0.5d0,dd2=5.00d0
!----v0 is adiabatic excitation energy 
   real(DP)  :: v0=0.1d0
   real(DP)  :: abs(nbinmax)
   real(DP)  :: emin=-0.9d0,emax=0.9d0
   integer :: nbinen=400
   save
   contains
!---- Absorption spectrum via reflection principle
   subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
      use mod_general
      use mod_system, ONLY: names
      use mod_harmon
      real(DP) x(:,:),y(:,:),z(:,:)
      real(DP) vx(:,:),vy(:,:),vz(:,:)
      real(DP) amt(:,:)
      integer :: ibinn,iw
      real(DP) dbin,anorm,dev,abs_mean,v1,v2,xx1,yy1

       dbin=(emax-emin)/nbinen

       do iw=1,nwalk
        xx1=x(1,iw)
        yy1=y(1,iw)
        v1=0.5*(k1*xx1**2+k2*yy1**2)
        v2=v0+0.5*(k1p*(xx1-dd1)**2+k2p*(yy1-dd2)**2)
        ibinn=(v2-v1-emin)/dbin
!        write(*,*)'ibinn',ibinn,v2,v1,nwalk
        abs(ibinn)=abs(ibinn)+1
       enddo



      if(modulo(it,nwrite).eq.0)then
!--print output every nwrite steps              
       open(180,file='abs.dat')
       open(181,file='deviation.dat')

       abs_mean=0.0d0
       anorm=0.0d0
       do ibinn=1,nbinen
        write(180,*)emin+ibinn*dbin,abs(ibinn)
        abs_mean=abs_mean+(emin+ibinn*dbin)*abs(ibinn)
        anorm=anorm+abs(ibinn)
       enddo
       abs_mean=abs_mean/anorm
       dev=0.0d0
       do ibinn=1,nbinen
        dev=dev+((emin+ibinn*dbin-abs_mean)**2)*abs(ibinn)
       enddo
       dev=dev/anorm
       write(181,*)it,dev,1/(2.0d0*dev)
       
       close(180)
       close(181)

      endif
   
      return
      end

      end module
                                        



