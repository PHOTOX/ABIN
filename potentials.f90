!-----Various analytical potentials
!---- Now includided Harmonic for diatomic molecule, 3D-harmonic and Morse + harmonic and morse hessians
!---- Parameters should be set in module mod_harmon or in input section system
!------------------------------------------------------
!------------------------------------------------------
!some constanst for analytical potentials      
      module mod_harmon
      use mod_array_size
      implicit none
!---constants for 3DHO
      real*8 :: k1=0.0d0,k2=0.0d0,k3=0.0d0
!---constants for 1D 2-particle harmonic oscillator
      real*8 :: k=0.000, r0=0.0d0
!----   CONSTANTS for morse potential ccc
!            V=De*(1-exp(-a(r-r0)))^2
      real*8 :: De=0.059167,a=-1.0d0
      real*8,allocatable :: hess(:,:,:)
      save
      contains

!------3D Harmonic Oscillator---only 1 particle!!
      subroutine force_2dho(x,y,z,fxab,fyab,fzab,eclas)
      use mod_array_size
      use mod_general
      implicit none
      real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(out) :: fxab(npartmax,nwalkmax),fyab(npartmax,nwalkmax),fzab(npartmax,nwalkmax)
      real*8,intent(out) :: eclas
      real*8             :: energy
      integer            :: iw

      eclas=0.0d0
      energy=0.0d0
      do iw=1,nwalk
        fxab(1,iw)=-k1*x(1,iw)
        fyab(1,iw)=-k2*y(1,iw)
        fzab(1,iw)=-k3*z(1,iw)
        energy=energy+0.5*k1*x(1,iw)**2+0.5*k2*y(1,iw)**2
        energy=energy+0.5*k3*z(1,iw)**2
      enddo

      eclas=energy/nwalk

      return
      end
      

!ccccccccHARMONIC OSCILLATOR--diatomic molecules--ccccccccccccccccccccc
      subroutine force_harmon(x,y,z,fxab,fyab,fzab,eclas)
      use mod_array_size
      use mod_general
      implicit none
      real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(out) :: fxab(npartmax,nwalkmax),fyab(npartmax,nwalkmax),fzab(npartmax,nwalkmax)
      real*8,intent(out) :: eclas
      real*8             :: dx,dy,dz,r,fac
      integer            :: i

        eclas=0.0d0

        do i=1,nwalk
         dx=x(2,i)-x(1,i)
         dy=y(2,i)-y(1,i)
         dz=z(2,i)-z(1,i)
         r=dx**2+dy**2+dz**2
         r=dsqrt(r)
         fac=k*(r-r0)/r
         fxab(1,i)=fac*dx
         fxab(2,i)=-fxab(1,i)
         fyab(1,i)=fac*dy
         fyab(2,i)=-fyab(1,i)
         fzab(1,i)=fac*dz
         fzab(2,i)=-fzab(1,i)
         eclas=eclas+0.5*k*(r-r0)**2/nwalk
        enddo

        end



      subroutine hess_harmon(x,y,z)
      use mod_array_size
      use mod_general
      implicit none
      real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8            :: dx,dy,dz,r,fac
      integer           :: i,ipom1,ipom2

        do i=1,nwalk

         dx=x(2,i)-x(1,i)
         dy=y(2,i)-y(1,i)
         dz=z(2,i)-z(1,i)
         r=dx**2+dy**2+dz**2
         r=dsqrt(r)
         fac=k*(r-r0)/r
         hess(1,1,i)=(k*dx**2/r**2-fac*dx**2/r**2+fac)/nwalk
         hess(2,2,i)=(k*dy**2/r**2-fac*dy**2/r**2+fac)/nwalk
         hess(3,3,i)=(k*dz**2/r**2-fac*dz**2/r**2+fac)/nwalk
         hess(4,4,i)=hess(1,1,i)
         hess(5,5,i)=hess(2,2,i)
         hess(6,6,i)=hess(3,3,i)
         hess(2,1,i)=(k*dx*dy/r**2-fac*dx*dy/r**2)/nwalk
         hess(3,1,i)=(k*dz*dx/r**2-fac*dz*dx/r**2)/nwalk
         hess(3,2,i)=(k*dz*dy/r**2-fac*dz*dy/r**2)/nwalk
         hess(4,1,i)=-hess(1,1,i)
         hess(4,2,i)=-hess(2,1,i)
         hess(4,3,i)=-hess(3,1,i)
         hess(5,1,i)=-hess(2,1,i)
         hess(5,2,i)=-hess(2,2,i)
         hess(5,3,i)=-hess(3,2,i)
         hess(5,4,i)=hess(2,1,i)
         hess(6,1,i)=-hess(3,1,i)
         hess(6,2,i)=-hess(3,2,i)
         hess(6,3,i)=-hess(3,3,i)
         hess(6,4,i)=hess(3,1,i)
         hess(6,5,i)=hess(3,2,i)
         do ipom1=1,5
          do ipom2=ipom1+1,6
          hess(ipom1,ipom2,i)=hess(ipom2,ipom1,i)
          enddo
         enddo
!      nwalk enddo
       enddo


      end


      subroutine force_morse(x,y,z,fxab,fyab,fzab,eclas)
      use mod_array_size
      use mod_general
      implicit none
      real*8,intent(in)  ::  x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(out) ::  fxab(npartmax,nwalkmax),fyab(npartmax,nwalkmax),fzab(npartmax,nwalkmax)
      real*8,intent(out) :: eclas
      real*8 :: dx,dy,dz,r,fac,ex
      integer :: i

!cccccccc  V=De*(1-exp(-a(r-r0)))^2
!NOT REALLY SURE about a
!if it is not set from input, we determine it from k(normaly used for
!harmon osciallator)
      if(a.le.0) a=dsqrt(k/2/De)
        eclas=0.0d0

        do i=1,nwalk
         dx=x(2,i)-x(1,i)
         dy=y(2,i)-y(1,i)
         dz=z(2,i)-z(1,i)
         r=dx**2+dy**2+dz**2
         r=dsqrt(r)
         ex=exp(-a*(r-r0))
         fac=2*a*ex*De*(1-ex)/r
         fxab(1,i)=fac*dx
         fxab(2,i)=-fxab(1,i)
         fyab(1,i)=fac*dy
         fyab(2,i)=-fyab(1,i)
         fzab(1,i)=fac*dz
         fzab(2,i)=-fzab(1,i)
         eclas=eclas+De*(1-ex)**2/nwalk
        enddo

        end


      subroutine hess_morse(x,y,z)
      use mod_array_size
      use mod_general
      implicit none
      real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8             ::  dx,dy,dz,r,fac,ex,fac2
      integer            :: i,ipom1,ipom2

!NOT REALLY SURE about a
       a=dsqrt(k/2/De)

       do i=1,nwalk

        dx=x(2,i)-x(1,i)
        dy=y(2,i)-y(1,i)
        dz=z(2,i)-z(1,i)
        r=dx**2+dy**2+dz**2
        r=dsqrt(r)
        ex=exp(-a*(r-r0))
        fac=2*a*ex*De*(1-ex)/r
        fac2=2*De*a**2*ex**2/r**2
        hess(1,1,i)=fac2*dx**2-(fac*dx**2)/r**2+fac-fac*a*dx**2/r
        hess(2,2,i)=fac2*dy**2-(fac*dy**2)/r**2+fac-fac*a*dy**2/r
        hess(3,3,i)=fac2*dz**2-(fac*dz**2)/r**2+fac-fac*a*dz**2/r
        hess(1,1,i)=hess(1,1,i)/nwalk
        hess(2,2,i)=hess(2,2,i)/nwalk
        hess(3,3,i)=hess(3,3,i)/nwalk
        hess(4,4,i)=hess(1,1,i)
        hess(5,5,i)=hess(2,2,i)
        hess(6,6,i)=hess(3,3,i)
        hess(2,1,i)=fac2*dx*dy-(fac*dx*dy)/r**2-fac*a*dx*dy/r
        hess(3,1,i)=fac2*dx*dz-(fac*dx*dz)/r**2-fac*a*dx*dz/r
        hess(3,2,i)=fac2*dz*dy-(fac*dz*dy)/r**2-fac*a*dz*dy/r
        hess(1,1,i)=hess(2,1,i)/nwalk
        hess(2,2,i)=hess(3,1,i)/nwalk
        hess(3,3,i)=hess(3,2,i)/nwalk
        hess(4,1,i)=-hess(1,1,i)
        hess(4,2,i)=-hess(2,1,i)
        hess(4,3,i)=-hess(3,1,i)
        hess(5,1,i)=-hess(2,1,i)
        hess(5,2,i)=-hess(2,2,i)
        hess(5,3,i)=-hess(3,2,i)
        hess(5,4,i)=hess(2,1,i)
        hess(6,1,i)=-hess(3,1,i)
        hess(6,2,i)=-hess(3,2,i)
        hess(6,3,i)=-hess(3,3,i)
        hess(6,4,i)=hess(3,1,i)
        hess(6,5,i)=hess(3,2,i)
        do ipom1=1,5
         do ipom2=ipom1+1,6
          hess(ipom1,ipom2,i)=hess(ipom2,ipom1,i)
         enddo
        enddo

!     nwalk enddo
      enddo


      end

   !theoretically only needs to be called once, but nah
   subroutine hess_2dho()
   use mod_array_size
   use mod_general
   implicit none
   integer :: i

   do i=1,nwalk
      hess(1,1,i)=k1/nwalk
      hess(2,2,i)=k2/nwalk
      hess(3,3,i)=k3/nwalk
      hess(2,1,i)=0.0d0
      hess(3,1,i)=0.0d0
      hess(3,2,i)=0.0d0
      hess(1,2,i)=0.0d0
      hess(1,3,i)=0.0d0
      hess(2,3,i)=0.0d0
   enddo

   end

end module mod_harmon
