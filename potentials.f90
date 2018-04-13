!  Various analytical potentials
!  Now including harmonic potential for diatomic molecule, 3D-harmonic and Morse
!  For all these, hessian is also available.
!  Now also includes double-well potential, originaly used for testing PLUMED.
!  
!  Parameters should be set in input section system

!------------------------------------------------------
!     some constansts for analytical potentials      
module mod_harmon
   use mod_const, only: DP
   implicit none
!--constants for 3DHO
   real(DP) :: k1=0.0d0,k2=0.0d0,k3=0.0d0
!--constants for 1D 2-particle harmonic oscillator
   real(DP) :: k=0.000d0, r0=0.0d0
!--constants for double well 
   real(DP) :: lambda_dw=0.0d0, D0_dw=0.0d0, k_dw, r0_dw
!--CONSTANTS for morse potential ccc
!  V=De*(1-exp(-a(r-r0)))^2
   real(DP) :: De=0.059167d0,a=-1.0d0
   real(DP),allocatable :: hess(:,:,:)
   save
   CONTAINS

!------3D Harmonic Oscillator---only 1 particle!!
   subroutine force_2dho(x,y,z,fxab,fyab,fzab,eclas)
      use mod_general, only: nwalk
      real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out) :: fxab(:,:),fyab(:,:),fzab(:,:)
      real(DP),intent(out) :: eclas
      real(DP)             :: energy
      integer              :: iw

      eclas = 0.0d0
      energy = 0.0d0
      do iw=1,nwalk
        fxab(1,iw) = -k1 * x(1,iw)
        fyab(1,iw) = -k2 * y(1,iw)
        fzab(1,iw) = -k3 * z(1,iw)
        energy = energy + 0.5d0*k1*x(1,iw)**2 + 0.5d0*k2*y(1,iw)**2
        energy = energy + 0.5d0*k3*z(1,iw)**2
      enddo

      eclas = energy / nwalk

      return
   end subroutine force_2dho
      
!oooooooo DOUBLE WELL potential -- ref. Tuckerman, Statistical Mechanics,p350 --oooooooooooooooo
  subroutine  force_doublewell(x,y,z,fxab,fyab,fzab,eclas)
      use mod_general, only: nwalk
      real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out) :: fxab(:,:),fyab(:,:),fzab(:,:)
      real(DP),intent(out) :: eclas
      integer              :: i

      eclas=0.0d0
 
      do i=1,nwalk
         fxab(1,i) = -4*D0_dw*x(1,i) * (x(1,i)**2-r0_dw**2) - lambda_dw*y(1,i)
         fyab(1,i) = -k_dw * y(1,i) - lambda_dw*x(1,i)
         fzab(1,i) = 0
         eclas = eclas + D0_dw*(x(1,i)**2-r0_dw**2)**2 + &
            0.5*k_dw*y(1,i)**2 + lambda_dw*x(1,i)*y(1,i)
      enddo

      eclas = eclas / nwalk

  end subroutine force_doublewell


!ccccccccHARMONIC OSCILLATOR--diatomic molecules--ccccccccccccccccccccc
   subroutine force_harmon(x,y,z,fxab,fyab,fzab,eclas)
      use mod_general, only: nwalk
      real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out) :: fxab(:,:),fyab(:,:),fzab(:,:)
      real(DP),intent(out) :: eclas
      real(DP)             :: dx,dy,dz,r,fac
      integer            :: i

      eclas = 0.0d0

      do i=1,nwalk
         dx = x(2,i) - x(1,i)
         dy = y(2,i) - y(1,i)
         dz = z(2,i) - z(1,i)
         r = dx**2 + dy**2 + dz**2
         r = sqrt(r)
         fac = k*(r-r0) / r
         fxab(1,i) = fac * dx
         fxab(2,i) = -fxab(1,i)
         fyab(1,i) = fac * dy
         fyab(2,i) = -fyab(1,i)
         fzab(1,i) = fac * dz
         fzab(2,i) = -fzab(1,i)
         eclas = eclas + 0.5d0*k*(r-r0)**2 / nwalk
      enddo

   end subroutine force_harmon



   subroutine hess_harmon(x,y,z)
      use mod_general, only: nwalk
      real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
      real(DP)            :: dx,dy,dz,r,fac
      integer           :: i,ipom1,ipom2

      do i=1,nwalk

         dx = x(2,i) - x(1,i)
         dy = y(2,i) - y(1,i)
         dz = z(2,i) - z(1,i)
         r = dx**2 + dy**2 + dz**2
         r = sqrt(r)
         fac = k*(r-r0) / r
         hess(1,1,i) =  (k*dx**2/r**2-fac*dx**2/r**2+fac)/nwalk
         hess(2,2,i) =  (k*dy**2/r**2-fac*dy**2/r**2+fac)/nwalk
         hess(3,3,i) =  (k*dz**2/r**2-fac*dz**2/r**2+fac)/nwalk
         hess(4,4,i) =  hess(1,1,i)
         hess(5,5,i) =  hess(2,2,i)
         hess(6,6,i) =  hess(3,3,i)
         hess(2,1,i) =  (k*dx*dy/r**2-fac*dx*dy/r**2)/nwalk
         hess(3,1,i) =  (k*dz*dx/r**2-fac*dz*dx/r**2)/nwalk
         hess(3,2,i) =  (k*dz*dy/r**2-fac*dz*dy/r**2)/nwalk
         hess(4,1,i) = -hess(1,1,i)
         hess(4,2,i) = -hess(2,1,i)
         hess(4,3,i) = -hess(3,1,i)
         hess(5,1,i) = -hess(2,1,i)
         hess(5,2,i) = -hess(2,2,i)
         hess(5,3,i) = -hess(3,2,i)
         hess(5,4,i) =  hess(2,1,i)
         hess(6,1,i) = -hess(3,1,i)
         hess(6,2,i) = -hess(3,2,i)
         hess(6,3,i) = -hess(3,3,i)
         hess(6,4,i) =  hess(3,1,i)
         hess(6,5,i) =  hess(3,2,i)

         do ipom1=1,5
            do ipom2=ipom1+1,6
               hess(ipom1,ipom2,i) = hess(ipom2,ipom1,i)
            enddo
         enddo

!     nwalk enddo
      enddo


    end subroutine hess_harmon


   subroutine force_morse(x,y,z,fxab,fyab,fzab,eclas)
      use mod_general, only: nwalk
      real(DP),intent(in)  ::  x(:,:),y(:,:),z(:,:)
      real(DP),intent(out) ::  fxab(:,:),fyab(:,:),fzab(:,:)
      real(DP),intent(out) :: eclas
      real(DP) :: dx,dy,dz,r,fac,ex
      integer :: i

!cccccccc  V=De*(1-exp(-a(r-r0)))^2
!NOT REALLY SURE about a
!if it is not set from input, we determine it from k(normaly used for
!harmon osciallator)
      if(a.le.0) a=sqrt(k/2/De)

      eclas = 0.0d0

      do i=1,nwalk
         dx = x(2,i) - x(1,i)
         dy = y(2,i) - y(1,i)
         dz = z(2,i) - z(1,i)
         r = dx**2 + dy**2 + dz**2
         r = sqrt(r)
         ex = exp(-a*(r-r0))
         fac = 2*a*ex*De*(1-ex) / r
         fxab(1,i) = fac*dx
         fxab(2,i) = -fxab(1,i)
         fyab(1,i) = fac*dy
         fyab(2,i) = -fyab(1,i)
         fzab(1,i) = fac*dz
         fzab(2,i) = -fzab(1,i)
         eclas = eclas + De*(1-ex)**2 / nwalk
      enddo

   end subroutine force_morse

   subroutine hess_morse(x,y,z)
      use mod_general, only: nwalk
      real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
      real(DP)             :: dx,dy,dz,r,fac,ex,fac2
      integer              :: i,ipom1,ipom2

!NOT REALLY SURE about a
      a=sqrt(k/2/De)

      do i=1,nwalk

         dx=x(2,i)-x(1,i)
         dy=y(2,i)-y(1,i)
         dz=z(2,i)-z(1,i)
         r=dx**2+dy**2+dz**2
         r=sqrt(r)
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


   end subroutine hess_morse

   !theoretically only needs to be called once, but nah
   subroutine hess_2dho()
   use mod_general, only: nwalk
   integer :: i

   do i=1,nwalk
      hess(1,1,i) = k1 / nwalk
      hess(2,2,i) = k2 / nwalk
      hess(3,3,i) = k3 / nwalk
      hess(2,1,i) = 0.0d0
      hess(3,1,i) = 0.0d0
      hess(3,2,i) = 0.0d0
      hess(1,2,i) = 0.0d0
      hess(1,3,i) = 0.0d0
      hess(2,3,i) = 0.0d0
   enddo

   end subroutine hess_2dho

end module mod_harmon


module mod_splined_grid
   use mod_const, only: DP
   implicit none
   integer, parameter :: MAX_GRID_SIZE = 10000
   integer  :: grid_size
   real(DP), allocatable :: second_derivatives(:)
   real(DP) :: x_grid(MAX_GRID_SIZE)
   real(DP) :: y_grid(MAX_GRID_SIZE)
   real(DP) :: x_max, x_min
   private
   public :: force_splined_grid, initialize_spline
   save
CONTAINS

   real(DP) function  potential_cubic_spline(X)
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: Y
      call splint(x_grid, y_grid, second_derivatives, grid_size, X, Y)
      potential_cubic_spline = Y
      return
   end function potential_cubic_spline

   subroutine force_splined_grid(x, y, z, fx, fy, fz, eclas)
      use mod_general,  only: nwalk
      use mod_utils,    only: abinerror
      real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out) :: fx(:,:),fy(:,:),fz(:,:)
      real(DP),intent(out) :: eclas
      real(DP) :: en_1, en_2, dx = 0.0001d0
!      real(DP) :: potential_cubic_spline
      integer  :: iw

      fy = 0.0d0; fz = 0.0d0; eclas = 0.0d0
      do iw = 1, nwalk

         if(x(1, iw).lt.x_min.or.x(1, iw).gt.x_max)then
            write(*,*)"ERROR: Particle got out of the grid!"
            call abinerror("force_splined_grid")
         end if

         eclas = eclas + potential_cubic_spline(x(1,iw))
            
         ! Let's just do numerical forces for now!
         ! use 3-point numerical derivative
         en_1 = potential_cubic_spline(x(1,iw)-dx)
         en_2 = potential_cubic_spline(x(1,iw)+dx)
         fx(1, iw) = (en_1-en_2) / 2 / dx

      end do

   end subroutine force_splined_grid

   subroutine initialize_spline()
      use mod_general, only: idebug
      implicit none
      character(len=100) :: chfile_pot='potential.dat'
      real(DP) :: yp1, ypn ! first derivatives at the grid edges
      integer  :: IUNIT=600, i
      real(DP) :: x, dx
    
      grid_size = 0

      ! TODO: make x_grid et al allocatable 
      open(IUNIT,file=chfile_pot, status='OLD', action='read')
      do
         read(IUNIT, *, end = 51)x_grid(grid_size+1), y_grid(grid_size+1)
         grid_size = grid_size + 1
      enddo
51    close(IUNIT)

      if(grid_size.lt.3)then
         write(*,*)"ERROR: Could not find grid in "//chfile_pot
         stop 1
      end if

      allocate(second_derivatives(grid_size))
      x_min = x_grid(1)
      x_max = x_grid(grid_size)

!     Settings for natural cubic spline
      yp1 = 1e31
      ypn = 1e31
!     On second thought, let's just set the first derivatives the same
!     as at the boundaries
!      yp1 = (xa(2)-xa(1)) / (ya(1)-ya(2))
!      ypn = 1e31
      call spline(x_grid, y_grid, grid_size, yp1, ypn, second_derivatives)

      ! Print out the splined potential just to be sure
      if (idebug.eq.1)then
         x = x_min
         dx = (x_max - x_min) / grid_size / 5
         open(IUNIT, file="potential_splined.dat", action="write")
!         write(IUNIT, *)x_min, x_max, grid_size
         do i = 1, grid_size * 5
            write(IUNIT,*)x, potential_cubic_spline(x)
            x = x + dx
         end do
         close(IUNIT)
      end if
   end subroutine initialize_spline



   ! From numerical recipies, sliglty modified
   SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit real*8(a-h,o-z)
      INTEGER, INTENT(IN)  :: N
      REAL(DP), INTENT(IN) :: X, XA(N), YA(N), Y2A(N)
      REAL(DP), INTENT(OUT) :: Y
      INTEGER  :: KLO, KHI, K
      KLO=1
      KHI=N
      DO WHILE (KHI-KLO.GT.1)
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      END DO
      H = XA(KHI) - XA(KLO)
      IF (H.EQ.0.) THEN
         write(*,*) 'Bad XA input.'
         stop 1 
      END IF
!     DH: when extrapolating out of boundaries, just do linear
!     TODO: do this also to the END boundary
!      if(x.lt.xa(1))then
!         Y2A(KHI) = 0
!         Y2A(KLO) = 0
!      end if
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO) + B*YA(KHI) + &
           ((A**3-A) * Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.

      ! TODO implement first derivative
      RETURN
   END

   SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER  :: n, NMAX
      REAL(DP)   :: yp1, ypn, x(n), y(n), y2(n)
      PARAMETER (NMAX=500)
      INTEGER i, k
      REAL(DP)  p, qn, sig, un, u(NMAX)

      if (yp1.gt..99e30) then
         y2(1) = 0.
         u(1) = 0.
      else
         y2(1) = -0.5
         u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.
         y2(i) = (sig-1.)/p
         u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
              / (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then
         qn = 0.
         un = 0.
      else
         qn = 0.5
         un = (3./(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1)) / (x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k = n-1, 1, -1
         y2(k) = y2(k) * y2(k+1) + u(k)
      enddo
      return
    END SUBROUTINE SPLINE

end module
