! Various analytical potentials:
! - one particle in harmonic potential well
! - two-particle harmonic rotor
! - two particles bound by Morse potential
! - one particle in a 2D double well potential
!   (originally used for testing PLUMED interface,
!    see MSc thesis of Jirka Suchan)

! Parameters are set in input section 'system'

module mod_potentials
   use mod_const, only: DP
   implicit none
   ! some constansts for analytical potentials
   ! constants particle in 3D harmonic oscillator
   real(DP) :: kx = 0.0D0, ky = 0.0D0, kz = 0.0D0
   ! constants for 2-particle harmonic oscillator
   real(DP) :: k = 0.000D0, r0 = 0.0D0
   ! constants for double well
   real(DP) :: lambda_dw = 0.0D0, D0_dw = 0.0D0, k_dw, r0_dw
   ! Dissociation energy constant for Morse potential
   real(DP) :: dissociation_energy = 1.0D0
   ! Array holding the hessian matrix
   ! TODO: Should probably be elsewhere with other global arrays
   real(DP), allocatable :: hess(:, :, :)
   save
contains

   ! 1 particle in a 3D harmonic potential
   subroutine force_harmonic_oscillator(x, y, z, fxab, fyab, fzab, eclas)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fxab(:, :), fyab(:, :), fzab(:, :)
      real(DP), intent(out) :: eclas
      integer :: iw

      eclas = 0.0D0
      do iw = 1, nwalk
         fxab(1, iw) = -kx * x(1, iw)
         fyab(1, iw) = -ky * y(1, iw)
         fzab(1, iw) = -kz * z(1, iw)
         eclas = eclas + 0.5D0 * (kx * x(1, iw)**2 + &
                                & ky * y(1, iw)**2 + &
                                & kz * z(1, iw)**2)
      end do

      eclas = eclas / nwalk
   end subroutine force_harmonic_oscillator

   ! 2D DOUBLE WELL potential
   ! According to: Tuckerman, Statistical Mechanics, page 350
   subroutine force_doublewell(x, y, fxab, fyab, eclas)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :)
      real(DP), intent(out) :: fxab(:, :), fyab(:, :)
      real(DP), intent(out) :: eclas
      integer :: i

      eclas = 0.0D0

      do i = 1, nwalk
         fxab(1, i) = -4 * D0_dw * x(1, i) * (x(1, i)**2 - r0_dw**2) - lambda_dw * y(1, i)
         fyab(1, i) = -k_dw * y(1, i) - lambda_dw * x(1, i)
         eclas = eclas + D0_dw * (x(1, i)**2 - r0_dw**2)**2 + &
                 0.5 * k_dw * y(1, i)**2 + lambda_dw * x(1, i) * y(1, i)
      end do

      eclas = eclas / nwalk

   end subroutine force_doublewell

   ! Diatomic molecule bound by harmonic potential
   subroutine force_harmonic_rotor(x, y, z, fx, fy, fz, eclas)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      real(DP) :: dx, dy, dz, r, fac
      integer :: i

      eclas = 0.0D0

      do i = 1, nwalk
         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         fac = k * (r - r0) / r
         fx(1, i) = fac * dx
         fx(2, i) = -fx(1, i)
         fy(1, i) = fac * dy
         fy(2, i) = -fy(1, i)
         fz(1, i) = fac * dz
         fz(2, i) = -fz(1, i)
         eclas = eclas + 0.5D0 * k * (r - r0)**2
      end do

      eclas = eclas / nwalk
   end subroutine force_harmonic_rotor

   subroutine hessian_harmonic_rotor(x, y, z)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: dx, dy, dz, r, fac
      integer :: i, ipom1, ipom2

      do i = 1, nwalk

         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         fac = k * (r - r0) / r
         hess(1, 1, i) = (k * dx**2 / r**2 - fac * dx**2 / r**2 + fac) / nwalk
         hess(2, 2, i) = (k * dy**2 / r**2 - fac * dy**2 / r**2 + fac) / nwalk
         hess(3, 3, i) = (k * dz**2 / r**2 - fac * dz**2 / r**2 + fac) / nwalk
         hess(4, 4, i) = hess(1, 1, i)
         hess(5, 5, i) = hess(2, 2, i)
         hess(6, 6, i) = hess(3, 3, i)
         hess(2, 1, i) = (k * dx * dy / r**2 - fac * dx * dy / r**2) / nwalk
         hess(3, 1, i) = (k * dz * dx / r**2 - fac * dz * dx / r**2) / nwalk
         hess(3, 2, i) = (k * dz * dy / r**2 - fac * dz * dy / r**2) / nwalk
         hess(4, 1, i) = -hess(1, 1, i)
         hess(4, 2, i) = -hess(2, 1, i)
         hess(4, 3, i) = -hess(3, 1, i)
         hess(5, 1, i) = -hess(2, 1, i)
         hess(5, 2, i) = -hess(2, 2, i)
         hess(5, 3, i) = -hess(3, 2, i)
         hess(5, 4, i) = hess(2, 1, i)
         hess(6, 1, i) = -hess(3, 1, i)
         hess(6, 2, i) = -hess(3, 2, i)
         hess(6, 3, i) = -hess(3, 3, i)
         hess(6, 4, i) = hess(3, 1, i)
         hess(6, 5, i) = hess(3, 2, i)

         do ipom1 = 1, 5
            do ipom2 = ipom1 + 1, 6
               hess(ipom1, ipom2, i) = hess(ipom2, ipom1, i)
            end do
         end do

      end do
   end subroutine hessian_harmonic_rotor

   ! Morse potential between two particles in 3D
   ! V = DE * (1 - exp(-a*(r-r0)) )^2
   ! where DE is dissociation energy and
   ! parameter 'a' is linked to the harmonic constant as
   ! a = sqrt( k / 2 / de)
   ! Here we actually use the harmonic constant as an input parameter.
   subroutine force_morse(x, y, z, fx, fy, fz, eclas)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      real(DP) :: dx, dy, dz, r, fac, ex
      real(DP) :: a
      integer :: i

      a = dsqrt(k / 2.0D0 / dissociation_energy)

      eclas = 0.0D0

      do i = 1, nwalk
         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         ex = exp(-a * (r - r0))
         fac = 2 * a * ex * dissociation_energy * (1 - ex) / r
         fx(1, i) = fac * dx
         fx(2, i) = -fx(1, i)
         fy(1, i) = fac * dy
         fy(2, i) = -fy(1, i)
         fz(1, i) = fac * dz
         fz(2, i) = -fz(1, i)
         eclas = eclas + dissociation_energy * (1 - ex)**2
      end do

      eclas = eclas / nwalk

   end subroutine force_morse

   subroutine hessian_morse(x, y, z)
      use mod_general, only: nwalk
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: dx, dy, dz, r, fac, ex, fac2
      real(DP) :: a
      integer :: i, ipom1, ipom2

      a = dsqrt(k / 2 / dissociation_energy)

      do i = 1, nwalk

         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         ex = exp(-a * (r - r0))
         fac = 2 * a * ex * dissociation_energy * (1 - ex) / r
         fac2 = 2 * dissociation_energy * a**2 * ex**2 / r**2
         hess(1, 1, i) = fac2 * dx**2 - (fac * dx**2) / r**2 + fac - fac * a * dx**2 / r
         hess(2, 2, i) = fac2 * dy**2 - (fac * dy**2) / r**2 + fac - fac * a * dy**2 / r
         hess(3, 3, i) = fac2 * dz**2 - (fac * dz**2) / r**2 + fac - fac * a * dz**2 / r
         hess(1, 1, i) = hess(1, 1, i) / nwalk
         hess(2, 2, i) = hess(2, 2, i) / nwalk
         hess(3, 3, i) = hess(3, 3, i) / nwalk
         hess(4, 4, i) = hess(1, 1, i)
         hess(5, 5, i) = hess(2, 2, i)
         hess(6, 6, i) = hess(3, 3, i)
         hess(2, 1, i) = fac2 * dx * dy - (fac * dx * dy) / r**2 - fac * a * dx * dy / r
         hess(3, 1, i) = fac2 * dx * dz - (fac * dx * dz) / r**2 - fac * a * dx * dz / r
         hess(3, 2, i) = fac2 * dz * dy - (fac * dz * dy) / r**2 - fac * a * dz * dy / r
         hess(1, 1, i) = hess(2, 1, i) / nwalk
         hess(2, 2, i) = hess(3, 1, i) / nwalk
         hess(3, 3, i) = hess(3, 2, i) / nwalk
         hess(4, 1, i) = -hess(1, 1, i)
         hess(4, 2, i) = -hess(2, 1, i)
         hess(4, 3, i) = -hess(3, 1, i)
         hess(5, 1, i) = -hess(2, 1, i)
         hess(5, 2, i) = -hess(2, 2, i)
         hess(5, 3, i) = -hess(3, 2, i)
         hess(5, 4, i) = hess(2, 1, i)
         hess(6, 1, i) = -hess(3, 1, i)
         hess(6, 2, i) = -hess(3, 2, i)
         hess(6, 3, i) = -hess(3, 3, i)
         hess(6, 4, i) = hess(3, 1, i)
         hess(6, 5, i) = hess(3, 2, i)
         do ipom1 = 1, 5
            do ipom2 = ipom1 + 1, 6
               hess(ipom1, ipom2, i) = hess(ipom2, ipom1, i)
            end do
         end do

      end do

   end subroutine hessian_morse

   ! This hessian is constant, so theoretically only needs
   ! to be called once, but we're callling it each step anyway for now.
   subroutine hessian_harmonic_oscillator()
      use mod_general, only: nwalk
      integer :: i

      do i = 1, nwalk
         hess(1, 1, i) = kx / nwalk
         hess(2, 2, i) = ky / nwalk
         hess(3, 3, i) = kz / nwalk
         hess(2, 1, i) = 0.0D0
         hess(3, 1, i) = 0.0D0
         hess(3, 2, i) = 0.0D0
         hess(1, 2, i) = 0.0D0
         hess(1, 3, i) = 0.0D0
         hess(2, 3, i) = 0.0D0
      end do
   end subroutine hessian_harmonic_oscillator

end module mod_potentials

! 1D numerical potential with cubic splines
module mod_splined_grid
   use mod_const, only: DP
   implicit none
   integer, parameter :: MAX_GRID_SIZE = 10000
   integer :: grid_size
   real(DP), allocatable :: second_derivatives(:)
   real(DP) :: x_grid(MAX_GRID_SIZE)
   real(DP) :: y_grid(MAX_GRID_SIZE)
   real(DP) :: x_max, x_min
   private
   public :: force_splined_grid, initialize_spline
   save
contains

   real(DP) function potential_cubic_spline(X)
      real(DP), intent(IN) :: X
      real(DP) :: Y
      call splint(x_grid, y_grid, second_derivatives, grid_size, X, Y)
      potential_cubic_spline = Y
      return
   end function potential_cubic_spline

   ! TODO: Remove y and z components from arguments, since we
   ! only support 1D spline at the moment
   subroutine force_splined_grid(x, fx, eclas)
      use mod_general, only: nwalk
      use mod_utils, only: abinerror
      real(DP), intent(in) :: x(:, :)
      real(DP), intent(out) :: fx(:, :)
      real(DP), intent(out) :: eclas
      real(DP) :: en_1, en_2, dx = 0.0001D0
      integer :: iw

      eclas = 0.0D0
      do iw = 1, nwalk

         if (x(1, iw) < x_min .or. x(1, iw) > x_max) then
            write (*, *) "ERROR: Particle got out of the grid!"
            call abinerror("force_splined_grid")
         end if

         eclas = eclas + potential_cubic_spline(x(1, iw))

         ! Let's just do numerical forces for now!
         ! use 3-point numerical derivative
         en_1 = potential_cubic_spline(x(1, iw) - dx)
         en_2 = potential_cubic_spline(x(1, iw) + dx)
         fx(1, iw) = (en_1 - en_2) / 2 / dx

      end do

   end subroutine force_splined_grid

   subroutine initialize_spline()
      use mod_general, only: idebug
      implicit none
      character(len=100) :: chfile_pot = 'potential.dat'
      real(DP) :: yp1, ypn ! first derivatives at the grid edges
      integer :: IUNIT = 600, i
      real(DP) :: x, dx

      grid_size = 0

      ! TODO: make x_grid et al allocatable
      open (IUNIT, file=chfile_pot, status='OLD', action='read')
      do
         read (IUNIT, *, end=51) x_grid(grid_size + 1), y_grid(grid_size + 1)
         grid_size = grid_size + 1
      end do
51    close (IUNIT)

      if (grid_size < 3) then
         write (*, *) "ERROR: Could not find grid in "//chfile_pot
         stop 1
      end if

      allocate (second_derivatives(grid_size))
      x_min = x_grid(1)
      x_max = x_grid(grid_size)

!     Settings for natural cubic spline
      yp1 = 1E31
      ypn = 1E31
!     On second thought, let's just set the first derivatives the same
!     as at the boundaries
!      yp1 = (xa(2)-xa(1)) / (ya(1)-ya(2))
!      ypn = 1e31
      call spline(x_grid, y_grid, grid_size, yp1, ypn, second_derivatives)

      ! Print out the splined potential just to be sure
      if (idebug == 1) then
         x = x_min
         dx = (x_max - x_min) / grid_size / 5
         open (IUNIT, file="potential_splined.dat", action="write")
!         write(IUNIT, *)x_min, x_max, grid_size
         do i = 1, grid_size * 5
            write (IUNIT, *) x, potential_cubic_spline(x)
            x = x + dx
         end do
         close (IUNIT)
      end if
   end subroutine initialize_spline

   ! From numerical recipies, sliglty modified
   subroutine SPLINT(XA, YA, Y2A, N, X, Y)
      implicit real * 8(a - h, o - z)
      integer, intent(IN) :: N
      real(DP), intent(IN) :: X, XA(N), YA(N), Y2A(N)
      real(DP), intent(OUT) :: Y
      integer :: KLO, KHI, K
      KLO = 1
      KHI = N
      do while (KHI - KLO > 1)
         K = (KHI + KLO) / 2
         if (XA(K) > X) then
            KHI = K
         else
            KLO = K
         end if
      end do
      H = XA(KHI) - XA(KLO)
      if (H == 0.) then
         write (*, *) 'Bad XA input.'
         stop 1
      end if
!     DH: when extrapolating out of boundaries, just do linear
!     TODO: do this also to the END boundary
!      if(x.lt.xa(1))then
!         Y2A(KHI) = 0
!         Y2A(KLO) = 0
!      end if
      A = (XA(KHI) - X) / H
      B = (X - XA(KLO)) / H
      Y = A * YA(KLO) + B * YA(KHI) + &
          ((A**3 - A) * Y2A(KLO) + (B**3 - B) * Y2A(KHI)) * (H**2) / 6.

      ! TODO implement first derivative
      return
   end subroutine SPLINT

   subroutine spline(x, y, n, yp1, ypn, y2)
      integer :: n, NMAX
      real(DP) :: yp1, ypn, x(n), y(n), y2(n)
      parameter(NMAX=500)
      integer i, k
      real(DP) p, qn, sig, un, u(NMAX)

      if (yp1 > .99E30) then
         y2(1) = 0.
         u(1) = 0.
      else
         y2(1) = -0.5
         u(1) = (3./(x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
      end if
      do i = 2, n - 1
         sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
         p = sig * y2(i - 1) + 2.
         y2(i) = (sig - 1.) / p
         u(i) = (6.*((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1)) &
                     / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
      end do
      if (ypn > .99E30) then
         qn = 0.
         un = 0.
      else
         qn = 0.5
         un = (3./(x(n) - x(n - 1))) * (ypn - (y(n) - y(n - 1)) / (x(n) - x(n - 1)))
      end if
      y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1.)
      do k = n - 1, 1, -1
         y2(k) = y2(k) * y2(k + 1) + u(k)
      end do
      return
   end subroutine SPLINE

end module
