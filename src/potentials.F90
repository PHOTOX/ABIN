! Various analytical potentials:
! - one particle in harmonic potential well
! - two-particle harmonic rotor
! - two particles bound by Morse potential
! - one particle in a double well potential
!   linearly coupled to harmonic oscillator
!   (originally used for testing PLUMED interface,
!    see MSc thesis of Jirka Suchan)

! User-definer parameters are set input section 'system'
module mod_potentials
   use mod_const, only: DP
   implicit none
   public
   private :: dw, morse, ho, hrot

   ! We use derived types to encapsulate potential parameters
   ! https://fortran-lang.org/learn/quickstart/derived_types

   ! Parameters for a double-well potential in x dimension
   ! linearly coupled to harmonic potential in y dimension
   ! See Mark E. Tuckerman, Statistical Mechanics:
   ! Theory and Molecular simulation,
   ! Chapter 8: Free energy calculations,
   ! Section 8.10, Adiabatic dynamics, p. 353, eq. 8.10.21
   ! V(x,y) = D0*(x^2 - r0^2)^2 + 1/2*k*y^2 + lambda*x*y
   type :: dw_params
      real(DP) :: lambda
      real(DP) :: d0
      real(DP) :: k
      real(DP) :: r0
   end type
   type(dw_params) :: dw

   ! Parameters for two particles bound by Morse potential
   ! V = D0 * (1 - exp(-a*(r-r0)) )^2
   type :: morse_params
      ! Dissociation energy constant
      real(DP) :: d0
      real(DP) :: a
      real(DP) :: r0
   end type
   type(morse_params) :: morse

   ! 3D harmonic oscillator
   ! E = 1/2 * (kx * x^2 + ky * y^2 + kz * z^2)
   type :: harm_osc_params
      real(DP) :: kx, ky, kz
   end type
   type(harm_osc_params) :: ho

   ! Two particles bound by harmonic potential
   ! E = 1/2 * k * (r - r0)^2
   type :: harm_rotor_params
      real(DP) :: k, r0
   end type
   type(harm_rotor_params) :: hrot
   save
contains

   ! 1 particle in a 3D harmonic potential
   subroutine harmonic_oscillator_init(natom, kx, ky, kz, vx, vy, vz)
      use mod_error, only: fatal_error
      use mod_general, only: irest
      use mod_utils, only: real_nonnegative
      use mod_system, only: dime, f
      integer, intent(in) :: natom
      real(DP), intent(in) :: kx, ky, kz
      real(DP), dimension(:, :), intent(inout) :: vx, vy, vz

      if (natom /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Harmonic potential is only for 1 particle')
      end if

      call real_nonnegative(kx, 'kx')
      call real_nonnegative(ky, 'ky')
      call real_nonnegative(kz, 'kz')
      ho = harm_osc_params(kx=kx, ky=ky, kz=kz)

      ! Zero-out velocity in inactive dimensions
      if (kx == 0.0D0 .and. irest == 0) then
         vx = 0.0D0
      end if
      if (ky == 0.0D0 .and. irest == 0) then
         vy = 0.0D0
      end if
      if (kz == 0.0D0 .and. irest == 0) then
         vz = 0.0D0
      end if

      f = 0
      dime = 0
      ! Determine number of dimensions
      if (kx /= 0.0D0) then
         dime = dime + 1
      end if
      if (ky /= 0.0D0) then
         dime = dime + 1
      end if
      if (kz /= 0.0D0) then
         dime = dime + 1
      end if
      if (dime == 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'At least one of kx, ky or kz must be non-zero')
      end if
   end subroutine harmonic_oscillator_init

   subroutine force_harmonic_oscillator(x, y, z, fx, fy, fz, eclas, nwalk)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: nwalk
      integer :: iw

      eclas = 0.0D0
      do iw = 1, nwalk
         fx(1, iw) = -ho%kx * x(1, iw)
         fy(1, iw) = -ho%ky * y(1, iw)
         fz(1, iw) = -ho%kz * z(1, iw)
         eclas = eclas + 0.5D0 * (ho%kx * x(1, iw)**2 + &
                                & ho%ky * y(1, iw)**2 + &
                                & ho%kz * z(1, iw)**2)
      end do

      eclas = eclas / nwalk
   end subroutine force_harmonic_oscillator

   ! This hessian is constant, so theoretically only needs
   ! to be called once, but we're callling it each step anyway for now.
   subroutine hessian_harmonic_oscillator(nwalk, hess)
      integer, intent(in) :: nwalk
      real(DP), dimension(:, :, :), intent(out) :: hess
      integer :: i

      do i = 1, nwalk
         hess(1, 1, i) = ho%kx / nwalk
         hess(2, 2, i) = ho%ky / nwalk
         hess(3, 3, i) = ho%kz / nwalk
         hess(2, 1, i) = 0.0D0
         hess(3, 1, i) = 0.0D0
         hess(3, 2, i) = 0.0D0
         hess(1, 2, i) = 0.0D0
         hess(1, 3, i) = 0.0D0
         hess(2, 3, i) = 0.0D0
      end do
   end subroutine hessian_harmonic_oscillator

   subroutine doublewell_init(natom, lambda, d0, k, r0, vy, vz)
      use mod_error, only: fatal_error
      use mod_utils, only: real_nonnegative
      use mod_system, only: dime, f
      integer, intent(in) :: natom
      real(DP), intent(in) :: lambda
      real(DP), intent(in) :: d0
      real(DP), intent(in) :: k
      real(DP), intent(in) :: r0
      real(DP), dimension(:, :), intent(inout) :: vy, vz

      if (natom /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Double-well potential is only for 1 particle')
      end if

      call real_nonnegative(k, 'k_dw')
      call real_nonnegative(d0, 'd0_dw')

      dw = dw_params(lambda=lambda, k=k, d0=d0, r0=r0)

      ! This particular double-well is 2D, unless we make it 1D
      dime = 2
      vz = 0.0D0
      if (dw%lambda == 0.0D0 .and. dw%k == 0.0D0) then
         dime = 1
         vy = 0.0D0
      end if
      f = 0
   end subroutine doublewell_init

   subroutine force_doublewell(x, y, fx, fy, eclas, nwalk)
      real(DP), intent(in) :: x(:, :), y(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: nwalk
      integer :: i

      eclas = 0.0D0

      do i = 1, nwalk
         fx(1, i) = -4 * dw%d0 * x(1, i) * (x(1, i)**2 - dw%r0**2) - dw%lambda * y(1, i)
         fy(1, i) = -dw%k * y(1, i) - dw%lambda * x(1, i)
         eclas = eclas + dw%d0 * (x(1, i)**2 - dw%r0**2)**2 + &
                 0.5D0 * dw%k * y(1, i)**2 + dw%lambda * x(1, i) * y(1, i)
      end do

      eclas = eclas / nwalk

   end subroutine force_doublewell

   ! Diatomic molecule bound by harmonic potential
   subroutine harmonic_rotor_init(natom, k, r0)
      use mod_error, only: fatal_error
      use mod_utils, only: real_positive
      integer, intent(in) :: natom
      real(DP), intent(in) :: k, r0

      if (natom /= 2) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Harmonic rotor is only for 2 particles')
      end if

      call real_positive(k, 'k')
      call real_positive(r0, 'r0')

      hrot = harm_rotor_params(k=k, r0=r0)
   end subroutine harmonic_rotor_init

   subroutine force_harmonic_rotor(x, y, z, fx, fy, fz, eclas, nwalk)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: nwalk
      real(DP) :: dx, dy, dz, r, fac
      integer :: i

      eclas = 0.0D0

      do i = 1, nwalk
         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         fac = hrot%k * (r - hrot%r0) / r
         fx(1, i) = fac * dx
         fx(2, i) = -fx(1, i)
         fy(1, i) = fac * dy
         fy(2, i) = -fy(1, i)
         fz(1, i) = fac * dz
         fz(2, i) = -fz(1, i)
         eclas = eclas + 0.5D0 * hrot%k * (r - hrot%r0)**2
      end do

      eclas = eclas / nwalk
   end subroutine force_harmonic_rotor

   subroutine hessian_harmonic_rotor(x, y, z, nwalk, hess)
      real(DP), dimension(:, :), intent(in) :: x, y, z
      integer, intent(in) :: nwalk
      real(DP), intent(out) :: hess(:, :, :)
      real(DP) :: dx, dy, dz, r, fac
      integer :: i, ipom1, ipom2

      do i = 1, nwalk

         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         fac = hrot%k * (r - hrot%r0) / r
         hess(1, 1, i) = (hrot%k * dx**2 / r**2 - fac * dx**2 / r**2 + fac) / nwalk
         hess(2, 2, i) = (hrot%k * dy**2 / r**2 - fac * dy**2 / r**2 + fac) / nwalk
         hess(3, 3, i) = (hrot%k * dz**2 / r**2 - fac * dz**2 / r**2 + fac) / nwalk
         hess(4, 4, i) = hess(1, 1, i)
         hess(5, 5, i) = hess(2, 2, i)
         hess(6, 6, i) = hess(3, 3, i)
         hess(2, 1, i) = (hrot%k * dx * dy / r**2 - fac * dx * dy / r**2) / nwalk
         hess(3, 1, i) = (hrot%k * dz * dx / r**2 - fac * dz * dx / r**2) / nwalk
         hess(3, 2, i) = (hrot%k * dz * dy / r**2 - fac * dz * dy / r**2) / nwalk
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
   subroutine morse_init(natom, k_morse, r0, d0)
      use mod_error, only: fatal_error
      use mod_utils, only: real_positive
      integer, intent(in) :: natom
      real(DP), intent(in) :: k_morse, r0, d0
      real(DP) :: a

      if (natom /= 2) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Morse potential is only for 2 particles')
      end if

      call real_positive(k_morse, 'k_morse')
      call real_positive(d0, 'd0_morse')
      call real_positive(r0, 'r0_morse')

      a = dsqrt(k_morse / 2.0D0 / d0)
      morse = morse_params(d0=d0, a=a, r0=r0)
   end subroutine morse_init

   subroutine force_morse(x, y, z, fx, fy, fz, eclas, nwalk)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: nwalk
      real(DP) :: dx, dy, dz, r, fac, ex
      integer :: i

      eclas = 0.0D0

      do i = 1, nwalk
         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         ex = exp(-morse%a * (r - morse%r0))
         fac = 2 * morse%a * ex * morse%d0 * (1 - ex) / r
         fx(1, i) = fac * dx
         fx(2, i) = -fx(1, i)
         fy(1, i) = fac * dy
         fy(2, i) = -fy(1, i)
         fz(1, i) = fac * dz
         fz(2, i) = -fz(1, i)
         eclas = eclas + morse%d0 * (1 - ex)**2
      end do

      eclas = eclas / nwalk

   end subroutine force_morse

   subroutine hessian_morse(x, y, z, nwalk, hess)
      real(DP), dimension(:, :), intent(in) :: x, y, z
      integer, intent(in) :: nwalk
      real(DP), dimension(:, :, :), intent(out) :: hess
      real(DP) :: dx, dy, dz, r, fac, ex, fac2
      real(DP) :: a
      integer :: i, ipom1, ipom2

      a = morse%a

      do i = 1, nwalk

         dx = x(2, i) - x(1, i)
         dy = y(2, i) - y(1, i)
         dz = z(2, i) - z(1, i)
         r = dx**2 + dy**2 + dz**2
         r = dsqrt(r)
         ex = exp(-a * (r - morse%r0))
         fac = 2 * a * ex * morse%d0 * (1 - ex) / r
         fac2 = 2 * morse%d0 * a**2 * ex**2 / r**2
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
         hess(2, 1, i) = hess(2, 1, i) / nwalk
         hess(3, 1, i) = hess(3, 1, i) / nwalk
         hess(3, 2, i) = hess(3, 2, i) / nwalk
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

end module mod_potentials
