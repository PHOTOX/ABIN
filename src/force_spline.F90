! 1D user-defined numerical potential with cubic splines
module mod_splined_grid
   use mod_const, only: DP
   implicit none
   private
   integer, parameter :: MAX_GRID_SIZE = 10000
   integer :: grid_size
   real(DP), allocatable :: second_derivatives(:)
   real(DP) :: x_grid(MAX_GRID_SIZE)
   real(DP) :: y_grid(MAX_GRID_SIZE)
   real(DP) :: x_max, x_min
   character(len=256) :: potential_file = 'potential.dat'

   public :: potential_file
   public :: force_splined_grid, initialize_spline
   save
contains

   real(DP) function potential_cubic_spline(x) result(y)
      real(DP), intent(IN) :: X
      call splint(x_grid, y_grid, second_derivatives, grid_size, x, y)
   end function potential_cubic_spline

   subroutine force_splined_grid(x, fx, eclas, walkmax)
      use mod_error, only: fatal_error
      real(DP), intent(in) :: x(:, :)
      real(DP), intent(out) :: fx(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: walkmax
      real(DP) :: en_1, en_2, dx = 0.0001D0
      integer :: iw

      eclas = 0.0D0
      do iw = 1, walkmax

         if (x(1, iw) < x_min .or. x(1, iw) > x_max) then
            call fatal_error(__FILE__, __LINE__, "Particle got out of the grid!")
         end if

         eclas = eclas + potential_cubic_spline(x(1, iw))

         ! Let's just do numerical forces for now!
         ! use 3-point numerical derivative
         en_1 = potential_cubic_spline(x(1, iw) - dx)
         en_2 = potential_cubic_spline(x(1, iw) + dx)
         fx(1, iw) = (en_1 - en_2) / 2 / dx

      end do

      eclas = eclas / walkmax

   end subroutine force_splined_grid

   subroutine initialize_spline()
      use mod_general, only: idebug
      real(DP) :: yp1, ypn ! first derivatives at the grid edges
      integer :: IUNIT, i
      real(DP) :: x, dx

      grid_size = 0

      open (newunit=IUNIT, file=potential_file, status='old', action='read')
      do
         read (IUNIT, *, end=51) x_grid(grid_size + 1), y_grid(grid_size + 1)
         grid_size = grid_size + 1
      end do
51    close (IUNIT)

      if (grid_size < 3) then
         write (*, *) "ERROR: Could not find grid in "//potential_file
         stop 1
      end if

      allocate (second_derivatives(grid_size))
      x_min = x_grid(1)
      x_max = x_grid(grid_size)

      ! Settings for natural cubic spline
      yp1 = 1E31
      ypn = 1E31
      ! On second thought, let's just set the first derivatives the same
      ! as at the boundaries
      !  yp1 = (xa(2)-xa(1)) / (ya(1)-ya(2))
      !  ypn = 1e31
      call spline(x_grid, y_grid, grid_size, yp1, ypn, second_derivatives)

      ! Print out the splined potential just to be sure
      if (idebug == 1) then
         x = x_min
         dx = (x_max - x_min) / grid_size / 5
         open (newunit=IUNIT, file="potential_splined.dat", action="write")
         do i = 1, grid_size * 5
            write (IUNIT, *) x, potential_cubic_spline(x)
            x = x + dx
         end do
         close (IUNIT)
      end if
   end subroutine initialize_spline

   ! From numerical recipies, sliglty modified
   ! TODO: implement first derivative
   subroutine splint(xa, ya, y2a, n, x, y)
      use mod_error, only: fatal_error
      integer, intent(in) :: n
      real(DP), intent(in) :: x, xa(n), ya(n), y2a(n)
      real(DP), intent(out) :: y
      real(DP) :: a, b, h
      integer :: klo, khi, k
      klo = 1
      khi = n
      do while (khi - klo > 1)
         k = (khi + klo) / 2
         if (xa(k) > x) then
            khi = k
         else
            klo = k
         end if
      end do
      h = xa(khi) - xa(klo)
      if (h == 0.0D0) then
         call fatal_error(__FILE__, __LINE__, "Bad xa input")
      end if

      a = (xa(khi) - x) / h
      b = (x - xa(klo)) / h
      y = a * ya(klo) + b * ya(khi) + &
          ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.
   end subroutine splint

   subroutine spline(x, y, n, yp1, ypn, y2)
      integer, intent(in) :: n
      integer, parameter :: NMAX = 500
      real(DP) :: yp1, ypn, x(n), y(n), y2(n)
      real(DP) :: p, qn, sig, un, u(NMAX)
      integer :: i, k

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
   end subroutine spline

end module
