! 1D user-defined numerical potential with cubic splines
! Original code taken from QDYN
! https://github.com/PHOTOX/qdyn
module mod_splined_grid
   use mod_const, only: DP
   implicit none
   private
   integer, parameter :: MAX_GRID_SIZE = 1000
   integer :: grid_size
   real(DP), allocatable :: second_derivatives(:)
   real(DP) :: x_grid(MAX_GRID_SIZE)
   real(DP) :: y_grid(MAX_GRID_SIZE)

   character(len=256) :: potential_file = 'potential.dat'

   public :: potential_file
   public :: force_splined_grid
   public :: initialize_spline, finalize_spline
   save
contains

   real(DP) function potential_cubic_spline(x) result(y)
      real(DP), intent(in) :: x
      call splint(x_grid, y_grid, second_derivatives, grid_size, x, y)
   end function potential_cubic_spline

   subroutine force_splined_grid(x, fx, eclas, walkmax)
      use mod_error, only: fatal_error
      real(DP), intent(in) :: x(:, :)
      real(DP), intent(out) :: fx(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: walkmax
      ! Displacement for numerical forces
      real(DP), parameter :: DX = 0.0001D0
      real(DP) :: en_1, en_2
      integer :: iw

      eclas = 0.0D0
      do iw = 1, walkmax

         if (x(1, iw) < x_grid(1) .or. x(1, iw) > x_grid(grid_size)) then
            call fatal_error(__FILE__, __LINE__, "Particle got out of the grid!")
         end if

         eclas = eclas + potential_cubic_spline(x(1, iw))

         ! Let's just do numerical forces for now!
         ! use 3-point numerical derivative
         en_1 = potential_cubic_spline(x(1, iw) - DX)
         en_2 = potential_cubic_spline(x(1, iw) + DX)
         fx(1, iw) = (en_1 - en_2) / 2 / DX

      end do

      eclas = eclas / walkmax

   end subroutine force_splined_grid

   subroutine initialize_spline(natom)
      use mod_system, only: dime, f
      use mod_error, only: fatal_error
      integer, intent(in) :: natom
      real(DP) :: yp1, ypn ! first derivatives at the grid edges

      if (natom /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Splined grid potential is only for 1 particle')
      end if
      dime = 1
      f = 0

      call read_grid(potential_file, x_grid, y_grid, grid_size)

      allocate (second_derivatives(grid_size))

      yp1 = 1E31
      ypn = 1E31
      call spline(x_grid, y_grid, grid_size, yp1, ypn, second_derivatives)

      call print_splined_potential("potential_splined.dat", x_grid, grid_size)
   end subroutine initialize_spline

   subroutine finalize_spline()
      if (allocated(second_derivatives)) then
         deallocate (second_derivatives)
      end if
   end subroutine

   subroutine read_grid(fname, x_grid, y_grid, grid_size)
      use mod_error, only: fatal_error
      use mod_general, only: my_rank
      character(len=*), intent(in) :: fname
      real(DP), dimension(:), intent(out) :: x_grid, y_grid
      integer, intent(out) :: grid_size
      integer :: u, iost

      grid_size = 0
      if (my_rank == 0) then
         print*,''
         print*,'Reading numerical potential grid from file '//trim(fname)
         print*,'First column: x-coordinate / bohrs'
         print*,'Second column: potential energy / atomic units'
      end if

      open (newunit=u, file=fname, status='old', action='read')
      do
         read (u, *, iostat=iost) x_grid(grid_size + 1), y_grid(grid_size + 1)
         if (iost < 0) then
            exit
         end if
         if (iost > 0) then
            call fatal_error(__FILE__, __LINE__, "Invalid line in file "//fname)
         end if
         grid_size = grid_size + 1
         if (grid_size == MAX_GRID_SIZE) then
            call fatal_error(__FILE__, __LINE__, "Maximum grid points exceeded")
         end if
      end do
      close (u)

      if (grid_size < 3) then
         call fatal_error(__FILE__, __LINE__, "Invalid potential grid in file "//fname)
      end if
   end subroutine

   subroutine print_splined_potential(fname, x_grid, grid_size)
      character(len=*), intent(in) :: fname
      real(DP), intent(in) :: x_grid(:)
      integer, intent(in) :: grid_size
      real(DP) :: x, dx
      integer :: u, i

      x = x_grid(1)
      dx = (x_grid(grid_size) - x_grid(1)) / grid_size / 5
      open (newunit=u, file=fname, action="write")
      do i = 1, grid_size * 5
         write (u, *) x, potential_cubic_spline(x)
         x = x + dx
      end do
      close (u)
   end subroutine

   ! From numerical recipies, slightly modified
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
