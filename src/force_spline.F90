! 1D user-defined numerical potential with cubic splines
! Original code taken from QDYN
! https://github.com/PHOTOX/qdyn
! TODO: This module is currently disabled until the cubic spline routines are re-implemented
module mod_splined_grid
   use mod_const, only: DP
   use mod_error, only: fatal_error
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
      call fatal_error(__FILE__, __LINE__, "Spline potential not implemented")
      ! To squash compiler warnings
      y = x
   end function potential_cubic_spline

   subroutine force_splined_grid(x, fx, eclas, walkmax)
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
      integer, intent(in) :: natom

      call fatal_error(__FILE__, __LINE__, "Spline potential not implemented")

      if (natom /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Splined grid potential is only for 1 particle')
      end if
      dime = 1
      f = 0

      call read_grid(potential_file, x_grid, y_grid, grid_size)

      call validate_grid(x_grid, grid_size)

      allocate (second_derivatives(grid_size))

      ! TODO: Initialize spline potential

      call print_splined_potential("potential_splined.dat", x_grid, grid_size)
   end subroutine initialize_spline

   subroutine finalize_spline()
      if (allocated(second_derivatives)) then
         deallocate (second_derivatives)
      end if
   end subroutine finalize_spline

   subroutine read_grid(fname, x_grid, y_grid, grid_size)
      use mod_files, only: stdout
      character(len=*), intent(in) :: fname
      real(DP), dimension(:), intent(out) :: x_grid, y_grid
      integer, intent(out) :: grid_size
      integer :: u, iost

      grid_size = 0

      write (stdout, *) ''
      write (stdout, *) 'Reading numerical potential grid from file '//trim(fname)
      write (stdout, *) 'First column: x-coordinate / bohrs'
      write (stdout, *) 'Second column: potential energy / atomic units'

      open (newunit=u, file=fname, status='old', action='read')
      do
         read (u, *, iostat=iost) x_grid(grid_size + 1), y_grid(grid_size + 1)
         if (iost < 0) then
            exit
         end if
         if (iost > 0) then
            close (u)
            call fatal_error(__FILE__, __LINE__, "Invalid line in file "//fname)
            return
         end if
         grid_size = grid_size + 1
         if (grid_size == MAX_GRID_SIZE) then
            close (u)
            call fatal_error(__FILE__, __LINE__, "Maximum grid points exceeded")
         end if
      end do
      close (u)

   end subroutine read_grid

   subroutine validate_grid(x_grid, ngrid)
      real(DP), dimension(ngrid), intent(in) :: x_grid
      integer, intent(in) :: ngrid
      integer :: i

      if (ngrid < 3) then
         call fatal_error(__FILE__, __LINE__, "Grid must contain at least 3 points")
         return
      end if

      do i = 1, ngrid - 1
         if (x_grid(i + 1) <= x_grid(i)) then
            call fatal_error(__FILE__, __LINE__, &
               & "grid x values do not increase monotonically")
            return
         end if
      end do
   end subroutine validate_grid

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
         write (u, '(G0,X,G0)') x, potential_cubic_spline(x)
         x = x + dx
      end do
      close (u)
   end subroutine print_splined_potential

end module mod_splined_grid
