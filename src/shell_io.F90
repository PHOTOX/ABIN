!module handeling IO
module mod_shell_io
   use mod_const, only: DP, ANG
   use mod_error, only: fatal_error
   use mod_files, only: stderr
   implicit none
   private
   public :: write_geom, get_shellscript, open_engrad_file, read_energy, read_forces
contains

   subroutine write_geom(fname, at_names, x, y, z, natom, iw)
      character(len=*), intent(in) :: fname
      character(len=2), dimension(:), intent(in) :: at_names
      real(DP), dimension(:, :), intent(in) :: x, y, z
      integer, intent(in) :: natom, iw
      ! Format for geom.dat; needed, so that Molpro can read it
      character(len=*), parameter :: fgeom = '(A2,3E25.17E2)'
      integer :: u, iat

      ! Delete the last geometry
      open (newunit=u, file=fname)
      close (u, status='delete')

      open (newunit=u, file=fname, status='new', action='write', access='sequential')
      do iat = 1, natom
         write (u, fgeom) at_names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
      end do
      close (u)
   end subroutine write_geom

   function get_shellscript(potential) result(shellscript)
      use mod_utils, only: toupper
      character(len=*), intent(in) :: potential
      character(:), allocatable :: shellscript
      logical :: exists

      shellscript = './'//trim(toupper(potential))//'/r.'//potential
      inquire (file=shellscript, exist=exists)
      if (.not. exists) then
         call fatal_error(__FILE__, __LINE__, 'Shell executable '//trim(shellscript)//' does not exist')
      end if
   end function get_shellscript

   integer function open_engrad_file(fname, abort) result(uengrad)
      character(len=*), intent(in) :: fname
      logical, intent(inout) :: abort
      logical :: exists
      integer :: iost

      inquire (file=fname, exist=exists)
      if (.not. exists) then
         write (stderr, '(A)') 'WARNING: File '//trim(fname)//' does not exist. Waiting...'
         ! Should flush HDD buffer
         call execute_command_line('sync')
         call execute_command_line('sleep 0.5')
      end if

      open (newunit=uengrad, file=fname, status='old', action='read', iostat=iost)
      if (iost /= 0) then
         write (stderr, '(A)') 'Could not open file '//trim(fname)//' for reading energy and gradients.'
         abort = .true.
!$OMP FLUSH(abort)
         return
      end if
   end function open_engrad_file

   real(DP) function read_energy(engrad_unit, abort) result(energy)
      integer, intent(in) :: engrad_unit
      logical, intent(inout) :: abort
      integer :: iost
#if __GNUC__ != 7 || __GNUC_MINOR__ >= 5
      character(len=300) :: fname
      logical :: lopened
#endif

      ! Read electronic energy from engrad.dat
      read (engrad_unit, *, iostat=iost) energy
      if (iost /= 0) then
         ! Working around a compiler bug in gfortran 7.3
#if __GNUC__ == 7 && __GNUC_MINOR__ < 5
         write (stderr, '(A)') 'ERROR: Could not read energy'
#else
         inquire (unit=engrad_unit, opened=lopened, name=fname)
         write (stderr, '(A)') 'Could not read energy from file '//trim(fname)
#endif
         abort = .true.
         energy = 0.0D0
!$OMP FLUSH(abort)
         return
      end if
   end function read_energy

   subroutine read_forces(fx, fy, fz, num_atom, iw, engrad_unit, abort)
      use mod_files, only: stderr
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      integer, intent(in) :: iw, engrad_unit, num_atom
      logical, intent(inout) :: abort
#if __GNUC__ != 7 || __GNUC_MINOR__ >= 5
      character(len=300) :: fname
      logical :: lopened
#endif
      integer :: iat, iost

      ! WARNING: The engrad file contains energy gradients, we need to convert to forces.
      do iat = 1, num_atom
         read (engrad_unit, *, iostat=iost) fx(iat, iw), fy(iat, iw), fz(iat, iw)
         if (iost /= 0) then
            ! Working around a compiler bug in gfortran 7.3
#if __GNUC__ == 7 && __GNUC_MINOR__ < 5
            write (stderr, '(A)') 'ERROR: Could not read gradients'
#else
            inquire (unit=engrad_unit, opened=lopened, name=fname)
            write (stderr, '(A)') 'ERROR: Could not read gradients from file '//trim(fname)
#endif
            abort = .true.
!$OMP FLUSH(abort)
            return
         end if
         ! Convert gradients to forces
         fx(iat, iw) = -fx(iat, iw)
         fy(iat, iw) = -fy(iat, iw)
         fz(iat, iw) = -fz(iat, iw)
      end do
   end subroutine read_forces
end module mod_shell_io
