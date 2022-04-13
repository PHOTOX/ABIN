! Various helper functions and subroutines that may be used throughout the program.
module mod_utils
   use mod_const, only: DP
   use mod_interfaces
   use mod_error, only: fatal_error
   implicit none
   public
contains

   real(DP) function get_distance(x, y, z, at1, at2, iw) result(r)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: at1, at2, iw
      character(len=*), parameter :: error_msg = 'Atom indices in get_distance() must be unique'

      if (at1 == at2) then
         r = 0
         call fatal_error(__FILE__, __LINE__, error_msg)
         return
      end if

      r = (x(at1, iw) - x(at2, iw))**2
      r = r + (y(at1, iw) - y(at2, iw))**2
      r = r + (z(at1, iw) - z(at2, iw))**2
      r = dsqrt(r)
   end function get_distance

   real(DP) function get_angle(x, y, z, at1, at2, at3, iw) result(angle)
      use mod_const, only: PI
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw
      real(DP) :: vec1x, vec1y, vec1z
      real(DP) :: vec2x, vec2y, vec2z
      integer :: at1, at2, at3
      character(len=*), parameter :: error_msg = 'Atom indices in get_angle() must be unique'

      if (at1 == at2 .or. at1 == at3 .or. at2 == at3) then
         angle = 0
         call fatal_error(__FILE__, __LINE__, error_msg)
         return
      end if

      vec1x = x(at1, iw) - x(at2, iw)
      vec1y = y(at1, iw) - y(at2, iw)
      vec1z = z(at1, iw) - z(at2, iw)
      vec2x = x(at3, iw) - x(at2, iw)
      vec2y = y(at3, iw) - y(at2, iw)
      vec2z = z(at3, iw) - z(at2, iw)
      angle = 180 / pi * acos((vec1x * vec2x + vec1y * vec2y + vec1z * vec2z) / &
               & (dsqrt(vec1x**2 + vec1y**2 + vec1z**2) * dsqrt(vec2x**2 + vec2y**2 + vec2z**2)))
   end function get_angle

   real(DP) function get_dihedral(x, y, z, at1, at2, at3, at4, iw, shiftdih)
      use mod_const, only: PI
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :), shiftdih
      integer, intent(in) :: iw
      real(DP) :: vec1x, vec1y, vec1z
      real(DP) :: vec2x, vec2y, vec2z
      real(DP) :: vec3x, vec3y, vec3z, sign
      real(DP) :: norm1x, norm1y, norm1z, norm2x, norm2y, norm2z
      integer :: at1, at2, at3, at4

      vec1x = x(at1, iw) - x(at2, iw)
      vec1y = y(at1, iw) - y(at2, iw)
      vec1z = z(at1, iw) - z(at2, iw)
      vec2x = x(at3, iw) - x(at2, iw)
      vec2y = y(at3, iw) - y(at2, iw)
      vec2z = z(at3, iw) - z(at2, iw)
      vec3x = x(at4, iw) - x(at3, iw)
      vec3y = y(at4, iw) - y(at3, iw)
      vec3z = z(at4, iw) - z(at3, iw)

      norm1x = vec1y * vec2z - vec1z * vec2y
      norm1y = vec1z * vec2x - vec1x * vec2z
      norm1z = vec1x * vec2y - vec1y * vec2x
      norm2x = vec3y * vec2z - vec3z * vec2y
      norm2y = vec3z * vec2x - vec3x * vec2z
      norm2z = vec3x * vec2y - vec3y * vec2x

      sign = norm1x * vec3x + norm1y * vec3y + norm1z * vec3z
      ! TODO: Refactor, make more intermediate results, e.g.
      ! norms of the normal vectors.
      ! TODO: Add error handling for malformed dihedral angles to prevent division by zero
      get_dihedral = 180 / pi * acos( &
             & (norm1x * norm2x + norm1y * norm2y + norm1z * norm2z) / &
             & (dsqrt(norm1x**2 + norm1y**2 + norm1z**2) * dsqrt(norm2x**2 + norm2y**2 + norm2z**2)) &
             & )

      if (sign > 0) get_dihedral = shiftdih - get_dihedral

      return
   end function get_dihedral

   function sanitize_string(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      character(len=len(string) + 200) :: error_msg
      character(len=1) :: ch
      integer :: c, i

      return_string = adjustl(string)
      do i = 1, len(trim(return_string))
         ch = return_string(i:i)
         c = iachar(ch)
         ! check for almost all nonalphabetical chars from ASCII table
         ! allow dash, slash and dot -/.
         if (c < 44 .or. (c > 57 .and. c < 65) .or. c > 172) then
            error_msg = 'Suspicious character "'//ch//'" in input string "'//string//'"'
            call fatal_error(__FILE__, __LINE__, error_msg)
         end if
      end do

   end function sanitize_string

   ! Convert FORTRAN string to zero-terminated C string
   ! and remove any leading and trailing spaces.
   function c_string(string)
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      character(kind=C_CHAR, len=*), intent(in) :: string
      character(kind=C_CHAR, len=len(string) + 1) :: c_string
      c_string = trim(adjustl(string))//C_NULL_CHAR
   end function c_string

   function validate_atom_name(atom_name) result(adjustl_name)
      character(len=*), intent(in) :: atom_name
      character(len=len(atom_name)) :: adjustl_name
      character(len=30) :: error_msg
      ! Note that sanitize_string does adjustl() on the string.
      ! This is quite confusing, maybe we should change that.
      adjustl_name = sanitize_string(atom_name)
      if (len_trim(adjustl_name) > 2) then
         error_msg = 'Incorrect atom name: "'//atom_name//'"'
         call fatal_error(__FILE__, __LINE__, error_msg)
      end if
   end function validate_atom_name

   character(len=2) function normalize_atom_name(atom_name) result(atname)
      character(len=*), intent(in) :: atom_name
      atname = validate_atom_name(atom_name)
      atname(1:1) = toupper(atname(1:1))
      ! Lower the second character of atom name,
      ! this looks niced in the output, and is also needed by TeraChem.
      if (len_trim(atname) == 2) then
         atname(2:2) = tolower(atname(2:2))
      else
         ! Not sure if this is neccessary?
         atname(2:2) = ' '
      end if
   end function normalize_atom_name

   integer function count_atoms_by_name(names, atom_name, natom) result(atom_count)
      character(len=2), intent(in) :: names(:)
      character(len=*), intent(in) :: atom_name
      integer, intent(in) :: natom
      character(len=2) :: norm_name
      integer :: iat
      atom_count = 0
      norm_name = normalize_atom_name(atom_name)
      ! We suppose that the names array is already normalized!
      do iat = 1, natom
         if (names(iat) == norm_name) then
            atom_count = atom_count + 1
         end if
      end do
   end function count_atoms_by_name

   function tolower(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      integer, parameter :: UPPER_A = iachar('A'), UPPER_Z = iachar('Z')
      integer, parameter :: DELTA_LOWER_UPPER = iachar('a') - iachar('A')
      integer :: c, i

      return_string = sanitize_string(string)
      do i = 1, len(trim(return_string))
         c = iachar(return_string(i:i))
         if (c >= UPPER_A .and. c <= UPPER_Z) c = c + DELTA_LOWER_UPPER
         return_string(i:i) = achar(c)
      end do
   end function tolower

   function toupper(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      integer, parameter :: LOWER_A = iachar('a'), LOWER_Z = iachar('z')
      integer, parameter :: DELTA_UPPER_LOWER = iachar('A') - iachar('a')
      integer :: c, i

      return_string = sanitize_string(string)

      do i = 1, len(trim(return_string))
         c = iachar(return_string(i:i))

         if (c >= LOWER_A .and. c <= LOWER_Z) c = c + DELTA_UPPER_LOWER
         return_string(i:i) = achar(c)
      end do
   end function toupper

   ! TODO: Remove this in favour of fatal_error
   subroutine abinerror(chcaller)
      use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
      use mod_interfaces, only: finish
      use mod_files, only: stdout
      character(len=*), intent(in) :: chcaller
      integer, dimension(8) :: time_end
      integer :: u
      ! When the file ERROR exists (perhaps as a remnant of a previous run),
      ! we append it. This also helps us with testing multiple failure check within end-to-end tests.
      open (newunit=u, file='ERROR', action='write', access='append')
      write (u, *) 'FATAL ERROR encountered in subroutine: ', chcaller
      write (u, *) 'Check standard output for further information.'
      close (unit=u)
      call date_and_time(VALUES=time_end)
      write (stdout, *) ''
      write (stdout, *) 'Ended with ERROR at:'
      write (stdout, "(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)") time_end(5), ':', &
         time_end(6), ':', time_end(7), '  ', time_end(3), '.', time_end(2), '.', &
         time_end(1)
      call flush (OUTPUT_UNIT)
      call finish(1)
      stop 1
   end subroutine abinerror

   subroutine print_xyz_arrays(fx, fy, fz, natom, nwalk)
      use mod_files, only: stdout
      real(DP), dimension(:, :), intent(in) :: fx, fy, fz
      integer, intent(in) :: natom, nwalk
      integer :: iat, iw

      do iw = 1, nwalk
         do iat = 1, natom
            write (stdout, *) fx(iat, iw), fy(iat, iw), fz(iat, iw)
         end do
      end do
   end subroutine print_xyz_arrays

   subroutine archive_file(chfile, time_step)
      use mod_general, only: iremd
      use mod_mpi, only: get_mpi_rank
      integer, intent(in) :: time_step
      character(len=*), intent(in) :: chfile
      character(len=200) :: chsystem
      character(len=50) :: chit, charch
      integer :: my_rank

      my_rank = get_mpi_rank()
      if (my_rank == 0 .or. iremd == 1) then
         write (chit, *) time_step
         if (iremd == 1) then
            write (charch, '(A,I2.2)') trim(chfile)//".", my_rank
         else
            charch = trim(chfile)
         end if
         chsystem = 'cp '//trim(charch)//'  '//trim(charch)//'.'//adjustl(chit)
         write (*, *) 'Archiving file ', trim(charch)
         write (*, *) trim(chsystem)
         call system(chsystem)
      end if
   end subroutine archive_file

   subroutine file_exists_or_exit(fname)
      character(len=*), intent(in) :: fname
      character(len=len(fname) + 30) :: error_msg
      logical :: exists
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
         error_msg = 'Could not find file "'//trim(fname)//'"'
         call fatal_error(__FILE__, __LINE__, error_msg)
      end if
   end subroutine file_exists_or_exit

   function get_formatted_date_and_time(time_data) result(formatted_string)
      character(len=25) :: formatted_string
      integer, dimension(8), intent(in) :: time_data
      formatted_string = ''
      ! time_data must be get from date_and_time() intrinsic
      ! e.g. 1:48:39   3.11.2020
      write (formatted_string, "(I2.2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)") time_data(5), ':', &
         time_data(6), ':', time_data(7), '  ', time_data(3), '.', time_data(2), '.', &
         time_data(1)
   end function get_formatted_date_and_time

   real(DP) function ekin_p(px, py, pz, mass, natom, nwalk)
      implicit none
      real(DP), dimension(:, :), intent(in) :: px, py, pz
      real(DP), dimension(:, :), intent(in) :: mass
      integer, intent(in) :: natom, nwalk
      real(DP) :: tmp
      integer :: iw, iat

      ekin_p = 0.0D0

      do iw = 1, nwalk
         do iat = 1, natom
            tmp = px(iat, iw)**2 + py(iat, iw)**2 + pz(iat, iw)**2
            ekin_p = ekin_p + 0.5D0 * tmp / mass(iat, iw)
         end do
      end do
   end function ekin_p

   subroutine real_nonnegative(a, varname)
      real(DP), intent(in) :: a
      character(len=*), intent(in) :: varname
      if (a < 0.0D0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'parameter '//varname//' must be >= 0.0')
      end if
   end subroutine real_nonnegative

   subroutine real_positive(a, varname)
      real(DP), intent(in) :: a
      character(len=*), intent(in) :: varname
      if (a <= 0.0D0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'parameter '//varname//' must be > 0.0')
      end if
   end subroutine real_positive

   subroutine int_nonnegative(i, varname)
      integer, intent(in) :: i
      character(len=*), intent(in) :: varname

      if (i < 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'parameter '//varname//' must be >= 0')
      end if
   end subroutine int_nonnegative

   subroutine int_positive(i, varname)
      integer, intent(in) :: i
      character(len=*), intent(in) :: varname

      if (i <= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'parameter '//varname//' must be > 0')
      end if
   end subroutine int_positive

   subroutine int_switch(i, varname)
      integer, intent(in) :: i
      character(len=*), intent(in) :: varname

      if (i /= 0 .and. i /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'parameter '//varname//' must be 0 or 1')
      end if
   end subroutine int_switch

end module mod_utils
