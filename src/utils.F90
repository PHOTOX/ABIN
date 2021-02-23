! Various helper functions and subroutines that may be used throughout the program.
module mod_utils
   use mod_const, only: DP
   use mod_interfaces
   implicit none
   public
contains

   real(DP) function get_distance(x, y, z, at1, at2, iw)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: at1, at2, iw
      real(DP) :: r
      if (at1 == at2) then
         write (*, *) 'ERROR: Atom indices in function get_distance are not unique!'
         call abinerror('get_distance')
      end if

      r = (x(at1, iw) - x(at2, iw))**2
      r = r + (y(at1, iw) - y(at2, iw))**2
      r = r + (z(at1, iw) - z(at2, iw))**2
      r = sqrt(r)
      get_distance = r
   end function get_distance

   real(DP) function get_angle(x, y, z, at1, at2, at3, iw)
      use mod_const, only: PI
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw
      real(DP) :: vec1x, vec1y, vec1z
      real(DP) :: vec2x, vec2y, vec2z
      integer :: at1, at2, at3

      if (at1 == at2 .or. at1 == at3 .or. at2 == at3) then
         write (*, *) 'ERROR: Atom indices in function get_angle are not unique!'
         call abinerror('get_angle')
      end if

      vec1x = x(at1, iw) - x(at2, iw)
      vec1y = y(at1, iw) - y(at2, iw)
      vec1z = z(at1, iw) - z(at2, iw)
      vec2x = x(at3, iw) - x(at2, iw)
      vec2y = y(at3, iw) - y(at2, iw)
      vec2z = z(at3, iw) - z(at2, iw)
      get_angle = 180 / pi * acos((vec1x * vec2x + vec1y * vec2y + vec1z * vec2z) / &
               & (sqrt(vec1x**2 + vec1y**2 + vec1z**2) * sqrt(vec2x**2 + vec2y**2 + vec2z**2)))

      return
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
      get_dihedral = 180 / pi * acos( &
             & (norm1x * norm2x + norm1y * norm2y + norm1z * norm2z) / &
             & (sqrt(norm1x**2 + norm1y**2 + norm1z**2) * sqrt(norm2x**2 + norm2y**2 + norm2z**2)) &
             & )

      if (sign > 0) get_dihedral = shiftdih - get_dihedral

      return
   end function get_dihedral

   function SanitizeString(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      integer :: c, i

      ! TODO: Not sure whether this is handling non-ASCII cases correctly
      return_string = adjustl(string)
      do i = 1, len(trim(return_string))
         c = iachar(return_string(i:i))
         ! check for almost all nonalphabetical chars from ASCII table
         ! allow dash, slash and dot -/.
         if (c < 44 .or. (c > 57 .and. c < 65) .or. c > 172) then
            write (*, *) 'Suspicious character found in one of the input strings: '//string
            write (*, *) 'ASCII index:', c, ' Position in the string:', i
            write (*, *) 'Please check your input files. Exiting...'
            call abinerror('SanitizeString')
         end if
      end do

   end function SanitizeString

   ! Convert FORTRAN string to zero-terminated C string
   ! and remove any leading and trailing spaces.
   function c_string(string)
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      character(kind=C_CHAR, len=*), intent(in) :: string
      character(kind=C_CHAR, len=len(string) + 1) :: c_string
      c_string = trim(adjustl(string))//C_NULL_CHAR
   end function c_string

   function UpperToLower(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      integer, parameter :: UPPER_A = iachar('A'), UPPER_Z = iachar('Z')
      integer, parameter :: DELTA_LOWER_UPPER = iachar('a') - iachar('A')
      integer :: c, i

      return_string = SanitizeString(string)
      do i = 1, len(trim(return_string))
         c = iachar(return_string(i:i))
         if (c >= UPPER_A .and. c <= UPPER_Z) c = c + DELTA_LOWER_UPPER
         return_string(i:i) = achar(c)
      end do
   end function UpperToLower

   function LowerToUpper(string) result(return_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: return_string
      integer, parameter :: LOWER_A = iachar('a'), LOWER_Z = iachar('z')
      integer, parameter :: DELTA_UPPER_LOWER = iachar('A') - iachar('a')
      integer :: c, i

      return_string = SanitizeString(string)

      do i = 1, len(trim(return_string))
         c = iachar(return_string(i:i))

         if (c >= LOWER_A .and. c <= LOWER_Z) c = c + DELTA_UPPER_LOWER
         return_string(i:i) = achar(c)
      end do
   end function LowerToUpper

   subroutine not_compiled_with(feature, caller)
      character(len=*), intent(in) :: feature, caller
      write (*, *) 'ERROR: ABIN was not compiled with '//feature
      call abinerror(caller)
   end subroutine not_compiled_with

   ! TODO: Maybe move this into a separate error handling module,
   ! together with finish(), and move it to a separate file
   ! Though need to figure out how to do it without having circular dependencies :-(
   ! abinerror is needed everywhere, so shouldn't depend on much,
   ! but finish() is basically dependent on everything. :-(
   !
   ! In any case, mod_utils should not depend on anything outside of mod_const and mod_general!
   subroutine abinerror(chcaller)
      use mod_general, only: my_rank
      character(len=*), intent(in) :: chcaller
      integer, dimension(8) :: time_end
      ! When the file ERROR exists (perhaps as a remnant of a previous run),
      ! we append it. This also helps us with testing multiple failure check within end-to-end tests.
      open (unit=500, file='ERROR', action='write', access='append')
      write (500, *) 'FATAL ERROR encountered in subroutine: ', chcaller
      write (500, *) 'Check standard output for further information.'
      close (unit=500)
      if (my_rank == 0) then
         call date_and_time(VALUES=time_end)
         write (*, *) ''
         write (*, *) 'Ended with ERROR at:'
         write (*, "(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)") time_end(5), ':', &
            time_end(6), ':', time_end(7), '  ', time_end(3), '.', time_end(2), '.', &
            time_end(1)
      end if
      call flush (6)
      call finish(1)
      stop 1
   end subroutine abinerror

   subroutine print_xyz_arrays(fx, fy, fz)
      use mod_general, only: nwalk, natom
      use mod_const, only: DP
      implicit none
      real(KIND=DP), intent(in) :: fx(:, :), fy(:, :), fz(:, :)
      integer :: iat, iw

      do iw = 1, nwalk
         do iat = 1, natom
            write (*, *) fx(iat, iw), fy(iat, iw), fz(iat, iw)
         end do
      end do

   end subroutine print_xyz_arrays

   subroutine archive_file(chfile, time_step)
      use mod_general, only: iremd, my_rank
      integer, intent(in) :: time_step
      character(len=*), intent(in) :: chfile
      character(len=200) :: chsystem
      character(len=50) :: chit, charch
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

   subroutine debug_output(msg)
      ! See SO about fortran std units
      ! https://stackoverflow.com/questions/8508590/standard-input-and-output-units-in-fortran-90#8508757
      use iso_fortran_env, only: ERROR_UNIT
      character(len=*), intent(in) :: msg
      write (ERROR_UNIT, '(A)') msg
      call flush (ERROR_UNIT)
   end subroutine debug_output

   subroutine file_exists_or_exit(fname)
      character(len=*), intent(in) :: fname
      logical :: exists
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
         write (*, '(2A)') 'ERROR: Could not find file '//trim(fname)
         call abinerror('file_exists_or_exit')
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

end module mod_utils
