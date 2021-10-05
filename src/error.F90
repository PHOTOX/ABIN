! New unit test friendly error handling.
!
! The code in this file has beed adapted from pFUnit_demos repository -
! https://github.com/Goddard-Fortran-Ecosystem/pFUnit_demos
! - and is distributed under the Apache 2.0 license.
!
! Background: code that stops/aborts is unfriendly to testing
! frameworks that ideally want to attempt running later tests. At the
! same time, an explicit dependency on pFUnit forces users of your
! application to install pFUnit even if they don't care about testing
! your code
!
! Solution:
!
!   During initialization the testing framework can override the
!   default behavior with a call to set_error_method(). This logic
!   can comfortably live in your test code, and thus does not
!   introduce any undesirable dependencies (just a bit of
!   obscurity).
module mod_error
   implicit none
   private

   ! This method replaces abinerror()
   public :: fatal_error
   ! This method is used by pFUnit to set it's own error handler,
   ! see unit_tests/throw_with_pfunit.F90
   public :: set_error_method

   abstract interface
      subroutine error(filename, line_number, message)
         character(len=*), intent(in) :: filename
         integer, intent(in) :: line_number
         character(len=*), intent(in) :: message
      end subroutine error
   end interface

   ! This is the defailt error handle, which gets
   ! overriden in pFUnit unit tests.
   procedure(error), pointer :: error_method => print_error_and_stop

contains

   subroutine set_error_method(method)
      procedure(error) :: method
      error_method => method
   end subroutine set_error_method

   ! filename and line_number parameters should be passed using the preprocessor
   ! defined constants __FILE__ and __LINE__
   subroutine fatal_error(filename, line_number, message)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: line_number
      character(len=*), intent(in) :: message

      call error_method(filename, line_number, message=message)
   end subroutine fatal_error

   subroutine print_error_and_stop(filename, line, message)
      character(*), intent(in) :: filename
      integer, intent(in) :: line
      character(*), intent(in) :: message
      ! In case of an error, ABIN will return this to the shell.
      integer, parameter :: ERROR_CODE = 1

      call print_error(filename, line, message)

      ! Try to finalize various modules gracefully.
      call finish(ERROR_CODE)

      stop ERROR_CODE
   end subroutine print_error_and_stop

   ! We print the error both to stderr and to file 'ERROR'.
   ! Presence of this file straightforwardly indicates a problem,
   ! and scripts checking for errors don't need to grep ABIN output.
   subroutine print_error(filename, line, message)
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
      character(*), intent(in) :: filename
      integer, intent(in) :: line
      character(*), intent(in) :: message
      character(len=*), parameter :: ERROR_FILE = 'ERROR'
      character(len=:), allocatable :: base_name
      integer :: iunit

      base_name = get_base_name(filename)

      write (ERROR_UNIT, '(a9, a, a1, i0, a)') &
           & 'ERROR at ', base_name, ':', line, &
           & ' '//adjustl(trim(message))

      ! We do not print line number in the ERROR file,
      ! because we check the ERROR file in tests, and any change
      ! in a line number of the error would broke the test.
      ! TODO: We should probably make the test suite more clever to
      ! ignore the line number when comparing the ERROR to ERROR.ref
      open (newunit=iunit, file=ERROR_FILE, action='write')
      write (iunit, '(a)') 'ERROR in '//base_name//': '//adjustl(trim(message))
      close (unit=iunit)

      call flush (OUTPUT_UNIT)
      call flush (ERROR_UNIT)
   end subroutine print_error

   ! Get get base filename from the full path.
   function get_base_name(filename) result(base_name)
      character(:), allocatable :: base_name
      character(*), intent(in) :: filename
      integer :: idx

      idx = scan(filename, '/', back=.true.)
      base_name = filename(idx + 1:)
   end function get_base_name

end module mod_error
