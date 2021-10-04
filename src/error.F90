module mod_error
   implicit none
   private

   public :: fatal_error
   ! This method is used by pFUnit to set it's own error handler.
   ! Therefore, this method and its name cannot be changed!
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

   ! This should superseed abinerror()
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
      ! We write out the error message to file 'ERROR'
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
      write (iunit, '(a)') 'ERROR in '//base_name//': '//message
      close (unit=iunit)

      call flush (OUTPUT_UNIT)
      call flush (ERROR_UNIT)
   end subroutine print_error

   function get_base_name(filename) result(base_name)
      character(:), allocatable :: base_name
      character(*), intent(in) :: filename

      integer :: idx

      idx = scan(filename, '/', back=.true.)

      base_name = filename(idx + 1:)

   end function get_base_name

end module mod_error
