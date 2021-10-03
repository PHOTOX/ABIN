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

   procedure(error), pointer :: error_method => Fail

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

   ! Part of abinerror should go here.
   subroutine Fail(filename, line, message)
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
      character(*), intent(in) :: filename
      integer, intent(in) :: line
      character(*), intent(in) :: message
      ! In case of any error, ABIN will return this
      integer, parameter :: ERROR_CODE = 1
      ! We write out the error message to file 'ERROR'
      ! Presence of this file straightforwardly indicates a problem.
      character(len=*), parameter :: ERROR_FILE = 'ERROR'
      character(len=:), allocatable :: base_name, full_message
      integer :: iunit, strlen

      base_name = get_base_name(filename)

      strlen = len_trim(message) + len(base_name) + 20
      allocate (character(len=strlen) :: full_message)

      write (full_message, '(a9, a, a1, i0, a)') &
           & 'ERROR at ', base_name, ':', line, &
           & ' '//adjustl(trim(message))

      write (ERROR_UNIT, '(A)') full_message

      open (newunit=iunit, file=ERROR_FILE, action='write')
      write (iunit, '(a)') full_message
      close (unit=iunit)

      call flush (OUTPUT_UNIT)
      call flush (ERROR_UNIT)

      ! Depending on the state of the program, this call might itself fail
      call finish(ERROR_CODE)

      stop ERROR_CODE
   end subroutine Fail

   function get_base_name(filename) result(base_name)
      character(:), allocatable :: base_name
      character(*), intent(in) :: filename

      integer :: idx

      idx = scan(filename, '/', back=.true.)

      base_name = filename(idx + 1:)

   end function get_base_name

end module mod_error
