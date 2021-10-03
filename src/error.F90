module mod_error
   implicit none
   private

   public :: fatal_error
   ! This method is used by pFUnit to set it's own error handler. 
   ! Therefore, this method and its name cannot be changed!
   public :: set_throw_method

   abstract interface
      subroutine error(filename, line_number, message)
         character(len=*), intent(in) :: filename
         integer, intent(in) :: line_number
         character(len=*), optional, intent(in) :: message
      end subroutine error
   end interface

   procedure (error), pointer :: throw_method => Fail

contains

   subroutine set_throw_method(method)
      procedure (error) :: method
      throw_method => method
   end subroutine set_throw_method

   ! This should superseed abinerror()
   subroutine fatal_error(filename, line_number, message)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: line_number
      character(len=*), optional, intent(in) :: message

      call throw_method(filename, line_number, message=message)
   end subroutine fatal_error


   ! Part of abinerror should go here.
   subroutine Fail(filename, line, message)
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
      character(*), intent(in) :: filename
      integer, intent(in) :: line
      character(*), optional, intent(in) :: message
      ! In case of any error, ABIN will return this
      integer, parameter :: ERROR_CODE = 1
      character(:), allocatable :: base_name
      
      base_name = get_base_name(filename)

      write (ERROR_UNIT,'(a, a, a1, i0, a)') &
           & 'ERROR at ', base_name, ':', line, &
           & ' '//adjustl(trim(message))//' '

      ! TODO: Call finalize here.
      stop ERROR_CODE
   end subroutine Fail

   function get_base_name(filename) result(base_name)
      character(:), allocatable :: base_name
      character(*), intent(in) :: filename

      integer :: idx

      idx = scan(filename, '/', back=.true.)

      base_name = filename(idx+1:)

   end function get_base_name

end module mod_error
