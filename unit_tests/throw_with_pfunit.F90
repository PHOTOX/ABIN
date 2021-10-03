module throw_with_pfunit_mod
   ! This is a method from ABIN, see src/error.F90
   use mod_error, only: set_error_method
   implicit none
   private

   public :: throw
   public :: initialize_throw

contains

   subroutine throw(file_name, line_number, message)
      use funit, only: SourceLocation
      use funit, only: pFUnit_throw => throw
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: line_number
      character(len=*), intent(in) :: message

      call pFUnit_throw(message, SourceLocation(file_name, line_number))

   end subroutine throw

   subroutine initialize_throw()
      call set_error_method(throw)
   end subroutine initialize_throw

   
end module throw_with_pfunit_mod
