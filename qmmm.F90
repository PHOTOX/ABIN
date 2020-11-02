! This is only a stub QMMM module for now,
! although natmm and natqm are initialized in init.F90
! and used throughout the code.
module mod_qmmm
   implicit none
   private
   public :: natqm, natmm
   integer :: natqm, natmm
   ! NOT used anywhere at the moment
   character(len=10)    :: qmmmtype='NA'
end module mod_qmmm
