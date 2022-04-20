! This is a template file for user-defined analysis function.
! This function will be called if anal_ext=1 in input (section general).
module mod_analyze_ext
   use mod_const, only: DP
   implicit none
   private
   public :: analyze_ext
   save
contains

   ! Modify this function and recompile ABIN
   ! If you need to do some custom analysis along trajectories.
   subroutine analyze_ext()
      use mod_general, only: it
      ! If you need other input arrays, import them from modules
      ! use mod_arrays, only: x, y, z
      character(len=100) :: chsystem
      character(len=50) :: chit
      external system

      ! Launch external BASH script.
      ! Can be used e.g. to analyze wavefunction on-the-fly
      write (chit, *) it
      chsystem = './analyze_ext.sh '//trim(chit)
      call system(chsystem)
   end subroutine analyze_ext

end module mod_analyze_ext
