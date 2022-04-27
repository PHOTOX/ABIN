! This is a template file for user-defined analysis function.
! This function will be called if anal_ext=1 in input (section general).
module mod_analyze_ext
   use mod_const, only: DP
   use mod_error, only: fatal_error
   implicit none
   private
   public :: analyze_ext
   save
contains

   ! Modify this function and recompile ABIN
   ! If you need to do some custom analysis along trajectories.
   subroutine analyze_ext(time_step)
      integer, intent(in) :: time_step
      ! If you need other input arrays, import them from modules
      ! use mod_arrays, only: x, y, z
      character(len=100) :: chsystem
      integer :: icmd, istat

      ! Launch external BASH script.
      ! Can be used e.g. to analyze wavefunction on-the-fly
      write (chsystem, '(A,I0)') './analyze_ext.sh ', time_step
      call execute_command_line(chsystem, exitstat=istat, cmdstat=icmd)
      if (icmd /= 0 .or. istat /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Failed to execute command '//trim(chsystem))
      end if
   end subroutine analyze_ext

end module mod_analyze_ext
