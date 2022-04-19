! -------------------------------------
! A program to generate 32-bit random integers
! based on fortran intrinsic random_number(),
! as implemented in ABIN.
!
! Can be used to deterministically generate
! random seeds for a set of trajectories.
!
! Compile as:
! $ gfortran -O2 -fopenmp -I"src/" utils/abin-randomint.f90 -Lsrc/ -labin -Lwater_potentials/ -lwater -lstdc++
!
! See subroutine print_help() for more information.
! -------------------------------------
program abin_randomint
   use mod_files, only: stdout_to_devnull
   use mod_prng_init, only: initialize_fortran_prng, random_ints, get_random_seed
   implicit none
   integer :: seed = -1, nran = 1
   integer, allocatable :: irans(:)
   integer :: i, u

   call get_cmdline(seed, nran)

   call stdout_to_devnull()

   if (seed <= 0) then
      seed = get_random_seed()
   end if

   open (newunit=u, file="iseed0", action="write")
   write (u, '(I0)') seed
   close (u)

   allocate (irans(nran))

   call initialize_fortran_prng(seed)

   call random_ints(irans, nran, testing_mode=.false.)

   do i = 1, nran
      write (*, '(I0)') irans(i)
   end do
end program

subroutine get_cmdline(iseed, nran)
   use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
   implicit none
   integer, intent(inout) :: iseed, nran
   character(len=300) :: arg
   integer :: i

   i = 0
   do while (i < command_argument_count())

      i = i + 1
      call get_command_argument(i, arg)

      select case (arg)
      case ('-h', '--help')
         call print_help()
         stop 0
      case ('-s', '--seed')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg, '(I10)') iseed
      case ('-n', '--num')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg, '(I10)') nran
      case default
         call print_help()
         print*,''
         write (ERROR_UNIT, *) 'ERROR: Invalid command line option "'//trim(arg)//'"'
         stop 1
      end select

   end do
end subroutine get_cmdline

subroutine print_help()
   implicit none
   print '(a)', 'Program for generating 32-bit random integers from interval [0,2147483647).'
   print '(a)', ''
   print '(a)', 'USAGE: ./abin-randomint [-s <seed>] [-n <n_ran>]'
   print '(a)', ''
   print '(a)', 'optional command line options:'
   print '(a)', ''
   print '(a)', '  -h, --help         print help and exit'
   print '(a)', '  -s/--seed <seed>   random seed to initialize the PRNG'
   print '(a)', '                     If negative, seed is taken from the OS (default)'
   print '(a)', '  -n/--num <n_ran>   how many random numbers (default=1)'
end subroutine print_help
