!-------------------------------------
!- Small utility to generate random integers.
!- Daniel Hollas, 21.8.2014

! Compile as: gfortran MyIRandom.f90 random_lewerence.f -o MyIRandom
! See subroutine  PrintHelp for more information.

program MyIRandom
implicit none
integer              :: iseed0,iseed1,nran,i
real                 :: x=0.0
integer              :: iout=6,iran
real*8,allocatable   :: ranv(:)
integer              :: imax=2147483647 !maximum default integer on most systems
! seed should be set to a large odd integer according to the manual
iseed1=7654321
! secnds(x) gives number of seconds-x elapsed since midnight

call Get_cmdline(iseed0, nran, imax)


! When iseed < 0, take iseed randomly based on current time in seconds
if ( iseed0 < 0 )then
   ! the 2*int(secnds(x)) is always even (int=gives integer) so seed is always odd
   iseed1=iseed1+2*int(secnds(x))
else
   iseed1=iseed0
end if

open(100,file="iseed0",action="write")
write(100,*)iseed1
close(100)

allocate( ranv(nran) )

call vranf(ranv,nran,iseed1,iout)

do i=1,nran
   ranv(i)=ranv(i)*imax
   iran=nint( ranv(i) )
   write(*,*)iran
end do

deallocate (ranv)

end

subroutine Get_cmdline(iseed0, nran, imax)
implicit none
integer,intent(out)  :: iseed0, nran
integer, intent(inout) :: imax
character(len=100)   :: arg
integer              :: i, cac

i=0
cac=command_argument_count()
if ( cac.lt.2.or.cac.gt.3)then
   call PrintHelp()
   call PrintInputError()
end if

call get_command_argument(1, arg)
read(arg,*)iseed0

call get_command_argument(2, arg)
read(arg,*)nran

if (cac.eq.3)then
   call get_command_argument(3, arg)
   read(arg,*)imax
end if

end subroutine Get_cmdline

subroutine PrintHelp()
implicit none
    print '(a)', 'Program for generating random integers from interval [0,imax).'
    print '(a)', 'It is based on substractive lagged fibonacci PRNG implemented by M. Lewerence.'
    print '(a)', 'USAGE: ./MyIrandom <seed> <n_ran> [<imax>]'
    print '(a)', ''
    print '(a)', 'Without cmdline options, this help is printed.'
    print '(a)', ''
    print '(a)', 'cmdline options: ( They must be in correct order!)'
    print '(a)', ''
    print '(a)', '  <seed>    random seed to initialize the PRNG.' 
    print '(a)', '            If negative, we use seed based on current time.' 
    print '(a)', '  <n_ran>   How many random numbers to generate?' 
    print '(a)', '  <imax>    (OPTIONAL) Maximum integer. (default is 2147483647)' 
end subroutine PrintHelp

subroutine PrintInputError()
  write(*,*)'Error during reading command line options. Exiting...'
  stop 1
end subroutine PrintInputError
