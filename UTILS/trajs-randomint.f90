!-------------------------------------
!- Small utility to generate random integers.
! Uses simple PRNG from numerical recipes,
! which was used in older version of ABIN.
! This PRNG is intended to generate only initial random seeds
! to start multiple ABIN trajecotries.
! For other usage, use abin-randomint, which should be much better

! TODO: needs testing

! Daniel Hollas, 15.2.2016

! See subroutine  PrintHelp for more information.

program abin_randomint
implicit none
integer              :: iseed0, iseed1, nran, i
real                 :: x = 0.0
integer              :: iout = 6, iran
real*8, allocatable  :: ranv(:)
integer              :: imax = 2147483647 !maximum default integer on most systems
real*8               :: ran1
   ! seed should be set to a large odd integer according to the manual
   iseed1 = 7654321
   ! secnds(x) gives number of seconds-x elapsed since midnight

   call Get_cmdline(iseed0, nran, imax)


! When iseed < 0, take iseed randomly based on current time in seconds
! TODO: we should make this better, use milliseconds or something....
! and probably base it on UNIX time or something
   if ( iseed0 < 0 )then
      ! the 2*int(secnds(x)) is always even (int=gives integer) so seed is always odd
      iseed1 = iseed1+2*int(secnds(x))
   else
      iseed1 = iseed0
   end if

   open(100,file="iseed0",action="write")
   write(100,*)iseed1
   close(100)

   allocate( ranv(nran) )

   ! Initialize the sequence
   ranv(1) = ran1(-iseed0)

   do i=1,nran
      ranv(i) = ran1(0)
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

   i = 0
   cac = command_argument_count()
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
    print '(a)', 'USAGE: ./abin-randomint <seed> <n_ran> [<imax>]'
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

! (PSEUDO) RANDOM NUMBER GENERATOR     
! ================================
!                            
! RETURNS A UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0. 
! SET IDUM TO ANY NEGATIVE INTEGER VALUE TO INITIALIZE 
! OR REINITIALIZE THE SEQUENCE.
!                                                                      
! FROM "NUMERICAL RECIPES"                                            
!
! SAVE statement for IFF, R, IXn      B. Schmidt, May  31, 1992
! RAN1 returned as REAL(DP)             B. Schmidt, May  31, 1992
!
!-------------------------------------------------------------------
REAL*8 FUNCTION RAN1(IDUM)

     real*8  :: r(97)
     integer ::  ix1,ix2,ix3
     integer  ::  idum
     integer  ::  j,iff

     real*8, parameter  :: RM1=3.8580247D-6, RM2=7.4373773D-6 
     integer, parameter :: M1=259200,IA1=7141,IC1=54773
     integer, parameter :: M2=134456,IA2=8121,IC2=28411
     integer, parameter :: M3=243000,IA3=4561,IC3=51349                    

     SAVE IFF, R, IX1,IX2,IX3    !?????????
     DATA IFF /0/                                               

     IF (IDUM.LT.0.OR.IFF.EQ.0) THEN                           
        IFF=1                                      
        IX1=MOD(IC1-IDUM,M1)                      
        IX1=MOD(IA1*IX1+IC1,M1)                  
        IX2=MOD(IX1,M2)                         
        IX1=MOD(IA1*IX1+IC1,M1)                
        IX3=MOD(IX1,M3)                       
        DO 11 J=1,97                         
           IX1=MOD(IA1*IX1+IC1,M1)           
           IX2=MOD(IA2*IX2+IC2,M2)          
           R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1
11      CONTINUE                              
        IDUM=1                               
     ENDIF                                 
     IX1=MOD(IA1*IX1+IC1,M1)              
     IX2=MOD(IA2*IX2+IC2,M2)             
     IX3=MOD(IA3*IX3+IC3,M3)            
     J=1+(97*IX3)/M3                   
     IF(J.GT.97.OR.J.LT.1)STOP      
     RAN1=  R(J) 
     R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1


     RETURN                             
END function ran1

