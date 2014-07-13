! Various small utilities, that may be used throughout the program
! Daniel Hollas         2014
module mod_utils
  implicit none
   private :: PrintHelp
   contains
   elemental function UpperToLower(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer, parameter          :: UPPER_A = iachar ('A'),UPPER_Z = iachar('Z')
      integer, parameter          :: DELTA_LOWER_UPPER = iachar('a')-iachar('A')
      integer                     :: c,i
   
      do i=1,len(string)
         c =iachar( string(i:i))
         if (c >= UPPER_A .and. c <= UPPER_Z) c = c + DELTA_LOWER_UPPER
         return_string(i:i) = achar(c)
      end do
   end function UpperToLower
   
   elemental function LowerToUpper(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer, parameter          :: LOWER_A = iachar ('a'),LOWER_Z = iachar('z')
      integer, parameter          :: DELTA_UPPER_LOWER = iachar('A')-iachar('a')
      integer                     :: c,i
   
      do i=1,len(string)
         c =iachar( string(i:i))
      if (c >= LOWER_A .and. c <= LOWER_Z) c = c + DELTA_UPPER_LOWER
      return_string(i:i) = achar(c)
      end do
   end function LowerToUpper

   subroutine abinerror(chcaller)
   character(len=*),intent(in)   :: chcaller
   open(unit=500,file='ERROR')
   write(500,*)'FATAL ERROR encountered in subroutine: ',chcaller
   write(500,*)'Check standard output for further information. Exiting now...'
   close(unit=500)
   call flush
   stop 1
   end subroutine abinerror

   subroutine printf(fx,fy,fz)
   use mod_general, only: nwalk,natom
   use mod_const, only: DP
   implicit none
   real(KIND=DP),intent(in) :: fx(:,:),fy(:,:),fz(:,:)
   integer :: iat,iw

   do iw=1,nwalk
      do iat=1,natom
         write(*,*)fx(iat,iw),fy(iat,iw),fz(iat,iw)
      enddo
   enddo

   end subroutine printf

   subroutine Get_cmdline(chinput, chcoords )
   character(len=*),intent(inout)   :: chinput, chcoords
   character(len=len(chinput))   :: arg
   integer            :: i
   
   i=0
   do while (i < command_argument_count())

     i=i+1
     call get_command_argument(i, arg)
   
      select case (arg)
      case ('-h', '--help')
         call PrintHelp()
         stop 0

      case ('-i')
         i=i+1
         call get_command_argument(i, arg)
         !-format specifier is needed here in case of slashes
         read(arg,'(A)')chinput
         chinput=trim(chinput)

      case ('-x')
         i=i+1
         call get_command_argument(i, arg)
         read(arg,'(A)')chcoords
         chcoords=trim(chcoords)

      case default
         write(*,*)'Invalid command line argument!'
         call abinerror('Get_cmdline')
      end select

   end do

   end subroutine Get_cmdline


   subroutine PrintHelp()
   implicit none
    print '(a)', 'ABIN: Multipurpose ab initio BOMD program.'
    print '(a)', 'Without cmdline options, the program use sensible default values.'
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help             print help and exit'
    print '(a)', '  -i <input_parameters>    default: input.in'
    print '(a)', '  -x <input_coordinates>   default: mini.dat'
   end subroutine PrintHelp

end module mod_utils
