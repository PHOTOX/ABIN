! Various small functions and subroutines, that may be used throughout the program.
! Daniel Hollas         2014
module mod_utils
   use mod_const, only: DP
   implicit none
   private :: PrintHelp
   contains

   real(DP) function get_distance(x,y,z,at1,at2,iw)
   real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
   integer,intent(in)  :: at1, at2, iw
   real(DP) :: r
   r=(x(at1,iw)-x(at2,iw))**2
   r=r+(y(at1,iw)-y(at2,iw))**2
   r=r+(z(at1,iw)-z(at2,iw))**2
   r=sqrt(r)
   get_distance=r
   end function get_distance

   function SanitizeString(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer                     :: c, i

      return_string = adjustl(string)
      do i=1,len(trim(return_string))
         c = iachar( return_string(i:i))
         ! check for almost all nonalphabetical chars from ASCII table
         ! allow dash, slash and dot -/.
         if (c < 44.or.(c>57.and.c<65).or.c>172)then
            write(*,*)'Suspicious character found in one of the input strings: '//string
            write(*,*)'ASCII index:', c,' Position in the string:', i
            write(*,*)'Please check your input files. Exiting...'
            call abinerror('UpperToLower')
         end if
      end do
      
   end function SanitizeString

   function UpperToLower(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer, parameter          :: UPPER_A = iachar ('A'), UPPER_Z = iachar('Z')
      integer, parameter          :: DELTA_LOWER_UPPER = iachar('a')-iachar('A')
      integer                     :: c,i
   
      return_string = SanitizeString(string)
      do i=1,len(trim(return_string))
         c = iachar( return_string(i:i))
         if (c >= UPPER_A .and. c <= UPPER_Z) c = c + DELTA_LOWER_UPPER
         return_string(i:i) = achar(c)
      end do
   end function UpperToLower
   
   function LowerToUpper(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer, parameter          :: LOWER_A = iachar ('a'), LOWER_Z = iachar('z')
      integer, parameter          :: DELTA_UPPER_LOWER = iachar('A')-iachar('a')
      integer                     :: c,i
   
      return_string = SanitizeString(string)

      do i=1,len(trim(return_string))
         c =iachar( return_string(i:i))

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

   subroutine Get_cmdline(chinput, chcoords, chveloc )
   character(len=*),intent(inout)   :: chinput, chcoords, chveloc
   character(len=len(chinput))   :: arg
   integer            :: i
   logical            :: lexist
   
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
      case ('-v')
         i=i+1
         call get_command_argument(i, arg)
         read(arg,'(A)')chveloc
         chveloc=trim(chveloc)
      case default
         write(*,*)'Invalid command line argument!'
         call abinerror('Get_cmdline')
      end select

   end do
   !check for existence of input files

   inquire(file=chcoords,exist=lexist)
   if (.not.lexist)then
      write(*,*)'FATAL: Input file does not exists!'
      write(*,*)chcoords
      stop 1
   end if

   inquire(file=chinput,exist=lexist)
   if (.not.lexist)then
      write(*,*)'FATAL: The following input file does not exists!'
      write(*,*)chinput
      stop 1
   end if

   if (chveloc.ne.'')then
      inquire(file=chveloc,exist=lexist)
      if (.not.lexist)then
         write(*,*)'FATAL: The following input file does not exists!'
         write(*,*)chveloc
         stop 1
      end if
   end if

   end subroutine Get_cmdline


   subroutine PrintHelp()
   implicit none
    print '(a)', 'ABIN: Multipurpose ab initio MD program.'
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help               print help and exit'
    print '(a)', '  -i <input_parameters>    default: input.in'
    print '(a)', '  -x <input_coordinates>   default: mini.dat'
    print '(a)', '  -v <input_velocities>    no default'
    print '(a)', ''
   end subroutine PrintHelp

end module mod_utils
