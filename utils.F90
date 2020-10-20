! Various helper functions and subroutines that may be used throughout the program.
module mod_utils
   use mod_const, only: DP
   implicit none
   public
   contains

   real(DP) function get_distance(x, y, z, at1, at2, iw)
   real(DP),intent(in) :: x(:,:), y(:,:), z(:,:)
   integer,intent(in)  :: at1, at2, iw
   real(DP) :: r
   r =     ( x(at1, iw) - x(at2, iw) )**2
   r = r + ( y(at1, iw) - y(at2, iw) )**2
   r = r + ( z(at1, iw) - z(at2, iw) )**2
   r = sqrt(r)
   get_distance = r
   end function get_distance

   function SanitizeString(string) result (return_string)
      character(len=*),intent(in) :: string
      character(len=len(string))  :: return_string
      integer                     :: c, i

      ! TODO: Not sure whether this is handling non-ASCII cases correctly
      return_string = adjustl(string)
      do i=1,len(trim(return_string))
         c = iachar( return_string(i:i))
         ! check for almost all nonalphabetical chars from ASCII table
         ! allow dash, slash and dot -/.
         if (c < 44.or.(c>57.and.c<65).or.c>172)then
            write(*,*)'Suspicious character found in one of the input strings: '//string
            write(*,*)'ASCII index:', c,' Position in the string:', i
            write(*,*)'Please check your input files. Exiting...'
            call abinerror('SanitizeString')
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

   ! TODO: Maybe move this into a separate error handling module, 
   ! together with finish(), and move it to a separate file
   ! Though need to figure out how to do it without having circular dependencies :-(
   ! abinerror is needed everywhere, so shouldn't depend on much,
   ! but finish() is basically dependent on everything. :-(
   !
   ! In any case, mod_utils should not depend on anything outside of mod_const and mod_general!
   subroutine abinerror(chcaller)
      use mod_general,  only: my_rank
      character(len=*),intent(in)   :: chcaller
      integer,dimension(8) :: time_end
      open(unit=500, file='ERROR', action='write', access='sequential')
      write(500,*)'FATAL ERROR encountered in subroutine: ',chcaller
      write(500,*)'Check standard output for further information.'
      close(unit=500)
      if (my_rank.eq.0)then
         call date_and_time(VALUES=time_end)
         write(*,*)''
         write(*,*)'Ended with ERROR at:'
         write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")time_end(5),':',&
           time_end(6),':',time_end(7),'  ',time_end(3),'.',time_end(2),'.',&
           time_end(1)
      end if
      call flush(6)
      call finish(1)
      stop 1
   end subroutine abinerror

   subroutine print_xyz_arrays(fx,fy,fz)
   use mod_general, only: nwalk,natom
   use mod_const, only: DP
   implicit none
   real(KIND=DP),intent(in) :: fx(:,:),fy(:,:),fz(:,:)
   integer :: iat, iw

   do iw=1,nwalk
      do iat=1,natom
         write(*,*)fx(iat,iw),fy(iat,iw),fz(iat,iw)
      enddo
   enddo

   end subroutine print_xyz_arrays

   subroutine archive_file(chfile, time_step)
   use mod_general,  ONLY: iremd, my_rank
   integer, intent(in)  :: time_step
   character(len=*), intent(in)  :: chfile
   character(len=200)   :: chsystem
   character(len=50)    :: chit, charch
   if(my_rank.eq.0.or.iremd.eq.1)then
      write (chit,*)time_step
      if(iremd.eq.1)then
         write(charch,'(A,I2.2)')trim(chfile)//".",my_rank
      else
         charch=trim(chfile)
      end if
      chsystem='cp '//trim(charch)//'  '//trim(charch)//'.'//adjustl(chit)
      write(*,*)'Archiving file ', trim(charch)
      write(*,*)trim(chsystem)
      call system(chsystem)
   end if
   end subroutine archive_file

   subroutine debug_output(msg)
   ! See SO about fortran std units
   ! https://stackoverflow.com/questions/8508590/standard-input-and-output-units-in-fortran-90#8508757
   use iso_fortran_env, only: error_unit
   character(len=*), intent(in) :: msg
   write(error_unit, '(A)') msg
   call flush(error_unit)
   end subroutine debug_output

   subroutine file_exists_or_exit(fname)
   character(len=*), intent(in) :: fname
   logical :: exists
   inquire(file=trim(fname), exist = exists)
   if (.not.exists)then
      write(*,'(2A)')'ERROR: Could not find file '//trim(fname)
      call abinerror('file_exists_or_exit')
   end if
   end subroutine file_exists_or_exit

end module mod_utils
