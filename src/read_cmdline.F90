module mod_cmdline
   private
   public :: get_cmdline

contains

   subroutine print_help()
      use mod_interfaces, only: print_compile_info
      implicit none
      integer, dimension(8) :: time_data

      call date_and_time(values=time_data)
      print '(a)', ''
      print '(a)', 'ABIN: Multipurpose ab initio MD program.'
      print '(a)', ''
      call print_compile_info()
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -h, --help               print help and exit'
      print '(a)', '  -i <input_parameters>    default: input.in'
      print '(a)', '  -x <input_coordinates>   default: mini.dat'
      print '(a)', '  -v <input_velocities>    no default'
      print '(a)', ''
   end subroutine print_help

   subroutine get_cmdline(chinput, chcoords, chveloc, tc_port)
      use mod_utils, only: abinerror, file_exists_or_exit
      character(len=*), intent(inout) :: chinput, chcoords, chveloc, tc_port
      character(len=len(chinput)) :: arg
      integer :: i

      i = 0
      do while (i < command_argument_count())

         i = i + 1
         call get_command_argument(i, arg)

         select case (arg)
         case ('-h', '--help')
            call print_help()
            stop 0
         case ('-i')
            i = i + 1
            call get_command_argument(i, arg)
            !-format specifier is needed here in case of slashes
            read (arg, '(A)') chinput
            chinput = trim(chinput)
         case ('-x')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') chcoords
            chcoords = trim(chcoords)
         case ('-v')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') chveloc
            chveloc = trim(chveloc)
         case ('-M')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') tc_port
            tc_port = trim(tc_port)
         case default
            write (*, '(A)') 'Invalid command line argument '//arg
            call print_help()
            call abinerror('get_cmdline')
         end select

      end do

      call file_exists_or_exit(chinput)
      ! Input velocities and coordinates are checked in init,
      ! because of REMD
   end subroutine get_cmdline

end module mod_cmdline
