module mod_cmdline
   use mod_interfaces, only: print_compile_info
   private
   public :: get_cmdline

contains

   subroutine print_help()
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
      print '(a)', '  -V, --version            print version and compiler options'
      print '(a)', '  -i <input_parameters>    default: input.in'
      print '(a)', '  -x <input_coordinates>   default: mini.dat'
      print '(a)', '  -v <input_velocities>    no default'
      print '(a)', ''
   end subroutine print_help

   subroutine get_cmdline(chinput, chcoords, chveloc, &
                        & tc_server_name, tcpb_input_file, tcpb_host, tcpb_port)
      use mod_error, only: fatal_error
      use mod_utils, only: file_exists_or_exit
      character(len=*), intent(inout) :: chinput, chcoords, chveloc
      character(len=*), intent(inout) :: tc_server_name   ! TeraChem MPI interface
      character(len=*), intent(inout) :: tcpb_input_file  ! TeraChem protobuf interface
      character(len=*), intent(inout) :: tcpb_host
      integer, intent(inout) :: tcpb_port
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
         case ('-V', '--version')
            call print_compile_info()
            stop 0
         case ('-i')
            i = i + 1
            call get_command_argument(i, arg)
            ! Format specifier is needed here in case of filenames with slashes/dashes
            read (arg, '(A)') chinput
            if (trim(chinput) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty -i argument. Provide filename with input parameters')
            end if
         case ('-x')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') chcoords
            if (trim(chcoords) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty -x argument. Provide filename with coordinates')
            end if
         case ('-v')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') chveloc
            if (trim(chveloc) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty -v argument. Provide filename with velocities')
            end if
         case ('-M')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') tc_server_name
            if (trim(tc_server_name) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty -M argument. Provide name of the TeraChem MPI server')
            end if
         case ('--tcpb-input-file')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') tcpb_input_file
            if (trim(tcpb_input_file) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty --tcpb-input-file argument. Provide name of the TeraChem input file')
            end if
         case ('--tcpb-host')
            i = i + 1
            call get_command_argument(i, arg)
            read (arg, '(A)') tcpb_host
            if (trim(tcpb_host) == '') then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty --tcpb-host argument. Provide hostname of the TeraChem TCPB server')
            end if
         case ('--tcpb-port')
            i = i + 1
            call get_command_argument(i, arg)
            ! NOTE: If the port number is too large, this will trim it!
            ! This is okay since in the launch script we determine port number automatically anyway.
            read (arg, '(I5)') tcpb_port
            ! NOTE: When the parameter is missing, at least GCC compilers assign 0 as default
            if (tcpb_port == 0) then
               call fatal_error(__FILE__, __LINE__, &
                  & 'Empty --tcpb-port argument. Provide port of the TeraChem TCPB server')
            end if
         case default
            call print_help()
            call fatal_error(__FILE__, __LINE__, &
               & 'Invalid command line argument "'//trim(arg)//'"')
         end select

      end do

      call file_exists_or_exit(chinput)
      ! Input velocities and coordinates are checked in init,
      ! because of REMD
   end subroutine get_cmdline

end module mod_cmdline
