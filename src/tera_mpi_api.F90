! Common subroutines for MPI interface with TeraChem
module mod_terampi
   use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
#ifdef USE_MPI
   use mpi
   use mod_mpi, only: handle_mpi_error, get_mpi_error_string
#endif
   implicit none
   private
   integer, parameter :: MPI_TAG_EXIT = 0
   integer, parameter :: MPI_TAG_ERROR = 13
   integer, parameter :: MPI_TAG_GRADIENT = 2
   ! Default tag that we send to TC
   integer, parameter :: TC_TAG = MPI_TAG_GRADIENT
   ! By default, take port name from a file
   character(len=*), parameter :: TC_PORT_FILE_NAME = 'port.txt.'

   integer :: nteraservers = 1
   ! How long do we wait for TC port [seconds]
   real(DP) :: max_wait_time = 30
   ! Sleep interval (microseconds) while waiting for TC calculation to finish.
   integer :: mpi_milisleep = 50

   public :: nteraservers
   public :: mpi_milisleep, max_wait_time
   public :: TC_TAG
#ifdef USE_MPI
   integer, allocatable :: tc_comms(:)
   logical, allocatable :: communication_established(:)

   public :: get_tc_communicator
   public :: wait_for_terachem
   public :: send_natom, send_atom_types_and_scrdir, send_coordinates
#endif
   public :: initialize_tc_servers, initialize_terachem_interface, finalize_terachem
   save

contains

#ifdef USE_MPI

   subroutine initialize_terachem_interface(tc_server_name)
      use mod_general, only: nwalk, iremd
      use mod_mpi, only: get_mpi_rank
      character(len=*) :: tc_server_name
      integer :: i, ierr

      if (nwalk > 1) then
         write (stdout, '(A)') 'WARNING: You are using PIMD with direct TeraChem interface.'
         write (stdout, '(A)') 'You should have "integrator regular" in the TeraChem input file'
      end if

      write (stdout, '(A,I0)') 'Number of TeraChem servers: ', nteraservers

      if (mpi_milisleep < 0) then
         mpi_milisleep = 0
      end if

      allocate (tc_comms(nteraservers))
      allocate (communication_established(nteraservers))
      tc_comms = MPI_COMM_NULL
      communication_established = .false.

      ! Setting MPI_ERRORS_RETURN error handler allows us to retry
      ! failed MPI_LOOKUP_NAME() call. It also allows us
      ! to send the exit signal to TeraChem upon encountering an error.
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      call handle_mpi_error(ierr)

      if (iremd == 1) then
         ! TODO: Make this more general so that the number of
         ! TC servers can be lower than number of replicas.
         ! Right now, we hardcode one TC server per replica
         ! and we ignore the nteraservers setting in the ABIN input.
         nteraservers = 1
         call connect_tc_server(trim(tc_server_name), get_mpi_rank() + 1)
      else
         ! Parallel calls to MPI_Comm_connect often lead to segfault,
         ! maybe a bug in MPICH. Commenting out until we debug further.
         !!$OMP PARALLEL DO PRIVATE(i)
         do i = 1, nteraservers
            call connect_tc_server(trim(tc_server_name), i)
         end do
         !!$OMP END PARALLEL DO
      end if
   end subroutine initialize_terachem_interface

   subroutine connect_tc_server(tc_server_name, itera)
      use mod_general, only: iremd
      character(len=*), intent(in) :: tc_server_name
      integer, intent(in) :: itera
      character(len=MPI_MAX_PORT_NAME) :: port_name
      integer :: ierr, newcomm
      character(len=1024) :: server_name
      character(len=1024) :: portfile

      if (tc_server_name /= '') then

         write (server_name, '(A,I0)') trim(adjustl(tc_server_name))//'.', itera
         call lookup_port_via_nameserver(trim(server_name), port_name)

      else

         write (portfile, '(A,I0)') TC_PORT_FILE_NAME, itera
         call read_tc_port_from_file(trim(portfile), port_name)

      end if

      write (stdout, '(2A)') 'Found TeraChem port: ', trim(port_name)
      write (stdout, '(A)') 'Establishing connection...'
      call flush (OUTPUT_UNIT)

      ! Establish new communicator via port name
      call MPI_Comm_connect(trim(port_name), MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
      call handle_mpi_error(ierr)
      write (stdout, '(A)') 'Connection established!'

      ! This is a horrible hack :-(
      if (iremd == 1) then
         tc_comms(1) = newcomm
         communication_established(1) = .true.
      else
         tc_comms(itera) = newcomm
         communication_established(itera) = .true.
      end if
   end subroutine connect_tc_server

   integer function get_tc_communicator(itera) result(tc_comm)
      integer, intent(in) :: itera
      tc_comm = tc_comms(itera)
   end function get_tc_communicator

   ! Look for server_name via MPI nameserver, get port name
   subroutine lookup_port_via_nameserver(server_name, port_name)
      use mod_general, only: idebug
      character(len=*), intent(in) :: server_name
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      real(DP) :: timer
      integer :: ierr

      port_name = ''

      write (stdout, '(2A)') 'Looking up TeraChem server under name:', server_name
      call flush (OUTPUT_UNIT)

      timer = MPI_WTIME()

      do

         call MPI_LOOKUP_NAME(server_name, MPI_INFO_NULL, port_name, ierr)
         if (ierr == MPI_SUCCESS) then
            ! Workaround for a bug in hydra_nameserver for MPICH versions < 3.3
            if (len_trim(port_name) == 0) then
               write (stdout, '(A)') 'Found empty port, retrying...'
            else
               exit
            end if
         end if

         ! Timeout after max_wait_time seconds
         if ((MPI_WTIME() - timer) > max_wait_time) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Server name '//server_name//' not found.')
         end if

         ! Let's wait a bit since too many calls
         ! to MPI_LOOKUP_NAME() can crash the hydra_nameserver process
         if (idebug > 1) then
            write (stdout, '(A)') 'Waiting for TC port'
         end if
         call milisleep(500)

      end do

   end subroutine lookup_port_via_nameserver

   ! Read TeraChem port from a file.
   ! TeraChem prints it's port to STDOUT, where the launch script
   ! can grep it into a file, which is read here.
   ! This is more portable then the nameserver approach,
   ! but have to rely on HDD.
   subroutine read_tc_port_from_file(portfile, port_name)
      character(len=*), intent(in) :: portfile
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      integer :: iunit, iost
      real(DP) :: timer

      write (stdout, '(A)') 'Reading TeraChem port name from file '//portfile
      port_name = ''
      timer = MPI_WTIME()

      do
         open (newunit=iunit, file=portfile, action="read", status="old", iostat=iost)
         if (iost == 0) then
            exit
         end if

         if ((MPI_WTIME() - timer) > max_wait_time) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Could not open file '//portfile)
         end if

         write (stdout, '(A)') 'WARNING: Cannot open file '//portfile
         call milisleep(500)

      end do

      read (iunit, '(A)', iostat=iost) port_name
      if (iost /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Could not read file '//portfile)
      end if

      close (iunit, status='delete')
   end subroutine read_tc_port_from_file

   subroutine initialize_tc_servers()
      use mod_qmmm, only: natqm
      use mod_system, only: names
      integer :: itera

!$OMP PARALLEL DO
      do itera = 1, nteraservers
         call send_natom(natqm, tc_comms(itera))

         call send_atom_types_and_scrdir(names, natqm, 0, tc_comms(itera), .false.)
      end do
!$OMP END PARALLEL DO
   end subroutine initialize_tc_servers

   subroutine finalize_terachem(abin_error_code)
      integer, intent(in) :: abin_error_code
      integer :: itera, ierr, empty, mpi_tag

      ! To make sure we attempt to send MPI_TAG_EXIT to all servers,
      ! we will ignore MPI errors.
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)

      mpi_tag = MPI_TAG_EXIT
      if (abin_error_code /= 0) then
         mpi_tag = MPI_TAG_ERROR
      end if

      do itera = 1, nteraservers

         if (.not.communication_established(itera)) cycle

         write (stdout, '(A,I0)') 'Shutting down TeraChem server id = ', itera

         call MPI_Send(empty, 0, MPI_INTEGER, 0, mpi_tag, tc_comms(itera), ierr)

         if (ierr /= MPI_SUCCESS) then
            write (stderr, '(A,I0)') 'MPI ERROR during shutdown of TeraChem server id = ', itera
            write (stderr, '(A)') 'Verify manually that the TeraChem server was terminated.'
            write (stderr, *) get_mpi_error_string(ierr)
         end if

         call MPI_Comm_free(tc_comms(itera), ierr)
         if (ierr /= MPI_SUCCESS) then
            write (stderr, *) get_mpi_error_string(ierr)
         end if

      end do

      deallocate (tc_comms, communication_established)
   end subroutine finalize_terachem

   subroutine wait_for_terachem(tc_comm)
      integer, intent(in) :: tc_comm
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      logical :: ready
      ! TODO: we need to somehow make sure that
      ! we don't wait forever if TeraChem crashes.
      ! At this moment, this is ensured at the BASH script level.

      if (mpi_milisleep <= 0) return

      ! The idea here is to reduce the CPU usage of MPI_Recv() by taking a brief nap.
      ! In most MPI implementations, MPI_Recv() is actively polling the other end
      ! (in this case TeraChem) and consumes a whole CPU core. That's clearly wasteful,
      ! since we're typically waiting for a long time for the ab initio result.

      ! Some implementation provide an option to change this behaviour,
      ! but I didn't figure out any for MPICH.
      ! Based according to an answer here:
      ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu
      ready = .false.
      do while (.not. ready)
         call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, ready, status, ierr)
         call handle_mpi_error(ierr)
         call milisleep(mpi_milisleep)
      end do
   end subroutine wait_for_terachem

   subroutine milisleep(milisec)
      use, intrinsic :: iso_c_binding, only: c_int, c_int32_t
      use mod_interfaces, only: usleep
      integer :: milisec
      integer(kind=c_int32_t) :: usec
      integer(kind=c_int) :: c_err

      usec = int(milisec * 1000, c_int)
      ! TODO: See usleep manpage, we probably should not sleep more than a second
      c_err = usleep(usec)
      if (c_err /= 0) then
         write (stderr, *) "usleep returned an error: ", c_err
         call fatal_error(__FILE__, __LINE__, "usleep failed!")
      end if
   end subroutine milisleep

   subroutine append_scrdir_name(buffer, offset, iw, remd_replica)
      use mod_general, only: iremd
      character(len=*), intent(inout) :: buffer
      integer, intent(in) :: offset
      integer, intent(in) :: iw
      integer, intent(in) :: remd_replica
      integer :: i

      i = offset
      write (buffer(i + 1:i + 12), '(A8, I4.4)') '++scrdir', iw
      if (iremd == 1) then
         write (buffer(i + 13:i + 23), '(A4, I4.4, A2)') 'rank', remd_replica, '++'
      else
         write (buffer(i + 13:i + 14), '(A2)') '++'
      end if
   end subroutine append_scrdir_name

   subroutine send_natom(num_atom, tc_comm)
      use mod_general, only: idebug
      integer, intent(in) :: num_atom
      integer, intent(in) :: tc_comm
      integer :: ierr

      ! Send natqm and the type of each qmatom
      if (idebug > 1) then
         write (stdout, '(A, I0)') 'Sending number of atoms  = ', num_atom
         call flush (OUTPUT_UNIT)
      end if
      call MPI_Send(num_atom, 1, MPI_INTEGER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_natom

   subroutine send_atom_types_and_scrdir(at_names, num_atom, iw, tc_comm, send_scrdir)
      use mod_general, only: idebug
      use mod_mpi, only: get_mpi_rank
      character(len=2), intent(in) :: at_names(:)
      logical, intent(in) :: send_scrdir
      integer, intent(in) :: num_atom
      integer, intent(in) :: iw
      integer, intent(in) :: tc_comm
      integer, parameter :: MAX_SCRDIR_LEN = 30
      character(len=2*num_atom + MAX_SCRDIR_LEN) :: buffer
      integer :: ierr, offset, iat
      integer :: num_char

      buffer = ''
      offset = 1
      do iat = 1, num_atom
         write (buffer(offset:offset + 1), '(A2)') at_names(iat)
         offset = offset + 2
      end do
      num_char = num_atom * 2

      if (send_scrdir) then
         call append_scrdir_name(buffer, num_atom * 2, iw, get_mpi_rank())
         num_char = len_trim(buffer)
      end if

      if (idebug > 1) then
         write (stdout, '(A)') 'Sending QM atom types: '
         write (stdout, '(A)') trim(buffer)
         call flush (OUTPUT_UNIT)
      end if

      call MPI_Send(buffer, num_char, MPI_CHARACTER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_atom_types_and_scrdir

   subroutine send_coordinates(x, y, z, num_atom, iw, tc_comm, units)
      use mod_general, only: idebug
      use mod_const, only: ANG
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: num_atom
      integer, intent(in) :: iw
      integer, intent(in) :: tc_comm
      character(len=*), intent(in) :: units
      real(DP), allocatable :: coords(:, :)
      integer :: ierr, iat

      allocate (coords(3, num_atom))

      ! Amber MPI interface
      if (units == 'angstrom') then
         do iat = 1, num_atom
            coords(1, iat) = x(iat, iw) / ANG
            coords(2, iat) = y(iat, iw) / ANG
            coords(3, iat) = z(iat, iw) / ANG
         end do
         ! FMS / Surface Hopping interface
      else if (units == 'bohr') then
         do iat = 1, num_atom
            coords(1, iat) = x(iat, iw)
            coords(2, iat) = y(iat, iw)
            coords(3, iat) = z(iat, iw)
         end do
      else
         call fatal_error(__FILE__, __LINE__, 'Incorrect units in send_coordinates')
      end if

      if (idebug > 1) then
         write (stdout, '(A)') 'Sending QM coords: '
         do iat = 1, num_atom
            write (stdout, *) 'Atom ', iat, ': ', coords(:, iat)
            call flush (OUTPUT_UNIT)
         end do
      end if
      call MPI_Send(coords, num_atom * 3, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_coordinates

#else

   subroutine initialize_tc_servers()
      use mod_error, only: not_compiled_with
      call not_compiled_with('MPI')
   end subroutine initialize_tc_servers

   subroutine initialize_terachem_interface(tc_server_name)
      use mod_error, only: not_compiled_with
      character(len=*), intent(in) :: tc_server_name
      character(len=len(tc_server_name)) :: dummy
      dummy = tc_server_name
      call not_compiled_with('MPI')
   end subroutine initialize_terachem_interface

   ! This must be a no-op, since it is called from finish()
   subroutine finalize_terachem(error_code)
      integer, intent(in) :: error_code
      integer :: i
      i = error_code
   end subroutine finalize_terachem

! USE_MPI
#endif

end module mod_terampi
