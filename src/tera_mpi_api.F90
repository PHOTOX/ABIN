! Common subroutines for MPI interface with TeraChem
module mod_terampi
   use mod_const, only: DP
#ifdef USE_MPI
   use mpi
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
   ! Sleep interval while waiting for TC calculation to finish.
   real(DP) :: mpi_sleep = 0.05

   public :: nteraservers
   public :: mpi_sleep, max_wait_time
   public :: TC_TAG
#ifdef USE_MPI
   integer, allocatable ::  tc_comms(:)

   ! TODO: Move handle_mpi_error to a dedicated MPI module
   public :: handle_mpi_error
   public :: get_tc_communicator
   public :: wait_for_terachem
   public :: finalize_terachem
   public :: send_natom, send_atom_types_and_scrdir, send_coordinates
#endif
   public :: initialize_tc_servers, initialize_terachem_interface
   save

CONTAINS

#ifdef USE_MPI

   subroutine initialize_terachem_interface(tc_server_name)
      use mod_general, only: nwalk
      character(len=*) :: tc_server_name
      integer :: i, ierr

      if (nwalk > 1) then
         write (*, '(A)') 'WARNING: You are using PIMD with direct TeraChem interface.'
         write (*, '(A)') 'You should have "integrator regular" in the TeraChem input file'
      end if
      write (*, '(A,I0)') 'Number of TeraChem servers: ', nteraservers

      if(mpi_sleep <= 0) then
         write(*,*)'WARNING: Parameter "mpi_sleep" must be positive!'
         write(*,*)'Setting it back to default value'
         mpi_sleep = 0.05D0
      end if

      allocate (tc_comms(nteraservers))
      tc_comms = MPI_COMM_NULL

      ! Setting MPI_ERRORS_RETURN error handler allows us to retry
      ! failed MPI_LOOKUP_NAME() call. It also allows us
      ! to send the exit signal to TeraChem upon encoutering an error.
 
      ! It might also be a good idea to write our own error handler by MPI_Errorhandler_Create()
      ! https://www.open-mpi.org/doc/current/man3/MPI_Comm_create_errhandler.3.php
      ! so that we don't have to call handle_mpi_error() after each MPI call.
      ! This error handler should call abinerror() and if possible should try send
      ! the error shutdown MPI_Send to TC (though we'd need to make sure we don't
      ! enter some weird endless loop!).
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      call handle_mpi_error(ierr)

      ! Connect to all TC servers concurrently.
      !$OMP PARALLEL DO
      do i = 1, nteraservers
         call connect_tc_server(tc_server_name, i)
      end do
      !$OMP END PARALLEL DO
   end subroutine initialize_terachem_interface


   subroutine connect_tc_server(tc_server_name, itera)
      use mod_utils, only: abinerror
      ! TODO: Figure out how to handle REMD
      ! use mod_general, only: iremd, my_rank
      character(len=*) :: tc_server_name
      integer, intent(in) :: itera
      character(len=MPI_MAX_PORT_NAME) :: port_name
      integer :: ierr, newcomm
      character(len=1024) :: server_name
      character(len=1024) :: portfile
 
      if (tc_server_name /= '') then
 
         write (server_name, '(A,I0)')trim(adjustl(tc_server_name))//'.', itera
         call lookup_port_via_nameserver(trim(server_name), port_name)
 
      else
 
         write (portfile, '(A,I0)') TC_PORT_FILE_NAME, itera
         call read_tc_port_from_file(trim(portfile), port_name)
 
      end if
 
      write (6, '(2A)') 'Found TeraChem port: ', trim(port_name)
      write (6, '(A)') 'Establishing connection...'
      call flush(6)

      ! Establish new communicator via port name
      call MPI_Comm_connect(trim(port_name), MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
      call handle_mpi_error(ierr)
      write(6, '(A)') 'Connection established!'
 
      tc_comms(itera) = newcomm
   end subroutine connect_tc_server


   integer function get_tc_communicator(itera) result(tc_comm)
      integer, intent(in) :: itera
      tc_comm = tc_comms(itera)
   end function get_tc_communicator


   ! Look for server_name via MPI nameserver, get port name
   subroutine lookup_port_via_nameserver(server_name, port_name)
      use mod_general, only: idebug
      use mod_utils, only: abinerror
      character(len=*), intent(in) :: server_name
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      real(DP) :: timer
      integer :: ierr

      port_name = ''

      write(*,'(2A)') 'Looking up TeraChem server under name:', server_name
      call flush(6)

      timer = MPI_WTIME()

      do

         call MPI_LOOKUP_NAME(server_name, MPI_INFO_NULL, port_name, ierr)
         if (ierr == MPI_SUCCESS) then
            ! Workaround for a bug in hydra_nameserver for MPICH versions < 3.3
            if (len_trim(port_name) == 0) then
               write(*,'(a)') 'Found empty port, retrying...'
            else
               exit
            end if
         end if

         ! Let's wait a bit since too many calls
         ! to MPI_LOOKUP_NAME() can crash the hydra_nameserver process
         if(idebug > 1)then
            write(*, '(A)')'Waiting for TC port'
         end if
         ! TODO: Try out how long should we sleep here.
         call system('sleep 0.5')
       
         ! Timeout after max_wait_time seconds
         if ( (MPI_WTIME()-timer) > max_wait_time ) then
            write (*, *) 'Server name '//server_name//' not found.'
            call abinerror("lookup_port_via_nameserver")
         end if

      end do

   end subroutine lookup_port_via_nameserver


   ! Read TeraChem port from a file.
   ! TeraChem prints it's port to STDOUT, where the launch script
   ! can grep it into a file, which is read here.
   ! This is more portable then the nameserver approach,
   ! but have to rely on HDD.
   subroutine read_tc_port_from_file(portfile, port_name)
      use mod_utils, only: abinerror
      character(len=*), intent(in) :: portfile
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      integer :: iunit, iost
      real(DP) :: timer

      write (*, '(A)') 'Reading TeraChem port name from file '//portfile
      port_name = ''
      timer = MPI_WTIME()

      do
         open (newunit=iunit, file=portfile, action="read", status="old", iostat=iost)
         if (iost == 0) then
            exit
         end if

         write (*, '(A)') 'WARNING: Cannot open file '//portfile
         call system('sleep 0.5')

         if ( (MPI_WTIME()-timer) > max_wait_time) then
            write (*, '(A)') 'ERROR: Could not open file '//portfile
            call abinerror('read_tc_port_from_file')
         end if

      end do

      read (iunit, '(A)', iostat=iost) port_name
      if (iost /= 0) then
          write (*, '(A)') 'ERROR reading file '//portfile
          call abinerror('read_tc_port_from_file')
      end if

      close(iunit, status='delete')
   end subroutine read_tc_port_from_file


   subroutine initialize_tc_servers()
      use mod_qmmm,  only: natqm
      use mod_system,only: names
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
      integer :: itera, ierr, empty
 
      ! Make sure we send MPI_TAG_EXIT to all servers.
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      do itera=1, nteraservers
 
         write (*, '(A,I0)') 'Shutting down TeraChem server id = ', itera
         if (abin_error_code == 0) then
            call MPI_Send(empty, 0, MPI_INTEGER, 0, MPI_TAG_EXIT, tc_comms(itera), ierr)
         else
            call MPI_Send(empty, 0, MPI_INTEGER, 0, MPI_TAG_ERROR, tc_comms(itera), ierr)
         end if
         if (ierr /= MPI_SUCCESS) then
            write(*,'(A,I0)')'I got a MPI Error when I tried to shutdown TeraChem server id = ', itera
            write(*,'(A)')'Verify manually that the TeraChem server was terminated.'
            call print_mpi_error(ierr)
         end if
 
         call MPI_Comm_free(tc_comms(itera), ierr)
         if (ierr /= MPI_SUCCESS) then
            call print_mpi_error(ierr)
         end if
 
      end do
 
      deallocate (tc_comms)
   end subroutine finalize_terachem

   subroutine print_mpi_error(mpi_err)
      character(len=MPI_MAX_ERROR_STRING) :: error_string
      integer, intent(in) :: mpi_err
      integer :: ierr, result_len

      call MPI_Error_string(mpi_err, error_string, result_len, ierr)
      if (ierr == MPI_SUCCESS) then
         write (*, '(A)') trim(error_string)
      end if
   end subroutine print_mpi_error

   subroutine handle_mpi_error(mpi_err)
      use mod_utils, only: abinerror
      integer, intent(in) :: mpi_err

      if (mpi_err /= MPI_SUCCESS) then
         call print_mpi_error(mpi_err)
         call abinerror('handle_mpi_error')
      end if
   end subroutine handle_mpi_error


   subroutine wait_for_terachem(tc_comm)
      integer, intent(in) :: tc_comm
      character(len=20) :: chsys_sleep
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      logical :: ready
      ! TODO: we need to somehow make sure that
      ! we don't wait forever if TeraChem crashes.
      ! At this moment, this is ensured at the BASH script level.
 
      ! The idea here is to reduce the CPU usage of MPI_Recv() via system call to 'sleep'.
      ! In most MPI implementations, MPI_Recv() is actively polling the other end
      ! (in this case TeraChem) and consumes a whole CPU core. That's clearly wasteful,
      ! since we're waiting for a long time for the ab initio result.

      ! Some implementation provide an option to change this behaviour,
      ! but I didn't figure out any for MPICH so we rather inelegantly call system sleep.
      ! Based according to an answer here:
      ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu
      ready = .false.
      write (chsys_sleep, '(A6, F10.4)') 'sleep ', mpi_sleep
      do while(.not.ready)
         call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, ready, status, ierr)
         call handle_mpi_error(ierr)
         call system(trim(chsys_sleep))
      end do
   end subroutine wait_for_terachem

   subroutine append_scrdir_name(buffer, offset, iw, remd_replica)
      use mod_general, only: iremd
      character(len=*), intent(inout) :: buffer
      integer, intent(in) :: offset
      integer, intent(in) :: iw
      integer, intent(in) :: remd_replica
      integer :: i

      i = offset
      write (buffer(i+1:i+12), '(A8, I4.4)') '++scrdir', iw
      if (iremd == 1) then
         write (buffer(i+13:i+23), '(A4, I4.4, A2)') 'rank', remd_replica, '++'
      else
         write (buffer(i+13:i+14), '(A2)') '++'
      end if
   end subroutine append_scrdir_name

   subroutine send_natom(num_atom, tc_comm)
      use mod_general, only: idebug
      integer, intent(in) :: num_atom
      integer, intent(in) :: tc_comm
      integer :: ierr

      ! Send natqm and the type of each qmatom
      if (idebug > 1) then
         write(6,'(A, I0)') 'Sending number of atoms  = ', num_atom
         call flush(6)
      end if
      call MPI_Send(num_atom, 1, MPI_INTEGER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_natom


   subroutine send_atom_types_and_scrdir(at_names, num_atom, iw, tc_comm, send_scrdir)
      use mod_general, only: my_rank, idebug
      character(len=2), intent(in) :: at_names(:)
      logical, intent(in) :: send_scrdir
      integer, intent(in) :: num_atom
      integer, intent(in) :: iw
      integer, intent(in) :: tc_comm
      integer, parameter :: MAX_SCRDIR_LEN = 30
      character(len=2 * num_atom + MAX_SCRDIR_LEN) :: buffer
      integer :: ierr, offset, iat
      integer :: num_char

      buffer = ''
      offset = 1
      do iat = 1, num_atom
         write(buffer(offset:offset+1), '(A2)') at_names(iat)
         offset = offset + 2
      end do
      num_char = num_atom * 2

      if (send_scrdir) then
         call append_scrdir_name(buffer, num_atom * 2, iw, my_rank)
         num_char = len_trim(buffer)
      end if
 
      if (idebug > 1) then
         write (6, '(A)') 'Sending QM atom types: '
         write (*, '(A)') trim(buffer)
         call flush(6)
      end if
 
      call MPI_Send(buffer, num_char, MPI_CHARACTER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_atom_types_and_scrdir


   subroutine send_coordinates(x, y, z, num_atom, iw, tc_comm)
      use mod_general, only: idebug
      use mod_const, only: ANG
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: num_atom
      integer, intent(in) :: iw
      integer, intent(in) :: tc_comm
      real(DP), allocatable :: coords(:, :)
      integer :: ierr, iat

      allocate (coords(3, num_atom))
      do iat = 1, num_atom
         coords(1, iat) = x(iat, iw) / ANG
         coords(2, iat) = y(iat, iw) / ANG
         coords(3, iat) = z(iat, iw) / ANG
      end do

      if (idebug > 1) then
         write(6, '(A)') 'Sending QM coords: '
         do iat = 1, num_atom
            write (6, *) 'Atom ', iat, ': ', coords(:, iat)
            call flush(6)
         end do 
      end if
      call MPI_Send(coords, num_atom * 3, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_coordinates

#else

   subroutine initialize_tc_servers()
      use mod_utils, only: not_compiled_with
      call not_compiled_with('MPI', 'initialize_tc_servers')
   end subroutine initialize_tc_servers
      
   subroutine initialize_terachem_interface(tc_server_name)
      use mod_utils, only: not_compiled_with
      character(len=*), intent(in) :: tc_server_name
      write (*, *) 'TC_SERVER_NAME=', tc_server_name
      call not_compiled_with('MPI', 'initialize_terachem_interface')
   end subroutine initialize_terachem_interface

! USE_MPI
#endif

end module mod_terampi
