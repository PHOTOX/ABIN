module mod_terampi
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD.
! Perform MPI communications with terachem.
!
! Currently supports:
! pure QM and ONIOM
! Adapted from Sander (Amber14)
!
! Original Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
! Modified by Daniel Hollas hollas@vscht.cz
! ----------------------------------------------------------------
   use mod_const, only: DP
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private
   integer, parameter :: MPI_TAG_ERROR = 13, MPI_TAG_EXIT = 0
   ! By default, take port name from a file
   character(len=*), parameter :: TC_PORT_FILE_NAME = 'port.txt.'
   character(len=1024) :: tc_server_name = ''
   integer, allocatable ::  tc_comms(:)

   integer :: nteraservers = 1
   ! How long do we wait for TC port [seconds]
   real(DP) :: max_wait_time = 30
   ! Sleep interval while waiting for TC calculation to finish.
   real(DP) :: mpi_sleep = 0.05

   public :: tc_server_name
   public :: nteraservers
   public :: mpi_sleep, max_wait_time
#ifdef USE_MPI
   ! TODO: Move handle_mpi_error to a dedicated MPI module
   public :: handle_mpi_error
   public :: get_tc_communicator
   public :: finalize_terachem, initialize_tc_servers, initialize_terachem_interface
#endif
   save

CONTAINS

#ifdef USE_MPI

   subroutine initialize_terachem_interface()
      use mod_general, only: nwalk
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
         call connect_tc_server(i)
      end do
      !$OMP END PARALLEL DO

   end subroutine initialize_terachem_interface

   subroutine connect_tc_server(itera)
   use mod_utils, only: abinerror
   ! TODO: Figure out how to handle REMD
!   use mod_general, only: iremd, my_rank
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
      use mod_general, only: idebug
      use mod_qmmm,  only: natqm
      use mod_system,only: names
      use mod_utils, only: abinerror
      integer :: ierr, itera
 
      !$OMP PARALLEL DO PRIVATE(ierr, itera)
      do itera = 1, nteraservers
         if (idebug > 0) then
            write (*, *) 'Sending initial number of QM atoms to TeraChem.'
         end if
         call MPI_Send(natqm, 1, MPI_INTEGER, 0, 2, tc_comms(itera), ierr)
         call handle_mpi_error(ierr)
 
         if (idebug > 0) then
            write (*, *) 'Sending initial QM atom names to TeraChem.'
         end if
         call MPI_Send(names, 2*natqm, MPI_CHARACTER, 0, 2, tc_comms(itera), ierr)
         call handle_mpi_error(ierr)
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

! USE_MPI
#endif

end module mod_terampi
