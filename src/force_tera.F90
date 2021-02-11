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
   public :: force_tera
#ifdef USE_MPI
   ! TODO: Move handle_mpi_error to a dedicated MPI module
   public :: handle_mpi_error
   public :: get_tc_communicator
   public :: finalize_terachem, initialize_tc_servers, initialize_terachem_interface
#endif
   save

CONTAINS

subroutine force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
   use mod_const,    only: DP, ANG
   use mod_utils,    only: abinerror
   use mod_general,  only: iqmmm
   use mod_interfaces, only: oniom
   real(DP),intent(in)     ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)  ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)  ::  eclas
   integer,intent(in)      ::  walkmax
   integer  :: iw, itera
   integer :: OMP_GET_THREAD_NUM

! DHnote: we cannot use Niklasson's propagator in TC if nwalk > 1
! This is a responsibility of the user

   if(modulo(walkmax, nteraservers).ne.0)then
      write(*,*)'ERROR: Parameter "nwalk" must be divisible by "nteraservers"!'
      call abinerror("force_tera")
   end if


   itera = 1

! NOTE: Parallelization accross TeraChem servers
!$OMP PARALLEL DO PRIVATE(itera)
   do iw=1, walkmax

      ! map OMP thread to TC server
!$    itera = OMP_GET_THREAD_NUM() + 1


#ifdef USE_MPI
      call send_tera(x, y, z, iw, tc_comms(itera))

      call receive_tera(fx, fy,fz, eclas, iw, walkmax, tc_comms(itera))
#endif

      ! ONIOM was not yet tested!!
      if (iqmmm.eq.1) then
         call oniom(x, y, z, fx, fy, fz, eclas, iw)
      end if

   end do
!$OMP END PARALLEL DO

end subroutine force_tera


#ifdef USE_MPI

subroutine send_tera(x, y, z, iw, tc_comm)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, iremd, my_rank
   use mod_system,only: names
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror
   real(DP),intent(in)     ::  x(:,:),y(:,:),z(:,:)
   integer,intent(in)      ::  iw, tc_comm
   real(DP) :: coords(3, size(x,1) )
   character(len=2) :: names_qm(size(x,1)+6)
   integer  :: ierr, iat

   do iat=1,natqm
      coords(1,iat) = x(iat,iw)/ANG
      coords(2,iat) = y(iat,iw)/ANG
      coords(3,iat) = z(iat,iw)/ANG
   end do

   ! We send these data to TC each step
   !call send_natom(natqm)
   !call send_atom_types_and_scrdir(names)
   !call send_coordinates(x, y, z)


   ! Send natqm and the type of each qmatom
   if (idebug.gt.1) then
      write(6,'(/, a, i0)') 'Sending natqm = ', natqm
      call flush(6)
   end if
   call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, tc_comm, ierr )
   call handle_mpi_error(ierr)

   do iat=1,natqm
      names_qm(iat) = names(iat)
   end do
   if ( idebug.gt.1 ) then
      write(6,'(/,a)') 'Sending QM atom types: '
      write(*,*)(names_qm(iat), iat=1,natqm)
      call flush(6)
   end if
   ! DH WARNING: this will not work for iw>199
   ! not really tested for iw>99
   ! TODO: refactor this mess
   write(names_qm(natqm+1),'(A2)')'++'
   write(names_qm(natqm+2),'(A2)')'sc'
   if(iw.gt.99)then
      write(names_qm(natqm+3),'(A1,I1)')'r',1
      write(names_qm(natqm+4),'(I2.2)')iw-100
   else
      write(names_qm(natqm+3),'(A1,I1)')'r',0
      write(names_qm(natqm+4),'(I2.2)')iw
   end if

   ! REMD HACK
   if (iremd.eq.1)then
      write(names_qm(natqm+5),'(I2.2)')my_rank
      write(names_qm(natqm+6),'(A2)')'++'
   else
      write(names_qm(natqm+5),'(A2)')'++'
   end if

   call MPI_Send( names_qm, 2*natqm+12, MPI_CHARACTER, 0, 2, tc_comm, ierr )
   call handle_mpi_error(ierr)

   ! Send QM coordinate array
   if ( idebug > 1 ) then
      write(6,'(a)') 'Sending QM coords: '
      do iat=1, natqm
         write(6,*) 'Atom ',iat,': ',coords(:,iat)
         call flush(6)
      end do 
   end if
   call MPI_Send( coords, natqm*3, MPI_DOUBLE_PRECISION, 0, 2, tc_comm, ierr ) 
   call handle_mpi_error(ierr)

   ! NOT IMPLEMENTED !
   !if (natmm_tera > 0) then
   !   call send_mm_data(x, y, z, iw, tc_comm)
   !end if 
end subroutine send_tera


! QM/MM via TC-MPI interface is currently not
! implemented so excluding this code from compilation.
#if 0
subroutine send_mm_data(x, y, z, iw, tc_comm)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug
   use mod_qmmm, only: natqm
   real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
   integer,intent(in) :: iw, comm
   real(DP) :: coords(3, size(x,1) )
   integer  :: ierr, iat
   real(DP),intent(in) :: coords(:,:)

   call send_natom(natmm_tera)

   if (idebug > 1) then
      write (6, '(a)') 'Sending charges: '
   end if
   call MPI_Send(mmcharges, natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, tc_comm, ierr)
   call handle_mpi_error(ierr)

   ! Send MM point charge coordinate array
   if (idebug > 1) then
      write (6, '(a)') 'Sending MM coordinates...'
   end if
   do iat = 1, natmm_tera
      coords(1,iat) = x(iat+natqm,iw) / ANG
      coords(2,iat) = y(iat+natqm,iw) / ANG
      coords(3,iat) = z(iat+natqm,iw) / ANG
   end do
   call MPI_Send(coords, 3 * natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, tc_comm, ierr)
   call handle_mpi_error(ierr)
end subroutine send_mm_data
#endif

subroutine receive_tera(fx, fy, fz, eclas, iw, walkmax, tc_comm)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, it, nwrite
   use mod_io, only: print_charges, print_dipoles
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror
   real(DP),intent(inout)  ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)  ::  eclas
   integer,intent(in)      ::  iw, walkmax, tc_comm
   real(DP) :: qmcharges( size(fx,1) )
   real(DP) :: dxyz_all(3, size(fx,1) )
   real(DP) :: escf                 ! SCF energy
   real(DP) :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iat
   logical  :: recv_ready
   character*50 :: chsys_sleep
   ! -----------------------------------
   ! Begin receiving data from terachem
   ! -----------------------------------

   ! TODO: we need to somehow make sure that
   ! we don't wait forever if terachem crashes
   ! At this moment, this is ensured at the BASH script level.

   ! DH reduce cpu usage comming from MPI_Recv() via system call to 'sleep'.
   ! Not elegant, but MPICH apparently does not currently provide better solution.
   ! Based according to an answer here:
   ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu

   recv_ready = .false.
   write(chsys_sleep,'(A6, F10.4)')'sleep ', mpi_sleep
   do while(.not.recv_ready)
      call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, recv_ready, status, ierr)
      call handle_mpi_error(ierr)
      call system(chsys_sleep)
   end do

   ! Energy
   if (idebug > 2) then
      write(6,'(a)') 'Waiting to receive scf energy from TeraChem...'
      call flush(6)
   end if
   call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
   call handle_mpi_error(ierr)
   ! Checking for TAG=1, which means that SCF did not converge
   if (status(MPI_TAG).eq.1)then
      write(*,*)'GOT TAG 1 from TeraChem: SCF probably did not converge.'
      call abinerror('force_tera')
   end if
   if ( idebug > 1 ) then
      write(6,'(a,es15.6)') 'Received scf energy from server:', escf
      call flush(6)
   end if

   ! Charges (Mulliken or other)
   if (idebug > 2) then
      write(6,'(a)') 'Waiting to receive charges...'
   end if
   call MPI_Recv(qmcharges(:), natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
   call handle_mpi_error(ierr)
   if (modulo(it, nwrite) == 0 .and. nteraservers == 1) then
      call print_charges(qmcharges, iw)
   end if

   ! Dipole moment
   if ( idebug > 2 ) then
      write(6,'(a)') 'Waiting to receive dipole moment...'
   end if
   ! QM dipole moment
   call MPI_Recv( dipmom(:,1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr )
   call handle_mpi_error(ierr)
   if ( idebug > 1 ) then
      write(6,'(a,4es15.6)') 'Received QM  dipole moment from server:', dipmom(:,1)
      call flush(6)
   end if
   ! TODO: Attach dipoles to global electronic structure type
   ! and print them elsewhere. Right now when we run concurrent
   ! TC servers, the printing is not deterministic.
   if (modulo(it, nwrite) == 0 .and. nteraservers == 1) then
      call print_dipoles(dipmom(:,1), iw, 1)
   end if

   ! MM dipole moment, disabled for now
!   call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr )
!   call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr )
   
   ! QM gradients
   if (idebug > 1) then
      write(*,'(A)') 'Waiting to receive gradients...'
   end if
   call MPI_Recv(dxyz_all, 3 * natqm, MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
   call handle_mpi_error(ierr)
   if ( idebug > 1 ) then
      write (6, '(A)') 'Received the following gradients from server:'
      do iat = 1, natqm
         write (6, *) 'Atom ',iat, ': ',dxyz_all(:,iat)
      end do
      call flush(6)
   end if

   do iat = 1, natqm
      fx(iat,iw) = -dxyz_all(1,iat)
      fy(iat,iw) = -dxyz_all(2,iat)
      fz(iat,iw) = -dxyz_all(3,iat)
   end do

   ! TODO: Divide by walkmax in forces.xyz
!$OMP ATOMIC
   eclas = eclas + escf / walkmax

end subroutine receive_tera

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
