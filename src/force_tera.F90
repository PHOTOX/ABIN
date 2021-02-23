module mod_terampi
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD.
! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
! So far, I was not able to make it work with OpenMPI.
! (but now that we use file based tera_port, it should work as well)
!
! Currently supports:
! pure QM and ONIOM
! Adapted from Sander (Amber14)
!
! Original Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
! Modified by Dan Hollas hollas@vscht.cz
! Date: September 2014 - September 2015
! ----------------------------------------------------------------
   use mod_const, only: DP
   implicit none
   private
   integer, parameter :: MAXTERASERVERS = 9
   integer, parameter :: MPI_TAG_ERROR = 13, MPI_TAG_EXIT = 0
   ! By default, take port name from a file
   character(len=1024) :: teraport = ''
   integer :: newcomms(MAXTERASERVERS) ! Communicator, initialized in mpi_init subroutine
!  DH WARNING, initial hack, we do not support TeraChem-based QM/MM yet
   integer :: natmm_tera = 0
   integer :: nteraservers = 1
   real(DP), allocatable :: mmcharges(:)
   real(DP) :: mpi_sleep = 0.05
   public :: teraport, newcomms, mpi_sleep, nteraservers
   public :: force_tera, natmm_tera
#ifdef USE_MPI
   ! TODO: Move handle_mpi_error to a dedicated MPI module
   public :: handle_mpi_error
   public :: finalize_terachem, initialize_terachem, connect_terachem
#endif
   save

contains

   subroutine force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_const, only: DP, ANG
      use mod_utils, only: abinerror
      use mod_general, only: iqmmm
      use mod_interfaces, only: oniom
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw, itera
      integer :: OMP_GET_THREAD_NUM

! DHnote: we cannot use Niklasson's propagator in TC if nwalk > 1
! This is a responsibility of the user

      if (modulo(walkmax, nteraservers) /= 0) then
         write (*, *) 'ERROR: Parameter "nwalk" must be divisible by "nteraservers"!'
         call abinerror("force_tera")
      end if

      itera = 1

! NOTE: Parallelization accross TeraChem servers
!$OMP PARALLEL DO PRIVATE(itera)
      do iw = 1, walkmax

         ! map OMP thread to TC server
!$       itera = OMP_GET_THREAD_NUM() + 1

#ifdef USE_MPI
         call send_tera(x, y, z, iw, newcomms(itera))

         call receive_tera(fx, fy, fz, eclas, iw, walkmax, newcomms(itera))
#endif

         ! ONIOM was not yet tested!!
         if (iqmmm == 1) call oniom(x, y, z, fx, fy, fz, eclas, iw)

      end do
!$OMP END PARALLEL DO

   end subroutine force_tera

#ifdef USE_MPI

   subroutine send_tera(x, y, z, iw, newcomm)
      use mpi
      use mod_const, only: DP, ANG
      use mod_general, only: idebug, iremd, my_rank
      use mod_system, only: names
      use mod_qmmm, only: natqm
      use mod_utils, only: abinerror
      use mod_interfaces, only: oniom
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw, newcomm
      real(DP) :: coords(3, size(x, 1))
      character(len=2) :: names_qm(size(x, 1) + 6)
      integer :: ierr, iat

      do iat = 1, natqm
         coords(1, iat) = x(iat, iw) / ANG
         coords(2, iat) = y(iat, iw) / ANG
         coords(3, iat) = z(iat, iw) / ANG
      end do

      ! -----------------------------------------
      ! Begin sending data each step to terachem
      ! -----------------------------------------

      ! Send natqm and the type of each qmatom
      if (idebug > 1) then
         write (6, '(/, a, i0)') 'Sending natqm = ', natqm
         call flush (6)
      end if
      call MPI_Send(natqm, 1, MPI_INTEGER, 0, 2, newcomm, ierr)
      call handle_mpi_error(ierr)

      do iat = 1, natqm
         names_qm(iat) = names(iat)
      end do
      if (idebug > 1) then
         write (6, '(/,a)') 'Sending QM atom types: '
         write (*, *) (names_qm(iat), iat=1, natqm)
         call flush (6)
      end if
      ! DH WARNING: this will not work for iw>199
      ! not really tested for iw>99
      ! TODO: refactor this mess
      write (names_qm(natqm + 1), '(A2)') '++'
      write (names_qm(natqm + 2), '(A2)') 'sc'
      if (iw > 99) then
         write (names_qm(natqm + 3), '(A1,I1)') 'r', 1
         write (names_qm(natqm + 4), '(I2.2)') iw - 100
      else
         write (names_qm(natqm + 3), '(A1,I1)') 'r', 0
         write (names_qm(natqm + 4), '(I2.2)') iw
      end if

      ! REMD HACK
      if (iremd == 1) then
         write (names_qm(natqm + 5), '(I2.2)') my_rank
         write (names_qm(natqm + 6), '(A2)') '++'
      else
         write (names_qm(natqm + 5), '(A2)') '++'
      end if

      call MPI_Send(names_qm, 2 * natqm + 12, MPI_CHARACTER, 0, 2, newcomm, ierr)
      call handle_mpi_error(ierr)

      ! Send QM coordinate array
      if (idebug > 1) then
         write (6, '(a)') 'Sending QM coords: '
         do iat = 1, natqm
            write (6, *) 'Atom ', iat, ': ', coords(:, iat)
            call flush (6)
         end do
      end if
      call MPI_Send(coords, natqm * 3, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
      call handle_mpi_error(ierr)

      if (natmm_tera > 0) then

         do iat = 1, natmm_tera
            coords(1, iat) = x(iat + natqm, iw) / ANG
            coords(2, iat) = y(iat + natqm, iw) / ANG
            coords(3, iat) = z(iat + natqm, iw) / ANG
         end do

         ! Send natmm and the charge of each atom
         if (idebug > 1) then
            write (6, '(a, i0)') 'Sending natmm = ', natmm_tera
            call flush (6)
         end if
         call MPI_Send(natmm_tera, 1, MPI_INTEGER, 0, 2, newcomm, ierr)
         call handle_mpi_error(ierr)

         if (idebug > 1) then
            write (6, '(a)') 'Sending charges: '
         end if
         call MPI_Send(mmcharges, natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
         call handle_mpi_error(ierr)

         ! Send MM point charge coordinate array
         if (idebug > 1) then
            write (6, '(a)') 'Sending charges coords: '
         end if

         call MPI_Send(coords, 3 * natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
         call handle_mpi_error(ierr)
      end if

   end subroutine send_tera

   subroutine receive_tera(fx, fy, fz, eclas, iw, walkmax, newcomm)
      use mpi
      use mod_const, only: DP, ANG
      use mod_general, only: idebug
      use mod_qmmm, only: natqm
      use mod_utils, only: abinerror
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: iw, walkmax, newcomm
      ! TODO: make qmcharges global variable
      real(DP) :: qmcharges(size(fx, 1))
      real(DP) :: dxyz_all(3, size(fx, 1))
      real(DP) :: escf ! SCF energy
      real(DP) :: dipmom(4, 3) ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr, iat
      logical :: ltest
      character*50 :: chsys_sleep
      ! -----------------------------------
      ! Begin receiving data from terachem
      ! -----------------------------------

      ! DH TODO: we need to somehow make sure that we don't wait forever if
      ! terachem crashes
      ! At this moment, this is ensured at the bash script level

      ! DH reduce cpu usage comming from MPI_Recv() via system call to 'sleep'.
      ! Not elegant, but MPICH apparently does not currently provide better solution.
      ! Based according to an answer here:
      ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu

      ltest = .false.
      write (chsys_sleep, '(A6, F10.4)') 'sleep ', mpi_sleep
      do while (.not. ltest)
         call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, ltest, status, ierr)
         call handle_mpi_error(ierr)
         call system(chsys_sleep)
      end do

      ! Energy
      if (idebug > 2) then
         write (6, '(a)') 'Waiting to receive scf energy from TeraChem...'
         call flush (6)
      end if
      call MPI_Recv(escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
      call handle_mpi_error(ierr)
      ! Checking for TAG=1, which means that SCF did not converge
      if (status(MPI_TAG) == 1) then
         write (*, *) 'GOT TAG 1 from TeraChem: SCF probably did not converge.'
         call abinerror('force_tera')
      end if
      if (idebug > 1) then
         write (6, '(a,es15.6)') 'Received scf energy from server:', escf
         call flush (6)
      end if

      ! Charges (Mulliken or other)
      if (idebug > 2) then
         write (6, '(a)') 'Waiting to receive charges...'
      end if
      call MPI_Recv(qmcharges(:), natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
      call handle_mpi_error(ierr)
      if (idebug > 2) then
         write (6, '(a)') 'Received the following charges from server:'
         do iat = 1, natqm
            write (6, *) 'Atom ', iat, ': ', qmcharges(iat)
         end do
         call flush (6)
      end if

      ! Dipole moment
      if (idebug > 2) then
         write (6, '(a)') 'Waiting to receive dipole moment...'
      end if
      ! QM dipole moment
      call MPI_Recv(dipmom(:, 1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
      call handle_mpi_error(ierr)
      if (idebug > 1) then
         write (6, '(a,4es15.6)') 'Received QM  dipole moment from server:', dipmom(:, 1)
         call flush (6)
      end if
      ! MM dipole moment, disabled for now
!   call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
!   if ( idebug > 1 ) then
!      write(6,'(a,4es15.6)') 'Received MM  dipole moment from server:', dipmom(:,2)
!      call flush(6)
!   end if
      ! TOT dipole moment
!   call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
!   if ( idebug > 1 ) then
!      write(6,'(a,4es15.6)') 'Received TOT dipole moment from server:', dipmom(:,3)
!      call flush(6)
!   end if

      ! QM gradients
      if (idebug > 1) then
         write (*, '(a)') 'Waiting to receive gradients...'
      end if
      call MPI_Recv(dxyz_all, 3 * (natqm + natmm_tera), MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
      call handle_mpi_error(ierr)
      if (idebug > 1) then
         write (6, '(a)') 'Received the following gradients from server:'
         do iat = 1, natqm + natmm_tera
            write (6, *) 'Atom ', iat, ': ', dxyz_all(:, iat)
         end do
         call flush (6)
      end if
!!$OMP CRITICAL
      do iat = 1, natqm + natmm_tera
         fx(iat, iw) = -dxyz_all(1, iat)
         fy(iat, iw) = -dxyz_all(2, iat)
         fz(iat, iw) = -dxyz_all(3, iat)
      end do
!!$OMP END CRITICAL

!$OMP ATOMIC
      eclas = eclas + escf / walkmax

   end subroutine receive_tera

   subroutine connect_terachem(itera)
      use mpi
      use mod_utils, only: abinerror
      use mod_general, only: iremd, my_rank, idebug
      integer, intent(in) :: itera
      character(len=MPI_MAX_PORT_NAME) :: port_name
      integer :: ierr, newcomm, iost
      real(DP) :: timer
      character(len=1024) :: server_name, portfile
      character(len=1) :: chtera

      ! -----------------------------------
      ! Look for server_name, get port name
      ! After 60 seconds, exit if not found
      ! -----------------------------------
      port_name = ''

      if (itera > MAXTERASERVERS) then
         write (*, *) 'ERROR: We currently support only ', MAXTERASERVERS, 'TC servers!'
         write (*, *) 'Shutting down...'
         write (*, *) 'Running TC servers might not be shutdown properly!'
         call abinerror('force_tera')
      end if

      timer = MPI_WTIME()

      ! This allows us to retry failed MPI_LOOKUP_NAME() call
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      call handle_mpi_error(ierr)

      write (chtera, '(I1)') itera

      if (iremd == 1) write (chtera, '(I1)') my_rank + 1
      if (teraport /= '') then
         server_name = trim(teraport)//'.'//trim(chtera)
         write (6, '(2a)') 'Looking up TeraChem server under name:', trim(server_name)
         call flush (6)

         do

            call MPI_LOOKUP_NAME(trim(server_name), MPI_INFO_NULL, port_name, ierr)
            if (ierr == MPI_SUCCESS) then
               ! This sometimes happens, I have no idea why.
               if (len_trim(port_name) == 0) then
                  write (6, '(a)') 'Found empty port, retrying...'
                  call system('sleep 1')
               else
                  write (6, '(2a)') 'Found port: ', trim(port_name)
                  call flush (6)
                  exit
               end if
            else
               ! Let's wait a bit
               ! Too many calls to MPI_LOOKUP_NAME can crash the hydra_nameserver process
               if (idebug > 1) write (6, '(A)') 'Waiting for TC port'
               call system('sleep 1')
            end if

            if ((MPI_WTIME() - timer) > 60) then ! Time out after 60 seconds
               write (*, *) 'Port "'//trim(server_name)//'" not found. Timed out after 60 seconds.'
               call abinerror("connect_to_terachem")
            end if

         end do

      else

         portfile = 'port.txt.'//chtera
         write (6, '(A)') 'Reading TeraChem port name from file '//trim(portfile)
         call system('sync') ! flush HDD buffer, not sure how portable this is
         open (500, file=portfile, action="read", status="old", iostat=iost)
         if (iost /= 0) then
            write (*, *) 'WARNING: Cannot open file '//trim(portfile)
            write (*, *) 'Will wait for 10s and try again...'
            call system('sleep 10')
            open (500, file=portfile, action="read", status="old")
         end if
         read (500, '(A)', iostat=iost) port_name
         if (iost /= 0) then
            write (*, *) 'ERROR when reading file '//trim(portfile)
            stop 1
         end if
         close (500)

      end if

      write (6, '(a)') 'Establishing connection to TeraChem...'
      ! ----------------------------------------
      ! Establish new communicator via port name
      ! ----------------------------------------
      call flush (6)
      call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
      call handle_mpi_error(ierr)
      write (6, '(a,i0)') 'Established a new communicator:', newcomm

      newcomms(itera) = newcomm

   end subroutine connect_terachem

   subroutine initialize_terachem()
      use mpi
      use mod_qmmm, only: natqm
      use mod_system, only: names
      use mod_utils, only: abinerror
      integer :: ierr, itera

      if (natmm_tera > 0) then
         allocate (mmcharges(natmm_tera))
      end if

      do itera = 1, nteraservers
         write (*, *) 'Sending initial number of QM atoms to TeraChem.'
         call MPI_Send(natqm, 1, MPI_INTEGER, 0, 2, newcomms(itera), ierr)
         call handle_mpi_error(ierr)
      end do

      do itera = 1, nteraservers
         write (*, *) 'Sending initial QM atom names to TeraChem.'
         call MPI_Send(names, 2 * natqm, MPI_CHARACTER, 0, 2, newcomms(itera), ierr)
         call handle_mpi_error(ierr)
      end do

      if (mpi_sleep <= 0) then
         write (*, *) 'WARNING: Parameter "mpi_sleep" must be positive!'
         write (*, *) 'Setting it back to default value'
         mpi_sleep = 0.05
      end if

   end subroutine initialize_terachem

   subroutine finalize_terachem(error_code)
      use mpi
      integer, intent(in) :: error_code
      integer :: ierr, itera
      integer :: empty = 0

      do itera = 1, nteraservers

         write (*, *) 'Shutting down TeraChem server; id=', itera
         if (error_code == 0) then
            call MPI_Send(empty, 1, MPI_INTEGER, 0, MPI_TAG_EXIT, newcomms(itera), ierr)
         else
            call MPI_Send(empty, 1, MPI_INTEGER, 0, MPI_TAG_ERROR, newcomms(itera), ierr)
         end if
         if (ierr /= 0) then
            write (*, *) 'I got a MPI Error when I tried to shutdown TeraChem server id =', itera
            write (*, *) 'Please, verify manually that the TeraChem server was terminated.'
            write (*, *) 'The error code was:', ierr
         end if

      end do
   end subroutine finalize_terachem

   ! TODO: call this after each MPI call
   subroutine handle_mpi_error(mpi_err, status, datatype, expected_count)
      use mpi
      use mod_utils, only: abinerror
      integer, intent(in) :: mpi_err
      integer, intent(in), optional :: status(MPI_STATUS_SIZE)
      integer, intent(in), optional :: datatype, expected_count
      integer :: received_count, ierr

      if (present(status) .and. status(MPI_TAG) == MPI_TAG_ERROR) then
         write (*, *) 'TeraChem sent an ERROR TAG. Exiting...'
         call abinerror('TeraChem ERROR')
      end if

      if (present(expected_count)) then
         ! Compare the length of received message and what we expected
         call MPI_GET_COUNT(status, datatype, received_count, ierr)
         if (received_count /= expected_count) then
            write (*, *) 'Received message of unexpected size'
            call abinerror('MPI ERROR')
         end if
      end if

      if (mpi_err /= 0) then
         write (*, *) 'Unspecified MPI Error, code:', ierr
         call abinerror('MPI ERROR')
      end if
   end subroutine handle_mpi_error
#endif

end module mod_terampi

