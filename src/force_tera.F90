module mod_force_tera
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
   use mod_error, only: fatal_error
   use mod_files, only: stderr
   use mod_terampi
#ifdef USE_MPI
   use mpi
   use mod_mpi, only: handle_mpi_error, check_recv_count
#endif
   implicit none
   private
   public :: force_tera
   save

contains

   subroutine force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: iqmmm
      use mod_shell_interface, only: oniom
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw, itc
      integer :: tc_comm
      integer :: OMP_GET_THREAD_NUM
      logical :: abort

      ! DHnote: we cannot use Niklasson's propagator in TC if nwalk > 1
      ! This is a responsibility of the user

      if (modulo(walkmax, nteraservers) /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Parameter "nwalk" must be divisible by "nteraservers"!')
      end if

      itc = 1
      abort = .false.

      ! Parallelization accross TeraChem servers
!$OMP PARALLEL DO PRIVATE(iw, itc, tc_comm)
      do iw = 1, walkmax

         ! See comment in force_abin() to understand this
!$OMP FLUSH(abort)
         if (.not. abort) then
            ! map OMP thread to TC server
!$          itc = OMP_GET_THREAD_NUM() + 1

#ifdef USE_MPI
            tc_comm = get_tc_communicator(itc)

            call send_tera(x, y, z, iw, tc_comm)

            call receive_tera(fx, fy, fz, eclas, iw, walkmax, tc_comm, abort)
            if (abort) cycle
#endif

            ! ONIOM was not yet tested!!
            if (iqmmm == 1) then
               call oniom(x, y, z, fx, fy, fz, eclas, iw, abort)
               if (abort) cycle
            end if

         end if

      end do
!$OMP END PARALLEL DO

      if (abort) then
         call fatal_error(__FILE__, __LINE__, &
            & 'External forces error')
      end if

   end subroutine force_tera

#ifdef USE_MPI

   subroutine send_tera(x, y, z, iw, tc_comm)
      use mod_system, only: names
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw, tc_comm

      ! We send these data to TC each step
      call send_natom(natqm, tc_comm)

      call send_atom_types_and_scrdir(names, natqm, iw, tc_comm, .true.)

      call send_coordinates(x, y, z, natqm, iw, tc_comm, 'angstrom')

      ! NOT IMPLEMENTED !
      !if (natmm_tera > 0) then
      !   call send_mm_data(x, y, z, iw, tc_comm)
      !end if
   end subroutine send_tera

! QM/MM via TC-MPI interface is currently not
! implemented so excluding this code from compilation.
#if 0
   subroutine send_mm_data(x, y, z, iw, tc_comm)
      use mod_const, only: ANG
      use mod_general, only: idebug
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw, comm
      real(DP) :: coords(3, size(x, 1))
      integer :: ierr, iat
      real(DP), intent(in) :: coords(:, :)

      call send_natom(natmm_tera)

      if (idebug > 1) then
         print '(a)', 'Sending charges: '
      end if
      call MPI_Send(mmcharges, natmm_tera, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      ! Send MM point charge coordinate array
      if (idebug > 1) then
         print '(a)', 'Sending MM coordinates...'
      end if
      do iat = 1, natmm_tera
         coords(1, iat) = x(iat + natqm, iw) / ANG
         coords(2, iat) = y(iat + natqm, iw) / ANG
         coords(3, iat) = z(iat + natqm, iw) / ANG
      end do
      call MPI_Send(coords, 3 * natmm_tera, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_mm_data
#endif

   subroutine receive_tera(fx, fy, fz, eclas, iw, walkmax, tc_comm, abort)
      use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
      use mod_const, only: ANG
      use mod_general, only: idebug, it, nwrite
      use mod_io, only: print_charges, print_dipoles
      use mod_qmmm, only: natqm
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: iw, walkmax
      integer, intent(in) :: tc_comm
      logical, intent(inout) :: abort
      real(DP) :: qmcharges(size(fx, 1))
      real(DP) :: dxyz_all(3, size(fx, 1))
      real(DP) :: escf ! SCF energy
      real(DP) :: dipmom(4, 3) ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr, iat

      call wait_for_terachem(tc_comm)

      ! Begin receiving data from TeraChem

      if (idebug > 2) then
         print '(a)', 'Waiting to receive energy from TeraChem...'
         call flush (OUTPUT_UNIT)
      end if
      call MPI_Recv(escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                  & MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 1, MPI_DOUBLE_PRECISION)
      ! Checking for TAG=1, which means that SCF did not converge
      if (status(MPI_TAG) == 1) then
         write (stderr, *) 'Got TAG 1 from TeraChem: SCF probably did not converge.'
         abort = .true.
!$OMP FLUSH(abort)
         return
      end if
      if (idebug > 1) then
         print '(A,ES15.6)', 'Received SCF energy from server:', escf
         call flush (OUTPUT_UNIT)
      end if

      if (idebug > 2) then
         print '(a)', 'Waiting to receive charges...'
      end if
      call MPI_Recv(qmcharges(:), natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                   & MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, natqm, MPI_DOUBLE_PRECISION)
      if (modulo(it, nwrite) == 0 .and. nteraservers == 1) then
         call print_charges(qmcharges, iw)
      end if

      if (idebug > 2) then
         print '(a)', 'Waiting to receive dipole moment...'
      end if
      ! QM dipole moment
      call MPI_Recv(dipmom(:, 1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                  & MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 4, MPI_DOUBLE_PRECISION)
      if (idebug > 1) then
         print '(a,4es15.6)', 'Received QM  dipole moment from server:', dipmom(:, 1)
         call flush (OUTPUT_UNIT)
      end if
      ! TODO: Attach dipoles to global electronic structure type
      ! and print them elsewhere. Right now when we run concurrent
      ! TC servers, the printing is not deterministic.
      if (modulo(it, nwrite) == 0 .and. nteraservers == 1) then
         call print_dipoles(dipmom(:, 1), iw, 1)
      end if

      ! MM dipole moment, disabled for now
!   call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr )
!   call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr )

      ! QM gradients
      if (idebug > 1) then
         write (*, '(A)') 'Waiting to receive gradients...'
      end if
      call MPI_Recv(dxyz_all, 3 * natqm, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 3 * natqm, MPI_DOUBLE_PRECISION)
      if (idebug > 1) then
         print '(A)', 'Received the following gradients from server:'
         do iat = 1, natqm
            print '(A)', 'Atom ', iat, ': ', dxyz_all(:, iat)
         end do
         call flush (OUTPUT_UNIT)
      end if

      do iat = 1, natqm
         fx(iat, iw) = -dxyz_all(1, iat)
         fy(iat, iw) = -dxyz_all(2, iat)
         fz(iat, iw) = -dxyz_all(3, iat)
      end do

!$OMP ATOMIC
      eclas = eclas + escf / walkmax

   end subroutine receive_tera

! USE_MPI
#endif

end module mod_force_tera
