module mod_force_mace
! ----------------------------------------------------------------
! Interface for MACE (Machine Learning Atomic Cluster Expansion)
! Perform MPI communications with the MACE Python server.
!
! Modeled after mod_force_tera for TeraChem.
! ----------------------------------------------------------------
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stderr
   use mod_mace_mpi
#ifdef USE_MPI
   use mpi
   use mod_mpi, only: handle_mpi_error, check_recv_count
#endif
   implicit none
   private
   public :: force_mace
   save

contains

   subroutine force_mace(x, y, z, fx, fy, fz, eclas, walkmax)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw
      integer :: mace_comm
      logical :: abort

      abort = .false.

!$OMP PARALLEL DO PRIVATE(iw, mace_comm)
      do iw = 1, walkmax

!$OMP FLUSH(abort)
         if (.not. abort) then

#ifdef USE_MPI
            mace_comm = get_mace_communicator()

            call send_mace(x, y, z, iw, mace_comm)

            call receive_mace(fx, fy, fz, eclas, iw, walkmax, mace_comm, abort)
            if (abort) cycle
#endif

         end if

      end do
!$OMP END PARALLEL DO

      if (abort) then
         call fatal_error(__FILE__, __LINE__, &
         & 'MACE external forces error')
      end if

   end subroutine force_mace

#ifdef USE_MPI

   subroutine send_mace(x, y, z, iw, comm)
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw, comm

      ! Send number of atoms each step (for protocol consistency)
      call send_mace_natom(natqm, comm)

      ! Send coordinates in Bohr
      call send_mace_coordinates(x, y, z, natqm, iw, comm)
   end subroutine send_mace

   subroutine receive_mace(fx, fy, fz, eclas, iw, walkmax, comm, abort)
      use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
      use mod_general, only: idebug
      use mod_qmmm, only: natqm
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: iw, walkmax
      integer, intent(in) :: comm
      logical, intent(inout) :: abort
      real(DP) :: energy
      real(DP) :: forces(3, size(fx, 1))
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr, iat

      call wait_for_mace(comm)

      ! Receive energy (in Hartree, already converted by server)
      if (idebug > 2) then
         print '(a)', 'MACE: Waiting to receive energy...'
         call flush (OUTPUT_UNIT)
      end if
      call MPI_Recv(energy, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
      & MPI_ANY_TAG, comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 1, MPI_DOUBLE_PRECISION)

      ! Check for error tag
      if (status(MPI_TAG) == 1) then
         write (stderr, *) 'Got error tag from MACE server. Evaluation failed.'
         abort = .true.
!$OMP FLUSH(abort)
         return
      end if

      if (idebug > 1) then
         print '(A,ES15.6)', 'MACE: Received energy [Hartree]:', energy
         call flush (OUTPUT_UNIT)
      end if

      ! Receive forces (in Hartree/Bohr, already converted by server)
      if (idebug > 2) then
         print '(a)', 'MACE: Waiting to receive forces...'
      end if
      call MPI_Recv(forces, 3 * natqm, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 3 * natqm, MPI_DOUBLE_PRECISION)

      if (idebug > 1) then
         print '(A)', 'MACE: Received forces [Hartree/Bohr]:'
         do iat = 1, natqm
            print*,'Atom ', iat, ': ', forces(:, iat)
         end do
         call flush (OUTPUT_UNIT)
      end if

      ! Forces are received as forces (not gradients), so no sign flip needed
      do iat = 1, natqm
         fx(iat, iw) = forces(1, iat)
         fy(iat, iw) = forces(2, iat)
         fz(iat, iw) = forces(3, iat)
      end do

!$OMP ATOMIC
      eclas = eclas + energy / walkmax

   end subroutine receive_mace

! USE_MPI
#endif

end module mod_force_mace
