module mod_force_mace
! ----------------------------------------------------------------
! Interface for MACE (Machine Learning Atomic Cluster Expansion)
! Perform MPI communications with the MACE Python server.
!
! Modeled after mod_force_tera for TeraChem.
! ----------------------------------------------------------------
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
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
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw
      integer :: mace_comm
      logical :: abort

      abort = .false.

      do iw = 1, walkmax

#ifdef USE_MPI
         mace_comm = get_mace_communicator()
         ! Send coordinates in Bohr
         call send_mace_coordinates(x, y, z, natqm, iw, mace_comm)

         call receive_mace(fx, fy, fz, eclas, iw, walkmax, mace_comm, abort)
#endif
         if (abort) then
            call fatal_error(__FILE__, __LINE__, 'MACE evaluation failed')
         end if

      end do

   end subroutine force_mace

#ifdef USE_MPI

   subroutine receive_mace(fx, fy, fz, eclas, iw, walkmax, comm, abort)
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

      ! Receive energy (in Hartree, already converted by server)
      if (idebug > 1) then
         write (stdout, '(a)') 'MACE: Waiting to receive energy...'
         call flush (stdout)
      end if
      call MPI_Recv(energy, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 1, MPI_DOUBLE_PRECISION)

      ! Check for error tag
      if (status(MPI_TAG) == MACE_TAG_ERROR) then
         write (stderr, *) 'Got error tag from MACE server. Evaluation failed.'
         abort = .true.
         return
      end if

      if (idebug > 1) then
         print '(A,ES15.6)', 'MACE: Received energy [Hartree]:', energy
         call flush (stdout)
      end if

      ! Receive forces (in Hartree/Bohr, already converted by server)
      if (idebug > 1) then
         print '(a)', 'MACE: Waiting to receive forces...'
      end if
      call MPI_Recv(forces, 3 * natqm, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, 3 * natqm, MPI_DOUBLE_PRECISION)

      if (idebug > 1) then
         print '(A)', 'MACE: Received forces [Hartree/Bohr]:'
         do iat = 1, natqm
            write (stdout, *) 'Atom ', iat, ': ', forces(:, iat)
         end do
         call flush (stdout)
      end if

      ! Forces are received as forces (not gradients), so no sign flip needed
      do iat = 1, natqm
         fx(iat, iw) = forces(1, iat)
         fy(iat, iw) = forces(2, iat)
         fz(iat, iw) = forces(3, iat)
      end do

      eclas = eclas + energy / walkmax

   end subroutine receive_mace

! USE_MPI
#endif

end module mod_force_mace
