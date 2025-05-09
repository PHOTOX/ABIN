! Temperature Replica-Exchange MD (also known as parallel tempering)
module mod_remd
#ifdef USE_MPI
   use mpi
   use mod_mpi, only: handle_mpi_error, get_mpi_size, get_mpi_rank
#endif
   use mod_const, only: DP
   use mod_files, only: stdout, UREMD
   use mod_error, only: fatal_error
   implicit none
   private

   integer, parameter :: MAX_REPLICA = 50

   ! Public subroutines
   public :: remd_init, remd_swap
   ! Public input parameters, read from namelist 'remd'
   public :: nreplica, nswap, deltaT, Tmax, temp_list

#ifdef USE_MPI
   ! Internal variables
   integer :: reps(MAX_REPLICA) = -1
   real(DP) :: ratios_cumul(MAX_REPLICA) = 0.0D0
#endif

   ! Number of replicas
   integer :: nreplica = 0
   ! Attempt to swap replicas every nwap step.
   integer :: nswap = -1
   ! Determine temperatures by specifying deltaT
   ! between replicas.
   real(DP) :: deltaT = -1
   ! Specify maximum temperature (and minimum temperature by temp),
   ! the other temperatures are determined by distributions in remd_init().
   real(DP) :: Tmax = -1
   ! Determine all temperatures manually.
   real(DP) :: temp_list(MAX_REPLICA) = -1
   save

contains

#ifdef USE_MPI
   subroutine remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
      use mod_general, only: it
      use mod_random, only: vranf
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: eclas
      real(DP) :: eclas_new, prob, probs(MAX_REPLICA), ran(MAX_REPLICA), rn
      integer, parameter :: tag_en = 10, tag_swap = 1
      integer :: source, ierr
      integer :: status(MPI_STATUS_SIZE), i
      integer :: my_rank
      logical :: lswap, lswaps(MAX_REPLICA)
      character(len=100) :: formt

      lswap = .false.

      ! Broadcast array of random numbers from rank 0
      my_rank = get_mpi_rank()
      if (my_rank == 0) call vranf(ran, nreplica - 1)

      ! TODO: This is probably no longer needed, each mpi rank has now
      ! "independent" prngs...
      call MPI_Bcast(ran, MAX_REPLICA, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call handle_mpi_error(ierr)
      if (my_rank /= 0) rn = ran(my_rank)

      ! Main Swapping loop, starting with the second replica and going up
      ! In GROMACS for some reason that I do not understand, they attempt to swap
      ! only even or odd replicas pres swapping event.
      ! Right now, we don't do that.
      do source = 1, nreplica - 1

         if (my_rank == source - 1) then

            ! TODO: Why MPI_REAL8 ??? eclas_new is DP as well!
            call MPI_Recv(eclas_new, 1, MPI_REAL8, source, tag_en, MPI_COMM_WORLD, status, ierr)
            call handle_mpi_error(ierr)
            call MPI_Send(eclas, 1, MPI_DOUBLE_PRECISION, source, tag_en, MPI_COMM_WORLD, ierr)
            call handle_mpi_error(ierr)
            call MPI_Recv(lswap, 1, MPI_LOGICAL, source, tag_swap, MPI_COMM_WORLD, status, ierr)
            call handle_mpi_error(ierr)

         else if (my_rank == source) then
            ! Replica with the higher temperature is the master here
            call MPI_Send(eclas, 1, MPI_REAL8, source - 1, tag_en, MPI_COMM_WORLD, ierr)
            call handle_mpi_error(ierr)
            call MPI_Recv(eclas_new, 1, MPI_DOUBLE_PRECISION, source - 1, tag_en, MPI_COMM_WORLD, status, ierr)
            call handle_mpi_error(ierr)
            prob = (eclas - eclas_new) * (1 / temp_list(my_rank + 1) - 1 / temp_list(my_rank))
            if (prob > 0.0D0) then
               prob = 1.0D0
            else
               prob = exp(prob)
            end if
            if (rn < prob) then
               lswap = .true.
            else
               lswap = .false.
            end if
            call MPI_Send(lswap, 1, MPI_LOGICAL, source - 1, tag_swap, MPI_COMM_WORLD, ierr)
            call handle_mpi_error(ierr)

         end if

         if (my_rank == source .or. my_rank == source - 1) then
            if (lswap) then
               eclas = eclas_new
               call swap_replicas(x, y, z, px, py, pz, fxc, fyc, fzc, source - 1, source)
            end if
         end if

      end do

      i = reps(my_rank + 1)
      ! TODO: we dont really need allgather
      ! We should really track replicas a posteriori via some script, all output is there
      call MPI_AllGather(i, 1, MPI_INTEGER, reps, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call handle_mpi_error(ierr)
      call MPI_Gather(prob, 1, MPI_DOUBLE_PRECISION, probs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call handle_mpi_error(ierr)
      call MPI_Gather(lswap, 1, MPI_LOGICAL, lswaps, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call handle_mpi_error(ierr)
      if (my_rank == 0) then
         do i = 1, nreplica - 1
            if (lswaps(i)) ratios_cumul(i) = ratios_cumul(i) + 1.0D0
         end do
         write (formt, '("(A,",I0,"F7.4)")') nreplica - 1
         write (UREMD, '(A)') '===================================='
         write (UREMD, formt) 'Swaping probabilities: ', (probs(i), i=2, nreplica)
         write (UREMD, *) 'Exchanges: ', (lswaps(i), i=1, nreplica - 1)
         write (UREMD, formt) 'Acceptance ratios', (ratios_cumul(i) * nswap / it, i=1, nreplica - 1)
         write (formt, '("(A,",I0,"I4)")') nreplica
         write (UREMD, formt) 'ReplicaTravel:', (reps(i), i=1, nreplica)
      end if

   end subroutine remd_swap

   subroutine swap_replicas(x, y, z, px, py, pz, fxc, fyc, fzc, rank1, rank2)
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      integer, intent(in) :: rank1, rank2
      real(DP), allocatable :: x_new(:, :), y_new(:, :), z_new(:, :), buff(:, :)
      real(DP), allocatable :: px_new(:, :), py_new(:, :), pz_new(:, :)
      real(DP), allocatable :: fxc_new(:, :), fyc_new(:, :), fzc_new(:, :)
      real(DP) :: scal
      integer :: status(MPI_STATUS_SIZE)
      integer :: my_rank

      integer, parameter :: tag_x = 11, tag_y = 12, tag_z = 13
      integer, parameter :: tag_px = 114, tag_py = 115, tag_pz = 116
      integer, parameter :: tag_fx = 17, tag_fy = 18, tag_fz = 19
      integer :: ierr, dest, size1, size2, irank

      my_rank = get_mpi_rank()

      ! First, rank1 sends data and rank2 receives
      size1 = size(x(:, 1))
      size2 = size(x(1, :))
      ! TODO: we need in principle ony 1 buffer
      allocate (buff(size1, size2))
      allocate (x_new(size1, size2))
      allocate (y_new(size1, size2))
      allocate (z_new(size1, size2))
      allocate (px_new(size1, size2))
      allocate (py_new(size1, size2))
      allocate (pz_new(size1, size2))
      allocate (fxc_new(size1, size2))
      allocate (fyc_new(size1, size2))
      allocate (fzc_new(size1, size2))

      ! Momenta are scaled according to the new temperature
      ! p_new = p_old * dsqrt(T_new/T_old)
      if (my_rank == rank1) then
         scal = dsqrt(temp_list(my_rank + 1) / temp_list(my_rank + 2))
      else if (my_rank == rank2) then
         scal = dsqrt(temp_list(my_rank + 1) / temp_list(my_rank))
      end if

      if (my_rank == rank1) then
         dest = rank2
         call MPI_Send(x, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, ierr)
         call MPI_Send(y, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, ierr)
         call MPI_Send(z, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, ierr)
         call MPI_Send(px, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, ierr)
         call MPI_Send(py, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, ierr)
         call MPI_Send(pz, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, ierr)
         call MPI_Send(fxc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, ierr)
         call MPI_Send(fyc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, ierr)
         call MPI_Send(fzc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, ierr)
         call MPI_Send(reps(my_rank + 1), 1, MPI_INTEGER, dest, tag_fz, MPI_COMM_WORLD, ierr)
         call MPI_Recv(x_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(y_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(z_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(px_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(py_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(pz_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fxc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fyc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fzc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(irank, 1, MPI_INTEGER, dest, tag_fz, MPI_COMM_WORLD, status, ierr)
         reps(my_rank + 1) = irank
      else if (my_rank == rank2) then
         dest = rank1
         call MPI_Recv(x_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(y_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(z_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(px_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(py_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(pz_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fxc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fyc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(fzc_new, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr)
         call MPI_recv(irank, 1, MPI_INTEGER, dest, tag_fz, MPI_COMM_WORLD, status, ierr)
         call MPI_Send(x, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, ierr)
         call MPI_Send(y, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, ierr)
         call MPI_Send(z, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, ierr)
         call MPI_Send(px, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, ierr)
         call MPI_Send(py, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, ierr)
         call MPI_Send(pz, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, ierr)
         call MPI_Send(fxc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, ierr)
         call MPI_Send(fyc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, ierr)
         call MPI_Send(fzc, size1 * size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, ierr)
         call MPI_Send(reps(my_rank + 1), 1, MPI_INTEGER, dest, tag_fz, MPI_COMM_WORLD, ierr)
         reps(my_rank + 1) = irank
      end if

!   Now all ranks have required data, so they can swap them
      x = x_new; y = y_new; z = z_new

      px = px_new * scal; py = py_new * scal; pz = pz_new * scal

      fxc = fxc_new; fyc = fyc_new; fzc = fzc_new

      deallocate (buff)
      deallocate (x_new); deallocate (y_new); deallocate (z_new)
      deallocate (px_new); deallocate (py_new); deallocate (pz_new)
      deallocate (fxc_new); deallocate (fyc_new); deallocate (fzc_new)

      if (my_rank == rank1) write (*, *) 'Succesfull swap between replicas ', rank1, rank2

   end subroutine swap_replicas

   subroutine check_remd_params()
      use mod_utils, only: int_positive
      use mod_general, only: pot, ipimd
      character(len=300) :: errormsg

      call int_positive(nswap, 'nswap')

      if (nreplica > MAX_REPLICA) then
         write (errormsg, '(A,I0,A)') 'Maximum number of replicas ', MAX_REPLICA, ' exceeded.'
         call fatal_error(__FILE__, __LINE__, errormsg)
      end if

      if (deltaT < 0 .and. Tmax < 0 .and. temp_list(1) < 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'You did not specify deltaT, Tmax, or temp_list for REMD.')
      else if (deltaT > 0 .and. Tmax > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'You can not specify both deltaT and  Tmax for REMD!')
      else if (deltaT > 0 .and. temp_list(1) > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'You can not specify both deltaT and  temp_list() for REMD!')
      else if (Tmax > 0 .and. temp_list(1) > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'You can not specify both Tmax and temp_list() for REMD!')
      end if

      if (pot == '_cp2k_') then
         call fatal_error(__FILE__, __LINE__, &
            & "REMD with internal CP2K interface not supported")
      end if

      if (ipimd > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'REMD is currently implemented only for classical dynamics (mdtype=="md")')
      end if

      if (nreplica < 2) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Number of REMD replicas nreplica must be at least 2')
      end if

      if (get_mpi_size() /= nreplica) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Number of MPI processes does not match number of REMD replicas.')
      end if
   end subroutine check_remd_params

   subroutine remd_init(temp, temp0)
      use mod_const, only: AUTOK
      real(DP), intent(inout) :: temp, temp0
      character(len=300) :: errormsg
      real(DP) :: koef
      integer :: ierr, i
      integer :: my_rank

      my_rank = get_mpi_rank()

      ! Sanity checks
      call check_remd_params()

      write (stdout, '(A,I0)') 'Number of REMD replicas: ', nreplica

      ! Set the temperatures
      if (deltaT > 0) then
         temp = temp + my_rank * deltaT
         if (temp0 > 0) temp0 = temp0 + my_rank * deltaT
      else if (Tmax > 0) then
         koef = (Tmax / temp)**(1.D0 / (nreplica - 1))
         temp = temp * koef**my_rank
         if (temp0 > 0) temp0 = temp0 * koef**my_rank
      end if
      if (temp_list(1) > 0) then
         temp = temp_list(my_rank + 1)
         if (temp0 > 0) temp0 = temp_list(my_rank + 1)
      else
         temp_list(my_rank + 1) = temp
      end if

      if (temp_list(my_rank + 1) <= 0) then
         write (errormsg, '(A,I0,A)') 'Temperature of replica ', my_rank, ' was not set or is negative.'
         call fatal_error(__FILE__, __LINE__, errormsg)
      end if

      call MPI_AllGather(temp, 1, MPI_DOUBLE_PRECISION, temp_list, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      call handle_mpi_error(ierr)
      do i = 1, nreplica
         write (stdout, '("Temperature of replica ",I0," is ",F6.2," K")') i - 1, temp_list(i)
      end do
      temp_list = temp_list / AUTOK

      reps(my_rank + 1) = my_rank
   end subroutine remd_init

#else

   subroutine remd_init(temp, temp0)
      use mod_error, only: not_compiled_with
      real(DP), intent(inout) :: temp, temp0
      temp = 0.0D0; temp0 = 0.0D0
      call not_compiled_with('MPI')
   end subroutine remd_init

   subroutine remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
      use mod_error, only: not_compiled_with
      real(DP), dimension(:, :), intent(inout) :: x, y, z, px, py, pz, fxc, fyc, fzc
      real(DP), intent(inout) :: eclas

      ! Assignments to squash compiler warnings
      x = 0.0D0; y = 0.0D0; z = 0.0D0
      px = 0.0D0; py = 0.0D0; pz = 0.0D0
      fxc = 0.0D0; fyc = 0.0D0; fzc = 0.0D0
      eclas = 0.0D0
      call not_compiled_with('MPI')
   end subroutine remd_swap

! USE_MPI
#endif

end module mod_remd
