! Temperature Replica-Exchange MD (also known as parallel tempering)
module mod_remd
   use mod_const, only: DP
   use mod_utils, only: abinerror
#ifdef USE_MPI
   implicit none
   private
   public :: remd_init, remd_swap
   public :: nreplica, nswap, deltaT, Tmax, temp_list
   integer, parameter :: MAX_REPLICA = 50
   integer :: nreplica, nswap = -1, reps(MAX_REPLICA) = -1
   real(DP) :: deltaT = -1, Tmax = -1
   real(DP) :: temp_list(MAX_REPLICA) = -1, ratios_cumul(MAX_REPLICA) = 0.0D0

contains

   !TODO: make general subroutine check_mpi_error in utils.f90 and check every ierr
   subroutine remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
      use mpi
      use mod_general, only: my_rank, it
      use mod_random, only: vranf
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: eclas
      real(DP) :: eclas_new, prob, probs(MAX_REPLICA), ran(MAX_REPLICA), rn
      integer :: source, ierr, tag_en = 10, tag_swap = 1
      integer :: status(MPI_STATUS_SIZE), i
      logical :: lswap = .false., lswaps(MAX_REPLICA)

      ! Broadcast array of random numbers from rank 0
      if (my_rank == 0) call vranf(ran, nreplica - 1)

      ! TODO: This is probably no longer needed, each mpi rank has now
      ! "independent" prngs...
      call MPI_Bcast(ran, MAX_REPLICA, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      if (my_rank /= 0) rn = ran(my_rank)

      ! Main Swapping loop, starting with the second replica and going up
      ! In GROMACS for some reason that I do not understand, they attempt to swap
      ! only even or odd replicas pres swapping event.
      ! Right now, we don't do that.
      do source = 1, nreplica - 1

         if (my_rank == source - 1) then

            call MPI_Recv(eclas_new, 1, MPI_REAL8, source, tag_en, MPI_COMM_WORLD, status, ierr)
            call MPI_Send(eclas, 1, MPI_DOUBLE_PRECISION, source, tag_en, MPI_COMM_WORLD, ierr)
            call MPI_Recv(lswap, 1, MPI_LOGICAL, source, tag_swap, MPI_COMM_WORLD, status, ierr)

         else if (my_rank == source) then
            ! Replica with the higher temperature is the master here
            call MPI_Send(eclas, 1, MPI_REAL8, source - 1, tag_en, MPI_COMM_WORLD, ierr)
            call MPI_Recv(eclas_new, 1, MPI_DOUBLE_PRECISION, source - 1, tag_en, MPI_COMM_WORLD, status, ierr)
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
            !      write(*,*)'Probability, ran. number, lswap', prob, rn, lswap
            call MPI_Send(lswap, 1, MPI_LOGICAL, source - 1, tag_swap, MPI_COMM_WORLD, ierr)

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
      call MPI_Gather(prob, 1, MPI_DOUBLE_PRECISION, probs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(lswap, 1, MPI_LOGICAL, lswaps, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if (my_rank == 0) then
         do i = 1, nreplica - 1
            if (lswaps(i)) ratios_cumul(i) = ratios_cumul(i) + 1.0D0
         end do
         write (20, *) '==============================================='
         write (20, *) 'Swaping probabilities: ', (probs(i), i=2, nreplica)
         write (20, *) 'Exchanges: ', (lswaps(i), i=1, nreplica - 1)
         write (20, *) 'Acceptance ratios', (ratios_cumul(i) * nswap / it, i=1, nreplica - 1)
         write (20, *) 'ReplicaTravel:', (reps(i), i=1, nreplica)
      end if

   end subroutine remd_swap

   subroutine swap_replicas(x, y, z, px, py, pz, fxc, fyc, fzc, rank1, rank2)
      use mpi
      use mod_general, only: my_rank
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      integer, intent(in) :: rank1, rank2
      real(DP), allocatable :: x_new(:, :), y_new(:, :), z_new(:, :), buff(:, :)
      real(DP), allocatable :: px_new(:, :), py_new(:, :), pz_new(:, :)
      real(DP), allocatable :: fxc_new(:, :), fyc_new(:, :), fzc_new(:, :)
      real(DP) :: scal
      integer :: status(MPI_STATUS_SIZE)

      integer :: tag_x = 11, tag_y = 12, tag_z = 13
      integer :: tag_px = 114, tag_py = 115, tag_pz = 116
      integer :: tag_fx = 17, tag_fy = 18, tag_fz = 19
      integer :: ierr, dest, size1, size2, irank

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
      ! p_new = p_old * sqrt(T_new/T_old)
      if (my_rank == rank1) then
         scal = dsqrt(temp_list(my_rank + 1) / temp_list(my_rank + 2))
      else if (my_rank == rank2) then
         scal = dsqrt(temp_list(my_rank + 1) / temp_list(my_rank))
      end if

!   write(*,*)'pz ', pz, my_rank
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

   subroutine remd_init(temp, temp0)
      use mpi
      use mod_const, only: AUTOK
      use mod_general, only: pot, my_rank, ipimd, irest
      real(DP), intent(inout) :: temp, temp0
      real(DP) :: koef
      integer :: ierr, i

!   call flush(6)
      if (irest == 0) then
         open (20, file='remd.out', action='write', access='sequential')
      else
         open (20, file='remd.out', action='write', access='append')
      end if
      ! First, do sanity check
      if (nswap < 0) then
         write (*, *) 'ERROR: nswap must be a positive integer!'
         call abinerror('remd_init')
      end if

      if (deltaT < 0 .and. Tmax < 0 .and. temp_list(1) < 0) then
         write (*, *) 'ERROR: You did not specify deltaT, Tmax, or temp_list for REMD.'
         call abinerror('remd_init')
      else if (deltaT > 0 .and. Tmax > 0) then
         write (*, *) 'ERROR: You can not specify both deltaT and  Tmax for REMD!'
         call abinerror('remd_init')
      else if (deltaT > 0 .and. temp_list(1) > 0) then
         write (*, *) 'ERROR: You can not specify both deltaT and  temp_list() for REMD!'
         call abinerror('remd_init')
      else if (Tmax > 0 .and. temp_list(1) > 0) then
         write (*, *) 'ERROR: You can not specify both Tmax and temp_list() for REMD!'
         call abinerror('remd_init')
      end if

      if (pot == '_cp2k_') then
         write (*, *) "ERROR: REMD with internal CP2K interface not supported!"
         call abinerror('remd_init')
      end if

      if (ipimd > 0) then
         write (*, *) 'REMD is currently implemented only for classical dynamics (ipimd.eq.0)!'
         call abinerror('remd_init')
      end if

      if (nreplica <= 1) then
         write (*, *) 'Incorrect number of REMD replicas!'
         call abinerror('remd_init')
      end if

      ! determine number of MPI processess and check with user input
      call MPI_Comm_size(MPI_COMM_WORLD, i, ierr)
      if (i /= nreplica) then
         write (*, *) 'ERROR: Number of MPI processes does not match number of replicas.'
         write (*, *) 'Number of processors used must equal "nreplica" parameter in input file.'
         call abinerror('remd_init')
      end if

      if (my_rank == 0) write (*, *) 'Number of REMD replicas: ', nreplica

      ! Set the temperatures
      if (deltaT > 0) then
         temp = temp + my_rank * deltaT
         if (temp0 > 0) temp0 = temp0 + my_rank * deltaT
      else if (Tmax > 0) then
         koef = (Tmax / temp)**(1.D0 / (nreplica - 1))
         write (*, *) 'koeff', koef
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
         write (*, '(A,I2,A)') 'ERROR: Temperature of replica ', my_rank, ' was not set or is negative.'
      end if

      ! not sure about this
      deltaT = deltaT / AUTOK

      call MPI_AllGather(temp, 1, MPI_DOUBLE_PRECISION, temp_list, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      if (my_rank == 0) then
         do i = 1, nreplica
            write (20, '(A,I2,A,F8.2)') 'Temperature of replica ', i - 1, ' is ', temp_list(i)
         end do
      end if
      temp_list = temp_list / AUTOK

      reps(my_rank + 1) = my_rank

   end subroutine remd_init

! USE_MPI
#endif

end module mod_remd
