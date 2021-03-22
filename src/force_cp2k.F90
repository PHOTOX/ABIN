module mod_cp2k
   use iso_c_binding
   use mod_const, only: DP
#ifndef USE_CP2K
   use mod_utils, only: not_compiled_with
#endif
   implicit none

   ! We are actually connecting to C interface,
   ! even though CP2K is written in Fortran
   interface
      subroutine CP2K_SET_POSITIONS(env_id, new_pos, sz)
         import :: C_DOUBLE
         integer, value :: env_id, sz
         real(C_DOUBLE), dimension(1:sz), intent(in) :: new_pos
      end subroutine CP2K_SET_POSITIONS

      subroutine CP2K_CALC_ENERGY_FORCE(env_id)
         integer, value :: env_id
      end subroutine CP2K_CALC_ENERGY_FORCE

      subroutine CP2K_DESTROY_FORCE_ENV(env_id)
         integer, value :: env_id
      end subroutine CP2K_DESTROY_FORCE_ENV

      subroutine CP2K_GET_POTENTIAL_ENERGY(env_id, E)
         import :: C_DOUBLE
         integer, value :: env_id
         real(C_DOUBLE), intent(out) :: E
      end subroutine CP2K_GET_POTENTIAL_ENERGY

      subroutine cp2k_create_force_env_comm(env_id, inpath, outpath, mpicomm)
         import :: C_CHAR
         integer, intent(out) :: env_id
         integer, value :: mpicomm
         character(LEN=1, KIND=C_CHAR), intent(IN) :: inpath(*), outpath(*)
      end subroutine cp2k_create_force_env_comm

      subroutine CP2K_GET_FORCES(env_id, f, sz)
         import :: C_DOUBLE
         integer, value :: env_id, sz
         real(C_DOUBLE), dimension(1:sz), intent(out) :: f
      end subroutine CP2K_GET_FORCES

   end interface

   private
   public :: init_cp2k, finalize_cp2k, force_cp2k
   public :: cp2k_mpi_beads
   logical :: cp2k_mpi_beads = .false.
#ifdef USE_CP2K
   integer :: f_env_id
#ifdef USE_MPI
   integer :: cp2k_rank, cp2k_mastercomm
#endif
   real(C_DOUBLE), dimension(:), pointer :: pos, force
#endif
   save

contains

#ifdef USE_CP2K

   subroutine init_cp2k()
      use mod_general, only: natom, nwalk, my_rank, mpi_world_size
      use mod_utils, only: abinerror
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
#ifdef USE_MPI
      use mpi
      integer :: bead, new_size
      integer :: cp2k_mpicomm
      character(len=300) :: chsed
      character(len=4) :: chbead
#endif
      integer :: ierr
      character(len=200, KIND=C_CHAR) :: cp2k_input_file = 'cp2k.inp'
      character(len=200, KIND=C_CHAR) :: cp2k_output_file = 'cp2k.out'

#ifdef USE_MPI
      call cp2k_init()

      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)

      if (mpi_world_size == 1 .or. nwalk == 1) cp2k_mpi_beads = .false.

      if (cp2k_mpi_beads) then
         bead = modulo(my_rank, nwalk) + 1

         if (modulo(mpi_world_size, nwalk) /= 0 .and. mpi_world_size >= nwalk) then
            write (*, *) 'ERROR:Number of MPI processes must be a multiple of nwalk!'
            call abinerror('init_cp2k')
         end if

         if (mpi_world_size < nwalk .and. modulo(nwalk, mpi_world_size) /= 0) then
            write (*, *) 'ERROR:Number of MPI processes not compatible with number of beads!'
            call abinerror('init_cp2k')
         end if

         ! Create new communicators for different beads
         call MPI_Comm_split(MPI_COMM_WORLD, bead, 0, cp2k_mpicomm, ierr)
         if (ierr /= 0) stop "Could not cplit communicators"

         call MPI_Comm_size(cp2k_mpicomm, new_size, ierr)
         call MPI_Comm_rank(cp2k_mpicomm, cp2k_rank, ierr)

         write (chbead, '(A,I3.3)') '.', bead
         chsed = 'sed "s/PROJECT .*/PROJECT RANK'//chbead//'/i" '//cp2k_input_file
         cp2k_output_file = trim(cp2k_output_file)//chbead//C_NULL_CHAR
         cp2k_input_file = trim(cp2k_input_file)//chbead//C_NULL_CHAR

         ! Create separate input file for each bead
         ! Each input file must have a unique project name
         ! Hence, we are using sed
         chsed = trim(chsed)//'>'//trim(cp2k_input_file)

         call system(chsed)
         call system('sync')

         ! Create a new communicator based on cp2k_rank
         ! This is used for gathering forces and energies between different beads
         call MPI_Comm_split(MPI_COMM_WORLD, cp2k_rank, 0, cp2k_mastercomm, ierr)
         if (ierr /= 0) stop "Could not cplit communicators"

         call cp2k_create_force_env_comm(f_env_id, cp2k_input_file, cp2k_output_file, cp2k_mpicomm)

      else

         cp2k_output_file = trim(cp2k_output_file)//C_NULL_CHAR
         cp2k_input_file = trim(cp2k_input_file)//C_NULL_CHAR
         call cp2k_create_force_env(f_env_id, cp2k_input_file, cp2k_output_file)

      end if
! Without MPI
#else

      ! TODO: Use the new c_string function
      cp2k_output_file = trim(cp2k_output_file)//C_NULL_CHAR
      cp2k_input_file = trim(cp2k_input_file)//C_NULL_CHAR

      call cp2k_init_without_mpi()

      call cp2k_create_force_env(f_env_id, cp2k_input_file, cp2k_output_file)

#endif

      allocate (pos(natom * 3), force(natom * 3), stat=ierr)
      if (ierr /= 0) stop "Could not allocate memory"

   end subroutine init_cp2k

   subroutine finalize_cp2k()
      deallocate (pos)
      deallocate (force)
      call cp2k_destroy_force_env(f_env_id)

#ifdef USE_MPI
      call cp2k_finalize()
#else
      call cp2k_finalize_without_mpi()
#endif

   end subroutine finalize_cp2k

   subroutine force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: natom, iqmmm, idebug, nwalk ! , my_rank, mpi_world_size
      use mod_utils, only: abinerror
      use mod_interfaces, only: oniom
      use mod_utils, only: abinerror
#ifdef USE_MPI
      use mpi
      integer :: status(MPI_STATUS_SIZE)
#endif
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: walkmax
      real(DP) :: e0, eclas_mpi
      integer :: iat, iw, ind
      integer :: cp2k_mastersize, cp2k_masterrank
#ifdef USE_MPI
      integer :: ierr, irank
      real(DP) :: fx_mpi(size(fx, 1), size(fx, 2))
      real(DP) :: fy_mpi(size(fx, 1), size(fx, 2))
      real(DP) :: fz_mpi(size(fx, 1), size(fx, 2))
#endif

      eclas = 0.0D0
      eclas_mpi = 0.0D0
      ! Auxiliary variable, initializing here to squash compiler warning
      e0 = 0.0D0

!   bead = modulo(my_rank, walkmax)
      cp2k_masterrank = 0
      cp2k_mastersize = 1
#ifdef USE_MPI
      call MPI_Comm_rank(cp2k_mastercomm, cp2k_masterrank, ierr)
      call MPI_Comm_size(cp2k_mastercomm, cp2k_mastersize, ierr)

      if (walkmax /= nwalk .and. cp2k_mpi_beads) then
         write (*, *) 'This feature is not supported with CP2 MPI interface.'
         call abinerror('force_cp2k')
      end if
#endif

      do iw = 1, walkmax

         if (cp2k_mpi_beads) then
            if (modulo(iw - 1, cp2k_mastersize) /= cp2k_masterrank) cycle
            !if (modulo(iw-1,walkmax).ne.bead.and.cp2k_mpi_beads.and.mpi_world_size.ge.walkmax) cycle
            ! In case number of MPI processes is smaller than nwalk
            !if(modulo(iw-1,mpi_world_size).ne.bead.and.cp2k_mpi_beads.and.mpi_world_size.lt.walkmax) cycle
         end if

         ind = 1
         do iat = 1, natom
            pos(ind) = x(iat, iw)
            pos(ind + 1) = y(iat, iw)
            pos(ind + 2) = z(iat, iw)
            ind = ind + 3
         end do

         if (idebug > 0) write (*, *) 'Setting positions into CP2K engine.', f_env_id
         call cp2k_set_positions(f_env_id, pos, natom * 3)

         if (idebug > 0) write (*, *) 'CP2K engine now calculates forces and energies.'
         call cp2k_calc_energy_force(f_env_id)

         call cp2k_get_potential_energy(f_env_id, e0)

         call cp2k_get_forces(f_env_id, force, natom * 3)

         eclas = eclas + e0 / walkmax

         ind = 1
         do iat = 1, natom
            fx(iat, iw) = force(ind)
            fy(iat, iw) = force(ind + 1)
            fz(iat, iw) = force(ind + 2)
            ind = ind + 3
         end do

         if (iqmmm == 1 .and. cp2k_masterrank == 0) then
            call oniom(x, y, z, fx, fy, fz, eclas, iw)
         end if

      end do

#ifdef USE_MPI

      if (cp2k_mpi_beads) then
         call MPI_Allreduce(eclas, eclas_mpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cp2k_mastercomm, ierr)
         eclas = eclas_mpi

         ! Now we collect forces
         ! First, collect all forces to rank0, than broadcast
         if (cp2k_masterrank == 0) then

            do irank = 1, cp2k_mastersize - 1
               call MPI_Recv(fx_mpi, natom * nwalk, MPI_DOUBLE_PRECISION, irank, MPI_ANY_TAG, cp2k_mastercomm, status, ierr)
               call MPI_Recv(fy_mpi, natom * nwalk, MPI_DOUBLE_PRECISION, irank, MPI_ANY_TAG, cp2k_mastercomm, status, ierr)
               call MPI_Recv(fz_mpi, natom * nwalk, MPI_DOUBLE_PRECISION, irank, MPI_ANY_TAG, cp2k_mastercomm, status, ierr)

               do iw = 1, nwalk
                  if (modulo(iw - 1, cp2k_mastersize) /= irank) cycle
                  fx(:, iw) = fx_mpi(:, iw)
                  fy(:, iw) = fy_mpi(:, iw)
                  fz(:, iw) = fz_mpi(:, iw)
               end do

            end do

         else

            call MPI_Send(fx, natom * nwalk, MPI_DOUBLE_PRECISION, 0, 2, cp2k_mastercomm, ierr)
            call MPI_Send(fy, natom * nwalk, MPI_DOUBLE_PRECISION, 0, 2, cp2k_mastercomm, ierr)
            call MPI_Send(fz, natom * nwalk, MPI_DOUBLE_PRECISION, 0, 2, cp2k_mastercomm, ierr)

         end if

         call MPI_Bcast(fx, nwalk * natom, MPI_DOUBLE_PRECISION, 0, cp2k_mastercomm, ierr)
         call MPI_Bcast(fy, nwalk * natom, MPI_DOUBLE_PRECISION, 0, cp2k_mastercomm, ierr)
         call MPI_Bcast(fz, nwalk * natom, MPI_DOUBLE_PRECISION, 0, cp2k_mastercomm, ierr)

      end if
#endif

   end subroutine force_cp2k

! CP2K endif
#else
   ! Dummy functions for compilation without CP2K
   subroutine init_cp2k()
      use mod_utils, only: abinerror
      write (*, *) 'ERROR: ABIN was not compiled with CP2K interface.'
      call not_compiled_with('internal CP2K interface', 'init_cp2k')
   end subroutine init_cp2k

   subroutine finalize_cp2k()
   end subroutine finalize_cp2k

   subroutine force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_utils, only: abinerror
      real(DP) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP) :: eclas
      integer :: walkmax
      ! Just to squash compiler warnings :-(
      x = fx; y = fy; z = fz
      eclas = 0.0D0
      walkmax = 0

      call not_compiled_with('internal CP2K interface', 'force_cp2k')
   end subroutine force_cp2k

#endif

end module mod_cp2k
