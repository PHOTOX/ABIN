! -----------------------------------------------------------------
!  ABIN: Multipurpose ab initio MD program.
!  The potential is calculated on-the-fly by an external program.
!------------------------------------------------------------------
!  Copyright (C) 2014             D.Hollas, M.Oncak, O.Svoboda and P.Slavicek
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program in the file LICENSE. If not, see <http://www.gnu.org/licenses/>.
program abin
   use mod_const, only: DP, AUtoFS
   use mod_arrays
   use mod_general
   use mod_sh, only: surfacehop, ntraj, sh_init, get_nacm, move_vars
   use mod_lz, only: lz_hop, en_array_lz, lz_rewind
   use mod_kinetic, only: temperature
   use mod_utils, only: abinerror, archive_file, get_formatted_date_and_time
   use mod_transform
   use mod_mdstep
   use mod_minimize, only: minimize
   use mod_analysis, only: analysis, restout
   use mod_interfaces
   use mod_en_restraint
   use mod_plumed
   use mod_terampi_sh, only: move_new2old_terash
   use mod_remd
#ifdef USE_MPI
   use mpi, only: MPI_COMM_WORLD, MPI_Barrier
#endif
   implicit none
   ! TODO: These should probably be defined and stored in some module, not here
   real(DP) :: dt = 20.0D0, eclas = 0.0D0, equant = 0.0D0
   integer :: itrj
   logical :: file_exists
   integer, dimension(8) :: time_start, time_end
   real(DP) :: total_cpu_time
   integer :: ierr
!$ integer :: nthreads, omp_get_max_threads

   call date_and_time(VALUES=time_start)

   ! INPUT AND INITIALIZATION SECTION
   call init(dt)

   ! This cannot be in init because of the namelist 'system'
   if (my_rank == 0) call clean_temp_files()

   if (irest == 1 .and. (my_rank == 0 .or. iremd == 1)) then
      call archive_file('restart.xyz', it)
   end if

   ! Surface hopping initialization
   if (ipimd == 2 .or. ipimd == 4) then
      call sh_init(x, y, z, vx, vy, vz)
   else if (ipimd == 5 .and. pot == '_tera_') then
      call sh_init(x, y, z, vx, vy, vz)
      call lz_rewind(en_array_lz)
   end if

!$ nthreads = omp_get_max_threads()
   if (my_rank == 0) then
!$    write (*, *) 'Number of OpenMP threads used = ', nthreads
      write (*, '(A)') 'Job started at: '//trim(get_formatted_date_and_time(time_start))
      write (*, *) ''
   end if

   ! Transform coordinates and velocities Path Integral MD
   ! (staging or normal modes)
   if (istage == 1 .or. inormalmodes > 0) then
      call initialize_pi_transforms(x, y, z, vx, vy, vz)
   end if

   ! Note that 'amt' equals 'am' for non-PI simulations
   px = amt * vx
   py = amt * vy
   pz = amt * vz

   if (ipimd == 3) then

      call minimize(x, y, z, fxc, fyc, fzc, eclas)

   else

      if (my_rank == 0) then
         write (*, *)
         write (*, *) '#      Step     Time [fs]'
      end if
!---------------- PROPAGATION-----------------------------------

#ifdef USE_MPI
      ! Without this Barrier, ranks > 0 do not write geom.dat in force_clas
      ! I don't know why the hell not.
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      ! getting initial forces and energies
      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
      if (ipimd == 1) then
         call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)
      end if

      ! Correct energy history for LZ
      if (ipimd == 5) call lz_rewind(en_array_lz)

      ! if we use reference potential with RESPA
      if (pot_ref /= 'none') then
         call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot_ref)
         call force_clas(fxc_diff, fyc_diff, fzc_diff, x, y, z, eclas, pot)
      end if

      ! setting initial values for SURFACE HOPPING
      if (ipimd == 2 .or. ipimd == 4) then
         do itrj = 1, ntraj
            if (irest /= 1) call get_nacm(itrj)
            call move_vars(vx, vy, vz, vx_old, vy_old, vz_old, itrj)
            if (pot == '_tera_' .or. restrain_pot == '_tera_') call move_new2old_terash()
         end do
      else if (ipimd == 5 .and. pot == '_tera_') then
         call move_new2old_terash()
      end if

      ! LOOP OVER TIME STEPS
      ! "it" variable is set to 0 or read from restart.xyz in subroutine init
      it = it + 1
      do it = (it), nstep

         ! TODO: Move this bit to a helper function check_for_exit()
#ifdef USE_MPI
         ! This is needed, because all ranks need to see fle EXIT
         ! Maybe we should get rid of it for performance reasons
         ! We could still call Abort from rank0, but we would not be sure that we are on
         ! the same timestep. Does it matter??
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
         inquire (FILE="EXIT", EXIST=file_exists)
         if (file_exists) then
            write (*, *) 'Found file EXIT. Exiting...'
            if (istage == 1) then
               call QtoX(vx, vy, vz, transxv, transyv, transzv)
               call QtoX(x, y, z, transx, transy, transz)
               call restout(transx, transy, transz, transxv, transyv, transzv, it - 1)
            else if (inormalmodes > 0) then
               call UtoX(x, y, z, transx, transy, transz)
               call UtoX(vx, vy, vz, transxv, transyv, transzv)
               call restout(transx, transy, transz, transxv, transyv, transzv, it - 1)
            else
               call restout(x, y, z, vx, vy, vz, it - 1)
            end if

#ifdef USE_MPI
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
            if (my_rank == 0) then
               call system('rm EXIT')
            end if

            exit ! break from time loop

         end if

         ! CALL the integrator, propagate through one time step
         select case (md)
         case (1)
            call respastep(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, fxc, fyc, fzc, fxq, fyq, fzq)
         case (2)
            call verletstep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
            ! include entire Ehrenfest step, in first step, we start from pure initial state so at first step we dont
            ! need NAMCE and take forces just as a grad E
         case (3)
            call respashake(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, fxc, fyc, fzc, fxq, fyq, fzq)
         case (4)
            call doublerespastep(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, fxc, fyc, fzc, fxq, fyq, fzq)
         end select

         sim_time = sim_time + dt

         vx = px / amt
         vy = py / amt
         vz = pz / amt

         ! SURFACE HOPPING SECTION
         ! SH is called here, Ehrenfest inside the Verlet MD step
         if (ipimd == 2) then

            call surfacehop(x, y, z, vx, vy, vz, vx_old, vy_old, vz_old, dt, eclas)

            px = amt * vx
            py = amt * vy
            pz = amt * vz

            ! TODO: this should be in the surfacehop routine
            if (pot == '_tera_') then
               call move_new2old_terash()
            end if

         end if

         ! LANDAU ZENER HOPPING
         if (ipimd == 5) then
            call lz_hop(x, y, z, vx, vy, vz, fxc, fyc, fzc, amt, dt, eclas, pot)
            px = amt * vx
            py = amt * vy
            pz = amt * vz

            if (pot == '_tera_') call move_new2old_terash()
         end if

#ifdef USE_MPI
         ! SWAP REMD REPLICAS
         if (iremd == 1 .and. modulo(it, nswap) == 0) then
            call remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
         end if
#endif

!--------------------SECTION of trajectory ANALYSIS
! In order to analyze the output, we have to perform the back transformation
! Transformed (cartesian) coordinates are stored in trans matrices.
! DHmod-21.12.2012 enter this section only every ncalc step

         ! maybe we should visit this section only when my_rank.eq.0
         ! this would not be the case for REMD

         if (modulo(it, ncalc) /= 0) then
            cycle
         end if

         ! TODO: Move this call inside analysis() subroutine
         call temperature(px, py, pz, amt, eclas)

         if (istage == 1) then

            call QtoX(vx, vy, vz, transxv, transyv, transzv)
            call QtoX(x, y, z, transx, transy, transz)
            call FQtoFX(fxc, fyc, fzc, transfxc, transfyc, transfzc)
            call analysis(transx, transy, transz, transxv, transyv, transzv,  &
                 &       transfxc, transfyc, transfzc, eclas, equant)

         else if (inormalmodes > 0) then

            call UtoX(x, y, z, transx, transy, transz)
            call UtoX(vx, vy, vz, transxv, transyv, transzv)
            call UtoX(fxc, fyc, fzc, transfxc, transfyc, transfzc)
            call analysis(transx, transy, transz, transxv, transyv, transzv,  &
                 &        transfxc, transfyc, transfzc, eclas, equant)
         else

            call analysis(x, y, z, vx, vy, vz, fxc, fyc, fzc, eclas, equant)

         end if

         if (modulo(it, nwrite) == 0 .and. my_rank == 0) then
            write (*, '(I20,F15.2)') it, sim_time * AUtoFS
            call flush (6)
         end if

         ! Time step loop
      end do

      ! DUMP restart file at the end of a run
      ! Because NCALC might be >1, we have to perform transformation to get the most
      ! recent coordinates and velocities
      it = it - 1

      if (istage == 1) then
         call QtoX(vx, vy, vz, transxv, transyv, transzv)
         call QtoX(x, y, z, transx, transy, transz)
         call restout(transx, transy, transz, transxv, transyv, transzv, it)
      else if (inormalmodes > 0) then
         call UtoX(x, y, z, transx, transy, transz)
         call UtoX(vx, vy, vz, transxv, transyv, transzv)
         call restout(transx, transy, transz, transxv, transyv, transzv, it)
      else
         call restout(x, y, z, vx, vy, vz, it)
      end if

      ! MINIMIZATION endif
   end if

   call finish(0)

   ! FINAL TIMING
   ! TODO: Maybe we should print some statistics more often
   ! i.e. hours per picosecond or sth like that
   if (my_rank == 0) then
      call cpu_time(total_cpu_time)
      write (*, '(A)') 'Total cpu time [s] (does not include ab initio calculations)'
      write (*, *) total_cpu_time
      write (*, '(A)') 'Total cpu time [hours] (does not include ab initio calculations)'
      write (*, *) total_cpu_time / 3600.

      write (*, '(A)') 'Job started at:  '//trim(get_formatted_date_and_time(time_start))
      call date_and_time(VALUES=time_end)
      write (*, '(A)') 'Job finished at: '//trim(get_formatted_date_and_time(time_end))
   end if

contains

   subroutine clean_temp_files()
      ! TODO: Implement "clean" bash function in abin interfaces
      ! that should be called here (and only if irest=0)
      call system('rm -f ERROR engrad*.dat.* nacm.dat hessian.dat.* geom.dat.*')
   end subroutine clean_temp_files

end program abin
