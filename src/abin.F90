! -----------------------------------------------------------------
!  ABIN: Multipurpose ab initio MD program.
!  The potential is calculated on-the-fly by an external program.
!------------------------------------------------------------------
!  Copyright (C) 2014    D.Hollas, J.Suchan, M.Oncak, O.Svoboda and P.Slavicek
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
   use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
   use mod_const, only: DP, AUtoFS
   use mod_arrays
   use mod_files, only: stdout
   use mod_general, only: sim_time, pot, pot_ref, iremd, ipimd, &
      & nwrite, nstep, ncalc, it, inormalmodes, istage, irest
   use mod_init, only: init
   use mod_sh, only: surfacehop, sh_init, get_nacm, move_vars
   use mod_lz, only: lz_hop, en_array_lz, lz_rewind
   use mod_kinetic, only: temperature
   use mod_utils, only: del_file, archive_file
   use mod_transform, only: initialize_pi_transforms, &
                           & qtox, utox, fqtofx
   use mod_mdstep, only: mdstep
   use mod_minimize, only: minimize
   use mod_analysis, only: analysis, restout
   use mod_interfaces
   use mod_en_restraint, only: restrain_pot
   use mod_terampi_sh, only: move_new2old_terash
   use mod_mpi, only: get_mpi_rank, mpi_barrier_wrapper
   use mod_remd, only: remd_swap, nswap
   implicit none
   real(DP) :: dt = 20.0D0, eclas = 0.0D0, equant = 0.0D0
   logical :: file_exists
   integer, dimension(8) :: time_start
   integer :: my_rank

   call date_and_time(values=time_start)

   ! INPUT AND INITIALIZATION SECTION
   call init(dt)

   my_rank = get_mpi_rank()

   ! This cannot be in init because of the namelist 'system'
   if (my_rank == 0) then
      call clean_temp_files()
   end if

   if (irest == 1 .and. (my_rank == 0 .or. iremd == 1)) then
      call archive_file('restart.xyz', it)
   end if

   ! Surface hopping initialization
   if (ipimd == 2) then
      call sh_init(x, y, z, vx, vy, vz)
   else if (ipimd == 5 .and. pot == '_tera_') then
      call sh_init(x, y, z, vx, vy, vz)
      call lz_rewind(en_array_lz)
   end if

   write (stdout, '(A)') 'Job started at: '//trim(get_formatted_date_and_time(time_start))
   write (stdout, *) ''

   ! Transform coordinates and velocities for Path Integral MD
   ! (staging or normal modes)
   if (istage == 1 .or. inormalmodes > 0) then
      call initialize_pi_transforms(x, y, z, vx, vy, vz)
   end if

   ! Note that 'amt' equals physical atomic masses in non-PI simulations
   px = amt * vx
   py = amt * vy
   pz = amt * vz

   if (ipimd == 3) then
      call minimize(x, y, z, fxc, fyc, fzc, eclas)
      call print_footer(time_start)
      call finish(0)
      stop 0
   end if

   write (stdout, *)
   write (stdout, '("#",10X,A,11X,A)') 'Step', 'Time [fs]'

   ! ---------------- PROPAGATION-------------

   ! Without this Barrier, ranks > 0 do not write geom.dat in force_clas
   ! I don't know why the hell not.
   call mpi_barrier_wrapper()

   ! Get initial forces and energies
   call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
   if (ipimd == 1) then
      call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)
   end if

   ! Correct energy history for LZ
   if (ipimd == 5) then
      call lz_rewind(en_array_lz)
   end if

   ! Get reference forces and energies for ab initio  MTS RESPA
   if (pot_ref /= '_none_') then
      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot_ref)
      call force_clas(fxc_diff, fyc_diff, fzc_diff, x, y, z, eclas, pot)
   end if

   ! Set initial values for surface hopping
   if (ipimd == 2) then
      if (irest /= 1) then
         call get_nacm(pot)
      end if
      call move_vars(vx, vy, vz, vx_old, vy_old, vz_old)
      if (pot == '_tera_' .or. restrain_pot == '_tera_') then
         call move_new2old_terash()
      end if
   else if (ipimd == 5 .and. pot == '_tera_') then
      call move_new2old_terash()
   end if

   ! LOOP OVER TIME STEPS
   ! "it" variable is set to 0 or read from restart.xyz in subroutine init
   it = it + 1
   do it = (it), nstep

      ! This barrier is needed so that all MPI processes see the file 'EXIT'
      ! if it is present (we delete it below before we stop the program).
      call mpi_barrier_wrapper()

      inquire (file="EXIT", exist=file_exists)
      if (file_exists) then
         write (stdout, *) 'Found file EXIT. Writing restart file and exiting.'
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

         call mpi_barrier_wrapper()

         if (my_rank == 0) then
            call del_file('EXIT')
         end if

         exit ! break from time loop

      end if

      ! PROPAGATE through one time step
      call mdstep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)

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

         if (pot == '_tera_') then
            call move_new2old_terash()
         end if
      end if

      ! SWAP REMD REPLICAS
      if (iremd == 1 .and. modulo(it, nswap) == 0) then
         call remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
      end if

      ! --- Trajectory analysis ---
      ! In order to analyze the output, we have to perform the back transformation
      ! Transformed (cartesian) coordinates are stored in trans matrices.

      ! Enter this section only every ncalc step
      if (modulo(it, ncalc) /= 0) then
         cycle
      end if

      call temperature(px, py, pz, amt, eclas)

      if (istage == 1) then

         call QtoX(vx, vy, vz, transxv, transyv, transzv)
         call QtoX(x, y, z, transx, transy, transz)
         call FQtoFX(fxc, fyc, fzc, transfxc, transfyc, transfzc)
         call analysis(transx, transy, transz, transxv, transyv, transzv,  &
              &       transfxc, transfyc, transfzc, eclas)

      else if (inormalmodes > 0) then

         call UtoX(x, y, z, transx, transy, transz)
         call UtoX(vx, vy, vz, transxv, transyv, transzv)
         call UtoX(fxc, fyc, fzc, transfxc, transfyc, transfzc)
         call analysis(transx, transy, transz, transxv, transyv, transzv,  &
              &        transfxc, transfyc, transfzc, eclas)
      else

         call analysis(x, y, z, vx, vy, vz, fxc, fyc, fzc, eclas)

      end if

      if (modulo(it, nwrite) == 0) then
         write (stdout, '(I15,F15.2)') it, sim_time * AUtoFS
         call flush (OUTPUT_UNIT)
      end if

      ! Time step loop
   end do

   ! Write restart file at the end of a run
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

   call print_footer(time_start)
   call finish(0)

contains

   subroutine clean_temp_files()
      ! TODO: Implement "clean" bash function in abin interfaces
      ! that should be called here (and only if irest=0)
      call execute_command_line('rm -f ERROR engrad*.dat.* nacm.dat hessian.dat.* geom.dat.*')
   end subroutine clean_temp_files

   function get_formatted_date_and_time(time_data) result(formatted_string)
      character(len=25) :: formatted_string
      integer, dimension(8), intent(in) :: time_data
      formatted_string = ''
      ! time_data must be get from date_and_time() intrinsic
      ! e.g. 1:48:39   3.11.2020
      write (formatted_string, "(I2.2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)") time_data(5), ':', &
         time_data(6), ':', time_data(7), '  ', time_data(3), '.', time_data(2), '.', &
         time_data(1)
   end function get_formatted_date_and_time

   subroutine print_footer(time_start)
      integer, dimension(8), intent(in) :: time_start
      integer, dimension(8) :: time_end

      call date_and_time(values=time_end)
      write (stdout, *) ''
      write (stdout, '(A)') 'Job finished successfully!'
      write (stdout, '(A)') 'Job started at:  '//trim(get_formatted_date_and_time(time_start))
      write (stdout, '(A)') 'Job finished at: '//trim(get_formatted_date_and_time(time_end))
   end subroutine print_footer

end program abin
