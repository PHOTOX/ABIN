! This module contains some routines that do analyses and I/O operations.
! It also contains routines performing restart.
! TODO: Separate restart to its own module
module mod_analysis
   use mod_const, only: DP
   use mod_files, only: stderr, stdout, UMOVIE, UFORCE
   use mod_utils, only: append_rank
   implicit none
   private
   public :: trajout, restout, analysis, restin
   ! These are the character string that we check for during restart.
   character(len=*), parameter :: chnose = 'NHC momenta', chQT = 'Quantum Thermostat', &
                                & chLT = 'Langevin Thermostat', &
                                & chSH = 'Coefficients for SH', &
                                & chAVG = 'Cumulative averages of various estimators', &
                                & chcoords = 'Cartesian Coordinates [au]', &
                                & chvel = 'Cartesian Velocities [au]', &
                                  chLZ = 'Energy difference history for LZ'

contains

   ! Contains all analysis stuff
   subroutine analysis(x, y, z, vx, vy, vz, fxc, fyc, fzc, eclas)
      use mod_analyze_ext, only: analyze_ext
      use mod_estimators, only: estimators
      use mod_general, only: it, ipimd, icv, nwrite, nwritef, nwritev, &
                             nrest, nwritex, nstep, anal_ext, idebug
      use mod_analyze_geometry
      use mod_io
      use mod_vinit, only: remove_rotations
      use mod_system, only: am
      implicit none
      !intent inout because of estimators, writing to nwalk+1
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: eclas

      if (modulo(it, nwrite) == 0 .and. idebug > 0) then
         call remove_rotations(x, y, z, vx, vy, vz, am, .false.)
      end if

      if (ipimd == 1 .or. icv == 1) then
         call estimators(x, y, z, fxc, fyc, fzc, eclas)
      end if

      if (modulo(it, nwritex) == 0) then
         call trajout(x, y, z, it)
      end if

      if (nwritev > 0) then
         if (modulo(it, nwritev) == 0) then
            call velout(vx, vy, vz)
         end if
      end if

      if (nwritef > 0) then
         if (modulo(it, nwritef) == 0) then
            call forceout(x, y, z, fxc, fyc, fzc, UFORCE)
         end if
      end if

      if (ndist >= 1 .and. modulo(it, nwrite) == 0) then
         call print_distances(x, y, z)
      end if

      if (nang >= 1 .and. modulo(it, nwrite) == 0) then
         call print_angles(x, y, z)
      end if

      if (ndih >= 1 .and. modulo(it, nwrite) == 0) then
         call print_dihedrals(x, y, z)
      end if

      if ((modulo(it, nrest) == 0) .or. it == nstep) then
         call restout(x, y, z, vx, vy, vz, it)
      end if

      if (anal_ext == 1) then
         ! Custom analysis function, see analyze_ext_template.F90
         call analyze_ext()
      end if

   end subroutine analysis

   ! TODO: move these subroutines to io.F90
   subroutine trajout(x, y, z, time_step)
      use mod_const, only: ANG
      use mod_files, only: UMOVIE
      use mod_general, only: nwalk, natom, sim_time
      use mod_mpi, only: get_mpi_rank
      use mod_system, only: names
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: time_step
      integer :: iat, iw
      character(len=20) :: chout
      logical :: lopened

      inquire (UMOVIE, opened=lopened)
      if (.not. lopened) then
         chout = append_rank('movie.xyz')
         open (UMOVIE, file=chout, access='append', action="write")
      end if

      do iw = 1, nwalk
         write (UMOVIE, *) natom
         ! In the future, we should get rid of the time step?
         write (UMOVIE, '(A10,I20,A15,F15.2)') 'Time step:', time_step, ' Sim. Time [au]', sim_time
         do iat = 1, natom
            ! printing with slightly lower precision for saving space
            write (UMOVIE, '(A2,3E18.8E2)') names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
         end do
      end do
   end subroutine trajout

   subroutine forceout(x, y, z, fx, fy, fz, fUNIT)
      use mod_general, only: nwalk, natom !, it
      use mod_system, only: names
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: fx(:, :), fy(:, :), fz(:, :)
      integer, intent(in) :: funit
      real(DP) :: fx_tot, fy_tot, fz_tot
      real(DP) :: fx_rot, fy_rot, fz_rot
      integer :: iat, iw

      ! Calculate net translational and rotational gradient (should be close to zero)
      fx_tot = 0.0D0; fy_tot = 0.0D0; fz_tot = 0.0D0
      fx_rot = 0.0D0; fy_rot = 0.0D0; fz_rot = 0.0D0

      do iw = 1, nwalk
         do iat = 1, natom
            fx_tot = fx_tot + fx(iat, iw)
            fy_tot = fy_tot + fy(iat, iw)
            fz_tot = fz_tot + fz(iat, iw)
            fx_rot = fx_rot + fy(iat, iw) * z(iat, iw) - fz(iat, iw) * y(iat, iw)
            fy_rot = fy_rot + fz(iat, iw) * x(iat, iw) - fx(iat, iw) * z(iat, iw)
            fz_rot = fz_rot + fx(iat, iw) * y(iat, iw) - fy(iat, iw) * x(iat, iw)
         end do
      end do

      ! TODO: Include the timestep in the output
      ! Either here:
      ! write(funit, *)natom, it
      ! or maybe better here
      ! write(funit,fkom)'time step', it, &
      !                & 'net force:',fx_tot,fy_tot,fz_tot, &
      !                & 'torque force:',fx_rot,fy_rot,fz_rot
      write (funit, *) natom
      write (funit, '(A10,3E13.5E2,A14,3E13.5E2)') &
         & 'net force:', fx_tot, fy_tot, fz_tot, &
         & ' torque force:', fx_rot, fy_rot, fz_rot

      do iw = 1, nwalk
         do iat = 1, natom
            ! Printing in atomic units
            write (funit, '(A2,3E18.10E2)') names(iat), fx(iat, iw), fy(iat, iw), fz(iat, iw)
         end do
      end do
   end subroutine forceout

   subroutine velout(vx, vy, vz)
      use mod_general, only: nwalk, natom, it
      use mod_system, only: names
      use mod_files, only: UVELOC
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer :: iat, iw

      write (UVELOC, '(I0)') natom
      write (UVELOC, '(A,I0)') 'Time step: ', it

      do iw = 1, nwalk
         do iat = 1, natom
            ! Printing in atomic units
            write (UVELOC, '(A2,3E18.10E2)') names(iat), vx(iat, iw), vy(iat, iw), vz(iat, iw)
         end do
      end do
   end subroutine velout

   subroutine restout(x, y, z, vx, vy, vz, time_step)
      use mod_general, only: icv, ihess, nwalk, ipimd, natom, &
                             iremd, pot, narchive, sim_time
      use mod_mpi, only: get_mpi_rank
      use mod_utils, only: archive_file, append_rank
      use mod_nhc, only: inose, nhc_restout
      use mod_estimators
      use mod_kinetic, only: entot_cumul, est_temp_cumul
      use mod_sh_integ, only: sh_write_wf
      use mod_sh, only: write_nacmrest, istate
      use mod_lz, only: lz_restout
      use mod_gle, only: gle_restout, pile_restout
      use mod_random, only: write_prng_state
      use mod_terampi_sh, only: write_wfn
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer, intent(in) :: time_step
      integer :: iw, my_rank, urest
      logical :: file_exists
      character(len=200) :: chout, chsystem

      my_rank = get_mpi_rank()

      if (pot == '_tera_' .and. &
         & (ipimd == 2 .or. ipimd == 5)) then
         call write_wfn()
      end if

      chout = append_rank('restart.xyz')

      inquire (file=chout, exist=file_exists)
      if (file_exists) then
         chsystem = 'cp '//trim(chout)//'  '//trim(chout)//'.old'
         if (iremd == 1) then
            call execute_command_line(chsystem)
         else if (my_rank == 0) then
            call execute_command_line(chsystem)
         end if
      end if

      open (newunit=urest, file=chout, action='write')

      write (urest, '(I0,X,ES24.16E3)') time_step, sim_time

      write (urest, *) chcoords
      call write_xyz(x, y, z, natom, nwalk, urest)

      write (urest, *) chvel
      call write_xyz(vx, vy, vz, natom, nwalk, urest)

      if (ipimd == 2) then
         call write_nacmrest()
         write (urest, *) chSH
         write (urest, *) istate
         call sh_write_wf(urest)
      end if

      if (ipimd == 5) then
         write (urest, *) chLZ
         call lz_restout(urest)
      end if

      if (inose == 1) then
         write (urest, *) chnose
         call nhc_restout(urest)
      end if

      if (inose == 2 .or. inose == 4) then
         write (urest, *) chQT
         call gle_restout(urest)
      end if
      if (inose == 3) then
         write (urest, *) chLT
         call pile_restout(urest)
      end if

      write (urest, *) chAVG
      write (urest, '(ES24.16E3)') est_temp_cumul
      write (urest, '(2ES25.16E3)') est_prim_cumul, est_vir_cumul
      write (urest, '(ES24.16E3)') entot_cumul

      if (icv == 1) then
         write (urest, '(3ES25.16E3)') est_prim2_cumul, est_prim_vir, est_vir2_cumul
         write (urest, '(2ES25.16E3)') cv_prim_cumul, cv_vir_cumul
         if (ihess == 1) then
            write (urest, '(ES24.16E3)') cv_dcv_cumul
            do iw = 1, nwalk
               write (urest, '(ES24.16E3)') cvhess_cumul(iw)
            end do
         end if
      end if

      ! write current state of PRNG
      call write_prng_state(urest)

      close (urest)

      if (modulo(time_step, narchive) == 0) then
         call archive_file('restart.xyz', time_step)
      end if

   contains

      subroutine write_xyz(x, y, z, natom, nwalk, urest)
         real(DP), dimension(:, :), intent(in) :: x, y, z
         integer, intent(in) :: natom, nwalk
         integer, intent(in) :: urest
         integer :: iat, iw

         do iw = 1, nwalk
            do iat = 1, natom
               write (urest, '(3ES25.16E3)') x(iat, iw), y(iat, iw), z(iat, iw)
            end do
         end do
      end subroutine write_xyz

   end subroutine restout

   ! Subroutine that reads from restart.xyz during restart
   ! It is called from subroutine init.
   subroutine restin(x, y, z, vx, vy, vz, it)
      use mod_general, only: icv, ihess, nwalk, ipimd, natom, &
                             pot, update_simtime
      use mod_nhc, only: readNHC, inose, nhc_restin
      use mod_mpi, only: get_mpi_rank
      use mod_estimators
      use mod_kinetic, only: entot_cumul, est_temp_cumul
      use mod_sh_integ, only: sh_read_wf
      use mod_sh, only: write_nacmrest, istate
      use mod_lz, only: lz_restin
      use mod_gle, only: gle_restin, pile_restin
      use mod_random, only: read_prng_state
      use mod_terampi_sh, only: read_wfn
      real(DP), intent(out) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: vx(:, :), vy(:, :), vz(:, :)
      integer, intent(out) :: it
      real(DP) :: sim_time
      integer :: iw, urest
      character(len=100) :: chtemp
      character(len=50) :: chin

      chin = append_rank('restart.xyz')

      write (stdout, *) 'irest=1: Reading geometry, velocities and other information from '//chin

      open (newunit=urest, file=chin, status="old", action="read")
      read (urest, *) it, sim_time
      call update_simtime(sim_time)

      read (urest, '(A)') chtemp
      call checkchar(chtemp, chcoords)
      call read_xyz(x, y, z, natom, nwalk, urest)

      read (urest, '(A)') chtemp
      call checkchar(chtemp, chvel)
      call read_xyz(vx, vy, vz, natom, nwalk, urest)

      if (ipimd == 2) then
         read (urest, '(A)') chtemp
         call checkchar(chtemp, chsh)
         read (urest, *) istate
         call sh_read_wf(urest)
      end if

      if (ipimd == 5) then
         read (urest, '(A)') chtemp
         call checkchar(chtemp, chlz)
         call lz_restin(urest, x, y, z, vx, vy, vz)
      end if

      if (inose == 1 .and. readNHC) then
         read (urest, '(A)') chtemp
         call checkchar(chtemp, chnose)
         call nhc_restin(urest)
      end if

      if (inose == 2 .or. inose == 4) then
         read (urest, '(A)') chtemp
         call checkchar(chtemp, chqt)
         call gle_restin(urest)
      end if

      if (inose == 3) then
         read (urest, '(A)') chtemp
         call checkchar(chtemp, chLT)
         call pile_restin(urest)
      end if

      ! reading cumulative averages of various estimators
      read (urest, '(A)') chtemp
      call checkchar(chtemp, chavg)
      read (urest, *) est_temp_cumul
      read (urest, *) est_prim_cumul, est_vir_cumul
      read (urest, *) entot_cumul

      if (icv == 1) then
         read (urest, *) est_prim2_cumul, est_prim_vir, est_vir2_cumul
         read (urest, *) cv_prim_cumul, cv_vir_cumul
         if (ihess == 1) then
            read (urest, *) cv_dcv_cumul
            do iw = 1, nwalk
               read (urest, *) cvhess_cumul(iw)
            end do
         end if
      end if

      ! Read PRNG state (optional)
      ! By now, PRNG is already initialized, so if the state
      ! is not found in the restart file, we simply march on.
      call read_prng_state(urest)

      close (urest)

      if (pot == '_tera_' .and. ipimd == 2) then
         call read_wfn()
      end if
      if (pot == '_tera_' .and. ipimd == 5) then
         call read_wfn()
      end if

   contains

      subroutine checkchar(chin, chref)
         use mod_error, only: fatal_error
         character(len=*), intent(in) :: chin, chref

         if (trim(adjustl(chin)) /= trim(chref)) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Invalid restart file restart.xyz'//new_line('a')//&
               & 'I read "'//trim(adjustl(chin))//'" but expected "'//trim(chref)//'"')
         end if
      end subroutine checkchar

      subroutine read_xyz(x, y, z, natom, nwalk, urest)
         real(DP), dimension(:, :), intent(out) :: x, y, z
         integer, intent(in) :: natom, nwalk
         integer, intent(in) :: urest
         integer :: iat, iw

         do iw = 1, nwalk
            do iat = 1, natom
               read (urest, '(3ES25.16E3)') x(iat, iw), y(iat, iw), z(iat, iw)
            end do
         end do
      end subroutine read_xyz

   end subroutine restin

end module mod_analysis
