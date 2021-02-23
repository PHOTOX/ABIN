! This module contains some routines that do analyses and I/O operations.
! It also contains routines performing restart.
! TODO: Separate restart to its own module
module mod_analysis
   use mod_const, only: DP
   use mod_files
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
   subroutine analysis(x, y, z, vx, vy, vz, fxc, fyc, fzc, eclas, equant)
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
      real(DP), intent(in) :: eclas, equant
      real(DP) :: energy

      ! eclas is the ab initio energy averaged per bead,
      ! equant is additional harmonic energy between PI beads (from force_quantum)
      ! TODO: Print equant or energy somewhere
      energy = eclas + equant

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
      use mod_general, only: nwalk, natom, iremd, my_rank, sim_time
      use mod_system, only: names
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: time_step
      integer :: iat, iw
      character(len=20) :: chout
      logical :: lopened

      inquire (UMOVIE, opened=lopened)

      if (.not. lopened) then
         chout = 'movie.xyz'
         if (iremd == 1) write (chout, '(A,I2.2)') trim(chout)//'.', my_rank
         open (UMOVIE, file=chout, access='append', action="write")
      end if

      ! printing with slightly lower precision for saving space
      ! could be probably much lower
10    format(A2, 3E18.8E2)

      do iw = 1, nwalk
         write (UMOVIE, *) natom
         ! In the future, we should get rid of the time step?
         write (UMOVIE, '(A10,I20,A15,F15.2)') 'Time step:', time_step, ' Sim. Time [au]', sim_time
         do iat = 1, natom
            write (UMOVIE, 10) names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
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
      character(len=40) :: fgeom, fkom

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

      fgeom = '(A2,3E18.10E2)'
      fkom = '(A10,3E13.5E2,A14,3E13.5E2)'

      ! TODO: Include somehow the timestep in the output
      ! Either here:
      ! write(funit, *)natom, it
      ! or maybe better here
      ! write(funit,fkom)'time step', it, &
      !                & 'net force:',fx_tot,fy_tot,fz_tot, &
      !                & 'torque force:',fx_rot,fy_rot,fz_rot
      write (funit, *) natom
      write (funit, fkom) 'net force:', fx_tot, fy_tot, fz_tot, &
                        & ' torque force:', fx_rot, fy_rot, fz_rot
      do iw = 1, nwalk
         do iat = 1, natom
            ! Printing in atomic units
            write (funit, fgeom) names(iat), fx(iat, iw), fy(iat, iw), fz(iat, iw)
         end do
      end do

   end subroutine forceout

   subroutine velout(vx, vy, vz)
      use mod_general, only: nwalk, natom, it
      use mod_system, only: names
      use mod_files, only: UVELOC
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer :: iat, iw
      character(len=20) :: fgeom

      fgeom = '(A2,3E18.10E2)'
      write (UVELOC, *) natom
      write (UVELOC, *) 'Time step:', it

      do iw = 1, nwalk
         do iat = 1, natom
            ! Printing in atomic units
            write (UVELOC, fgeom) names(iat), vx(iat, iw), vy(iat, iw), vz(iat, iw)
         end do
      end do

   end subroutine velout

   subroutine restout(x, y, z, vx, vy, vz, time_step)
      use mod_general, only: icv, ihess, nwalk, ipimd, natom, &
                             iremd, my_rank, pot, narchive, sim_time
      use mod_utils, only: archive_file
      use mod_nhc, only: inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
      use mod_estimators
      use mod_kinetic, only: entot_cumul, est_temp_cumul
      use mod_sh_integ, only: sh_write_wf
      use mod_sh, only: write_nacmrest, ntraj, istate
      use mod_lz, only: lz_restout
      use mod_gle
      use mod_random
      use mod_terampi_sh, only: write_wfn
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer, intent(in) :: time_step
      integer :: iat, iw, inh, itrj, is
      logical :: file_exists
      character(len=200) :: chout, chsystem, chformat

      if (pot == '_tera_' .and. ipimd == 2) call write_wfn()
      if (pot == '_tera_' .and. ipimd == 5) call write_wfn()

      if (iremd == 1) then
         write (chout, '(A,I2.2)') 'restart.xyz.', my_rank
      else
         chout = 'restart.xyz'
      end if

      inquire (FILE=chout, EXIST=file_exists)
      if (file_exists) then
         chsystem = 'cp '//trim(chout)//'  '//trim(chout)//'.old'
         if (iremd == 1) then
            call system(chsystem)
         else if (my_rank == 0) then
            call system(chsystem)
         end if
      end if

      ! open(102, file=chout, action='WRITE',recl=250)
      ! intel compilers don't write too many columns on single line
      open (102, file=chout, action='WRITE')

      write (102, *) time_step, sim_time

      write (102, *) chcoords
      do iw = 1, nwalk
         do iat = 1, natom
            write (102, *) x(iat, iw), y(iat, iw), z(iat, iw)
         end do
      end do

      write (102, *) chvel

      do iw = 1, nwalk
         do iat = 1, natom
            write (102, *) vx(iat, iw), vy(iat, iw), vz(iat, iw)
         end do
      end do

      if (ipimd == 2) then
         call write_nacmrest()
         ! hack, only one trajectory supported at this point
         itrj = 1
         write (102, *) chSH
         write (102, *) istate(itrj)
         call sh_write_wf(102, itrj)
      end if

      if (ipimd == 5) then
         write (102, *) chLZ
         call lz_restout(102)
      end if

      if (inose == 1) then
         write (102, *) chnose

         if (imasst == 1) then
            do inh = 1, nchain
               do iw = 1, nwalk
                  do iat = 1, natom
                     write (102, *) pnhx(iat, iw, inh), pnhy(iat, iw, inh), pnhz(iat, iw, inh)
                  end do
               end do
            end do

         else

            do inh = 1, nchain
               do iw = 1, nwalk
                  do iat = 1, nmolt
                     write (102, *) pnhx(iat, iw, inh)
                  end do
               end do
            end do
         end if

      end if

      if (inose == 2) then
         write (102, *) chQT
         write (chformat, '(A1,I1,A7)') '(', ns, 'E25.16)'
         do iw = 1, nwalk
            do iat = 1, natom * 3
               write (102, fmt=chformat) (ps(iat, is, iw), is=1, ns)
            end do
         end do
         write (102, *) langham
      end if
      if (inose == 3) then
         write (102, *) chLT
         write (102, *) langham
      end if

      write (102, *) chAVG
      write (102, *) est_temp_cumul
      write (102, *) est_prim_cumul, est_vir_cumul
      write (102, *) entot_cumul

      if (icv == 1) then
         write (102, '(3E25.16)') est_prim2_cumul, est_prim_vir, est_vir2_cumul
         write (102, '(2E25.16)') cv_prim_cumul, cv_vir_cumul
         if (ihess == 1) then
            write (102, *) cv_dcv_cumul
            do iw = 1, nwalk
               write (102, *) cvhess_cumul(iw)
            end do
         end if
      end if

      ! write current state of PRNG
      call rsavef(102)

      close (102)

      if (modulo(time_step, narchive) == 0) then
         call archive_file('restart.xyz', time_step)
      end if

   end subroutine restout

   ! Subroutine that reads from restart.xyz during restart
   ! It is called from subroutine init.
   subroutine restin(x, y, z, vx, vy, vz, it)
      use mod_general, only: icv, ihess, nwalk, ipimd, natom, &
                             iremd, my_rank, pot, sim_time
      use mod_nhc, only: readNHC, inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
      use mod_estimators
      use mod_kinetic, only: entot_cumul, est_temp_cumul
      use mod_sh_integ, only: sh_read_wf
      use mod_sh, only: write_nacmrest, ntraj, istate
      use mod_lz, only: lz_restin
      use mod_gle
      use mod_random
      use mod_terampi_sh, only: read_wfn
      real(DP), intent(out) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: vx(:, :), vy(:, :), vz(:, :)
      integer, intent(out) :: it
      integer :: iat, iw, inh, itrj, is
      character(len=100) :: chtemp
      character(len=200) :: chformat
      logical :: prngread
      character(len=20) :: chin

      prngread = .false.

      if (iremd == 1) then
         write (chin, '(A,I2.2)') 'restart.xyz.', my_rank
      else
         chin = 'restart.xyz'
      end if

      if (my_rank == 0) then
         write (*, *) 'irest=1, Reading geometry, velocities and other information from restart.xyz'
      end if

      ! open(111, file=chin, status = "OLD", action = "READ", recl=1500)
      open (111, file=chin, status="OLD", action="READ")
      read (111, *) it, sim_time
      read (111, '(A)') chtemp
      call checkchar(chtemp, chcoords)
      do iw = 1, nwalk
         do iat = 1, natom
            read (111, *) x(iat, iw), y(iat, iw), z(iat, iw)
         end do
      end do

      read (111, '(A)') chtemp
      call checkchar(chtemp, chvel)
      do iw = 1, nwalk
         do iat = 1, natom
            read (111, *) vx(iat, iw), vy(iat, iw), vz(iat, iw)
         end do
      end do

      if (ipimd == 2) then
         read (111, '(A)') chtemp
         call checkchar(chtemp, chsh)
         ! only 1 trajectory supported at this point
         itrj = 1
         read (111, *) istate(itrj)
         call sh_read_wf(111, itrj)
      end if

      if (ipimd == 5) then
         read (111, '(A)') chtemp
         call checkchar(chtemp, chlz)
         call lz_restin(111, x, y, z, vx, vy, vz)
      end if

      if (inose == 1 .and. readNHC == 1) then
         read (111, '(A)') chtemp
         call checkchar(chtemp, chnose)
         if (imasst == 1) then
            do inh = 1, nchain
               do iw = 1, nwalk
                  do iat = 1, natom
                     read (111, *) pnhx(iat, iw, inh), pnhy(iat, iw, inh), pnhz(iat, iw, inh)
                  end do
               end do
            end do

         else

            do inh = 1, nchain
               do iw = 1, nwalk
                  do iat = 1, nmolt
                     read (111, *) pnhx(iat, iw, inh)
                  end do
               end do
            end do

         end if

      end if

      if (inose == 2 .and. readQT == 1) then
         read (111, '(A)') chtemp
         call checkchar(chtemp, chqt)
         write (chformat, '(A1,I1,A7)') '(', ns, 'E25.16)'
         do iw = 1, nwalk
            do iat = 1, natom * 3
               do is = 1, ns - 1
                  !read(111, fmt=chformat)(ps(iat, is, iw), is = 1, ns)
                  read (111, '(1E25.16)', advance="no") ps(iat, is, iw)
               end do
               read (111, '(1E25.16)') ps(iat, ns, iw)
            end do
         end do
         read (111, *) langham
      end if

      if (inose == 3 .and. readQT == 1) then
         read (111, '(A)') chtemp
         call checkchar(chtemp, chLT)
         read (111, *) langham
      end if

      ! reading cumulative averages of various estimators
      read (111, '(A)') chtemp
      call checkchar(chtemp, chavg)
      read (111, *) est_temp_cumul
      read (111, *) est_prim_cumul, est_vir_cumul
      read (111, *) entot_cumul

      if (icv == 1) then
         read (111, *) est_prim2_cumul, est_prim_vir, est_vir2_cumul
         read (111, *) cv_prim_cumul, cv_vir_cumul
         if (ihess == 1) then
            read (111, *) cv_dcv_cumul
            do iw = 1, nwalk
               read (111, *) cvhess_cumul(iw)
            end do
         end if
      end if

      ! Trying to restart PRNG
      ! prngread is optional argument determining, whether we write or read
      ! currently,prngread is not used, since vranf is initialize BEFORE restart
      ! and is possibly rewritten here
      call rsavef(111, prngread)

      close (111)

      if (pot == '_tera_' .and. ipimd == 2) call read_wfn()
      if (pot == '_tera_' .and. ipimd == 5) call read_wfn()

   contains

      subroutine checkchar(chin, chref)
         use mod_utils, only: abinerror
         character(len=*) :: chin, chref

         if (trim(adjustl(chin)) /= trim(chref)) then
            write (*, *) 'ERROR while reading from restart.xyz.'
            write (*, *) 'I read: ', chin
            write (*, *) 'but expected: ', chref
            call abinerror('restin')
         end if
      end subroutine checkchar

   end subroutine restin

end module mod_analysis
