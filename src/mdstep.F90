! Implementations of various MD integrators.
! - Velocity Verlet with thermostatting for classical MD.
! - RESPA for PIMD
! - RESPA for multiple time-step MD with reference potential.
module mod_mdstep
   use mod_const, only: DP
   use mod_general, only: update_simtime
   use mod_nhc, only: inose
   use mod_transform
   implicit none
   private
   public :: initialize_integrator

   ! Denotes number of internal steps in multiple-time-step RESPA algorithms
   ! nabin - internal steps in PIMD integration,
   integer, public :: nabin = 50
   ! nstep_ref - internal steps with reference potential pot_ref
   integer, public :: nstep_ref = 1

   abstract interface
      subroutine integrator(x, y, z, px, py, pz, amt, dt, E_pot, fx, fy, fz)
         import :: DP
         implicit none
         real(DP), dimension(:, :), intent(inout) :: x, y, z, px, py, pz
         real(DP), dimension(:, :), intent(in) :: amt
         real(DP), intent(in) :: dt
         real(DP), intent(inout) :: E_pot
         real(DP), dimension(:, :), intent(inout) :: fx, fy, fz
      end subroutine integrator
   end interface

   procedure(integrator), pointer, public :: mdstep => null()

contains

   subroutine check_input()
      use mod_utils, only: int_positive
      call int_positive(nabin, 'nabin')
      call int_positive(nstep_ref, 'nstep_ref')
   end subroutine check_input

   subroutine initialize_integrator(dt, ipimd, inormalmodes, nshake, pot, pot_ref)
      use mod_const, only: AUtoFS
      use mod_error, only: fatal_error
      use mod_files, only: stdout
      real(DP), intent(in) :: dt
      integer, intent(in) :: ipimd, inormalmodes, nshake
      character(len=*), intent(in) :: pot, pot_ref

      call check_input()

      if (ipimd == 0) then
         write (stdout, *) 'Using velocity Verlet integrator'
         mdstep => verletstep
         nabin = 1
      end if

      ! Selecting proper integrator for a given MD type
      if (ipimd == 2 .or. ipimd == 5) then
         mdstep => verletstep
         nabin = 1
      else if (ipimd == 1 .and. inormalmodes /= 1) then
         write (stdout, *) 'Using RESPA integrator for PIMD.'
         mdstep => respastep
      else if (ipimd == 1 .and. inormalmodes == 1) then
         write (stdout, *) 'Using velocity Verlet integrator with analytical PI normal mode propagation.'
         mdstep => verletstep
         nabin = 1
      end if

      if (nshake /= 0) then
         mdstep => respashakestep
      end if

      if (pot_ref /= '_none_') then
         mdstep => doublerespastep
         write (stdout, '(A)') 'Using Multiple Time-Step RESPA integrator'
         write (stdout, '(A)') "Reference (cheap) potential is "//trim(pot_ref)
         write (stdout, '(A, F6.2)') "with timestep [fs] ", dt / nstep_ref * AUtoFS
         write (stdout, '(A)') "Full potential is "//trim(pot)
         write (stdout, '(A, F6.2)') "with timestep [fs] ", dt * AUtoFS
         if (ipimd /= 0) then
            call fatal_error(__FILE__, __LINE__, &
               & 'ab initio MTS is implemented only for classical MD')
         end if
         if (nshake > 0) then
            call fatal_error(__FILE__, __LINE__, &
               & 'ab initio MTS is not implemented with SHAKE, apologies.')
         end if
      end if

      ! NOTE: This checks for a programming error, we sanitize user input earlier
      if (.not. associated(mdstep) .and. ipimd /= 3) then
         call fatal_error(__FILE__, __LINE__, 'invalid integrator')
      end if
   end subroutine initialize_integrator

   subroutine shiftX(rx, ry, rz, px, py, pz, mass, dt)
      real(DP), intent(inout) :: rx(:, :), ry(:, :), rz(:, :)
      real(DP), intent(in) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: mass(:, :)
      real(DP), intent(in) :: dt

      RX = RX + dt * PX / MASS
      RY = RY + dt * PY / MASS
      RZ = RZ + dt * PZ / MASS

   end subroutine shiftX

   subroutine shiftP(px, py, pz, fx, fy, fz, dt)
      use mod_system, only: conatom
      use mod_vinit, only: constrainP
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(in) :: dt

      PX = PX + dt * FX
      PY = PY + dt * FY
      PZ = PZ + dt * FZ

      if (conatom > 0) call constrainP(px, py, pz, conatom)

   end subroutine shiftP

   ! Based on:
   ! Efficient stochastic thermostatting of path integral molecular dynamics
   ! J. Chem. Phys. 133, 124104 ?2010?
   subroutine propagate_nm(x, y, z, px, py, pz, m, dt)
      use mod_const, only: DP, PI
      use mod_general, only: nwalk, natom
      use mod_nhc, only: temp
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: m(:, :)
      real(DP), intent(in) :: dt
      real(DP), dimension(natom, nwalk) :: x_old, y_old, z_old
      real(DP), dimension(natom, nwalk) :: px_old, py_old, pz_old
      real(DP), dimension(nwalk) :: omega(nwalk)
      real(DP) :: omega_n, om, omt, c, s
      integer :: iat, iw, k

      ! expecting m = am = physical masses

      x_old = x
      y_old = y
      z_old = z
      px_old = px
      py_old = py
      pz_old = pz

      omega_n = temp * nwalk

      do iw = 2, nwalk
         k = iw - 1
         omega(iw) = 2 * omega_n * sin(k * PI / nwalk)
      end do

      ! First, propagate centroid
      iw = 1
      do iat = 1, natom
         x(iat, iw) = x(iat, iw) + dt * px(iat, iw) / m(iat, iw)
         y(iat, iw) = y(iat, iw) + dt * py(iat, iw) / m(iat, iw)
         z(iat, iw) = z(iat, iw) + dt * pz(iat, iw) / m(iat, iw)
      end do

      ! eq 23 from J. Chem. Phys. 133, 124104 2010
      ! exact propagation of a free ring polymer in normal mode coordinates
      do iw = 2, nwalk
         om = omega(iw)
         omt = omega(iw) * dt
         c = cos(omt)
         s = sin(omt)
         do iat = 1, natom
            ! Propagate positions
            x(iat, iw) = x_old(iat, iw) * c &
                         + px_old(iat, iw) * s / m(iat, iw) / om
            y(iat, iw) = y_old(iat, iw) * c &
                         + py_old(iat, iw) * s / m(iat, iw) / om
            z(iat, iw) = z_old(iat, iw) * c &
                         + pz_old(iat, iw) * s / m(iat, iw) / om

            ! propagate momenta
            px(iat, iw) = px_old(iat, iw) * c &
                          - x_old(iat, iw) * s * m(iat, iw) * om
            py(iat, iw) = py_old(iat, iw) * c &
                          - y_old(iat, iw) * s * m(iat, iw) * om
            pz(iat, iw) = pz_old(iat, iw) * c &
                          - z_old(iat, iw) * s * m(iat, iw) * om

         end do

      end do

   end subroutine propagate_nm

   subroutine thermostat(px, py, pz, amt, dt)
      use mod_nhc, only: inose, shiftNHC_yosh, shiftNHC_yosh_mass, imasst
      use mod_gle, only: gle_step, pile_step
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt

      if (inose == 1) then

         if (imasst == 1) then
            call shiftNHC_yosh_mass(px, py, pz, amt, dt)
         else
            call shiftNHC_yosh(px, py, pz, amt, dt)
         end if

      else if (inose == 2 .or. inose == 4) then

         ! colored-noise thermostats
         call gle_step(px, py, pz, amt)

      else if (inose == 3) then

         ! white-noise thermostat (Langevin, PILE)
         call pile_step(px, py, pz, amt)

      end if

   end subroutine thermostat

   ! VELOCITY VERLET INTEGRATOR WITH THERMOSTATS
   ! At this moment it does not contain shake routines,
   ! which are only in respashake function.
   ! Works with NHC, GLE and PIGLET thermostats.
   ! Contains propagation of normal modes according to:
   ! eq 23 from J. Chem. Phys. 133, 124104 2010
   subroutine verletstep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
      use mod_general, only: pot, inormalmodes, en_restraint
      use mod_interfaces, only: force_clas
      use mod_en_restraint
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas

      if (inose > 0) call thermostat(px, py, pz, amt, dt / 2)

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)
      if (en_restraint >= 1) call energy_restraint(x, y, z, px, py, pz, eclas)

      if (inormalmodes == 1) then
         ! Warning, initial hack, passing amt here
         call propagate_nm(x, y, z, px, py, pz, amt, dt)
      else
         call shiftX(x, y, z, px, py, pz, amt, dt)
      end if

      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      if (inose > 0) call thermostat(px, py, pz, amt, dt / 2)

      call update_simtime(dt)
   end subroutine verletstep

   ! RESPA ALGORITHM  10.12.2012
   ! General algorithm of the propagation:
   ! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)
   ! further info is before subroutine respa_shake
   subroutine respastep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
      use mod_general, only: pot
      use mod_interfaces, only: force_clas, force_quantum
      use mod_arrays, only: fxq, fyq, fzq, amg
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas
      real(DP) :: equant
      integer :: iabin

      call thermostat(px, py, pz, amt, dt / (2 * nabin))

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      do iabin = 1, nabin

         if (inose > 0 .and. iabin /= 1) call thermostat(px, py, pz, amt, dt / (2 * nabin))

         call shiftP(px, py, pz, fxq, fyq, fzq, dt / (2 * nabin))

         call shiftX(x, y, z, px, py, pz, amt, dt / nabin)

         call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)

         call shiftP(px, py, pz, fxq, fyq, fzq, dt / (2 * nabin))

         if (inose > 0 .and. iabin /= nabin) call thermostat(px, py, pz, amt, dt / (2 * nabin))

      end do

      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      call thermostat(px, py, pz, amt, dt / (2 * nabin))

      call update_simtime(dt)
   end subroutine respastep

   ! RESPA ALGORITHM WITH RATTLE     10.12.2012
   ! General algorithm of the propagation:
   ! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)
   !
   ! Implementation of shake+RESPA according to www.chim.unifi.it/orac/MAN/node9.html
   !
   ! Some discussion concerning the Shake+nose-hoover may be found in:
   ! p.199,M.Tuckermann, Statistical mechanics.. (2010)
   ! Basically, it says that global NHC should not compromise constraints.
   ! However, this doesn't apply to massive thermostatting that we use for PIMD!
   !
   ! Therefore, we use different thermostating scheme here (subroutine shiftNHC_yosh).
   ! User must define "molecules", whose atoms share constraints.
   ! Atoms which are not part of any constraint can be molecules themselves.
   ! We then append each "molecule" with its own thermostat.
   ! Definition of molecules are controlled by variables nmolt,natmolt nshakemol.
   ! Note that molecules must be in sequential order.
   ! For global thermostattting (i.e. dynamics is not disturbed that much), specify the system as one molecule.
   ! But don't do this with PIMD!
   !
   ! For example, for system of chloride and 5 rigid water molecules,
   ! there will be 6 molecules and 6 NHC chains.
   subroutine respashakestep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
      use mod_general, only: pot
      use mod_nhc, only: shiftNHC_yosh
      use mod_shake, only: shake, nshake
      use mod_interfaces, only: force_clas, force_quantum
      use mod_arrays, only: amg, fxq, fyq, fzq
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP) :: equant
      integer :: iabin

      if (inose == 1) then
         call shiftNHC_yosh(px, py, pz, amt, dt / (2 * nabin))
      end if

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      ! RATTLE
      if (nshake >= 1) call shake(x, y, z, px, py, pz, amt, 0, 1)

      if (inose == 1) then
         call shiftNHC_yosh(px, py, pz, amt, -dt / (2 * nabin))
      end if

      ! RESPA LOOP
      do iabin = 1, nabin

         if (inose == 1) call shiftNHC_yosh(px, py, pz, amt, dt / (2 * nabin))

         call shiftP(px, py, pz, fxq, fyq, fzq, dt / (2 * nabin))

         ! RATTLE
         if (nshake >= 1) call shake(x, y, z, px, py, pz, amt, 0, 1)

         call shiftX(x, y, z, px, py, pz, amt, dt / nabin)

         ! SHAKE , iq=1, iv=0
         if (NShake > 0) call shake(x, y, z, px, py, pz, amt, 1, 0)

         call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)

         call shiftP(px, py, pz, fxq, fyq, fzq, dt / (2 * nabin))

         if (inose == 1) call shiftNHC_yosh(px, py, pz, amt, dt / (2 * nabin))

      end do
      ! END OF RESPA LOOP

      if (inose == 1) call shiftNHC_yosh(px, py, pz, amt, -dt / (2 * nabin))

      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      if (inose == 1) call shiftNHC_yosh(px, py, pz, amt, dt / (2 * nabin))

      call update_simtime(dt)
   end subroutine respashakestep

   ! Double RESPA algorithm using reference low-cost potential pot_ref with smaller time step
   subroutine doublerespastep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
      use mod_general, only: pot, pot_ref
      use mod_interfaces, only: force_clas, force_quantum
      use mod_arrays, only: fxc_diff, fyc_diff, fzc_diff, &
                          & fxq, fyq, fzq, &
                          & amg
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas
      real(DP) :: equant
      real(DP) :: dtsm
      integer :: iabin, iref

      dtsm = dt / (nabin * nstep_ref)

      if (inose > 0) call thermostat(px, py, pz, amt, dtsm / 2)

      call shiftP(px, py, pz, fxc_diff, fyc_diff, fzc_diff, dt / 2)

      do iref = 1, nstep_ref

         if (inose > 0 .and. iref /= 1) call thermostat(px, py, pz, amt, dtsm / 2)

         call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2 / nstep_ref)

         do iabin = 1, nabin

            if (inose > 0 .and. iabin /= 1) call thermostat(px, py, pz, amt, dtsm / 2)

            call shiftP(px, py, pz, fxq, fyq, fzq, dtsm / 2)

            call shiftX(x, y, z, px, py, pz, amt, dtsm)

            call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)

            call shiftP(px, py, pz, fxq, fyq, fzq, dtsm / 2)

            if (inose > 0 .and. iabin /= nabin) call thermostat(px, py, pz, amt, dtsm / 2)

         end do

         call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot_ref)

         call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2 / nstep_ref)

         if (inose > 0 .and. iref /= nstep_ref) call thermostat(px, py, pz, amt, dtsm / 2)

      end do

      call force_clas(fxc_diff, fyc_diff, fzc_diff, x, y, z, eclas, pot)

      call shiftP(px, py, pz, fxc_diff, fyc_diff, fzc_diff, dt / 2)

      if (inose > 0) call thermostat(px, py, pz, amt, dtsm / 2)

      call update_simtime(dt)
   end subroutine doublerespastep

end module mod_mdstep
