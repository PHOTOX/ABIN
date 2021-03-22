! Implementations of various MD integrators.
! - Velocity Verlet with thermostatting for classical MD.
! - RESPA for PIMD
! - RESPA for multiple time-step MD with reference potentials.
module mod_mdstep
   use mod_const, only: DP
   use mod_kinetic, only: ekin_p
   use mod_utils, only: abinerror
   use mod_transform
   implicit none
   private
   public :: verletstep, respastep, respashake, doublerespastep

contains

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

   subroutine thermostat(px, py, pz, amt, dt)
      use mod_nhc, only: inose, shiftNHC_yosh, shiftNHC_yosh_mass, imasst
      use mod_gle, only: gle_step, pile_step, langham
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: dt

      if (inose == 1) then

         if (imasst == 1) then
            call shiftNHC_yosh_mass(px, py, pz, amt, dt)
         else
            call shiftNHC_yosh(px, py, pz, amt, dt)
         end if

         ! colored-noise thermostats
      else if (inose == 2) then

         langham = langham + ekin_p(px, py, pz)
         call gle_step(px, py, pz, amt)
         langham = langham - ekin_p(px, py, pz)

         ! white-noise thermostat (PILE)
      else if (inose == 3) then

         langham = langham + ekin_p(px, py, pz)
         call pile_step(px, py, pz, amt)
         langham = langham - ekin_p(px, py, pz)

      end if

   end subroutine thermostat

!  VELOCITY VERLET INTEGRATOR WITH THERMOSTATS
!  At this moment it does not contain shake routines,
!  which are only in respashake function.
!  Works with NHC, GLE and PIGLET thermostats.
!  Contains propagation of normal modes according to:
!  eq 23 from J. Chem. Phys. 133, 124104 2010
   subroutine verletstep(x, y, z, px, py, pz, amt, dt, eclas, fxc, fyc, fzc)
      use mod_general, only: pot, ipimd, inormalmodes, en_restraint
      use mod_nhc, only: inose
      use mod_interfaces, only: force_clas, propagate_nm
      ! Not yet implemented
      !use mod_sh, only: ehrenfest_forces
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

      if (ipimd == 4) then
         ! This is only a stub for now
         write (*, *) 'ERROR: Ehrenfest MD not implemented yet!!'
         call abinerror('verletstep')
         !call ehrenfest_forces(x, y, z, fxc, fyc, fzc, px, py, pz, dt, eclas)
      end if

      call shiftP(px, py, pz, fxc, fyc, fzc, dt / 2)

      if (inose > 0) call thermostat(px, py, pz, amt, dt / 2)

   end subroutine verletstep

   ! RESPA ALGORITHM  10.12.2012
   ! General algorithm of the propagation:
   ! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)
   ! further info is before subroutine respa_shake
   subroutine respastep(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, &
                        fxc, fyc, fzc, fxq, fyq, fzq)
      use mod_general, only: nabin, pot
      use mod_nhc, only: inose
      use mod_interfaces, only: force_clas, force_quantum
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: fxq(:, :), fyq(:, :), fzq(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amg(:, :), amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas, equant
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

   end subroutine respastep

! RESPA ALGORITHM WITH RATTLE     10.12.2012
! General algorithm of the propagation:
! see eq. 64 in Mol. Phys. 1996, vol. 87, 1117 (Martyna et al.)

! Implementation of shake+RESPA according to www.chim.unifi.it/orac/MAN/node9.html

! Some discussion concerning the Shake+nose-hoover may be found in:
! p.199,M.Tuckermann, Statistical mechanics.. (2010)
! Basically, it says that global NHC should not compromise constraints.
! However, this doesn't apply to massive thermostatting that we use for PIMD!

! Therefore, we use different thermostating scheme here (subroutine shiftNHC_yosh).
! User must define "molecules", whose atoms share constraints.
! Atoms which are not part of any constraint can be molecules themselves.
! We then append each "molecule" with its own thermostat.
! Definition of molecules are controlled by variables nmolt,natmolt nshakemol.
! Note that molecules must be in sequential order.
! For global thermostattting (i.e. dynamics is not disturbed that much), specify the system as one molecule.
! But don't do this with PIMD!

! For example, for system of chloride and 5 rigid water molecules,
! there will be 6 molecules and 6 NHC chains.

   subroutine respashake(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, &
                         fxc, fyc, fzc, fxq, fyq, fzq)
      use mod_general, only: nabin, pot
      use mod_nhc, only: inose, shiftNHC_yosh, shiftNHC_yosh_mass
      use mod_shake, only: shake, nshake
      use mod_interfaces, only: force_clas, force_quantum
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amg(:, :), amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas, equant
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: fxq(:, :), fyq(:, :), fzq(:, :)
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

   end subroutine respashake

   ! Double RESPA ALGORITHM using reference semiempirical potential (or any other potential)
   subroutine doublerespastep(x, y, z, px, py, pz, amt, amg, dt, equant, eclas, &
                              fxc, fyc, fzc, fxq, fyq, fzq)
      use mod_general, only: nabin, pot, pot_ref, nstep_ref
      use mod_nhc, only: inose
      use mod_interfaces, only: force_clas, force_quantum
      use mod_arrays, only: fxc_diff, fyc_diff, fzc_diff
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(inout) :: fxq(:, :), fyq(:, :), fzq(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amg(:, :), amt(:, :)
      real(DP), intent(in) :: dt
      real(DP), intent(inout) :: eclas, equant
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

   end subroutine doublerespastep

end module mod_mdstep
