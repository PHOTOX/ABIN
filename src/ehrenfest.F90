! Simple Ehrenfest dynamics in a real adiabatic electronic basis (atomic units).

module mod_eh
   use mod_const, only: DP, AUtoFS
   use mod_error, only: fatal_error
   use mod_files, only: stdout, UPOP, UPES
   use mod_eh_integ, only: nstate, phase, popsumthr, eh_set_initialwf, &
                           eh_select_integrator, eh_integrate_wf, get_elpop, &
                           get_coherence, check_popsum
   implicit none
   private

   public :: en_array, nacx, nacy, nacz, eh_phase_nacme, nac_accu1, nac_accu2, &
             read_eh_input, print_eh_input, eh_init, eh_advance, eh_mean_force, &
             state_fx, state_fy, state_fz, old_en, old_nacx, old_nacy, old_nacz, eh_output

   real(DP), allocatable :: en_array(:)
   real(DP), allocatable :: nacx(:, :, :), nacy(:, :, :), nacz(:, :, :)
   integer :: nac_accu1 = 5, nac_accu2 = 5
   integer :: istate_init = 1
   integer :: substep = 100
   character(len=50) :: integ = 'butcher'
   character(len=50) :: couplings = 'analytic'
   real(DP), allocatable :: state_fx(:, :, :), state_fy(:, :, :), state_fz(:, :, :)
   real(DP), allocatable :: old_en(:)
   real(DP), allocatable :: old_nacx(:, :, :), old_nacy(:, :, :), old_nacz(:, :, :)
   ! Initial mean-field potential plus nuclear kinetic energy, in atomic units.
   real(DP), public, protected :: entot0 = 0.0_DP

   namelist /eh/ nstate, istate_init, substep, integ, couplings, popsumthr, nac_accu1, nac_accu2
contains

   subroutine read_eh_input(param_unit)
      use mod_utils, only: tolower
      integer, intent(in) :: param_unit
      integer :: iost

      rewind (param_unit)
      read (param_unit, eh, iostat=iost)
      if (iost /= 0) call fatal_error(__FILE__, __LINE__, 'Missing or invalid &eh namelist')
      rewind (param_unit)
      integ = tolower(integ)
      couplings = tolower(couplings)
      call check_eh_parameters()
      phase = 0
   end subroutine read_eh_input

   subroutine check_eh_parameters()
      use mod_general, only: nwalk, irest, iqmmm, iremd, en_restraint, pot_ref, iknow
      use mod_nhc, only: inose
      use mod_shake, only: nshake
      use mod_system, only: conatom
      use mod_plumed, only: iplumed
      use mod_sbc, only: isbc

      ! If you know, you know. If you don't know, you don't know.
      if (iknow == 1) then
         write (stdout, '(A)') 'Ehrenfest dynamics: iknow = 1, so the check of parameters is skipped.'
         return
      end if

      ! Check that the Ehrenfest parameters are consistent with the dynamics.
      if (nwalk /= 1) &
         call fatal_error(__FILE__, __LINE__, 'nwalk must be 1 for Ehrenfest dynamics')
      if (irest /= 0) &
         call fatal_error(__FILE__, __LINE__, 'irest must be 0 for Ehrenfest dynamics')
      if (iqmmm /= 0) &
         call fatal_error(__FILE__, __LINE__, 'iqmmm must be 0 for Ehrenfest dynamics')
      if (iremd /= 0) &
         call fatal_error(__FILE__, __LINE__, 'iremd must be 0 for Ehrenfest dynamics')
      if (nshake /= 0) &
         call fatal_error(__FILE__, __LINE__, 'nshake must be 0 for Ehrenfest dynamics')
      if (en_restraint /= 0.0_DP) &
         call fatal_error(__FILE__, __LINE__, 'en_restraint must be 0.0 for Ehrenfest dynamics')
      if (pot_ref /= '_none_') &
         call fatal_error(__FILE__, __LINE__, 'pot_ref must be _none_ for Ehrenfest dynamics')
      if (conatom /= 0) &
         call fatal_error(__FILE__, __LINE__, 'conatom must be 0 for Ehrenfest dynamics')
      if (iplumed /= 0) &
         call fatal_error(__FILE__, __LINE__, 'iplumed must be 0 for Ehrenfest dynamics')
      if (isbc /= 0) &
         call fatal_error(__FILE__, __LINE__, 'isbc must be 0 for Ehrenfest dynamics')
      if (nstate < 2) &
         call fatal_error(__FILE__, __LINE__, 'nstate must be >= 2 for Ehrenfest dynamics since for'// &
                          'nstate=1the Ehrenfest dynamics reduces to classical dynamics')
      if (istate_init < 1 .or. istate_init > nstate) &
         call fatal_error(__FILE__, __LINE__, 'istate_init is out of bounds')
      if (substep < 1) &
         call fatal_error(__FILE__, __LINE__, 'substep must be >= 1 for Ehrenfest dynamics')
      if (nac_accu1 < 1 .or. nac_accu2 < 1) &
         call fatal_error(__FILE__, __LINE__, 'Invalid Ehrenfest NAC accuracy')
      if (couplings /= 'analytic') &
         call fatal_error(__FILE__, __LINE__, 'Ehrenfest requires analytic NAC vectors for the nuclear force')
      if (inose > 0) &
         call fatal_error(__FILE__, __LINE__, 'Ehrenfest dynamics is not compatible with thermostats in theory')
   end subroutine check_eh_parameters

   subroutine print_eh_input()
      write (stdout, nml=eh, delim='APOSTROPHE')
   end subroutine print_eh_input

   subroutine eh_init(x, y, z, vx, vy, vz)
      use mod_general, only: natom, pot
      use mod_interfaces, only: eh_force
      use mod_kinetic, only: ekin_v
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: initial_fx(size(x, 1), size(x, 2))
      real(DP) :: initial_fy(size(y, 1), size(y, 2))
      real(DP) :: initial_fz(size(z, 1), size(z, 2)), initial_energy

      call check_eh_parameters()

      !alocate variables
      allocate (nacx(natom, nstate, nstate), nacy(natom, nstate, nstate), nacz(natom, nstate, nstate))
      allocate (state_fx(natom, 1, nstate), state_fy(natom, 1, nstate), state_fz(natom, 1, nstate))
      allocate (en_array(nstate), old_en(nstate))
      allocate (old_nacx(natom, nstate, nstate), old_nacy(natom, nstate, nstate), old_nacz(natom, nstate, nstate))

      ! setup initial values
      en_array = 0.0_DP
      nacx = 0.0_DP; nacy = 0.0_DP; nacz = 0.0_DP
      call eh_set_initialwf(istate_init)
      call eh_select_integrator(integ)
      write (stdout, '(A)') 'Initializing Ehrenfest dynamics (adiabatic mean field)'
      call eh_force(x, y, z, initial_fx, initial_fy, initial_fz, initial_energy, pot)
      entot0 = initial_energy + ekin_v(vx, vy, vz)
   end subroutine eh_init

   ! Keep real adiabatic NAC vectors continuous between nuclear endpoints.
   subroutine eh_phase_nacme()
      integer :: i, j
      real(DP) :: overlap
      do i = 1, nstate - 1
         do j = i + 1, nstate
            overlap = sum(old_nacx(:, i, j) * nacx(:, i, j) + old_nacy(:, i, j) * nacy(:, i, j) + old_nacz(:, i, j) * nacz(:, i, j))
            if (overlap < 0.0_DP) then
               !flip the sign of the NAC vectors to keep them continuous
               nacx(:, i, j) = -nacx(:, i, j)
               nacy(:, i, j) = -nacy(:, i, j)
               nacz(:, i, j) = -nacz(:, i, j)
               !use antisimmetry
               nacx(:, j, i) = -nacx(:, i, j)
               nacy(:, j, i) = -nacy(:, i, j)
               nacz(:, j, i) = -nacz(:, i, j)
            end if
         end do
      end do
   end subroutine eh_phase_nacme

   subroutine eh_mean_force(fx, fy, fz, energy)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: energy
      real(DP) :: weight
      integer :: i, j

      fx = 0.0_DP
      fy = 0.0_DP
      fz = 0.0_DP
      energy = 0.0_DP

      do i = 1, nstate
         ! weight = |c_i|^2
         weight = get_elpop(i)
         ! as I got the populatin I also evaluate the mean energy
         energy = energy + weight * en_array(i)

         ! contribution from the diferent states to the mean force
         fx = fx + weight * state_fx(:, :, i)
         fy = fy + weight * state_fy(:, :, i)
         fz = fz + weight * state_fz(:, :, i)

         ! contribution from the NACs to the mean force
         ! utilising symetry of the NACs, only sum over i<j
         do j = i + 1, nstate
            ! F = sum_i |c_i|^2 F_i - 2 sum_{i<j} Re(c_i* c_j)(E_j-E_i) d_{ij}
            weight = -2 * get_coherence(i, j) * (en_array(j) - en_array(i))
            fx(:, 1) = fx(:, 1) + weight * nacx(:, i, j)
            fy(:, 1) = fy(:, 1) + weight * nacy(:, i, j)
            fz(:, 1) = fz(:, 1) + weight * nacz(:, i, j)
         end do
      end do
   end subroutine eh_mean_force

   ! evaluates the change in the electronic state populations
   subroutine eh_advance(vx, vy, vz, dt)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: dt
      ! later known as dotproduct
      real(DP) :: sigma0(nstate, nstate), sigma1(nstate, nstate)
      real(DP) :: e0(nstate), e1(nstate)
      real(DP) :: d0(nstate, nstate), d1(nstate, nstate)
      real(DP) :: fr0, fr1
      real(DP) :: norm
      integer :: i, j

      ! the weights of the the change in the electronic state populations are given by the NACs and the velocities
      do i = 1, nstate
         do j = 1, nstate
            sigma0(i, j) = sum(vx(:, 1) * old_nacx(:, i, j) + vy(:, 1) * old_nacy(:, i, j) + vz(:, 1) * old_nacz(:, i, j))
            sigma1(i, j) = sum(vx(:, 1) * nacx(:, i, j) + vy(:, 1) * nacy(:, i, j) + vz(:, 1) * nacz(:, i, j))
         end do
      end do

      ! substep integration of the electronic state populations
      do i = 1, substep
         fr0 = real(i - 1, DP) / substep
         fr1 = real(i, DP) / substep
         e0 = (1.0_DP - fr0) * old_en + fr0 * en_array
         e1 = (1.0_DP - fr1) * old_en + fr1 * en_array
         d0 = (1.0_DP - fr0) * sigma0 + fr0 * sigma1
         d1 = (1.0_DP - fr1) * sigma0 + fr1 * sigma1
         call eh_integrate_wf(e0, e1, d0, d1, dt / substep)
         norm = check_popsum()
      end do
   end subroutine eh_advance

   subroutine eh_output()
      use mod_general, only: it, nwrite, sim_time
      integer :: i
      real(DP) :: norm
      if (modulo(it, nwrite) /= 0) return
      norm = check_popsum()
      write (UPOP, '(F15.6,*(ES20.10))') sim_time * AUtoFS, (get_elpop(i), i=1, nstate), norm
      write (UPES, '(F15.6,*(ES20.10))') sim_time * AUtoFS, en_array
   end subroutine eh_output

end module mod_eh
