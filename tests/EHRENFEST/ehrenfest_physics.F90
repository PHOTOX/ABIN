! Standalone analytic checks, run with make ehrenfest-test.
program ehrenfest_physics
   use mod_const, only: DP, ANG, imaginary_unit
   use mod_general, only: set_natom, nwalk, pot
   use mod_system, only: am
   use mod_nhc, only: inose
   use mod_eh, only: eh_init, eh_mean_force, eh_advance, en_array, nacx, old_en, old_nacx, entot0
   use mod_interfaces, only: eh_force
   use mod_eh_integ, only: nstate, eh_read_wf, eh_write_wf, get_elpop, get_coherence, check_popsum, &
                           phase, eh_read_phase_bin, eh_select_integrator, eh_integrate_wf, eh_set_energy_shift
   use mod_sh_integ, only: sh_nstate => nstate, sh_set_initialwf, sh_get_elpop => get_elpop
   implicit none
   real(DP) :: x(2, 1), y(2, 1), z(2, 1), fx(2, 1), fy(2, 1), fz(2, 1)
   real(DP) :: f1(2, 1), f2(2, 1)
   real(DP) :: e1, e2, energy, norm
   real(DP) :: expected(2)
   real(DP) :: vx(2, 1), vy(2, 1), vz(2, 1)
   real(DP) :: sigma, omega, dt, expected_population
   integer :: u

   call set_natom(2)
   nwalk = 1; nstate = 2; inose = 0; pot = '_nai_'
   sh_nstate = 3
   call sh_set_initialwf(2)
   x(:, 1) = [0.0_DP, 6.93_DP * ANG]
   y = 0.0_DP; z = 0.0_DP
   allocate (am(2))
   am = [2.0_DP, 3.0_DP]
   vx(:, 1) = [0.01_DP, -0.02_DP]
   vy(:, 1) = [0.03_DP, -0.04_DP]
   vz(:, 1) = [0.05_DP, -0.06_DP]
   call eh_init(x, y, z, vx, vy, vz)
   call eh_mean_force(f1, fy, fz, e1)
   if (abs(e1 - en_array(1)) > 1.E-14_DP) error stop 'Initial electronic energy'
   if (abs(entot0 - e1 - 0.5_DP * sum(am * (vx(:, 1)**2 + vy(:, 1)**2 + vz(:, 1)**2))) > &
       1.E-14_DP) error stop 'Initial total energy'
   if (maxval(abs(old_en - en_array)) > 1.E-14_DP) error stop 'Initial energy history'
   if (maxval(abs(old_nacx - nacx)) > 1.E-14_DP) error stop 'Initial NAC history'
   call set_wf(0.0_DP, 1.0_DP)
   call eh_mean_force(f2, fy, fz, e2)

   ! A real equal superposition has an off-diagonal force. Changing its
   ! relative phase to pi/2 must remove that term, at the same populations.
   call set_wf(sqrt(0.5_DP), sqrt(0.5_DP))
   call eh_mean_force(fx, fy, fz, energy)
   expected = (f1(:, 1) + f2(:, 1)) / 2 - (e2 - e1) * nacx(:, 1, 2)
   if (maxval(abs(fx(:, 1) - expected)) > 1.E-13_DP) error stop 'Coherence force'
   if (abs(energy - (e1 + e2) / 2) > 1.E-14_DP) error stop 'Mean energy'
   call set_wf(sqrt(0.5_DP), sqrt(0.5_DP), imaginary=.true.)
   call eh_mean_force(fx, fy, fz, energy)
   if (maxval(abs(fx - (f1 + f2) / 2)) > 1.E-13_DP) error stop 'Imaginary coherence'
   if (abs(sum(fx)) > 1.E-13_DP) error stop 'Momentum conservation'

   ! At fixed geometry and constant velocity, the 2x2 electronic generator
   ! has an analytic exponential, including both the gap and NAC terms.
   call set_wf(1.0_DP, 0.0_DP)
   call eh_force(x, y, z, fx, fy, fz, energy, pot)
   vx(:, 1) = [-0.002_DP, 0.002_DP]
   vy = 0.0_DP; vz = 0.0_DP
   sigma = sum(vx(:, 1) * nacx(:, 1, 2))
   omega = sqrt(sigma**2 + ((e2 - e1) / 2)**2)
   dt = 100.0_DP
   expected_population = (sigma / omega)**2 * sin(omega * dt)**2
   call eh_advance(vx, vy, vz, dt)
   if (abs(get_elpop(2) - expected_population) > 1.E-10_DP) error stop 'Electronic propagation'
   norm = check_popsum()
   if (abs(norm - 1.0_DP) > 1.E-10_DP) error stop 'Electronic norm'
   if (sh_nstate /= 3) error stop 'SH state count was modified by Ehrenfest'
   if (sh_get_elpop(2) /= 1.0_DP) error stop 'SH wavefunction was modified by Ehrenfest'
   call check_electronic_integrators()
   print*,'Ehrenfest force, coherence, electronic propagation and norm checks passed.'

contains
   subroutine check_electronic_integrators()
      character(len=7), parameter :: methods(3) = [character(len=7) :: 'butcher', 'rk4', 'euler']
      real(DP) :: en(2), coupling(2, 2), tolerance, elapsed, step, frequency, half_gap
      real(DP) :: expected_coherence
      complex(2*DP) :: expected_c(2), actual_c(2), stored_c(2), global_phase
      integer :: mode, method, k, phase_unit, wf_unit

      en = [0.1_DP, 0.3_DP]
      coupling = 0.0_DP
      coupling(1, 2) = 0.04_DP
      coupling(2, 1) = -coupling(1, 2)
      step = 0.01_DP
      elapsed = 100 * step
      half_gap = (en(2) - en(1)) / 2
      frequency = sqrt(half_gap**2 + coupling(1, 2)**2)
      global_phase = exp(-imaginary_unit * sum(en) * elapsed / 2)
      expected_c(1) = global_phase * (cos(frequency * elapsed) + &
                                      imaginary_unit * half_gap / frequency * sin(frequency * elapsed))
      expected_c(2) = global_phase * coupling(1, 2) / frequency * sin(frequency * elapsed)
      expected_coherence = real(conjg(expected_c(1)) * expected_c(2), kind=DP)
      call eh_set_energy_shift(0.0_DP)

      do mode = 0, 1
         phase = mode
         do method = 1, size(methods)
            call set_wf(1.0_DP, 0.0_DP)
            open (newunit=phase_unit, status='scratch', form='unformatted')
            write (phase_unit) 0.0_DP, 0.0_DP
            write (phase_unit) 0.0_DP, 0.0_DP
            rewind (phase_unit)
            call eh_read_phase_bin(phase_unit)
            close (phase_unit)
            call eh_select_integrator(trim(methods(method)))
            do k = 1, 100
               call eh_integrate_wf(en, en, coupling, coupling, step)
            end do
            call read_coefficients(actual_c)
            stored_c = actual_c
            if (mode == 1) actual_c = actual_c * exp(-imaginary_unit * en * elapsed)
            tolerance = 1.E-11_DP
            if (methods(method) == 'euler') tolerance = 1.E-3_DP
            if (maxval(abs(actual_c - expected_c)) > tolerance) error stop 'Complex electronic amplitudes'
            if (abs(get_coherence(1, 2) - expected_coherence) > tolerance) error stop 'Complex coherence'

            ! Preserve the extended precision through the complex coefficient file format.
            open (newunit=wf_unit, status='scratch', action='readwrite')
            call eh_write_wf(wf_unit)
            rewind (wf_unit)
            call eh_read_wf(wf_unit)
            close (wf_unit)
            call read_coefficients(actual_c)
            if (maxval(abs(actual_c - stored_c)) > 1.E-30_DP) error stop 'Complex wavefunction I/O'
         end do
      end do
      phase = 0
      print*,'Complex amplitudes, coherence and I/O passed for all integrators and both phase modes.'
   end subroutine check_electronic_integrators

   subroutine read_coefficients(coefficients)
      complex(2*DP), intent(out) :: coefficients(:)
      integer :: wf_unit, state
      open (newunit=wf_unit, status='scratch', action='readwrite')
      call eh_write_wf(wf_unit)
      rewind (wf_unit)
      do state = 1, size(coefficients)
         read (wf_unit, *) coefficients(state)
      end do
      close (wf_unit)
   end subroutine read_coefficients

   subroutine set_wf(c1, c2, imaginary)
      real(DP), intent(in) :: c1, c2
      logical, intent(in), optional :: imaginary
      real(DP) :: re2, im2
      re2 = c2; im2 = 0.0_DP
      if (present(imaginary)) then
         if (imaginary) then
            re2 = 0.0_DP; im2 = c2
         end if
      end if
      open (newunit=u, status='scratch', action='readwrite')
      write (u, *) cmplx(c1, 0.0_DP, kind=2 * DP)
      write (u, *) cmplx(re2, im2, kind=2 * DP)
      rewind (u)
      call eh_read_wf(u)
      close (u)
   end subroutine set_wf
end program ehrenfest_physics
