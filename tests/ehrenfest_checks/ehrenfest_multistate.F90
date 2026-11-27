! Standalone electronic propagation check with prescribed energies and couplings.
! Compare high-state transfer with the analytic two-level result, then compare
! all five populations with SH and complex amplitudes with a power-series reference.
program ehrenfest_multistate
   use mod_const, only: DP, imaginary_unit
   use mod_eh_integ, only: nstate, phase, eh_set_initialwf, eh_set_energy_shift, &
                           eh_select_integrator, eh_integrate_wf, eh_read_wf, eh_write_wf, &
                           eh_read_phase_bin, get_elpop, check_popsum
   use mod_sh_integ, only: sh_nstate => nstate, sh_phase => phase, sh_set_initialwf, sh_set_energy_shift, &
                           sh_select_integrator, sh_integrate_wf, sh_read_wf, sh_read_phase_bin, sh_get_elpop => get_elpop
   implicit none
   integer, parameter :: ns = 5, steps = 2000, checkpoints = 4
   integer, parameter :: sample_every = steps / checkpoints
   character(len=7), parameter :: methods(3) = [character(len=7) :: 'butcher', 'rk4', 'euler']
   real(DP), parameter :: h = 0.01_DP, offset = -78.0_DP
   real(DP) :: en0(ns), en1(ns), d0(ns, ns), d1(ns, ns), expected, frequency, gap, tolerance, norm
   complex(2*DP) :: reference(ns, checkpoints)
   real(DP) :: diagonal_phase(ns, checkpoints), amplitude_error
   real(DP) :: parity_error, max_error, max_transfer
   integer :: method, mode, initial, k, i

   ! Configure both propagators with the same number of states and energy shift.
   nstate = ns
   sh_nstate = ns
   call eh_set_initialwf(ns)
   call sh_set_initialwf(ns)
   call eh_set_energy_shift(offset)
   call sh_set_energy_shift(offset)
   max_error = 0.0_DP
   parity_error = 0.0_DP
   amplitude_error = 0.0_DP
   call build_reference(reference, diagonal_phase)

   ! Repeat both scenarios with every integrator and phase convention.
   do mode = 0, 1
      phase = mode
      sh_phase = mode
      do method = 1, size(methods)
         call eh_select_integrator(trim(methods(method)))
         call sh_select_integrator(trim(methods(method)))
         ! Cover both ABIN state 4 -> 3 and spectroscopic S4 -> S3 (5 -> 4).
         do initial = 4, 5
            call reset(initial)
            ! Couple only the initial state and its lower neighbor to isolate a two-level system.
            en0 = offset + [0.0_DP, 0.01_DP, 0.02_DP, 0.03_DP, 0.04_DP]
            d0 = 0.0_DP
            d0(initial - 1, initial) = 0.04_DP
            d0(initial, initial - 1) = -0.04_DP
            do k = 1, steps
               call eh_integrate_wf(en0, en0, d0, d0, h)
               call sh_integrate_wf(en0, en0, d0, d0, h)
            end do
            ! The constant two-level generator gives an exact transfer probability.
            ! Allow a wider tolerance for the first-order Euler integrator.
            gap = (en0(initial) - en0(initial - 1)) / 2
            frequency = sqrt(gap**2 + 0.04_DP**2)
            expected = (0.04_DP / frequency * sin(frequency * steps * h))**2
            tolerance = 1.E-10_DP
            if (method == 3) tolerance = 1.E-3_DP
            max_error = max(max_error, abs(get_elpop(initial - 1) - expected))
            if (abs(get_elpop(initial - 1) - expected) > tolerance) error stop 'High-state transfer'
            if (abs(get_elpop(initial) - (1 - expected)) > tolerance) error stop 'Initial-state depletion'
            ! Uncoupled spectator states must remain empty and the total norm must stay near one.
            do i = 1, ns
               if (i /= initial .and. i /= initial - 1) then
                  if (get_elpop(i) /= 0.0_DP) error stop 'Spectator-state leakage'
               end if
            end do
            norm = check_popsum()
            if (abs(norm - 1) > tolerance) error stop 'High-state norm'
            call compare_populations()
         end do

         ! All pairs nonzero, varying independently, to exercise multistate interference.
         call reset(5)
         max_transfer = 0.0_DP
         call endpoint(0.0_DP, en0, d0)
         do k = 1, steps
            call endpoint(k * h, en1, d1)
            call eh_integrate_wf(en0, en1, d0, d1, h)
            call sh_integrate_wf(en0, en1, d0, d1, h)
            call compare_populations()
            if (mod(k, sample_every) == 0) call compare_amplitudes(k / sample_every)
            max_transfer = max(max_transfer, 1.0_DP - get_elpop(5))
            ! Use this endpoint as the start of the next interpolation interval.
            en0 = en1
            d0 = d1
         end do
         ! Require measurable transfer to catch a wavefunction that never evolves.
         if (max_transfer < 0.01_DP) error stop 'Multistate propagation did not transfer population'
      end do
   end do
   print*,'Five-state transfer passed for all integrators and both phase modes.'
   print*,'Maximum population difference between EH and SH:', parity_error
   print*,'Maximum analytic two-state error (including Euler):', max_error
   print*,'Maximum five-state amplitude error (including Euler):', amplitude_error

contains
   ! Seed both propagators with the same pure state and clear accumulated electronic phases.
   subroutine reset(state)
      integer, intent(in) :: state
      integer :: unit, j
      real(DP) :: amplitude
      ! EH reads complex values. SH reads separate real and imaginary columns.
      open (newunit=unit, status='scratch', form='formatted')
      do j = 1, ns
         amplitude = 0.0_DP
         if (j == state) amplitude = 1.0_DP
         write (unit, *) cmplx(amplitude, 0.0_DP, kind=2 * DP)
      end do
      rewind (unit)
      call eh_read_wf(unit)
      close (unit)
      open (newunit=unit, status='scratch', form='formatted')
      do j = 1, ns
         amplitude = 0.0_DP
         if (j == state) amplitude = 1.0_DP
         write (unit, *) amplitude, 0.0_DP
      end do
      rewind (unit)
      call sh_read_wf(unit)
      close (unit)
      ! Both phase readers accept the same binary matrix of pairwise phase differences.
      open (newunit=unit, status='scratch', form='unformatted')
      do j = 1, ns
         write (unit) spread(0.0_DP, 1, ns)
      end do
      rewind (unit)
      call eh_read_phase_bin(unit)
      rewind (unit)
      call sh_read_phase_bin(unit)
      close (unit)
   end subroutine reset

   ! Generate time-varying state energies and antisymmetric couplings at a step endpoint.
   subroutine endpoint(time, en, coupling)
      real(DP), intent(in) :: time
      real(DP), intent(out) :: en(ns), coupling(ns, ns)
      integer :: a, b
      coupling = 0.0_DP
      do a = 1, ns
         en(a) = offset + 0.01_DP * a + 0.003_DP * sin(0.1_DP * a * time)
         do b = a + 1, ns
            coupling(a, b) = 0.002_DP * (a + b) * cos(0.003_DP * (a * b) * time)
            coupling(b, a) = -coupling(a, b)
         end do
      end do
   end subroutine endpoint

   ! Solve in the adiabatic representation, using the same prescribed linear
   ! interpolation over each interval as the production integrators.
   subroutine build_reference(samples, phases)
      complex(2*DP), intent(out) :: samples(ns, checkpoints)
      real(DP), intent(out) :: phases(ns, checkpoints)
      complex(2*DP) :: coefficients(ns), start(ns, ns), finish(ns, ns)
      real(DP) :: en(ns), coupling(ns, ns), previous_en(ns), accumulated_phase(ns)
      integer :: step, sample

      coefficients = (0.0_DP, 0.0_DP)
      coefficients(ns) = (1.0_DP, 0.0_DP)
      accumulated_phase = 0.0_DP
      call endpoint(0.0_DP, en, coupling)
      start = generator(en, coupling)
      do step = 1, steps
         previous_en = en
         call endpoint(step * h, en, coupling)
         finish = generator(en, coupling)
         call series_step(start, finish - start, coefficients)
         accumulated_phase = accumulated_phase + h * ((previous_en - offset) + (en - offset)) / 2
         if (mod(step, sample_every) == 0) then
            sample = step / sample_every
            samples(:, sample) = coefficients
            phases(:, sample) = accumulated_phase
         end if
         start = finish
      end do
   end subroutine build_reference

   ! Combine the energy and coupling terms into a single complex matrix.
   function generator(en, coupling) result(matrix)
      real(DP), intent(in) :: en(ns), coupling(ns, ns)
      complex(2*DP) :: matrix(ns, ns)
      integer :: j

      matrix = cmplx(-coupling, kind=2 * DP)
      do j = 1, ns
         matrix(j, j) = -imaginary_unit * (en(j) - offset)
      end do
   end function generator

   ! For dc/ds = h*(start + s*delta)*c and c(s) = sum(a_m*s**m), the
   ! coefficients satisfy (m+1)*a_(m+1) = h*(start*a_m + delta*a_(m-1)).
   ! At s=1 the reference is the sum of these coefficients.
   subroutine series_step(start, delta, coefficients)
      complex(2*DP), intent(in) :: start(ns, ns), delta(ns, ns)
      complex(2*DP), intent(inout) :: coefficients(ns)
      complex(2*DP) :: previous(ns), term(ns), next(ns)
      integer :: order

      previous = (0.0_DP, 0.0_DP)
      term = coefficients
      do order = 1, 12
         next = (h * (matmul(start, term) + matmul(delta, previous))) / order
         coefficients = coefficients + next
         previous = term
         term = next
      end do
      if (maxval(abs(term)) > 1.E-20_DP) error stop 'Five-state reference series did not converge'
   end subroutine series_step

   ! Undo the diagonal phase transformation before comparing complex amplitudes.
   subroutine compare_amplitudes(sample)
      integer, intent(in) :: sample
      complex(2*DP) :: actual(ns)
      real(DP) :: error, limit
      integer :: unit, j

      open (newunit=unit, status='scratch', form='formatted')
      call eh_write_wf(unit)
      rewind (unit)
      do j = 1, ns
         read (unit, *) actual(j)
      end do
      close (unit)
      if (mode == 1) actual = actual * exp(-imaginary_unit * diagonal_phase(:, sample))
      limit = 1.E-10_DP
      ! Phase-mode stages interpolate the accumulated phases within each step.
      if (mode == 1) limit = 1.E-8_DP
      if (methods(method) == 'euler') limit = 1.E-3_DP
      error = maxval(abs(actual - reference(:, sample)))
      amplitude_error = max(amplitude_error, error)
      if (error > limit) then
         print*,'Integrator: ', trim(methods(method)), ', phase mode: ', mode, ', checkpoint: ', sample
         print*,'Maximum amplitude error: ', error, ', tolerance: ', limit
         error stop 'Five-state complex electronic amplitudes'
      end if
   end subroutine compare_amplitudes

   ! Check every state's EH/SH population and retain the largest difference for reporting.
   subroutine compare_populations()
      integer :: j
      real(DP) :: error
      do j = 1, ns
         error = abs(get_elpop(j) - sh_get_elpop(j))
         parity_error = max(parity_error, error)
         if (error > 1.E-11_DP) error stop 'EH/SH electronic population mismatch'
      end do
   end subroutine compare_populations
end program ehrenfest_multistate
