! Check high-state transfer and compare electronic propagation on identical data.
program ehrenfest_multistate
   use mod_const, only: DP
   use mod_eh_integ, only: nstate, phase, eh_set_initialwf, eh_set_energy_shift, &
                           eh_select_integrator, eh_integrate_wf, eh_read_wf, eh_read_phase_bin, get_elpop, check_popsum
   use mod_sh_integ, only: sh_nstate => nstate, sh_phase => phase, sh_set_initialwf, sh_set_energy_shift, &
                           sh_select_integrator, sh_integrate_wf, sh_read_wf, sh_read_phase_bin, sh_get_elpop => get_elpop
   implicit none
   integer, parameter :: ns = 5, steps = 2000
   character(len=7), parameter :: methods(3) = [character(len=7) :: 'butcher', 'rk4', 'euler']
   real(DP), parameter :: h = 0.01_DP, offset = -78.0_DP
   real(DP) :: en0(ns), en1(ns), d0(ns, ns), d1(ns, ns), expected, frequency, gap, tolerance, norm
   real(DP) :: parity_error, max_error, max_transfer
   integer :: method, mode, initial, k, i

   nstate = ns
   sh_nstate = ns
   call eh_set_initialwf(ns)
   call sh_set_initialwf(ns)
   call eh_set_energy_shift(offset)
   call sh_set_energy_shift(offset)
   max_error = 0.0_DP
   parity_error = 0.0_DP

   do mode = 0, 1
      phase = mode
      sh_phase = mode
      do method = 1, size(methods)
         call eh_select_integrator(trim(methods(method)))
         call sh_select_integrator(trim(methods(method)))
         ! Cover both ABIN state 4 -> 3 and spectroscopic S4 -> S3 (5 -> 4).
         do initial = 4, 5
            call reset(initial)
            en0 = offset + [0.0_DP, 0.01_DP, 0.02_DP, 0.03_DP, 0.04_DP]
            d0 = 0.0_DP
            d0(initial - 1, initial) = 0.04_DP
            d0(initial, initial - 1) = -0.04_DP
            do k = 1, steps
               call eh_integrate_wf(en0, en0, d0, d0, h)
               call sh_integrate_wf(en0, en0, d0, d0, h)
            end do
            gap = (en0(initial) - en0(initial - 1)) / 2
            frequency = sqrt(gap**2 + 0.04_DP**2)
            expected = (0.04_DP / frequency * sin(frequency * steps * h))**2
            tolerance = 1.E-10_DP
            if (method == 3) tolerance = 1.E-3_DP
            max_error = max(max_error, abs(get_elpop(initial - 1) - expected))
            if (abs(get_elpop(initial - 1) - expected) > tolerance) error stop 'High-state transfer'
            if (abs(get_elpop(initial) - (1 - expected)) > tolerance) error stop 'Initial-state depletion'
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
            max_transfer = max(max_transfer, 1.0_DP - get_elpop(5))
            en0 = en1
            d0 = d1
         end do
         if (max_transfer < 0.01_DP) error stop 'Multistate propagation did not transfer population'
      end do
   end do
   print*,'Five-state transfer passed for all integrators and both phase modes.'
   print*,'Maximum population difference between EH and SH:', parity_error
   print*,'Maximum analytic error (including Euler):', max_error

contains
   subroutine reset(state)
      integer, intent(in) :: state
      integer :: unit, j
      real(DP) :: amplitude
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
