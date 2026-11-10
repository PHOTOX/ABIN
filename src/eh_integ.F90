! Electronic wavefunction propagation for Ehrenfest dynamics.
! Adapted from sh_integ.F90 with independent wavefunction and integrator state.
! Butcher (5th order), RK4 and Euler use the same numerical schemes as SH.
module mod_eh_integ
   use mod_const, only: DP, imaginary_unit
   use mod_error, only: fatal_error
   implicit none
   private
   public :: eh_set_initialwf, eh_set_energy_shift
   public :: eh_select_integrator, eh_integrate_wf
   public :: check_popsum

   ! Restart functionality
   public :: eh_write_phase_bin, eh_read_phase_bin
   public :: eh_write_wf, eh_read_wf
   public :: get_coherence, get_elpop

   ! Input parameters
   public :: phase
   public :: nstate, popsumthr

   ! Complex amplitudes
   complex(2*DP), allocatable :: cel(:)
   ! phase for phase=1, i.e. the phase of the wavefunction in the adiabatic representation
   real(DP), allocatable, dimension(:, :) :: gama

   real(DP) :: Eshift, popsumthr = 0.001D0
   ! TODO: phase variable should be named differently (as in sh_integ.F90) to avoid confusion with the phase of the wavefunction (gama).
   integer :: phase = 0
   ! Number of electronic states
   integer :: nstate = 1

   ! The Ehrenfest integrator is selected at runtime and stored in a procedure pointer.
   abstract interface
      subroutine eh_integrator(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
         import :: DP
         implicit none
         ! Interpolated potential energies
         real(DP), intent(in), dimension(:) :: en_array_int, en_array_newint
         ! Interpolated dotproduct = nacm x velocity
         real(DP), intent(in), dimension(:, :) :: dotproduct_int, dotproduct_newint
         ! Electronic timestep
         real(DP), intent(in) :: dtp
      end subroutine eh_integrator
   end interface

   procedure(eh_integrator), pointer :: eh_integrate_wf => null()

contains

   subroutine eh_select_integrator(integrator)
      character(len=*), intent(in) :: integrator

      if (integrator == 'butcher') then
         eh_integrate_wf => butcherstep
      else if (integrator == 'rk4') then
         eh_integrate_wf => rk4step
      else if (integrator == 'euler') then
         eh_integrate_wf => eulerstep
      else
         call fatal_error(__FILE__, __LINE__, &
            & 'Invalid Ehrenfest integrator "'//trim(integrator)//'"')
      end if
   end subroutine eh_select_integrator

   ! Calculate electronic population of a given state
   real(DP) function get_elpop(state) result(el_pop)
      integer, intent(in) :: state
      el_pop = real(conjg(cel(state)) * cel(state))
   end function get_elpop

   ! Electronic coherence: real part of conjg(c_i)*c_j in the adiabatic representation.
   real(DP) function get_coherence(i, j) result(coherence)
      integer, intent(in) :: i, j
      complex(2*DP) :: density_element

      density_element = conjg(cel(i)) * cel(j)

      !contribution from the phase of the wavefunction in the adiabatic representation
      if (phase == 1) density_element = density_element * exp(imaginary_unit * gama(i, j))

      coherence = real(density_element)
   end function get_coherence

   ! Integrate the Ehrenfest wavefunction over one step
   subroutine integstep(k, en, y, dotprod, gam, dtp)
      complex(2*DP), intent(out) :: k(:)
      complex(2*DP), intent(in) :: y(:)
      real(DP), intent(in) :: en(:), dotprod(:, :), gam(:, :), dtp
      integer :: ist1, ist2

      ! contribution to the change from the state it self
      k = -imaginary_unit * en * y
      if (phase == 0) then
         do ist1 = 1, nstate
            do ist2 = 1, nstate
               ! contribution from the other states
               if (ist1 /= ist2) then
                  k(ist1) = k(ist1) - y(ist2) * dotprod(ist1, ist2)
               end if
            end do
         end do
         ! same as above but for the flow of the code done separately
         ! ifs are heavy
      else if (phase == 1) then
         do ist1 = 1, nstate
            do ist2 = 1, nstate
               if (ist1 /= ist2) then
                  k(ist1) = k(ist1) - y(ist2) * exp(imaginary_unit * gam(ist1, ist2)) * dotprod(ist1, ist2)
               end if
            end do
         end do
      end if

      k = dtp * k
   end subroutine integstep

   ! Integrate the phase (gama) of the Ehrenfest wavefunction over one step
   subroutine integ_gama(en_array_int, en_array_newint, dtp)
      real(DP), intent(in) :: en_array_int(:)
      real(DP), intent(in) :: en_array_newint(:)
      real(DP), intent(in) :: dtp
      real(DP) :: g
      integer :: ist1, ist2

      do ist1 = 1, nstate
         do ist2 = 1, nstate
            if (ist1 /= ist2) then
               g = (en_array_int(ist1) + en_array_newint(ist1)) / 2
               g = g - (en_array_int(ist2) + en_array_newint(ist2)) / 2

               gama(ist1, ist2) = gama(ist1, ist2) + g * dtp
            end if
         end do
      end do

   end subroutine integ_gama

   ! euler integrator for the Ehrenfest (not ideal but it is an option)
   subroutine eulerstep(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in) :: en_array_int(:), en_array_newint(:)
      real(DP), intent(in) :: dotproduct_int(:, :), dotproduct_newint(:, :), dtp
      complex(2*DP) :: k1(nstate)

      call integstep(k1, en_array_int + eshift, cel, dotproduct_int, gama, dtp)
      cel = cel + k1
      call integ_gama(en_array_int, en_array_newint, dtp)
   end subroutine eulerstep

   ! rk4 integrator for the Ehrenfest
   subroutine rk4step(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in) :: en_array_int(:), en_array_newint(:)
      real(DP), intent(in) :: dotproduct_int(:, :), dotproduct_newint(:, :), dtp
      real(DP) :: dotprod2(nstate, nstate), en0(nstate), en1(nstate), en2(nstate)
      real(DP) :: gam0(nstate, nstate), gam1(nstate, nstate), gam2(nstate, nstate)
      complex(2*DP) :: k1(nstate), k2(nstate), k3(nstate), k4(nstate)

      en0 = en_array_int + eshift
      en1 = en_array_newint + eshift
      en2 = (en0 + en1) / 2
      dotprod2 = (dotproduct_int + dotproduct_newint) / 2
      gam0 = gama
      if (phase == 1) call integ_gama(en_array_int, en_array_newint, dtp)
      gam1 = gama
      gam2 = (gam0 + gam1) / 2

      call integstep(k1, en0, cel, dotproduct_int, gam0, dtp)
      call integstep(k2, en2, cel + k1 / 2, dotprod2, gam2, dtp)
      call integstep(k3, en2, cel + k2 / 2, dotprod2, gam2, dtp)
      call integstep(k4, en1, cel + k3, dotproduct_newint, gam1, dtp)
      cel = cel + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
   end subroutine rk4step

   ! butcher integrator for the Ehrenfest
   subroutine butcherstep(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in) :: en_array_int(:), en_array_newint(:)
      real(DP), intent(in) :: dotproduct_int(:, :), dotproduct_newint(:, :), dtp
      real(DP) :: dotprod2(nstate, nstate), dotprod4(nstate, nstate), dotprod34(nstate, nstate)
      real(DP) :: en0(nstate), en1(nstate), en2(nstate), en4(nstate), en34(nstate)
      real(DP) :: gam0(nstate, nstate), gam1(nstate, nstate), gam2(nstate, nstate)
      real(DP) :: gam4(nstate, nstate), gam34(nstate, nstate)
      complex(2*DP) :: k1(nstate), k2(nstate), k3(nstate), k4(nstate), k5(nstate), k6(nstate)

      en0 = en_array_int + eshift
      en1 = en_array_newint + eshift
      en2 = (en0 + en1) / 2
      en4 = (en0 + en2) / 2
      en34 = (en2 + en1) / 2
      dotprod2 = (dotproduct_int + dotproduct_newint) / 2
      dotprod4 = (dotproduct_int + dotprod2) / 2
      dotprod34 = (dotprod2 + dotproduct_newint) / 2
      gam0 = gama
      if (phase == 1) call integ_gama(en_array_int, en_array_newint, dtp)
      gam1 = gama
      gam2 = (gam0 + gam1) / 2
      gam4 = (gam0 + gam2) / 2
      gam34 = (gam2 + gam1) / 2

      call integstep(k1, en0, cel, dotproduct_int, gam0, dtp)
      call integstep(k2, en4, cel + k1 / 4, dotprod4, gam4, dtp)
      call integstep(k3, en4, cel + k1 / 8 + k2 / 8, dotprod4, gam4, dtp)
      call integstep(k4, en2, cel - k2 / 2 + k3, dotprod2, gam2, dtp)
      call integstep(k5, en34, cel + 3 * k1 / 16 + 9 * k4 / 16, dotprod34, gam34, dtp)
      call integstep(k6, en1, cel - 3 * k1 / 7 + 2 * k2 / 7 + 12 * k3 / 7 - 12 * k4 / 7 + 8 * k5 / 7, dotproduct_newint, gam1, dtp)
      cel = cel + 7 * k1 / 90 + 32 * k3 / 90 + 12 * k4 / 90 + 32 * k5 / 90 + 7 * k6 / 90
   end subroutine butcherstep

   subroutine eh_write_phase_bin(iunit)
      integer, intent(in) :: iunit
      integer :: ist1, ist2
      do ist1 = 1, nstate
         write (iunit) (gama(ist1, ist2), ist2=1, nstate)
      end do
   end subroutine eh_write_phase_bin

   subroutine eh_read_phase_bin(iunit)
      integer, intent(in) :: iunit
      integer :: ist1, ist2

      if (.not. allocated(gama)) allocate (gama(nstate, nstate))

      do ist1 = 1, nstate
         read (iunit) (gama(ist1, ist2), ist2=1, nstate)
      end do
   end subroutine eh_read_phase_bin

   !TODO: expand to have multiple initial states
   subroutine eh_set_initialwf(initial_state)
      integer, intent(in) :: initial_state

      allocate (gama(nstate, nstate))
      gama = 0.0_DP
      allocate (cel(nstate))
      cel = (0.0_DP, 0.0_DP)
      cel(initial_state) = (1.0_DP, 0.0_DP)
   end subroutine eh_set_initialwf

   subroutine eh_set_energy_shift(potential_energy)
      real(DP), intent(in) :: potential_energy
      Eshift = -potential_energy
   end subroutine eh_set_energy_shift

   subroutine eh_write_wf(outunit)
      integer, intent(in) :: outunit
      integer :: ist1
      do ist1 = 1, nstate
         write (outunit, *) cel(ist1)
      end do
   end subroutine eh_write_wf

   subroutine eh_read_wf(inunit)
      integer, intent(in) :: inunit
      integer :: ist1
      complex(2*DP) :: coefficient

      if (.not. allocated(cel)) allocate (cel(nstate))
      ! One complex coefficient (real, imaginary) per state.
      do ist1 = 1, nstate
         read (inunit, *) coefficient
         cel(ist1) = coefficient
      end do
   end subroutine eh_read_wf

   !better to check sometimes
   real(DP) function check_popsum()
      use mod_files, only: stderr
      real(DP) :: popsum
      integer :: ist1

      popsum = 0.0D0
      do ist1 = 1, nstate
         popsum = popsum + get_elpop(ist1)
      end do

      if (abs(popsum - 1.0D0) > popsumthr) then
         write (stderr, *) 'ERROR: Sum of electronic populations = ', popsum
         write (stderr, *) 'which differs from 1.0 by more than popsumthr = ', popsumthr
         write (stderr, *) 'Increase the number of substeps or use more accurate integrator.'
         call fatal_error(__FILE__, __LINE__, 'Inaccurate Ehrenfest integration led to unphysical populations')
      end if

      check_popsum = popsum
   end function check_popsum

end module mod_eh_integ
