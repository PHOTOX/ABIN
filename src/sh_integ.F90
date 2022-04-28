! Routines for handling and propagation of electronic WF during Surface Hopping.
! Currently, the only production method is the Butcher 5th order integrator.
! Euler and 4th-order Runge-Kutta methods are only here for debugging purposes.
! This module is the parent of mod_sh which contains the driver SH routine.
module mod_sh_integ
   use mod_const, only: DP
   use mod_error, only: fatal_error
   private
   public :: sh_integrate_wf, sh_set_initialwf, sh_write_phase_bin, sh_read_phase_bin
   public :: sh_decoherence_correction, check_popsum, sh_TFS_transmat
   public :: sh_write_wf, sh_read_wf, sh_debug_wf
   public :: integ, phase, el_pop
   public :: nstate, popsumthr
   public :: correct_decoherence

   ! Electronic State coefficients
   ! TODO: Would be nice if we used the complex type
   real(DP), allocatable, dimension(:) :: cel_re, cel_im
   real(DP), allocatable, dimension(:, :) :: gama
   ! TODO: Get rid of this global, calculate elpop on the fly
   real(DP), allocatable, dimension(:) :: el_pop

   real(DP) :: eshift, popsumthr = 0.001D0
   ! TODO: phase variable should be named differently
   integer :: phase = 0
   ! Number of electronic states
   integer :: nstate = 1
   ! Numerical integrator (default is Butcher 5-th order)
   character(len=50) :: integ = 'butcher'
   ! NOTE: To exactly replicate older data, switch this to .false.
   ! See the comment in subroutine sh_decoherence_correction() for explanation.
   logical :: correct_decoherence = .true.

contains

   subroutine sh_integrate_wf(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in), dimension(:) :: en_array_int, en_array_newint
      real(DP), intent(in), dimension(:, :) :: dotproduct_int, dotproduct_newint
      real(DP), intent(in) :: dtp ! Time step

      if (integ == 'butcher') then

         call butcherstep(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)

      else if (integ == 'rk4') then

         call rk4step(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)

      else if (integ == 'euler') then

         call eulerstep(en_array_int, en_array_newint, dotproduct_int, dtp)

      end if

      call sh_calc_elpop()

   end subroutine sh_integrate_wf

   subroutine sh_calc_elpop()
      integer :: ist1
      do ist1 = 1, nstate
         el_pop(ist1) = cel_re(ist1)**2 + cel_im(ist1)**2
      end do
   end subroutine sh_calc_elpop

   ! Calculates transitions matrix according to the Tully's fewest switches algorithm
   subroutine sh_TFS_transmat(dotproduct_int, dotproduct_newint, ist, pop0, t, dtp)
      ! Interpolated dotproducts between velocities and NA couplings
      real(DP), intent(in) :: dotproduct_int(:, :)
      real(DP), intent(in) :: dotproduct_newint(:, :)
      ! Index of current state
      integer, intent(in) :: ist
      ! Population of the current state
      real(DP), intent(in) :: pop0
      ! Electronic time step
      real(DP), intent(in) :: dtp
      ! OUTPUT: Transition probability matrix
      real(DP), intent(out) :: t(:, :)
      real(DP) :: a_re, a_im
      integer :: ist2

      t = 0.00D0

      do ist2 = 1, nstate
         a_re = (cel_re(ist) * cel_re(ist2) + cel_im(ist) * cel_im(ist2))
         if (phase == 1) then
            a_im = (-cel_im(ist) * cel_re(ist2) + cel_re(ist) * cel_im(ist2))
            t(ist, ist2) = a_re * cos(gama(ist, ist2)) - sin(gama(ist, ist2)) * a_im
            t(ist, ist2) = t(ist, ist2) * (dotproduct_int(ist, ist2) + dotproduct_newint(ist, ist2))
         else
            ! t(ist,ist2)=2*a_re*dotproduct_int(ist,ist2)
            t(ist, ist2) = a_re * (dotproduct_int(ist, ist2) + dotproduct_newint(ist, ist2))
         end if
         t(ist, ist2) = t(ist, ist2) * dtp / (pop0 + 1D-20)
      end do

      do ist2 = 1, nstate
         if (t(ist, ist2) > 1.0_DP) then
            call fatal_error(__FILE__, __LINE__, &
               & 'ERROR: Hopping probability greater than 1.')
         end if
         if (t(ist, ist2) < 0.0D0) t(ist, ist2) = 0.0D0
      end do

   end subroutine sh_TFS_transmat

   subroutine integstep(k_re, k_im, en, y_re, y_im, dotprod, gam, dtp)
      real(DP), intent(out) :: k_re(:), k_im(:)
      real(DP), intent(in) :: dotprod(:, :), gam(:, :)
      real(DP), intent(in) :: en(:), y_im(:), y_re(:), dtp
      real(DP) :: g
      integer :: ist1, ist2

      if (phase == 0) then
         do ist1 = 1, nstate
            k_re(ist1) = en(ist1) * y_im(ist1)
            k_im(ist1) = -en(ist1) * y_re(ist1)
            do ist2 = 1, nstate
               if (ist1 /= ist2) then
                  k_re(ist1) = k_re(ist1) - y_re(ist2) * dotprod(ist1, ist2)
                  k_im(ist1) = k_im(ist1) - y_im(ist2) * dotprod(ist1, ist2)
               end if
            end do
            k_re(ist1) = dtp * k_re(ist1)
            k_im(ist1) = dtp * k_im(ist1)
         end do
      end if

      if (phase == 1) then
         do ist1 = 1, nstate
            k_re(ist1) = 0.0_DP
            k_im(ist1) = 0.0_DP
            do ist2 = 1, nstate
               if (ist1 /= ist2) then
                  g = gam(ist1, ist2)
                  k_re(ist1) = k_re(ist1) - (y_re(ist2) * cos(g) - y_im(ist2) * sin(g)) * dotprod(ist1, ist2)
                  k_im(ist1) = k_im(ist1) - (y_im(ist2) * cos(g) + y_re(ist2) * sin(g)) * dotprod(ist1, ist2)
               end if
            end do
            k_re(ist1) = dtp * k_re(ist1)
            k_im(ist1) = dtp * k_im(ist1)
         end do
      end if
   end subroutine integstep

   subroutine integ_gama(en_array_int, en_array_newint, dtp)
      real(DP), intent(in) :: en_array_int(:)
      real(DP), intent(in) :: en_array_newint(:), dtp
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

   subroutine eulerstep(en_array_int, en_array_newint, dotproduct_int, dtp)
      real(DP), intent(in), dimension(:) :: en_array_int, en_array_newint
      real(DP), intent(in), dimension(:, :) :: dotproduct_int
      real(DP) :: dtp
      real(DP) :: dotprod0(nstate, nstate), gam0(nstate, nstate)
      real(DP) :: k1_re(nstate), k1_im(nstate)
      real(DP) :: y_im(nstate), y_re(nstate)
      real(DP) :: en0(nstate)
      integer :: ist1, ist2

      do ist1 = 1, nstate
         en0(ist1) = en_array_int(ist1) + eshift
         y_re(ist1) = cel_re(ist1)
         y_im(ist1) = cel_im(ist1)
         do ist2 = 1, nstate
            dotprod0(ist1, ist2) = dotproduct_int(ist1, ist2)
            gam0(ist1, ist2) = gama(ist1, ist2)
         end do
      end do

      call integstep(k1_re, k1_im, en0, y_re, y_im, dotprod0, gam0, dtp)

      do ist1 = 1, nstate
         cel_re(ist1) = cel_re(ist1) + k1_re(ist1)
         cel_im(ist1) = cel_im(ist1) + k1_im(ist1)
      end do

      call integ_gama(en_array_int, en_array_newint, dtp)
   end subroutine eulerstep

   subroutine rk4step(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in) :: en_array_int(:)
      real(DP), intent(in) :: en_array_newint(:)
      real(DP), intent(in) :: dotproduct_int(:, :)
      real(DP), intent(in) :: dotproduct_newint(:, :), dtp
      real(DP) :: dotprod2(nstate, nstate)
      real(DP) :: dotprod0(nstate, nstate), dotprod1(nstate, nstate)
      real(DP) :: k1_re(nstate), k1_im(nstate)
      real(DP) :: k2_re(nstate), k2_im(nstate)
      real(DP) :: k3_re(nstate), k3_im(nstate)
      real(DP) :: k4_re(nstate), k4_im(nstate)
      real(DP) :: y_im(nstate), y_re(nstate)
      real(DP) :: en0(nstate), en1(nstate), en2(nstate)
      real(DP) :: gam0(nstate, nstate), gam1(nstate, nstate), gam2(nstate, nstate)
      integer :: ist1, ist2

!     initial interpolations
      do ist1 = 1, nstate
         en0(ist1) = en_array_int(ist1) + eshift
         en1(ist1) = en_array_newint(ist1) + eshift
         en2(ist1) = (en0(ist1) + en1(ist1)) / 2
         y_re(ist1) = cel_re(ist1)
         y_im(ist1) = cel_im(ist1)
         do ist2 = 1, nstate
            dotprod0(ist1, ist2) = dotproduct_int(ist1, ist2)
            dotprod1(ist1, ist2) = dotproduct_newint(ist1, ist2)
            dotprod2(ist1, ist2) = (dotprod0(ist1, ist2) + dotprod1(ist1, ist2)) / 2
            if (phase == 1) gam0(ist1, ist2) = gama(ist1, ist2)
         end do
      end do

      if (phase == 1) then
         call integ_gama(en_array_int, en_array_newint, dtp)
         do ist1 = 1, nstate
            do ist2 = 1, nstate
               gam1(ist1, ist2) = gama(ist1, ist2)
               gam2(ist1, ist2) = (gam0(ist1, ist2) + gam1(ist1, ist2)) / 2
            end do
         end do
      end if

      call integstep(k1_re, k1_im, en0, y_re, y_im, dotprod0, gam0, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + k1_re(ist1) / 2
         y_im(ist1) = cel_im(ist1) + k1_im(ist1) / 2
      end do
      call integstep(k2_re, k2_im, en2, y_re, y_im, dotprod2, gam2, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + k2_re(ist1) / 2
         y_im(ist1) = cel_im(ist1) + k2_im(ist1) / 2
      end do
      call integstep(k3_re, k3_im, en2, y_re, y_im, dotprod2, gam2, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + k3_re(ist1)
         y_im(ist1) = cel_im(ist1) + k3_im(ist1)
      end do
      call integstep(k4_re, k4_im, en1, y_re, y_im, dotprod1, gam1, dtp)

      do ist1 = 1, nstate
         cel_re(ist1) = cel_re(ist1) + k1_re(ist1) / 6 + k2_re(ist1) / 3 + k3_re(ist1) / 3 + k4_re(ist1) / 6
         cel_im(ist1) = cel_im(ist1) + k1_im(ist1) / 6 + k2_im(ist1) / 3 + k3_im(ist1) / 3 + k4_im(ist1) / 6
      end do

   end subroutine rk4step

   subroutine butcherstep(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)
      real(DP), intent(in), dimension(:) :: en_array_int, en_array_newint
      real(DP), intent(in), dimension(:, :) :: dotproduct_int, dotproduct_newint
      real(DP), intent(in) :: dtp
      real(DP) :: dotprod2(nstate, nstate), dotprod4(nstate, nstate), dotprod34(nstate, nstate)
      real(DP) :: dotprod0(nstate, nstate), dotprod1(nstate, nstate)
      real(DP) :: k1_re(nstate), k1_im(nstate)
      real(DP) :: k2_re(nstate), k2_im(nstate)
      real(DP) :: k3_re(nstate), k3_im(nstate)
      real(DP) :: k4_re(nstate), k4_im(nstate)
      real(DP) :: k5_re(nstate), k5_im(nstate)
      real(DP) :: k6_re(nstate), k6_im(nstate)
      real(DP) :: y_im(nstate), y_re(nstate)
      real(DP) :: en0(nstate), en1(nstate), en2(nstate), en4(nstate), en34(nstate)
      real(DP) :: gam0(nstate, nstate), gam1(nstate, nstate), gam2(nstate, nstate)
      real(DP) :: gam4(nstate, nstate), gam34(nstate, nstate)
      integer :: ist1, ist2

!     initial interpolations
      do ist1 = 1, nstate
         en0(ist1) = en_array_int(ist1) + eshift
         en1(ist1) = en_array_newint(ist1) + eshift
         en2(ist1) = (en0(ist1) + en1(ist1)) / 2
         en4(ist1) = (en0(ist1) + en2(ist1)) / 2
         en34(ist1) = (en2(ist1) + en1(ist1)) / 2
         y_re(ist1) = cel_re(ist1)
         y_im(ist1) = cel_im(ist1)
         do ist2 = 1, nstate
            dotprod0(ist1, ist2) = dotproduct_int(ist1, ist2)
            dotprod1(ist1, ist2) = dotproduct_newint(ist1, ist2)
            dotprod2(ist1, ist2) = (dotprod0(ist1, ist2) + dotprod1(ist1, ist2)) / 2
            dotprod4(ist1, ist2) = (dotprod0(ist1, ist2) + dotprod2(ist1, ist2)) / 2
            dotprod34(ist1, ist2) = (dotprod2(ist1, ist2) + dotprod1(ist1, ist2)) / 2
            if (phase == 1) gam0(ist1, ist2) = gama(ist1, ist2)
         end do
      end do

      if (phase == 1) then
         call integ_gama(en_array_int, en_array_newint, dtp)
         do ist1 = 1, nstate
            do ist2 = 1, nstate
               gam1(ist1, ist2) = gama(ist1, ist2)
               gam2(ist1, ist2) = (gam0(ist1, ist2) + gam1(ist1, ist2)) / 2
               gam4(ist1, ist2) = (gam0(ist1, ist2) + gam2(ist1, ist2)) / 2
               gam34(ist1, ist2) = (gam2(ist1, ist2) + gam1(ist1, ist2)) / 2
            end do
         end do
      end if

      call integstep(k1_re, k1_im, en0, y_re, y_im, dotprod0, gam0, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + k1_re(ist1) / 4
         y_im(ist1) = cel_im(ist1) + k1_im(ist1) / 4
      end do
      call integstep(k2_re, k2_im, en4, y_re, y_im, dotprod4, gam4, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + k1_re(ist1) / 8 + k2_re(ist1) / 8
         y_im(ist1) = cel_im(ist1) + k1_im(ist1) / 8 + k2_im(ist1) / 8
      end do
      call integstep(k3_re, k3_im, en4, y_re, y_im, dotprod4, gam4, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) - k2_re(ist1) / 2 + k3_re(ist1)
         y_im(ist1) = cel_im(ist1) - k2_im(ist1) / 2 + k3_im(ist1)
      end do
      call integstep(k4_re, k4_im, en2, y_re, y_im, dotprod2, gam2, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) + 3 * k1_re(ist1) / 16 + 9 * k4_re(ist1) / 16
         y_im(ist1) = cel_im(ist1) + 3 * k1_im(ist1) / 16 + 9 * k4_im(ist1) / 16
      end do
      call integstep(k5_re, k5_im, en34, y_re, y_im, dotprod34, gam34, dtp)

      do ist1 = 1, nstate
         y_re(ist1) = cel_re(ist1) - 3 * k1_re(ist1) / 7 + 2 * k2_re(ist1) / 7 + 12 * k3_re(ist1) / 7 &
                      - 12 * k4_re(ist1) / 7 + 8 * k5_re(ist1) / 7
         y_im(ist1) = cel_im(ist1) - 3 * k1_im(ist1) / 7 + 2 * k2_im(ist1) / 7 + 12 * k3_im(ist1) / 7 &
                      - 12 * k4_im(ist1) / 7 + 8 * k5_im(ist1) / 7
      end do
      call integstep(k6_re, k6_im, en1, y_re, y_im, dotprod1, gam1, dtp)

      do ist1 = 1, nstate
         cel_re(ist1) = cel_re(ist1) + &
                            & 7 * k1_re(ist1) / 90 + 32 * k3_re(ist1) / 90 + 12 * k4_re(ist1) / 90 + &
                            & 32 * k5_re(ist1) / 90 + 7 * k6_re(ist1) / 90
         cel_im(ist1) = cel_im(ist1) + &
                            & 7 * k1_im(ist1) / 90 + 32 * k3_im(ist1) / 90 + 12 * k4_im(ist1) / 90 + &
                            & 32 * k5_im(ist1) / 90 + 7 * k6_im(ist1) / 90
      end do

   end subroutine butcherstep

   ! Simple decoherence correction per Grannuci, Persico (2007)
   ! "Critical appraisal of the fewest switches algorithm for surface hopping"
   ! https://doi.org/10.1063/1.2715585
   subroutine sh_decoherence_correction(potential_energies, alpha, kinetic_energy, current_state, dtp)
      real(DP), intent(in) :: potential_energies(:)
      real(DP), intent(in) :: alpha, kinetic_energy, dtp
      integer, intent(in) :: current_state
      integer :: ist1
      real(DP) :: delta_e, tau, scaling_factor, renormalization_factor, sum_norm

      do ist1 = 1, nstate
         if (ist1 == current_state) cycle

         delta_e = abs(potential_energies(ist1) - potential_energies(current_state))
         tau = (1.0D0 + alpha / kinetic_energy) / delta_e
         scaling_factor = dexp(-dtp / tau)

         ! WARNING: The Eq. 17 in the Persico paper is wrong,
         ! the scaling factor should be damping populations, NOT the WF coefficients.
         ! In previous ABIN versions (and in some other SH packages), this was incorrect.
         ! To reproduce older data, the user can set correct_decoherence=.false. in namelist sh.
         if (correct_decoherence) then
            scaling_factor = dsqrt(scaling_factor)
         end if

         cel_re(ist1) = cel_re(ist1) * scaling_factor
         cel_im(ist1) = cel_im(ist1) * scaling_factor
      end do

      ! Renormalize the current state
      sum_norm = 1.0D0
      do ist1 = 1, nstate
         if (ist1 /= current_state) then
            sum_norm = sum_norm - cel_re(ist1)**2 - cel_im(ist1)**2
         end if
      end do

      renormalization_factor = sum_norm / (cel_re(current_state)**2 + cel_im(current_state)**2 + 1.0D-7)

      ! Following should never happen as we check for popsumthr later in this subroutine
      if (renormalization_factor < 0.0D0) then
         write (*, *) 'renormalization_factor=', renormalization_factor, ' but should be > 0'
         write (*, *) 'This usually means inaccurate integration of electronic SE.'
         write (*, *) 'Increase number of substeps or use more accurate integrator.'
         call fatal_error(__FILE__, __LINE__, 'Invalid decoherence renormalization.')
      end if

      renormalization_factor = dsqrt(renormalization_factor)

      cel_re(current_state) = cel_re(current_state) * renormalization_factor
      cel_im(current_state) = cel_im(current_state) * renormalization_factor

      call sh_calc_elpop()

   end subroutine sh_decoherence_correction

   subroutine sh_write_phase_bin(iunit)
      integer, intent(in) :: iunit
      integer :: ist1, ist2
      do ist1 = 1, nstate
         write (iunit) (gama(ist1, ist2), ist2=1, nstate)
      end do
   end subroutine sh_write_phase_bin

   subroutine sh_read_phase_bin(iunit)
      integer, intent(in) :: iunit
      integer :: ist1, ist2
      do ist1 = 1, nstate
         read (iunit) (gama(ist1, ist2), ist2=1, nstate)
      end do
   end subroutine sh_read_phase_bin

   subroutine sh_set_initialwf(initial_state, initial_poten, irest)
      real(DP), intent(in) :: initial_poten
      integer, intent(in) :: initial_state
      ! Did we restart or not?
      integer, intent(in) :: irest
      integer :: ist1

      allocate (gama(nstate, nstate))
      allocate (el_pop(nstate))

      gama = 0.0D0

      if (irest == 0) then
         allocate (cel_re(nstate), cel_im(nstate))
         do ist1 = 1, nstate
            cel_re(ist1) = 0.0D0
            cel_im(ist1) = 0.0D0
         end do

         cel_re(initial_state) = 1.0D0
      end if

      ! Eshift hard set as the potential energy of the ground state
      ! Probably should be handled better (and differently for different trajs)
      Eshift = -initial_poten

      call sh_calc_elpop()

   end subroutine sh_set_initialwf

   subroutine sh_write_wf(outunit)
      integer, intent(in) :: outunit
      integer :: ist1
      do ist1 = 1, nstate
         write (outunit, *) cel_re(ist1), cel_im(ist1)
      end do
   end subroutine sh_write_wf

   subroutine sh_read_wf(inunit)
      integer, intent(in) :: inunit
      integer :: ist1

      if (.not. allocated(cel_re)) then
         allocate (cel_re(nstate), cel_im(nstate))
      end if

      do ist1 = 1, nstate
         read (inunit, *) cel_re(ist1), cel_im(ist1)
      end do
   end subroutine sh_read_wf

   real(DP) function check_popsum()
      real(DP) :: popsum
      integer :: ist1

      popsum = 0.0D0
      do ist1 = 1, nstate
         popsum = popsum + el_pop(ist1)
      end do

      if (abs(popsum - 1.0D0) > popsumthr) then
         write (*, *) 'ERROR:Sum of electronic populations = ', popsum
         write (*, *) 'which differs from 1.0 by more than popsumthr = ', popsumthr
         write (*, *) 'Increase the number of substeps or use more accurate integrator.'
         call fatal_error(__FILE__, __LINE__, 'Invalid SH integration')
      end if

      check_popsum = popsum
   end function check_popsum

   ! This prints detailed debug info about WF during integration
   subroutine sh_debug_wf(ist, stepfs, t)
      use mod_files, only: UBKL, UWFCOEF, UPHASE
      integer, intent(in) :: ist ! current electronic state
      real(DP), intent(in) :: stepfs
      real(DP), intent(in) :: t(:, :)
      integer :: ist1, ist2
      character(len=500) :: formt

      write (formt, '(A7,I3,A7)') '(F15.4,', nstate, 'E20.10)'
      write (UBKL, fmt=formt) stepfs, (t(ist, ist1), ist1=1, nstate)

      write (UWFCOEF, fmt=formt, advance="no") stepfs, (cel_re(ist1), ist1=1, nstate)
      write (formt, '(A1,I3,A7)') '(', nstate, 'E20.10)'
      write (UWFCOEF, fmt=formt, advance="no") (cel_im(ist1), ist1=1, nstate)
      write (UWFCOEF, *) ''

      if (phase == 1) then
         write (UPHASE, '(F15.2,E20.10)', advance="no") stepfs, gama(2, 1)
         do ist1 = 3, nstate
            write (formt, '(A1,I3,A7)') '(', ist1 - 1, 'E20.10)'
            write (UPHASE, fmt=formt, advance="no") (gama(ist1, ist2), ist2=1, ist1 - 1)
         end do
         write (UPHASE, *) ''
      end if

   end subroutine sh_debug_wf

end module mod_sh_integ
