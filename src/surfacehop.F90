! Driver routines for Surface Hopping dynamics
! D. Hollas, O. Svoboda, P. Slavicek, M. Oncak
module mod_sh
   use mod_const, only: DP
   use mod_files, only: stdout, stderr
   use mod_error, only: fatal_error
   use mod_sh_integ
   implicit none
   private

   ! INPUT PARAMETERS
   public :: read_sh_input, print_sh_input
   public :: inac
   ! Main subroutines
   public :: surfacehop, sh_init
   ! Helper subroutines
   public :: move_vars, get_nacm, write_nacmrest, read_nacmrest
   public :: check_CIVector
   public :: set_current_state

   ! Variables holding SH state
   ! TODO: Make these protected and write Set methods
   public :: en_array, tocalc
   public :: nacx, nacy, nacz

   ! Initial electronic state
   integer :: istate_init = -1

   ! Numerical integrator for SH wavefunction
   ! (default is Butcher 5-th order)
   character(len=50) :: integ = 'butcher'
   ! Number of substeps for integrating electronic Schrodinger eq.
   integer :: substep = 100

   ! Controls calculations of Non-adiabatic Couplings (NAC)
   ! couplings = 'analytic'  (inac=0) - Analytical NAC (default)
   ! couplings = 'baeck-an'  (inac=1) - Baeck-An couplings
   ! couplings = 'none'      (inac=2) - Do not compute couplings
   integer :: inac = 0 ! for working within the code
   character(len=50) :: couplings = 'analytic' ! for reading the input file
   ! energy history array and time-derivate couplings (sigma_ba) necessary for Beack-An couplings
   real(DP), allocatable :: en_hist_array(:, :), sigma_ba(:, :) ! sigma_ba is actually the dotproduct

   ! 1 - Turn OFF hopping
   integer :: nohop = 0

   ! How to adjust velocity after hop:
   ! velocity_rescaling = 'nac_then_velocity' (adjmom=0) - Adjust velocity along the NAC vector, if not possible,
   ! try the velocity vector (default)
   ! velocity_rescaling = 'velocity'  (adjmom=1) - Rescale along the velocity vector
   integer :: adjmom = 0 ! for working within the code
   character(len=50) :: velocity_rescaling = 'nac_then_velocity' ! for reading the input file
   ! 1 - Reverse momentum direction after frustrated hop
   integer :: revmom = 0

   ! Numerical accuracy of MOLPRO NACME, 10^-nac_accu
   ! Default accuracy (equals MOLPRO defaults)
   integer :: nac_accu1 = 7
   ! If calculation fails with default, we retry with nac_accu2
   integer :: nac_accu2 = 5

   ! Decoherence correction parameter (a.u.)
   ! To turn off decoherence correction, set decoh_alpha = 0.0D0
   real(DP) :: decoh_alpha = 0.1D0

   ! Compute NACME only if states are less then deltaE apart (eV)
   ! The default value (5 eV) is very conservative.
   real(DP) :: deltaE = 5.0D0
   ! Compute NACME only if at least one of the two states has population > popthr
   real(DP) :: popthr = 0.001D0

   ! Thresholds for energy conservations (in eV)
   ! The default values are too permisive
   ! you should definitely tighten them in production simulations!
   ! Stop simulation if total energy difference in successive steps is > energydifthr
   real(DP) :: energydifthr = 1.0D0
   ! Stop simulation if total energy drift since the beginning exceeds energydriftthr
   real(DP) :: energydriftthr = 1.0D0

   ! Surface Hopping with Induced Transitions,
   ! typically used with methods that can't describe S1S0 intersections,
   ! such as TDDFT or ADC2.
   ! In this case, we force a hop to S0 when the energy difference
   ! between S1 and S0 drops below user defined threshold.
   logical :: SHwIT = .false.
   real(DP) :: dE_S0S1_thr = 0.0D0 !eV
   ! NA Couplings
   real(DP), allocatable :: nacx(:, :, :), nacy(:, :, :), nacz(:, :, :)
   ! *old variables holds data from the previous step
   real(DP), allocatable :: nacx_old(:, :, :), nacy_old(:, :, :), nacz_old(:, :, :)
   real(DP), allocatable :: en_array(:), en_array_old(:)
   ! Initial absolute electronic energy, needed for monitoring energy drift
   real(DP) :: entot0
   ! nstate x nstate matrix to determine which derivatives to compute
   ! off-diagonal elements correspond to NAC
   ! diagonal elements correspond to electronic gradients
   integer, allocatable :: tocalc(:, :)
   ! Current electronic state
   integer, public, protected :: istate
   ! Do not consider hopping into this state
   ! Useful e.g. for ethylene 2-state SA3 dynamics
   ! This is a bit of an ugly hack, it would be more general to have an array
   ! of states that are calculated but ignored.
   integer :: ignore_state = 0

   namelist /sh/ istate_init, nstate, substep, deltae, integ, couplings, nohop, phase, decoh_alpha, popthr, ignore_state, &
      nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, velocity_rescaling, revmom, &
      shwit, dE_S0S1_thr, correct_decoherence
   save

contains

   ! Initialization routine (ugly, but works)
   subroutine sh_init(x, y, z, vx, vy, vz)
      use mod_const, only: AUtoEV
      use mod_general, only: irest, natom, pot
      use mod_interfaces, only: force_clas
      use mod_kinetic, only: ekin_v
      use mod_files, only: nacmefile_init
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: dum_fx(size(x, 1), size(x, 2))
      real(DP) :: dum_fy(size(y, 1), size(y, 2))
      real(DP) :: dum_fz(size(z, 1), size(z, 2))
      real(DP) :: dum_eclas

      write (stdout, *) ''
      write (stdout, '(A)') 'Initializing Surface Hopping'
      call check_sh_parameters()
      ! TODO: Write a print_sh_header()
      ! with all major parameters included and explained
      if (ignore_state /= 0) then
         write (stdout, '(A,I0,A)') 'Ignoring state ', ignore_state, ' for hopping'
      end if

      ! Determining the initial state
      if (irest /= 1) then
         call set_current_state(istate_init)
      end if

      if (nohop == 1) then
         write (stdout, '(A,I0)') 'WARNING: nohop=1, Running adiabatic dynamics on state ', istate
      end if

      deltaE = deltaE / AUtoEV

      allocate (nacx(natom, nstate, nstate))
      allocate (nacy(natom, nstate, nstate))
      allocate (nacz(natom, nstate, nstate))
      allocate (nacx_old(natom, nstate, nstate))
      allocate (nacy_old(natom, nstate, nstate))
      allocate (nacz_old(natom, nstate, nstate))
      nacx = 0.0D0
      nacy = 0.0D0
      nacz = 0.0D0
      nacx_old = nacx
      nacy_old = nacx
      nacz_old = nacx

      allocate (en_array(nstate))
      allocate (en_array_old(nstate))
      en_array = 0.0D0
      en_array_old = en_array

      ! Initialize the history array we use to calculate the Baeck-An couplings
      if (inac == 1) then
         allocate (en_hist_array(nstate, 4)) !last 3 energies (1: current, 2: n-1, 3: n-2, 4: n-3)
         en_hist_array = 0.0D0
         allocate (sigma_ba(nstate, nstate)) !this is equivalent to dotproduct, but I will need to store old and new values
         sigma_ba = 0.0D0
      end if

      allocate (tocalc(nstate, nstate))
      tocalc = 0
      tocalc(istate, istate) = 1

      ! Compute initial wavefunction and energies,
      ! used for subsequent determination of TOCALC according to deltaE threshold
      ! NOTE: We're computing forces as well, since the force_abin() interface
      ! currently does not allow for computing only energies.
      dum_eclas = 0.0D0
      dum_fx = 0.0D0; dum_fy = 0.0D0; dum_fz = 0.0D0
      call force_clas(dum_fx, dum_fy, dum_fz, x, y, z, dum_eclas, pot)

      ! open nacme_all.dat if NACME is calculated
      if (inac == 0) call nacmefile_init()

      ! When restarting, initial SH WF was already read from the restart file
      if (irest == 0) then
         call sh_set_initialwf(istate)
      end if
      call sh_set_energy_shift(en_array(1))
      call sh_select_integrator(integ)

      entot0 = en_array(istate) + ekin_v(vx, vy, vz)
      entot0 = entot0 * AUtoEV

      call set_tocalc()

      if (irest == 1) then
         call read_nacmrest()
      end if

   end subroutine sh_init

   subroutine read_sh_input(param_unit)
      use mod_utils, only: tolower
      integer, intent(in) :: param_unit

      rewind (param_unit)
      read (param_unit, sh)
      rewind (param_unit)

      integ = tolower(integ)
   end subroutine read_sh_input

   subroutine print_sh_input()
      write (stdout, nml=sh, delim='APOSTROPHE')
   end subroutine print_sh_input

   subroutine set_current_state(current_state)
      integer, intent(in) :: current_state
      istate = current_state
   end subroutine set_current_state

   subroutine check_sh_parameters()
      use mod_utils, only: int_positive, int_nonnegative, int_switch, real_nonnegative
      use mod_general, only: iknow
      use mod_chars, only: chknow
      logical :: error

      error = .false.

      if (integ /= 'euler' .and. integ /= 'rk4' .and. integ /= 'butcher') then
         write (stderr, '(A)') 'variable integ must be "euler", "rk4" or "butcher".'
         error = .true.
      end if

      if (integ /= 'butcher') then
         write (stderr, '(A)') 'WARNING: variable integ is not "butcher", which is the default and most accurate.'
         write (stderr, '(A)') chknow
         if (iknow /= 1) error = .true.
      end if

      if (istate_init > nstate) then
         write (stderr, '(A)') 'Initial state > number of computed states.'
         error = .true.
      end if

      ! converting input 'couplings' into inac which is used in the code
      select case (couplings)
      case ('analytic')
         inac = 0
         write (stdout, '(A)') 'Using analytic ab initio couplings.'
      case ('baeck-an')
         inac = 1
         write (stdout, '(A)') 'Using approximate Baeck-An couplings.'
      case ('none')
         inac = 2
         write (stdout, '(A)') 'Ignoring nonadiabatic couplings.'
      case default
         write (stderr, '(A)') 'Parameter "couplings" must be "analytic", "baeck-an" or "none".'
         error = .true.
      end select

      ! converting input 'velocity_rescaling' into inac which is used in the code
      select case (velocity_rescaling)
      case ('nac_then_velocity')
         adjmom = 0
         write (stdout, '(A)') 'Rescaling velocity along the NAC vector after hop.'
         write (stdout, '(A)') 'If there is not enough energy, try rescaling along the velocity vector.'
      case ('velocity')
         adjmom = 1
         write (stdout, '(A)') 'Rescaling velocity along the momentum vector after hop.'
      case default
         write (stderr, '(A)') 'Parameter "velocity_rescaling" must be "nac_then_velocity" or "velocity".'
         error = .true.
      end select

      if (adjmom == 0 .and. inac == 1) then
         write (stderr, '(A)') 'Combination of velocity_rescaling="nac_then_velocity" and couplings="baeck-an" is not possible.'
         write (stderr, '(A)') 'Velocity cannot be rescaled along NAC when using Baeck-An.'
         write (stderr, '(A)') 'Change velocity_rescaling="velocity" to rescale along the velocity vector.'
         error = .true.
      end if

      if (inac == 2 .and. nohop == 0 .and. .not. shwit) then
         write (stdout, '(A)') 'WARNING: For simulations without couplings(="none") hopping probability cannot be determined.'
         nohop = 1
      end if

      if (nac_accu1 <= nac_accu2) then
         write (stderr, '(A)') 'WARNING: nac_accu1 < nac_accu2'
         write (stderr, '(A,I0)') 'Computing NACME only with default accuracy 10^-', nac_accu1
      end if

      if (decoh_alpha == 0.0D0) then
         write (stdout, '(A)') "Turning OFF decoherence correction"
      end if

      if (ignore_state == istate_init) then
         write (stderr, '(A)') 'ERROR: ignore_state == istate_init'
         write (stderr, '(A)') 'I cannot start on an ignored state'
         error = .true.
      end if
      if (ignore_state > nstate) then
         write (stderr, '(A)') 'Ignored state > number of computed states.'
         error = .true.
      end if

      call int_switch(nohop, 'nohop')
      call int_switch(adjmom, 'adjmom')
      call int_switch(adjmom, 'revmom')

      call int_positive(istate_init, 'istate_init')
      call int_positive(nstate, 'nstate')
      call int_positive(substep, 'substep')
      call int_positive(nac_accu1, 'nac_accu1')

      call int_nonnegative(nac_accu2, 'nac_accu2')
      call int_nonnegative(ignore_state, 'ignore_state')

      call real_nonnegative(decoh_alpha, 'decoh_alpha')
      call real_nonnegative(deltaE, 'deltaE')
      call real_nonnegative(popthr, 'popthr')
      call real_nonnegative(popsumthr, 'popsumthr')
      call real_nonnegative(energydifthr, 'energydifthr')
      call real_nonnegative(energydriftthr, 'energydriftthr')
      call real_nonnegative(dE_S0S1_thr, 'dE_S0S1_thr')

      if (error) then
         call fatal_error(__FILE__, __LINE__, 'Invalid Surface Hopping parameter(s)')
      end if

   end subroutine check_sh_parameters

   subroutine get_nacm(pot)
      character(len=*), intent(in) :: pot
      integer :: iost, ist1, ist2
      integer :: num_nacme

      ! In TeraChem SH interface, we already got NACME
      if (pot == '_tera_') return
      if (pot == '_nai_') return ! for NaI model, the couplings are calcualted together with forces

      ! Check whether we need to compute any NACME
      num_nacme = 0
      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate
            num_nacme = num_nacme + tocalc(ist1, ist2)
         end do
      end do

      if (inac > 0 .or. num_nacme == 0) return

      ! Calculate NACME using default accuracy
      call calc_nacm(pot, nac_accu1)

      iost = read_nacm()

      if (iost /= 0 .and. nac_accu1 > nac_accu2) then
         write (stderr, '(A)') 'WARNING: Some NACME not computed. Trying with decreased accuracy'
         write (stderr, '(A,I0)') 'Calling script r.'//trim(pot)//'.nacm with accuracy: 10^-', nac_accu2
         call calc_nacm(pot, nac_accu2)
         iost = read_nacm()
      end if

      if (iost /= 0) then
         call fatal_error(__FILE__, __LINE__, 'Some NACMEs could not be read')
      end if

      ! we always have to set tocalc because we change it in read_nacm()
      call set_tocalc()
   end subroutine get_nacm

   ! In this routine, we decide which gradients and which NACME we need to compute
   ! TOCALC is a symmetric matrix, upper-triangle defines NACME,
   ! the diagonal defines gradients
   subroutine set_tocalc()
      integer :: ist1, ist2
      real(DP) :: pop, pop2

      tocalc = 0

      ! The diagonal holds information about gradients that we need
      ! for SH, we just need the gradient of the current state
      tocalc(istate, istate) = 1

      ! Do not calculate NACME for ADIABATIC dynamics
      if (inac == 2) then
         return
      end if

      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate
            if (abs(en_array(ist1) - en_array(ist2)) < deltaE) then
               tocalc(ist1, ist2) = 1
               tocalc(ist2, ist1) = 1
            end if
         end do
      end do

      if (popthr > 0) then
         ! COMPUTE NACME only if population of the states is gt.popthr
         do ist1 = 1, nstate - 1
            pop = get_elpop(ist1)
            do ist2 = ist1 + 1, nstate
               pop2 = get_elpop(ist2)
               if (pop < popthr .and. pop2 < popthr .and. ist1 /= istate .and. ist2 /= istate) then
                  tocalc(ist1, ist2) = 0
               end if
            end do
         end do
      end if

      if (ignore_state > 0) then
         do ist1 = 1, nstate
            tocalc(ist1, ignore_state) = 0
            tocalc(ignore_state, ist1) = 0
         end do
      end if
   end subroutine set_tocalc

   subroutine write_nacmrest()
      use mod_general, only: narchive, it
      use mod_qmmm, only: natqm
      use mod_utils, only: rename_file, archive_file
      integer :: ist1, ist2, iat
      integer :: iunit1
      character(len=*), parameter :: restart_file = 'restart_sh.bin'

      call rename_file(restart_file, trim(restart_file)//'.old')

      open (newunit=iunit1, file=restart_file, action='write', status="new", access="sequential", form="unformatted")

      do ist1 = 1, nstate - 1

         do ist2 = ist1 + 1, nstate

            if (inac == 0) then

               do iat = 1, natqm
                  write (iunit1) nacx(iat, ist1, ist2), &
                               & nacy(iat, ist1, ist2), &
                               & nacz(iat, ist1, ist2)
               end do

            end if

         end do

      end do

      ! Printing total energy at t=0, so that we can safely restart
      ! and we do not break checking for energy drift
      ! For now, energy jump check is not handled well.
      write (iunit1) entot0

      if (phase == 1) then
         call sh_write_phase_bin(iunit1)
      end if

      close (iunit1)

      if (modulo(it, narchive) == 0) then
         call archive_file(restart_file, it)
      end if

   end subroutine write_nacmrest

   subroutine read_nacmrest()
      use mod_general, only: it
      use mod_qmmm, only: natqm
      use mod_utils, only: archive_file
      character(len=*), parameter :: restart_file = 'restart_sh.bin'
      integer :: iost, ist1, ist2, iat
      integer :: iunit1

      write (stdout, *) 'Reading Surface Hopping restart data from file '//trim(restart_file)

      open (newunit=iunit1, file=restart_file, action="read", status="old", access="sequential", form="unformatted")

      do ist1 = 1, nstate - 1

         do ist2 = ist1 + 1, nstate

            if (inac == 0) then

               do iat = 1, natqm
                  read (iunit1, iostat=iost) nacx(iat, ist1, ist2), &
                                                        & nacy(iat, ist1, ist2), &
                                                        & nacz(iat, ist1, ist2)
                  if (iost /= 0) then
                     call fatal_error(__FILE__, __LINE__, &
                        & 'Could not read NACME from restart file '//trim(restart_file))
                  end if
                  nacx(iat, ist2, ist1) = -nacx(iat, ist1, ist2)
                  nacy(iat, ist2, ist1) = -nacy(iat, ist1, ist2)
                  nacz(iat, ist2, ist1) = -nacz(iat, ist1, ist2)
               end do

            end if

         end do
      end do

      ! Reading total energy at t=0, so that we can safely restart
      ! and we do not break checking for energy drift
      ! For now, energy jump check is not handled well.
      read (iunit1) entot0

      if (phase == 1) then
         call sh_read_phase_bin(iunit1)
      end if

      close (iunit1)
      call archive_file(restart_file, it)

   end subroutine read_nacmrest

   integer function read_nacm()
      use mod_qmmm, only: natqm
      integer :: iost, ist1, ist2, iat, iunit

      iost = 0 ! needed if each tocalc=0
      open (newunit=iunit, file='nacm.dat')
      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate

            ! TODO: We should have some validation here
            ! that we're indeed reading NACME between correct states.
            if (tocalc(ist1, ist2) == 1) then

               do iat = 1, natqm ! reading only for QM atoms
                  read (iunit, *, IOSTAT=iost) nacx(iat, ist1, ist2), &
                                             & nacy(iat, ist1, ist2), &
                                             & nacz(iat, ist1, ist2)
                  if (iost == 0) then
                     tocalc(ist1, ist2) = 0 ! marking as read, useful if we do decreased accuracy
                     nacx(iat, ist2, ist1) = -nacx(iat, ist1, ist2)
                     nacy(iat, ist2, ist1) = -nacy(iat, ist1, ist2)
                     nacz(iat, ist2, ist1) = -nacz(iat, ist1, ist2)
                  else
                     close (iunit, status='delete')
                     write (stderr, '(A,I0,A,I0,A)') 'WARNING: NACME between states ', ist1, ' and ', ist2, ' not read.'
                     read_nacm = iost
                     return
                  end if
               end do

               ! if tocalc
            end if

         end do
      end do

      close (iunit, status='delete')
      read_nacm = iost
   end function read_nacm

   subroutine calc_nacm(pot, nac_accu)
      use mod_utils, only: toupper
      use mod_general, only: it
      character(len=*), intent(in) :: pot
      integer, intent(in) :: nac_accu
      integer :: ist1, ist2, u, itrj
      integer :: icmd, istat
      character(len=100) :: chsystem

      open (newunit=u, file='state.dat')
      write (u, '(I2)') nstate
      ! we print upper triangular part of tocalc matrix to file state.dat
      ! tocalc(i,j)==1 -> calculate NACME between states i and j
      ! tocalc(,)==0 -> don't calculate NACME
      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate
            write (u, '(I1,A1)', advance='no') tocalc(ist1, ist2), ' '
         end do
      end do
      close (u)

      chsystem = './'//trim(toupper(pot))//'/r.'//trim(pot)//'.nacm '
      itrj = 1
      write (chsystem, '(A,X,I0,X,I4.3,X,I0,X,A)') trim(chsystem), it, itrj, nac_accu, ' < state.dat'

      call execute_command_line(trim(chsystem), exitstat=istat, cmdstat=icmd)

      if (icmd /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Could not execute command "'//trim(chsystem)//'"')
      end if

      ! We continue on, since we can retry calculation with decreased accuracy
      if (istat /= 0) then
         write (stderr, *) 'WARNING: Command "'//trim(chsystem)//'" exited with non-zero status'
      end if
   end subroutine calc_nacm

   ! Calculate Baeck-An couplings
   ! Implementation was based on these two articles:
   ! original by Barbatti: https://doi.org/10.12688/openreseurope.13624.2
   ! another implementation by Truhlar: https://doi.org/10.1021/acs.jctc.1c01080
   ! In the numeric implementation, we follow Barbatti and use a higher order formula.
   subroutine calc_baeckan(dt)
      use mod_general, only: it
      integer :: ist1, ist2
      real(DP), intent(in) :: dt
      real(DP) :: de(4), de2dt2, argument

      sigma_ba = 0.0D0

      ! we don't have sufficient history until step 4
      if (it < 4) return

      do ist1 = 1, nstate
         do ist2 = ist1 + 1, nstate
            ! If ignore_state is set, we do not calculate sigma (dotproduct) for this state
            if (ignore_state == ist1 .or. ignore_state == ist2) cycle
            de = en_hist_array(ist2, :) - en_hist_array(ist1, :)
            ! Second derivative (de2dt2) comes from Eq. 16 from https://doi.org/10.12688/openreseurope.13624.2
            de2dt2 = (2.0D0 * de(1) - 5.0D0 * de(2) + 4.0D0 * de(3) - de(4)) / dt**2
            argument = de2dt2 / de(1)
            if (argument > 0.0D0) then
               sigma_ba(ist2, ist1) = dsqrt(argument) / 2.0D0
            end if
            sigma_ba(ist1, ist2) = -sigma_ba(ist2, ist1)
         end do
      end do

   end subroutine calc_baeckan

   ! move arrays from new step to old step
   subroutine move_vars(vx, vy, vz, vx_old, vy_old, vz_old)
      use mod_general, only: natom
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(out) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      integer :: ist1, ist2, iat

      do ist1 = 1, nstate
         en_array_old(ist1) = en_array(ist1)

         if (inac == 0) then ! nedelame, pokud nacitame rovnou dotprodukt
            do ist2 = 1, nstate
               do iat = 1, natom
                  nacx_old(iat, ist1, ist2) = nacx(iat, ist1, ist2)
                  nacy_old(iat, ist1, ist2) = nacy(iat, ist1, ist2)
                  nacz_old(iat, ist1, ist2) = nacz(iat, ist1, ist2)
               end do
            end do
         end if
      end do

      ! Shift the energy history for Baeck-An couplings
      if (inac == 1) then
         ! Move old energies by 1
         en_hist_array(:, 4) = en_hist_array(:, 3)
         en_hist_array(:, 3) = en_hist_array(:, 2)
         en_hist_array(:, 2) = en_hist_array(:, 1)
         ! new energy is stored before the couplings are calculated
         ! I avoided doing the same as with LZSH energy history tracking because then I need to modify force_abin, force_terash and
         ! every potential in potentials_sh (and all possible future ones). This way it is kept private and does not depend on the
         ! way energies are calculated.
         ! See commit: https://github.com/PHOTOX/ABIN/pull/186/commits/918f6837a76ec0f41b575d3ca948253eed2f30cc
      end if

      vx_old = vx
      vy_old = vy
      vz_old = vz
   end subroutine move_vars

   !******************************
   ! This is the main SH routine !
   !******************************
   subroutine surfacehop(x, y, z, vx, vy, vz, vx_old, vy_old, vz_old, dt, eclas)
      use mod_general, only: natom, nwrite, idebug, it, sim_time, pot
      use mod_const, only: AUTOFS
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP), intent(inout) :: eclas
      real(DP), intent(in) :: dt
      ! Interpolated quantities
      real(DP), dimension(natom, 1) :: vx_int, vy_int, vz_int
      real(DP), dimension(natom, 1) :: vx_newint, vy_newint, vz_newint
      real(DP), dimension(nstate) :: en_array_int, en_array_newint
      real(DP), dimension(natom, nstate, nstate) :: nacx_int, nacy_int, nacz_int
      real(DP), dimension(natom, nstate, nstate) :: nacx_newint, nacy_newint, nacz_newint
      real(DP), dimension(nstate, nstate) :: dotproduct_int, dotproduct_newint, sigma_ba_old
      ! Switching probabilities
      real(DP) :: t(nstate, nstate)
      ! Cumulative switching probabilities
      real(DP) :: t_tot(nstate, nstate)
      real(DP) :: popsum !populations
      integer :: ist1, itp
      integer :: ihop
      ! Shortened variable for current state (istate)
      ! TODO: Rename this!
      integer :: ist ! =istate
      real(DP) :: fr
      real(DP) :: Ekin, dtp
      real(DP) :: pop0
      ! Simulation time in femtoseconds
      real(DP) :: stepfs
      character(len=500) :: formt

      call check_energy(vx_old, vy_old, vz_old, vx, vy, vz)
      call check_energydrift(vx, vy, vz)

      t_tot = 1.0D0

      ! First, calculate NACME
      if (inac == 0) then ! Analytic ab initio couplings
         ! For TeraChem MPI / FMS interface, NAC are already computed!
         if (pot /= '_tera_' .and. pot /= '_nai_') then
            nacx = 0.0D0
            nacy = 0.0D0
            nacz = 0.0D0
            ! Compute and read NACME (MOLPRO-SH interface)
            call get_nacm(pot)
         end if
         ! TODO: Should we call this with TeraChem?
         ! I think TC already phases the couplings internally.
         call phase_nacme(nacx_old, nacy_old, nacz_old, nacx, nacy, nacz)
      else if (inac == 1) then ! Baeck-An couplings
         ! saving the current energy to the energy history (shifting was already done in previous step in move_vars)
         en_hist_array(:, 1) = en_array(:)
         sigma_ba_old = sigma_ba ! saving old sigma_ba
         call calc_baeckan(dt)
      end if

      ! smaller time step
      dtp = dt / substep

      ! MAIN LOOP
      ! Smaller time step for electronic population transfer
      do itp = 1, substep

         ist = istate

         ! pop0 is later used for Tully's fewest switches
         pop0 = get_elpop(ist)

         ! INTERPOLATION
         if ((inac == 0) .or. (inac == 2)) then
            fr = real(itp, DP) / real(substep, DP)
            call interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_newint, vy_newint, vz_newint, &
                             nacx_newint, nacy_newint, nacz_newint, en_array_newint, &
                             dotproduct_newint, fr)

            fr = real(itp - 1, DP) / real(substep, DP)
            call interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                             nacx_int, nacy_int, nacz_int, en_array_int, &
                             dotproduct_int, fr)
         else if (inac == 1) then
            fr = real(itp, DP) / real(substep, DP)
            call interpolate_ba(vx, vy, vz, vx_old, vy_old, vz_old, vx_newint, vy_newint, vz_newint, &
                                en_array_newint, dotproduct_newint, sigma_ba, sigma_ba_old, fr)

            fr = real(itp - 1, DP) / real(substep, DP)
            call interpolate_ba(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                                en_array_int, dotproduct_int, sigma_ba, sigma_ba_old, fr)

         end if

         ! Integrate electronic wavefunction for one dtp time step
         call sh_integrate_wf(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, dtp)

         ! Check whether total population is 1.0
         popsum = check_popsum()

         ! Calculate switching probabilities according to Tully's fewest switches
         ! Probabilities are stored in matrix t
         ! (although at this point we calculate only the relevant part of the matrix)
         call sh_TFS_transmat(dotproduct_int, dotproduct_newint, istate, pop0, t, dtp)

         if (idebug > 1) then
            stepfs = (sim_time + dtp * itp - substep * dtp) * AUtoFS
            call sh_debug_wf(ist, stepfs, t)
         end if

         do ist1 = 1, nstate
            ! Cumulative probability over whole big step
            ! This is only auxiliary variable used for output
            t_tot(ist, ist1) = t_tot(ist, ist1) * (1 - t(ist, ist1))
         end do

         ! Should we hop before decoherence?
         ! Every article says something different.
         ! Newton-X hops before decoherence.

         ! HOPPING SECTION
         call random_hop(transmat=t, old_state=ist, new_state=ihop)

         if (nohop == 1 .and. ihop /= 0) then
            write (stdout, '(A,I0)') 'WARNING: Ignoring hop to state ', ihop
         end if

         ! Did HOP occur?
         if (nohop /= 1 .and. ihop /= 0) then
            ! NOTE: Hop can still fail due to insufficient kinetic energy
            ! ("frustrated hop")
            if (adjmom == 0) then
               call try_hop_nacme_rescale(vx, vy, vz, ist, ihop, eclas)
            else if (adjmom == 1) then
               call try_hop_simple_rescale(vx, vy, vz, ist, ihop, eclas)
            else
               call fatal_error(__FILE__, __LINE__, 'Invalid adjmom value')
            end if

            ! NOTE: We're writing this geometry even for frustrated hops!
            ! Not sure if that is desired?
            call write_hopgeom(x, y, z, old_state=ist, new_state=ihop, timestep=it)
         end if

         ! Apply decoherence correction from Persico et al
         if (decoh_alpha > 0) then

            Ekin = ekin_v(vx_int, vy_int, vz_int)
            if (Ekin > 1.0D-4) then ! Decoherence diverges for zero velocities
               call sh_decoherence_correction(en_array_int, decoh_alpha, Ekin, istate, dtp)
            end if

         end if

         !itp loop
      end do

      if (SHwIT) then
         call shwit_check(en_array(1), en_array(2), dE_S0S1_thr, &
               vx, vy, vz, eclas)
      end if


      ! set tocalc array for the next step
      call set_tocalc()

      popsum = check_popsum()

      call move_vars(vx, vy, vz, vx_old, vy_old, vz_old)

      if (modulo(it, nwrite) == 0) then
         call write_sh_output()
      end if

   contains

      subroutine write_sh_output()
         use mod_general, only: sim_time
         use mod_const, only: ANG, AUTOFS
         use mod_files, only: UPOP, UPROB, UPES, UNACME, UDOTPROD
         integer :: ist1, ist2, iat
         real(DP) :: stepfs

         ! Simulation time in femtoseconds
         stepfs = sim_time * AUtoFS

         ! Write electronic populations
         write (formt, '(A10,I0,A13)') '(F15.2,I3,', nstate, 'F10.5,1F10.7)'
         write (UPOP, fmt=formt) stepfs, istate, (get_elpop(ist1), ist1=1, nstate), popsum

         ! Write hopping probabilities
         t_tot = 1 - t_tot ! up to know, t_tot was the probability of not hopping
         write (formt, '(A10,I0,A6)') '(F15.2,I3,', nstate, 'F10.5)'
         write (UPROB, fmt=formt) stepfs, istate, (t_tot(ist, ist1), ist1=1, nstate)

         ! Write potential energies to PES.dat
         write (formt, '(A7,I0,A7)') '(F15.2,', nstate, 'E20.10)'
         write (UPES, fmt=formt) stepfs, (en_array(ist1), ist1=1, nstate)

         if (inac == 0) write (UNACME, *) 'Time step:', it
         do ist1 = 1, nstate - 1
            do ist2 = ist1 + 1, nstate

               if (ist1 == 1 .and. ist2 == 2) then
                  write (UDOTPROD, '(F15.2,E20.10)', advance='no') stepfs, dotproduct_int(ist1, ist2)
               else
                  write (UDOTPROD, '(E20.10)', advance='no') dotproduct_int(ist1, ist2)
               end if

               if (inac == 0) then
                  write (UNACME, *) 'NACME between states:', ist1, ist2
                  do iat = 1, natom
                     write (UNACME, '(3E20.10)') nacx(iat, ist1, ist2), &
                                               & nacy(iat, ist1, ist2), &
                                               & nacz(iat, ist1, ist2)
                  end do
               end if

            end do
         end do
         write (UDOTPROD, *) ''
      end subroutine write_sh_output

   end subroutine surfacehop

   subroutine phase_nacme(nacx_old, nacy_old, nacz_old, nacx, nacy, nacz)
      real(DP), intent(in), dimension(:, :, :) :: nacx_old, nacy_old, nacz_old
      real(DP), intent(inout), dimension(:, :, :) :: nacx, nacy, nacz
      integer :: ist1, ist2, iat
      integer :: natom
      real(DP) :: vect_olap

      natom = size(nacx, 1)
      ! Calculating overlap between nacmes
      ! This is crucial, since NACME vectors can change
      ! orientation 180 degrees between time steps and we need to correct that
      do ist1 = 1, nstate
         do ist2 = 1, nstate

            vect_olap = 0.0D0
            do iat = 1, natom
               vect_olap = vect_olap + &
                  & nacx_old(iat, ist1, ist2) * nacx(iat, ist1, ist2)
               vect_olap = vect_olap + &
                  & nacy_old(iat, ist1, ist2) * nacy(iat, ist1, ist2)
               vect_olap = vect_olap + &
                  & nacz_old(iat, ist1, ist2) * nacz(iat, ist1, ist2)
            end do

            if (vect_olap < 0) then
               do iat = 1, natom
                  nacx(iat, ist1, ist2) = -nacx(iat, ist1, ist2)
                  nacy(iat, ist1, ist2) = -nacy(iat, ist1, ist2)
                  nacz(iat, ist1, ist2) = -nacz(iat, ist1, ist2)
               end do
            end if

         end do
      end do
   end subroutine phase_nacme

   subroutine write_hopgeom(x, y, z, old_state, new_state, timestep)
      use mod_system, only: names
      use mod_const, only: ANG
      real(DP), intent(in), dimension(:, :) :: x, y, z
      integer, intent(in) :: old_state, new_state, timestep
      character(len=100) :: formt
      integer :: natom
      integer :: iat, u

      natom = size(names)

      write (formt, '("hopgeom.",I0,".",I0,".",I0,".xyz")') old_state, new_state, timestep
      open (newunit=u, file=trim(formt), action='write')
      write (u, *) natom
      write (u, *) ''
      do iat = 1, natom
         write (u, '(A,3ES25.16E3)') names(iat), x(iat, 1) / ANG, y(iat, 1) / ANG, z(iat, 1) / ANG
      end do
      close (u)
   end subroutine write_hopgeom

   subroutine random_hop(transmat, old_state, new_state)
      use mod_random, only: vranf
      ! Transition matrix
      real(DP), intent(in), dimension(:, :) :: transmat
      ! Current state index
      integer, intent(in) :: old_state
      ! New state index (0 if no hop occurs)
      integer, intent(out) :: new_state
      real(DP) :: prob(nstate)
      ! Pseudorandom number
      real(DP) :: rdnum(1)
      integer :: ist1

      ! We return 0 if there's no hop
      new_state = 0

      ! Auxiliary calculations of probabilities on a number line
      prob = 0.0D0

      ! If we are in the ground state, we cannot jump into the ground state :-)
      if (old_state /= 1) then
         ! Probability of jumping from the current state to the ground state
         prob(1) = transmat(old_state, 1)
      end if

      do ist1 = 2, nstate
         if (ist1 == old_state) then
            prob(ist1) = prob(ist1 - 1)
         else
            prob(ist1) = prob(ist1 - 1) + transmat(old_state, ist1)
         end if
      end do

      ! Get one random number between 0 and 1
      call vranf(rdnum, 1)

      ! Determine, whether we hopped or not
      do ist1 = 1, nstate
         if (ist1 == old_state) cycle
         if (rdnum(1) < prob(ist1)) then
            new_state = ist1
            exit
         end if
      end do
   end subroutine random_hop

   subroutine try_hop_nacme_rescale(vx, vy, vz, instate, outstate, eclas)
      use mod_general, only: natom, pot
      use mod_system, only: am
      use mod_files, only: UPOP
      use mod_kinetic, only: ekin_v
      use mod_arrays, only: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: instate, outstate
      real(DP) :: a_temp, b_temp, c_temp, g_temp
      real(DP) :: ekin, ekin_new
      integer :: iat, iw

      iw = 1
      write (stdout, '(A,I0,A,I0)') 'Trying to hop from state ', instate, ' to state ', outstate
      write (stdout, *) "Checking kinetic energy in the direction of NACME vector"

      ! Checking for frustrated hop
      a_temp = 0.D0
      b_temp = 0.D0

      ekin = ekin_v(vx, vy, vz)

      do iat = 1, natom
         a_temp = a_temp + nacx(iat, instate, outstate)**2 / am(iat)
         a_temp = a_temp + nacy(iat, instate, outstate)**2 / am(iat)
         a_temp = a_temp + nacz(iat, instate, outstate)**2 / am(iat)
         b_temp = b_temp + nacx(iat, instate, outstate) * vx(iat, iw)
         b_temp = b_temp + nacy(iat, instate, outstate) * vy(iat, iw)
         b_temp = b_temp + nacz(iat, instate, outstate) * vz(iat, iw)
      end do
      a_temp = 0.5D0 * a_temp
      c_temp = b_temp**2 + 4 * a_temp * (en_array(instate) - en_array(outstate))

      if (a_temp <= 0.0D0) then
         write (stdout, *) 'WARNING: NACME vector is zero, rescaling velocities isotropically along the velocity vector'
         call try_hop_simple_rescale(vx, vy, vz, instate, outstate, eclas)
         return
      end if

      if (c_temp < 0) then
         write (stdout, *) 'WARNING:  Not enough kinetic energy in the direction of NAC vector.'
         write (stdout, *) 'Trying isotropic velocity rescaling instead'
         ! Try, whether there is enough total kinetic energy and scale velocities.
         call try_hop_simple_rescale(vx, vy, vz, instate, outstate, eclas)
         return
      end if

      call set_current_state(outstate)
      eclas = en_array(outstate)
      write (stdout, '(A,I0,A,I0)') '# Hop occured from state ', instate, ' to state ', outstate
      write (stdout, '(A, E20.10)') 'Potential energy difference / a.u. = ', en_array(outstate) - en_array(instate)
      write (stdout, '(A, E20.10)') 'Total kinetic energy before hop / a.u. = ', ekin

      ! Rescaling the velocities

      if (b_temp < 0) then
         g_temp = (b_temp + dsqrt(c_temp)) / 2.0D0 / a_temp
      end if

      if (b_temp >= 0) then
         g_temp = (b_temp - dsqrt(c_temp)) / 2.0D0 / a_temp
      end if

      iw = 1
      do iat = 1, natom
         vx(iat, iw) = vx(iat, iw) - g_temp * nacx(iat, instate, outstate) / am(iat)
         vy(iat, iw) = vy(iat, iw) - g_temp * nacy(iat, instate, outstate) / am(iat)
         vz(iat, iw) = vz(iat, iw) - g_temp * nacz(iat, instate, outstate) / am(iat)
      end do
      ekin_new = ekin_v(vx, vy, vz)

      write (stdout, '(A,2E20.10)') '# dE_pot     E_kin-total', &
                               & en_array(outstate) - en_array(instate), ekin

      call set_tocalc()
      write (stdout, *) '# Calculating forces for the new state.'
      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

   end subroutine try_hop_nacme_rescale

   subroutine interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                          nacx_int, nacy_int, nacz_int, en_array_int, &
                          dotproduct_int, fr)
      use mod_general, only: natom
      real(DP), intent(out) :: dotproduct_int(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP), intent(out) :: vx_int(:, :), vy_int(:, :), vz_int(:, :)
      real(DP), intent(out) :: nacx_int(:, :, :)
      real(DP), intent(out) :: nacy_int(:, :, :)
      real(DP), intent(out) :: nacz_int(:, :, :)
      real(DP), intent(out) :: en_array_int(:)
      ! How far are we interpolating?
      real(DP), intent(in) :: fr
      real(DP) :: frd
      integer :: iat, iw, ist1, ist2 !iteration counters

      iw = 1
      frd = 1.0D0 - fr

      do ist1 = 1, nstate

         en_array_int(ist1) = en_array(ist1) * fr + en_array_old(ist1) * frd
         do ist2 = 1, nstate
            dotproduct_int(ist1, ist2) = 0.0D0
            do iat = 1, natom
               nacx_int(iat, ist1, ist2) = nacx(iat, ist1, ist2) * fr + &
                                               & nacx_old(iat, ist1, ist2) * frd
               nacy_int(iat, ist1, ist2) = nacy(iat, ist1, ist2) * fr + &
                                               & nacy_old(iat, ist1, ist2) * frd
               nacz_int(iat, ist1, ist2) = nacz(iat, ist1, ist2) * fr + &
                                               & nacz_old(iat, ist1, ist2) * frd
               vx_int(iat, iw) = vx(iat, iw) * fr + vx_old(iat, iw) * frd
               vy_int(iat, iw) = vy(iat, iw) * fr + vy_old(iat, iw) * frd
               vz_int(iat, iw) = vz(iat, iw) * fr + vz_old(iat, iw) * frd
               dotproduct_int(ist1, ist2) = dotproduct_int(ist1, ist2) + &
                                               & vx_int(iat, iw) * nacx_int(iat, ist1, ist2)
               dotproduct_int(ist1, ist2) = dotproduct_int(ist1, ist2) + &
                                               & vy_int(iat, iw) * nacy_int(iat, ist1, ist2)
               dotproduct_int(ist1, ist2) = dotproduct_int(ist1, ist2) + &
                                               & vz_int(iat, iw) * nacz_int(iat, ist1, ist2)
            end do
         end do
      end do

   end subroutine interpolate

   ! interpolation of time-derivative coupling calculated via Baeck-An approximation
   ! this routine interpolates sigma_ba between integration steps
   subroutine interpolate_ba(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                             en_array_int, dotproduct_int, sigma_ba, sigma_ba_old, fr)
      use mod_general, only: natom
      real(DP), intent(in) :: sigma_ba(:, :), sigma_ba_old(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :) ! for velocity interpolation
      real(DP), intent(in) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP), intent(out) :: dotproduct_int(:, :)
      real(DP), intent(out) :: en_array_int(:)
      real(DP), intent(out) :: vx_int(:, :), vy_int(:, :), vz_int(:, :) ! interpolated velocities
      ! How far are we interpolating?
      real(DP), intent(in) :: fr
      real(DP) :: frd
      integer :: ist1, ist2, iw, iat !iteration counters

      frd = 1.0D0 - fr

      do ist1 = 1, nstate
         en_array_int(ist1) = en_array(ist1) * fr + en_array_old(ist1) * frd
         do ist2 = 1, nstate
            ! interpolating dot product
            dotproduct_int(ist1, ist2) = sigma_ba(ist1, ist2) * fr + sigma_ba_old(ist1, ist2) * frd
         end do
      end do

      ! interpolating velocity which is necessary for Ekin in the decoherence correction
      iw = 1
      do iat = 1, natom
         vx_int(iat, iw) = vx(iat, iw) * fr + vx_old(iat, iw) * frd
         vy_int(iat, iw) = vy(iat, iw) * fr + vy_old(iat, iw) * frd
         vz_int(iat, iw) = vz(iat, iw) * fr + vz_old(iat, iw) * frd
      end do

   end subroutine interpolate_ba

   subroutine try_hop_simple_rescale(vx, vy, vz, instate, outstate, eclas)
      use mod_general, only: pot
      use mod_kinetic, only: ekin_v
      use mod_arrays, only: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: instate, outstate
      real(DP) :: de, ekin, alfa, ekin_new

      ekin = 0.0D0
      ekin_new = 0.0D0

      dE = en_array(outstate) - en_array(instate)
      ekin = ekin_v(vx, vy, vz)

      if (ekin >= de) then

         alfa = dsqrt(1 - de / ekin)

         vx = alfa * vx
         vy = alfa * vy
         vz = alfa * vz

         ! Switch states
         call set_current_state(outstate)
         eclas = en_array(outstate)
         ekin_new = ekin_v(vx, vy, vz)

         write (*, '(A,I0,A,I0)') '# Hop occured from state ', instate, ' to state ', outstate
         write (*, '(A)') '# Adjusting velocities by isotropic scaling.'

         call set_tocalc()
         write (*, '(A)') '# Calculating forces for the new state.'
         call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

      else

         write (*, '(A,I0,A,I0)') '# Frustrated Hop occured from state ', &
                                    & instate, ' to state ', outstate
         if (revmom == 1) then
            write (*, '(A)') '# Reversing momentum direction.'
            vx = -vx
            vy = -vy
            vz = -vz
         end if

      end if

      write (*, '(A,E17.10,A,E17.10)') '# deltaE_pot / a.u. = ', dE, ' E_kin-total / a.u. = ', ekin

   end subroutine try_hop_simple_rescale

   subroutine check_energy(vx_old, vy_old, vz_old, vx, vy, vz)
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP) :: ekin, ekin_old, entot, entot_old

      ekin = ekin_v(vx, vy, vz)
      ekin_old = ekin_v(vx_old, vy_old, vz_old)

      entot = (ekin + en_array(istate)) * AUtoEV

      ! TODO: But what if we hopped to another state?
      ! en_array_old(istate) would not point to the correct energy,
      ! right?
      ! We need istate_old array ... but we should just refactor the whole thing...
      ! and make traj_data structure or something
      entot_old = (ekin_old + en_array_old(istate)) * AUtoEV

      if (abs(entot - entot_old) > energydifthr) then
         write (stderr, *) 'Total energy difference [eV] is:', entot - entot_old
         write (stderr, *) 'The threshold was:', energydifthr
         write (stderr, *) 'Ekin_old, Ekin, Epot_old, E_pot', &
            ekin_old, ekin, en_array_old(istate), en_array(istate)
         call fatal_error(__FILE__, __LINE__, 'Poor energy conservation')
      end if

   end subroutine check_energy

   subroutine shwit_check(en_S0, en_S1, threshold_ev, vx, vy, vz, eclas)
      use mod_const, only: AUtoEV
      real(DP), intent(in) :: en_S0, en_S1
      real(DP), intent(in) :: threshold_ev
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :), eclas
      real(DP) :: dE_S0S1

      ! We only hop from S1 to S0
      if (istate /= 2) return

      dE_S0S1 = (en_S1 - en_S0) * AUtoEV

      if (dE_S0S1 < threshold_ev) then
         write (stdout, *) 'SHwIT: S1 - S0 gap dropped below threshold!'
         write (stdout, *) dE_S0S1, ' < ', threshold_ev
         write (stdout, *) "Let's jump to S0 and continue!"
         ! This global flag is checked in the main loop in abin.F90
         ! STOP_SIMULATION = .true.
         call try_hop_simple_rescale(vx, vy, vz, 2, 1, eclas)
         ! Horrible hack to exchange c_el coefficients between S1 - S0
         call shwit_switch()
      end if
   end subroutine shwit_check

   subroutine check_energydrift(vx, vy, vz)
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: ekin, entot

      ekin = ekin_v(vx, vy, vz)

      entot = (ekin + en_array(istate)) * AUtoEV

      if (abs(entot - entot0) > energydriftthr) then
         write (stderr, *) 'Total energy drift [eV] is:', entot - entot0
         write (stderr, *) 'The threshold was:', energydriftthr
         call fatal_error(__FILE__, __LINE__, 'Energy drift exceeded threshold value')
      end if

   end subroutine check_energydrift

   integer function check_CIVector(CIvecs, CIvecs_old, ci_len, nstates)
      use mod_const, only: AUtoFS
      use mod_files, only: UDOTPRODCI
      use mod_general, only: sim_time
      real(DP), intent(in) :: CIvecs(:, :), CIvecs_old(:, :)
      integer, intent(in) :: ci_len, nstates
      real(DP) :: cidotprod(nstates)
      integer :: ist1, i
      character(len=20) :: formt

      do ist1 = 1, nstates
         cidotprod(ist1) = 0.0D0
         do i = 1, ci_len
            cidotprod(ist1) = cidotprod(ist1) + CIvecs(i, ist1) * CIvecs_old(i, ist1)
         end do
         if (cidotprod(ist1) < 1 / dsqrt(2.0D0)) then
            write (*, *) "Warning: cidotprod too low."
         end if
      end do

      write (formt, '(A7,I3,A7)') '(F15.2,', nstates, 'F15.10)'
      write (UDOTPRODCI, fmt=formt) sim_time * AUtoFS, (cidotprod(ist1), ist1=1, nstates)

      check_CIVector = 0
      return
   end function check_CIVector

end module mod_sh
