! Driver routines for Landau-Zener excited state dynamics
!    Based on PHYSICAL REVIEW A 84, 014701 (2011)
!    Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
!    Andrey K. Belyaev and Oleg V. Lebedev
! Implemented by: J. Suchan, (J. Chalabala)
! Modified by J. Janos (2023)

! No S/T TeraChem functionality YET
! No energy drift check
! No PIMD functionality

module mod_lz
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
   use mod_sh, only: set_current_state
   implicit none
   private
   ! Public subroutines
   public :: lz_init, print_lz_input, read_lz_input, write_lz_data
   public :: lz_hop, lz_rewind, lz_restin, lz_restout, lz_finalize
   ! Other variables that must be public (but perhaps shouldn't be!)
   ! Caveat: Every time we call force_clas en_array_lz is updated
   public :: nstate_lz
   public :: en_array_lz, tocalc_lz, istate_lz

   ! Energy difference (eV), up to which we consider possibility of LZ hops
   real(DP) :: deltaE_lz = 1.0D0
   ! Maximum energy difference (eV) between two consecutive steps, used in check_energy_lz()
   real(DP) :: energydifthr_lz = 0.5D0
   ! Initial electronic state
   integer :: initstate_lz = 1
   integer :: nstate_lz, nsinglet_lz = 0, ntriplet_lz = 0
   ! Do not consider hopping into this state
   ! Same meaning as in mod_sh
   integer, public :: ignore_state = 0

   ! Module variables
   integer :: istate_lz
   integer, allocatable, dimension(:) :: tocalc_lz
   real(DP), allocatable, dimension(:, :) :: en_array_lz, en_array_lz_backup
   real(DP), allocatable, dimension(:, :) :: fx_old, fy_old, fz_old
   real(DP), allocatable, dimension(:, :) :: px_temp, py_temp, pz_temp
   real(DP), allocatable, dimension(:, :) :: x_prev, y_prev, z_prev
   real(DP), allocatable, dimension(:, :) :: vx_prev, vy_prev, vz_prev

   namelist /lz/ initstate_lz, nstate_lz, nsinglet_lz, ntriplet_lz, ignore_state, deltaE_lz, energydifthr_lz
   save

contains

   subroutine read_lz_input(param_unit)
      integer, intent(in) :: param_unit

      rewind (param_unit)
      read (param_unit, lz)
      rewind (param_unit)
   end subroutine read_lz_input

   subroutine print_lz_input()
      write (stdout, nml=lz, delim='APOSTROPHE')
   end subroutine print_lz_input

   subroutine check_lz_parameters()
      use mod_utils, only: int_positive, int_nonnegative, real_nonnegative, real_positive

      if (initstate_lz > nstate_lz) then
         call fatal_error(__FILE__, __LINE__, &
            & '(LZ): Initial state > number of computed states.')
      end if

      if (nsinglet_lz == 0 .and. ntriplet_lz == 0 .and. nstate_lz > 0) then
         ! Assume singlet states by default
         nsinglet_lz = nstate_lz
      end if

      if ((nsinglet_lz + ntriplet_lz) /= nstate_lz) then
         call fatal_error(__FILE__, __LINE__, &
            & '(LZ): Sum of singlet and triplet states must give total number of states.')
      end if

      if (ignore_state == initstate_lz) then
         call fatal_error(__FILE__, __LINE__, &
            & 'ignore_state == initstate_lz, cannot start simulation on ignored state')
      end if

      if (ignore_state > nstate_lz) then
         call fatal_error(__FILE__, __LINE__, &
            & '(LZ): ignore_state > nstate_lz')
      end if

      if (ignore_state /= 0) then
         write (stdout, '(A,I0,A)') 'Ignoring state ', ignore_state, ' for hopping'
      end if

      if (nsinglet_lz > 2 .or. ntriplet_lz > 2) then
         write (*, *) 'WARNING: LZ was derived for a two-state problem. More states might cause unphysical behavior.'
         write (stdout, *) ''
      end if

      call int_positive(initstate_lz, 'initstate_lz')
      call int_positive(nstate_lz, 'nstate_lz')
      call int_nonnegative(ignore_state, 'ignore_state')
      call int_nonnegative(nsinglet_lz, 'nsinglet_lz')
      call int_nonnegative(ntriplet_lz, 'ntriplet_lz')
      call real_positive(deltaE_lz, 'deltaE_lz')
      call real_nonnegative(energydifthr_lz, 'energydifthr_lz')
   end subroutine check_lz_parameters

   !Initialization
   subroutine lz_init(pot)
      use mod_general, only: natom, nwalk
      character(len=*), intent(in) :: pot
      integer :: ist1

      call check_lz_parameters()

      ! Initial state
      istate_lz = initstate_lz

      allocate (tocalc_lz(nstate_lz))

      !Which gradients to compute
      do ist1 = 1, nstate_lz
         if (ist1 == istate_lz) then
            tocalc_lz(ist1) = 1
         else
            tocalc_lz(ist1) = 0
         end if
      end do

      !Allocate energy arrays
      allocate (en_array_lz(nstate_lz, 4), en_array_lz_backup(nstate_lz, 4)) !last 3 energies (1: current, 2: n-1, 3: n-3)
      allocate (fx_old(natom, nwalk), fy_old(natom, nwalk), fz_old(natom, nwalk))
      allocate (px_temp(natom, nwalk), py_temp(natom, nwalk), pz_temp(natom, nwalk))
      allocate (x_prev(natom, nwalk), y_prev(natom, nwalk), z_prev(natom, nwalk), &
                vx_prev(natom, nwalk), vy_prev(natom, nwalk), vz_prev(natom, nwalk))
      en_array_lz = 0.0D0

      ! TERA-MPI parameters
      if (pot == '_tera_') then
         call lz_init_terash()
      end if

   end subroutine lz_init

   ! TODO: Break coupling between LZ and SH modules
   ! The following function is needed for TERA-MPI interface
   subroutine lz_init_terash()
      use mod_general, only: natom
      use mod_sh, only: inac, en_array, tocalc, nacx, nacy, nacz
      use mod_sh_integ, only: nstate

      if (ntriplet_lz > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Landau-Zener with Singlet-Triplet transitions not implemented with TeraChem interface')
      end if

      nstate = nstate_lz !Needed in init_terash
      inac = 2 !Turns off couplings calculation

      !Based on sh_init() routine, sharing most of the functions
      !TODO: Break dependence - separation of MPI interface needed
      allocate (en_array(nstate_lz))
      allocate (nacx(natom, nstate_lz, nstate_lz))
      allocate (nacy(natom, nstate_lz, nstate_lz))
      allocate (nacz(natom, nstate_lz, nstate_lz))
      allocate (tocalc(nstate, nstate))
      en_array = 0.0D0
      tocalc = 0
      tocalc(istate_lz, istate_lz) = 1
      call set_current_state(istate_lz)
   end subroutine lz_init_terash

   ! This routine is used in the file interface defined in force_abin.F90
   subroutine write_lz_data()
      integer :: ist, u

      open (newunit=u, file='state.dat')
      write (u, '(I0)') nstate_lz

      ! Print for which state we need gradient
      ! First we have singlets, then triplets
      do ist = 1, nstate_lz
         if (tocalc_lz(ist) == 1) write (u, '(I0)') ist
      end do
      write (u, '(I0)') nsinglet_lz
      write (u, '(I0)') ntriplet_lz
      close (u)
   end subroutine write_lz_data

   !LZ singlets hop
   subroutine lz_hop(x, y, z, vx, vy, vz, fxc, fyc, fzc, amt, dt, eclas, chpot)
      use mod_const, only: ANG, AUTOFS, PI, AUTOEV, AUTOCM
      use mod_files, only: UPOP, UPROB, UPES
      use mod_general, only: natom, pot, nwrite, it, sim_time
      use mod_random, only: vranf
      use mod_kinetic, only: ekin_v
      use mod_interfaces, only: force_clas
      use mod_system, only: am
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(inout) :: eclas
      real(DP), intent(in) :: dt
      character(len=*), intent(in) :: chpot

      real(DP) :: prob(nstate_lz), prob2(nstate_lz), probc
      real(DP) :: ran(10)
      real(DP) :: en_diff(4), second_der, second_der_back, der_check, soc_matrix(nsinglet_lz, ntriplet_lz)
      integer :: ihop, icross, ihist, ist1, iat, ibeg, iend !, istatus
      integer :: S_to_T !itest, iost
      integer :: ist ! =istate_lz
      real(DP) :: Ekin, Epot, Epot2, dE, hop_rdnum, hop_rdnum2, vel_rescale, stepfs, molveloc, grad_diff
      real(DP), parameter :: one = 1.0D0
      character(len=100) :: formt, fmt_in, fmt_out

      !TODO: energy drift
      call check_energy_lz(en_array_lz, istate_lz, energydifthr_lz)
      !call check_energydrift(vx, vy, vz)

      !Current state
      ist = istate_lz

      !========ADIABATIC (same spin) LZ HOPPING========!
      !Compute probabilities for same spin as current state (s-s or t-t allowed only)
      !(But we should compute singlet+triplet probs together and decide in the END)

      prob = 0.0D0
      prob2 = 0.0D0
      probc = 0.0D0

      if (ist <= nsinglet_lz) then !Singlet 1->nsinglet_lz
         ibeg = 1
         iend = nsinglet_lz
      else !Triplet nsinglet_lz+1->nsinglet_lz+ntriplet_lz
         ibeg = nsinglet_lz + 1
         iend = nstate_lz
      end if

      do ist1 = ibeg, iend
         if (ist1 == ist) cycle
         if (ist1 == ignore_state) cycle
         ! only closest upper and lower states are considered for hopping
         if (ist1 > (ist + 1) .or. ist1 < (ist - 1)) cycle

         do ihist = 1, 4
            en_diff(ihist) = abs(en_array_lz(ist, ihist) - en_array_lz(ist1, ihist))
         end do

         ! the energy condition is considered here, otherwise it messes up probabilities
         if (abs(en_diff(2) * AUTOEV) > deltaE_lz) cycle

         ! Three point minima of adiabatic splitting Zjk
         if ((en_diff(1) > en_diff(2)) .and. (en_diff(2) < en_diff(3)) .and. (it > 2)) then
            second_der = ((en_diff(3) - 2 * en_diff(2) + en_diff(1)) / dt**2)
            prob(ist1) = exp(-PI / 2 * (dsqrt(en_diff(2)**3 / second_der)))
            write (fmt_in, '(I2.2)') ist
            write (fmt_out, '(I2.2)') ist1
            write (stdout, *) "Three-point minimum (", trim(fmt_in), "->", trim(fmt_out), &
               ") dE/a.u.", en_diff(2), "Probability:", prob(ist1)
            if (prob(ist1) > 1) then
               write (stderr, *) prob
               call fatal_error(__FILE__, __LINE__, 'LZ probability > 1')
            end if

            ! Adding check for false 3P-minima
            ! only for significant probabilities, for tiny probabilities this does not make sense
            ! e.g. for parallel states (the whole LZ does not make sense for this case)
            if (prob(ist1) > 0.01D0) then
               ! Calculating backward second derivative formula. We are not using the newest energy
               ! but the previous three. Discontinuity would not affect this formula.
               second_der_back = ((en_diff(2) - 2 * en_diff(3) + en_diff(4)) / dt**2)
               ! We have central and backward second derivative formulas and compare them.
               ! If the change is too large, we either almost hit the CI or we have discontinuity.
               der_check = abs((second_der - second_der_back) / second_der)
               ! If they differ by more then 130%, we have aa unphysical change of curvature --> certain discontinuity.
               if (der_check > 1.3D0) then
                  write (stdout, *) "ERROR: Change of curvature --> discontinuity in PES!"
                  write (stdout, *) "Probability set to 0!"
                  prob(ist1) = 0.0D0
                  ! 30% threshold was set empirically and should capture most discontinuities
                  ! yet it can also be a conical intersection. Thus, we just issue an warning
                  ! and let the user to evaluate on his own.
               else if (der_check > 0.3D0) then
                  write (stdout, *) "WARNING: Possible discontinuity in PES! Check PES.dat!"
               end if
            end if

         end if
      end do

      ! LZ warning
      if (sum(prob) > 1) then
         write (stdout, *) "WARNING: Sum of hopping probabilities > 1. Breakdown of LZ assumptions"
      end if

      !Hop?
      ihop = 0
      call vranf(ran, 1)
      hop_rdnum = ran(1)

      ! Determine, whether we hopped or not
      do ist1 = ibeg, iend
         if (ist1 == ist) cycle
         if (hop_rdnum < prob(ist1)) then
            prob2(ist1) = prob(ist1) !All probable hops
         else if (prob(ist1) > 0) then
            write (fmt_in, '(I2.2)') ist1
            write (stdout, *) "NO hop on state ", trim(fmt_in), ", Random n:", hop_rdnum
         end if
      end do

      ! Determine on which state (weighted sampling of all probable hops)
      call vranf(ran, 1)
      hop_rdnum2 = ran(1)
      if (sum(prob2) /= 0.0D0) then
         do ist1 = ibeg, iend
            if (ist1 == ist) cycle
            probc = probc + prob2(ist1) / sum(prob2) !1 if only one state considerable
            if (hop_rdnum2 < probc .and. ihop == 0) then
               ihop = ist1
            else if (prob2(ist1) > 0) then
               write (fmt_in, '(I2.2)') ist1
               write (stdout, *) "NO hop on state ", trim(fmt_in), ", Random n2:", hop_rdnum2
            end if
         end do
      end if

      if (ihop /= 0) then
         Ekin = ekin_v(vx, vy, vz)
         dE = (en_array_lz(ihop, 2) - en_array_lz(ist, 2))
         Epot = abs(en_array_lz(ist, 1) - en_array_lz(ihop, 1)) ! This is not used anywhere in the code
         Epot2 = abs(en_array_lz(ist, 2) - en_array_lz(ihop, 2))
         write (fmt_in, '(I2.2)') istate_lz
         write (fmt_out, '(I2.2)') ihop
         ! Energy conservation criteria
         if ((dE < Ekin) .and. (abs(Epot2 * AUTOEV) < deltaE_lz)) then
            !HOP
            write (stdout, *) "Adiabatic HOP! (", trim(fmt_in), "->", trim(fmt_out), ") dE/a.u.", Epot2, &
               "Random n:", hop_rdnum
            !We need to get to previous geometry, adjust its velocity according to
            !target state and do 1 step forward on the new state
            istate_lz = ihop
            if (pot == '_tera_') call set_current_state(ihop) !TERA-MPI TODO: refactor, separate lz from sh for _tera_ usage completely
            !a) Previous geometry
            x = x_prev; y = y_prev; z = z_prev; 
            !b) Previous velocities
            do iat = 1, natom
               vx(iat, 1) = vx_prev(iat, 1)
               vy(iat, 1) = vy_prev(iat, 1)
               vz(iat, 1) = vz_prev(iat, 1)
            end do
            !c) Gradient for new state
            do ist1 = 1, nstate_lz
               if (ist1 == istate_lz) then
                  tocalc_lz(ist1) = 1
               else
                  tocalc_lz(ist1) = 0
               end if
            end do
            en_array_lz_backup = en_array_lz
            call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

            !d) Simple velocity rescaling (https://doi.org/10.1063/1.4882073)
            vel_rescale = dsqrt(1 - (dE / Ekin))
            do iat = 1, natom
               vx(iat, 1) = vx(iat, 1) * vel_rescale
               vy(iat, 1) = vy(iat, 1) * vel_rescale
               vz(iat, 1) = vz(iat, 1) * vel_rescale
            end do

            !d) Do a verlet step forward from previous position
            !verletstep() function in mdstep.f90 is not declared yet, expressing explicitly
            px_temp = amt * vx; py_temp = amt * vy; pz_temp = amt * vz; 
            px_temp = px_temp + dt / 2 * fxc; py_temp = py_temp + dt / 2 * fyc; pz_temp = pz_temp + dt / 2 * fzc; 
            x = x + dt * px_temp / amt; y = y + dt * py_temp / amt; z = z + dt * pz_temp / amt; 
            call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
            px_temp = px_temp + dt / 2 * fxc; py_temp = py_temp + dt / 2 * fyc; pz_temp = pz_temp + dt / 2 * fzc; 
            vx = px_temp / amt; vy = py_temp / amt; vz = pz_temp / amt; 
            !e) Update only last energy history entry
            en_array_lz(:, 2) = en_array_lz_backup(:, 2)
            en_array_lz(:, 3) = en_array_lz_backup(:, 3)
            en_array_lz(:, 4) = en_array_lz_backup(:, 4)

         else
            write (stdout, *) "NO adiabatic HOP (", trim(fmt_in), "->", trim(fmt_out), ")dE/a.u.", Epot2, &
               "Probability:", prob(ihop), "Random n:", hop_rdnum
         end if

      end if

      !========DIABATIC (singlet x triplet) LZ HOPPING========!
      !formula from article: Spin-forbidden and Spin-allowed Cyclopropenone (c-H2C3O) Formation in Interstellar Medium
      if (ntriplet_lz > 0) then
         ist = istate_lz
         prob = 0.0D0

         Ekin = ekin_v(vx, vy, vz)
         molveloc = dsqrt(2.0D0 * Ekin / sum(am(1:natom)))
         call vranf(ran, 1)
         hop_rdnum = ran(1)
         icross = 0

         !Compute if the current state is crossing a different multiplicity state?
         if (ist > nsinglet_lz) then
            ibeg = 1
            iend = nsinglet_lz
            S_to_T = 0
         else
            ibeg = nsinglet_lz + 1
            iend = nstate_lz
            S_to_T = 1
         end if

         do ist1 = ibeg, iend
            do ihist = 1, 2
               en_diff(ihist) = (en_array_lz(ist, ihist) - en_array_lz(ist1, ihist))
            end do
            ! Crossing of Zjk - the enerrgy differrence changes sign
            if ((((en_diff(1) > 0) .and. (en_diff(2) < 0)) .or. ((en_diff(1) < 0) .and. (en_diff(2) > 0))) .and. it > 1) then
               !We need SOC and gradient of the other state
               call lz_getsoc(soc_matrix, chpot)
               !Which model to use?
               !1 - lets compute it from gradients
               !call lz_getgraddiff(grad_diff, ist, ist1, x, y, z, vx, vy, vz, fxc, fyc, fzc)
               !2 - predict it from last 2 energy differences
               grad_diff = abs((en_diff(1) - en_diff(2)) / dt) !grad diff will contain d(wT-wS)
               write (fmt_in, '(I2.2)') ist
               write (fmt_out, '(I2.2)') ist1
               if (S_to_T == 1) then
                  write (stdout, *) "S/T CROSS (", trim(fmt_in), "->", trim(fmt_out), &
                    &  ") SOC^2: ", soc_matrix(ist, ist1 - nsinglet_lz)
                  prob(ist1) = 1 - exp(-2 * PI * ((soc_matrix(ist, ist1 - nsinglet_lz) / (AUTOCM**2)) / grad_diff))
                  write (stdout, *) "Hop probability: ", prob(ist1)
               else if (S_to_T == 0) then
                  write (stdout, *) "T/S CROSS (", trim(fmt_in), "->", trim(fmt_out), &
                     & ") SOC^2: ", soc_matrix(ist1, ist - nsinglet_lz)
                  prob(ist1) = 1 - exp(-2 * PI * ((soc_matrix(ist1, ist - nsinglet_lz) / (AUTOCM**2)) / grad_diff))
                  write (stdout, *) "Hop probability: ", prob(ist1)
               end if
            end if
         end do

         ! Determine, whether we hopped or not
         do ist1 = ibeg, iend
            if (ist1 == ist) cycle
            if (hop_rdnum < prob(ist1)) then
               icross = ist1
               exit
            end if
         end do

         if (icross /= 0) then
            !We do not hop exactly in the crossing point, minor velocity adjustment
            !HOP
            Epot2 = abs(en_array_lz(istate_lz, 2) - en_array_lz(icross, 2))
            write (fmt_in, '(I2.2)') istate_lz
            write (fmt_out, '(I2.2)') icross
            write (stdout, *) "Nonadiabatic (S/T) HOP! (", trim(fmt_in), "->", trim(fmt_out), ")dGrad/a.u.", grad_diff, &
               "MolVeloc", molveloc, "Probability:", prob(icross), "Random n:", hop_rdnum
            istate_lz = icross
            !Get new forces
            do iat = 1, natom
               fx_old(iat, 1) = fxc(iat, 1)
               fy_old(iat, 1) = fyc(iat, 1)
               fz_old(iat, 1) = fzc(iat, 1)
            end do
            !Gradient to compute
            do ist1 = 1, nstate_lz
               if (ist1 == istate_lz) then
                  tocalc_lz(ist1) = 1
               else
                  tocalc_lz(ist1) = 0
               end if
            end do
            en_array_lz_backup = en_array_lz
            call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
            en_array_lz = en_array_lz_backup
            !Adjust velocities from previous step to new forces
            do iat = 1, natom
               vx(iat, 1) = vx(iat, 1) + (dt / 2.0D0) * (-fx_old(iat, 1) + fxc(iat, 1)) / amt(iat, 1)
               vy(iat, 1) = vy(iat, 1) + (dt / 2.0D0) * (-fy_old(iat, 1) + fyc(iat, 1)) / amt(iat, 1)
               vz(iat, 1) = vz(iat, 1) + (dt / 2.0D0) * (-fz_old(iat, 1) + fzc(iat, 1)) / amt(iat, 1)
            end do

            !Simple velocity scaling
            !vel_rescale = dsqrt(1-(Epot/Ekin))
            !do iat=1, natom
            !    vx(iat,1) = vx(iat,1) * vel_rescale
            !    vy(iat,1) = vy(iat,1) * vel_rescale
            !    vz(iat,1) = vz(iat,1) * vel_rescale
            !end do

         end if

      end if

      !Archive
      x_prev = x; y_prev = y; z_prev = z; vx_prev = vx; vy_prev = vy; vz_prev = vz; 
      !Write
      if (modulo(it, nwrite) == 0) then
         stepfs = sim_time * AUtoFS

         write (formt, '(A7,I3,A7)') '(F15.2,', nstate_lz, 'E20.10)'
         write (UPES, fmt=formt) stepfs, (en_array_lz(ist1, 1), ist1=1, nstate_lz)

         write (formt, '(A10,I3,A13)') '(F15.2,I3,', nstate_lz, 'F10.5,1F10.7)'
         write (UPOP, fmt=formt) stepfs, istate_lz, (real(tocalc_lz(ist1)), ist1=1, nstate_lz), one
      end if

   end subroutine lz_hop

   subroutine lz_getsoc(soc_matrix, chpot)
      use mod_utils, only: toupper

      real(DP), intent(inout) :: soc_matrix(:, :)
      character(len=*), intent(in) :: chpot

      integer :: iost, istatus, icmd, row, col
      integer :: soc_unit
      character(len=*), parameter :: chSOC = 'SOC.dat'
      character(len=100) :: chsystem
      logical :: file_exists

      !Perform SOC matrix computation
      chsystem = './'//trim(toupper(chpot))//'/r.'//trim(chpot)//'.soc' !r.pot.soc

      inquire (FILE=chsystem, EXIST=file_exists)
      if (.not. file_exists) then
         call fatal_error(__FILE__, __LINE__, &
            & 'File '//trim(chsystem)//' does not exist!')
      end if

      write (chsystem, '(A,I4.3)') trim(chsystem)//" ", 1 !iw
      chsystem = trim(chsystem)//' < ./state.dat'
      call execute_command_line(trim(chsystem), exitstat=istatus, cmdstat=icmd)
      if (icmd /= 0 .or. istatus /= 0) then
         write (stderr, *) 'ERROR: Something went wrong during the execution of the ab initio external program.'
         write (stderr, *) 'Inspect the output files in folder '//trim(toupper(chpot))//"/"
         write (stderr, *) 'CALL:', chsystem
         call fatal_error(__FILE__, __LINE__, 'Could not compute SOC matrix')
      end if

      ! make sure that the file exist and flush the disc buffer
      inquire (FILE=chSOC, EXIST=file_exists)
      if (.not. file_exists) then
         write (stderr, *) 'WARNING:File ', chSOC, ' does not exist. Waiting..'
         call execute_command_line('sync')
         call execute_command_line('sleep 0.5')
      end if

      open (newunit=soc_unit, file=chSOC, status='old', ACTION='READ', IOSTAT=iost)
      if (iost /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Could not open file '//trim(chSOC))
      end if
      !Read the SOC values H^2 [cm^-1]
      !     t1   t2
      ! s1  ..   ..
      ! s2  ..   ..
      do row = 1, nsinglet_lz
         read (soc_unit, *) (soc_matrix(row, col), col=1, ntriplet_lz)
      end do
      ! TODO: Should this file be deleted here?
      close (soc_unit)
   end subroutine lz_getsoc

   subroutine check_energy_lz(en_array_lz, istate_lz, energydifthr_lz)
      use mod_const, only: AUTOEV
      use mod_general, only: it

      real(DP), intent(in) :: en_array_lz(:, :)
      integer, intent(in) :: istate_lz
      real(DP), intent(in) :: energydifthr_lz

      ! Did the PES changed abruptly? LZ might induce unphysical hops.
      if ((abs(en_array_lz(istate_lz, 1) - en_array_lz(istate_lz, 2)) * AUTOEV) >= energydifthr_lz .and. it > 1) then
         write (stderr, *) "!!!QM energy of current state changed more than", energydifthr_lz, "eV."
         write (stderr, *) "State ", istate_lz, ", change: ", &
                    & (abs(en_array_lz(istate_lz, 1) - en_array_lz(istate_lz, 2))) * AUTOEV, " eV"
         write (stderr, *) "LZ unstable. Terminating."
         call fatal_error(__FILE__, __LINE__, 'Insufficient energy conservation in LZ')
      end if

   end subroutine check_energy_lz

   subroutine lz_rewind(en_array_lz)
      !Routine for deleting interim calculations from LZ energy history (e.g. initialization of forces)
      real(DP), intent(inout) :: en_array_lz(:, :)

      en_array_lz(:, 1) = en_array_lz(:, 2)
      en_array_lz(:, 2) = en_array_lz(:, 3)
   end subroutine lz_rewind

#if 0
   subroutine lz_getgraddiff(grad_diff, ist, ist1, x, y, z, vx, vy, vz, fxc, fyc, fzc)
      !Computes gradient difference * velocity to predict dEgap/dt
      use mod_general, only: natom, pot
      use mod_interfaces, only: force_clas

      real(DP), intent(inout) :: grad_diff
      integer, intent(in) :: ist, ist1
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP) :: eclas
      integer :: i, iat

      eclas = 0.0D0

      !Gradient to compute
      do i = 1, nstate_lz
         if (i == ist1) then
            tocalc_lz(i) = 1
         else
            tocalc_lz(i) = 0
         end if
      end do

      call force_clas(fx_old, fy_old, fz_old, x, y, z, eclas, pot) !state2

      !Gradient difference of two states
      !predicting dE change
      grad_diff = 0.0D0
      do iat = 1, natom
         grad_diff = grad_diff - vx(iat, 1) * (fx_old(iat, 1) - fxc(iat, 1))
         grad_diff = grad_diff - vy(iat, 1) * (fy_old(iat, 1) - fyc(iat, 1))
         grad_diff = grad_diff - vz(iat, 1) * (fz_old(iat, 1) - fzc(iat, 1))
      end do
      grad_diff = abs(grad_diff)

      !grad_diff = 0.0d0
      !do iat=1, natom
      !   grad_diff = grad_diff + abs(fx_old(iat,1)-fxc(iat,1))
      !   grad_diff = grad_diff + abs(fy_old(iat,1)-fyc(iat,1))
      !   grad_diff = grad_diff + abs(fz_old(iat,1)-fzc(iat,1))
      !end do

      !Back to original
      do i = 1, nstate_lz
         if (i == ist) then
            tocalc_lz(i) = 1
         else
            tocalc_lz(i) = 0
         end if
      end do

   end subroutine lz_getgraddiff
#endif

   subroutine lz_restout(fileunit)
      integer, intent(in) :: fileunit
      integer :: ist

      write (fileunit, *) istate_lz
      do ist = 1, nstate_lz
         write (fileunit, '(3ES25.16E3)') en_array_lz(ist, 1), en_array_lz(ist, 2), en_array_lz(ist, 3), en_array_lz(ist, 4)
      end do

   end subroutine lz_restout

   subroutine lz_restin(fileunit, x, y, z, vx, vy, vz)
      use mod_general, only: pot
      integer, intent(in) :: fileunit
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer :: ist

      read (fileunit, *) istate_lz
      !Reset which gradients to compute
      do ist = 1, nstate_lz
         if (ist == istate_lz) then
            tocalc_lz(ist) = 1
         else
            tocalc_lz(ist) = 0
         end if
      end do

      !TERA-MPI parameters
      if (pot == '_tera_') then
         call set_current_state(istate_lz)
      end if

      do ist = 1, nstate_lz
         read (fileunit, *) en_array_lz(ist, 1), en_array_lz(ist, 2), en_array_lz(ist, 3), en_array_lz(ist, 4)
      end do

      !Store previous step values
      x_prev = x; y_prev = y; z_prev = z; vx_prev = vx; vy_prev = vy; vz_prev = vz; 
   end subroutine lz_restin

   subroutine lz_finalize()
      ! Deallocate arrays
      if (allocated(tocalc_lz)) then
         deallocate (tocalc_lz, en_array_lz, en_array_lz_backup, fx_old, fy_old, fz_old)
         deallocate (px_temp, py_temp, pz_temp, x_prev, y_prev, z_prev, vx_prev, vy_prev, vz_prev)
      end if
   end subroutine lz_finalize

end module mod_lz
