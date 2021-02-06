! Driver routines for Landau-Zener excited state dynamics
!    Based on PHYSICAL REVIEW A 84, 014701 (2011)
!    Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
!    Andrey K. Belyaev and Oleg V. Lebedev
! Implemented by: J. Suchan, (J. Chalabala)

! No S/T TeraChem functionality YET
! No energy drift check
! No PIMD functionality

module mod_lz
   use mod_const, only: DP
   use mod_utils, only: abinerror
   use mod_general, only: nwalk, pot, irest
   use mod_array_size, only: NSTMAX, NTRAJMAX
   use mod_sh, only: istate_init, istate, inac !TERA-MPI interface
   use mod_sh_integ, only: nstate
   implicit none
   !private
   public :: lz_init, lz_hop, lz_rewind, lz_restin, lz_restout, lz_finalize !Routines
   public :: initstate_lz, nstate_lz, nsinglet_lz, ntriplet_lz, deltaE_lz, energydifthr_lz !User defined variables
   public :: en_array_lz, tocalc_lz, istate_lz !Routine variables
   !Caveat: Every time we call force_class en_array_lz is updated

   real(DP) :: deltaE_lz = 1.0D0 !Energy difference, up to which we consider possibility of LZ hops
   real(DP) :: energydifthr_lz = 0.5D0 !Maximum energy difference (eV) between two consecutive steps, used in check_energy_lz()
   ! Initial electronic state
   integer :: initstate_lz = 1
   integer :: tocalc_lz(NSTMAX)
   integer :: nstate_lz, nsinglet_lz = 0, ntriplet_lz = 0
   integer :: calcsoc_lz = 0

   !Module variables
   integer :: istate_lz
   real(DP), allocatable :: en_array_lz(:, :), en_array_lz_backup(:, :)
   real(DP), allocatable :: fx_old(:, :), fy_old(:, :), fz_old(:, :)
   real(DP), allocatable :: px_temp(:, :), py_temp(:, :), pz_temp(:, :)
   real(DP), allocatable :: x_prev(:, :), y_prev(:, :), z_prev(:, :), &
                            vx_prev(:, :), vy_prev(:, :), vz_prev(:, :)
   save

contains

   !Initialization
   subroutine lz_init()
      use mod_general, only: natom
      integer :: ist1

      !Initial state
      istate_lz = initstate_lz

      !Which gradients to compute
      do ist1 = 1, nstate_lz
         if (ist1 == istate_lz) then
            tocalc_lz(ist1) = 1
         else
            tocalc_lz(ist1) = 0
         end if
      end do

      !Allocate energy arrays
      allocate (en_array_lz(nstate_lz, 3), en_array_lz_backup(nstate_lz, 3)) !last 3 energies (1: current, 2: n-1, 3: n-3)
      allocate (fx_old(natom, nwalk + 1), fy_old(natom, nwalk + 1), fz_old(natom, nwalk + 1))
      allocate (px_temp(natom, nwalk + 1), py_temp(natom, nwalk + 1), pz_temp(natom, nwalk + 1))
      allocate (x_prev(natom, nwalk + 1), y_prev(natom, nwalk + 1), z_prev(natom, nwalk + 1), &
                vx_prev(natom, nwalk + 1), vy_prev(natom, nwalk + 1), vz_prev(natom, nwalk + 1))
      en_array_lz = 0.0D0

      !TERA-MPI parameters
      if (pot == '_tera_') then
         nstate = nstate_lz
         istate_init = initstate_lz
         inac = 2
      end if

   end subroutine lz_init

   !LZ singlets hop
   subroutine lz_hop(x, y, z, vx, vy, vz, fxc, fyc, fzc, amt, dt, eclas, chpot)
      use mod_const, only: ANG, AUTOFS, PI, AUTOEV, AUTOCM
      use mod_files, only: UPOP, UPROB, UPES
      use mod_general, only: natom, pot, nwrite, it, sim_time
      use mod_random, only: vranf
      use mod_kinetic, only: ekin_v
      use mod_interfaces, only: force_clas
      use mod_utils, only: abinerror
      !use mod_files,    ONLY: MAXUNITS
      use mod_system, only: am
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(inout) :: eclas
      real(DP), intent(in) :: dt
      character(len=*), intent(in) :: chpot

      real(DP) :: prob(NSTMAX), prob2(NSTMAX), probc
      real(DP) :: ran(10)
      real(DP) :: en_diff(3), second_der, soc_matrix(nsinglet_lz, ntriplet_lz)
      integer :: ihop, icross, ihist, ist1, iat, ibeg, iend !, istatus
      integer :: S_to_T !itest, iost
      integer :: ist ! =istate_lz
      real(DP) :: Ekin, Epot, Epot2, dE, hop_rdnum, hop_rdnum2, vel_rescale, stepfs, molveloc, grad_diff
      real(DP) :: one = 1.0D0
      character(len=100) :: formt, fmt_in, fmt_out !, chSOC='SOC.dat'

      !TODO: energy drift
      call check_energy_lz(en_array_lz, istate_lz, energydifthr_lz)
      !call check_energydrift(vx, vy, vz, itrj)

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
         do ihist = 1, 3
            en_diff(ihist) = abs(en_array_lz(ist, ihist) - en_array_lz(ist1, ihist))
         end do
         ! Three point minima of adiabatic splitting Zjk
         if ((en_diff(1) > en_diff(2)) .and. (en_diff(2) < en_diff(3)) .and. (it > 2)) then
            second_der = ((en_diff(3) - 2 * en_diff(2) + en_diff(1)) / dt**2)
            prob(ist1) = exp(-PI / 2 * (sqrt(en_diff(2)**3 / second_der)))
            write (fmt_in, '(I2.2)') ist
            write (fmt_out, '(I2.2)') ist1
            write (*, *) "Three-point minimum (", trim(fmt_in), "->", trim(fmt_out), &
               ") dE/a.u.", en_diff(2), "Probability:", prob(ist1)
            if (prob(ist1) > 1) then
               call abinerror('landau_zener_prob')
            end if

         end if
      end do
      !write(*,*) "diff1",en_diff(1),"diff2",en_diff(2),"diff3",en_diff(3)

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
            write (*, *) "NO hop on state ", trim(fmt_in), ", Random n:", hop_rdnum
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
               write (*, *) "NO hop on state ", trim(fmt_in), ", Random n2:", hop_rdnum2
            end if
         end do
      end if

      if (ihop /= 0) then
         Ekin = ekin_v(vx, vy, vz)
         dE = (en_array_lz(ihop, 2) - en_array_lz(ist, 2))
         Epot = abs(en_array_lz(ist, 1) - en_array_lz(ihop, 1))
         Epot2 = abs(en_array_lz(ist, 2) - en_array_lz(ihop, 2))
         write (fmt_in, '(I2.2)') istate_lz
         write (fmt_out, '(I2.2)') ihop
         ! Energy conservation criteria
         if ((dE < Ekin) .and. (abs(Epot2 * AUTOEV) < deltaE_lz)) then
            !HOP
            write (*, *) "Adiabatic HOP! (", trim(fmt_in), "->", trim(fmt_out), ") dE/a.u.", Epot2, &
               "Random n:", hop_rdnum
            !We need to get to previous geometry, adjust its velocity according to
            !target state and do 1 step forward on the new state
            istate_lz = ihop
            if (pot == '_tera_') istate = ihop !TERA-MPI TODO: refactor, separate lz from sh for _tera_ usage completely
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
            vel_rescale = sqrt(1 - (dE / Ekin))
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

         else
            write (*, *) "NO adiabatic HOP (", trim(fmt_in), "->", trim(fmt_out), ")dE/a.u.", Epot2, &
               "Probability:", prob(ihop), "Random n:", hop_rdnum
         end if

      end if

      !========DIABATIC (singlet x triplet) LZ HOPPING========!
      !formula from article: Spin-forbidden and Spin-allowed Cyclopropenone (c-H2C3O) Formation in Interstellar Medium
      if (ntriplet_lz > 0) then
         ist = istate_lz
         prob = 0.0D0

         Ekin = ekin_v(vx, vy, vz)
         molveloc = sqrt(2.0D0 * Ekin / sum(am(1:natom)))
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
                  write (*, *) "S/T CROSS (", trim(fmt_in), "->", trim(fmt_out), ") SOC^2: ", soc_matrix(ist, ist1 - nsinglet_lz)
                  prob(ist1) = 1 - exp(-2 * PI * ((soc_matrix(ist, ist1 - nsinglet_lz) / (AUTOCM**2)) / grad_diff))
                  write (*, *) "Hop probability: ", prob(ist1)
               else if (S_to_T == 0) then
                  write (*, *) "T/S CROSS (", trim(fmt_in), "->", trim(fmt_out), ") SOC^2: ", soc_matrix(ist1, ist - nsinglet_lz)
                  prob(ist1) = 1 - exp(-2 * PI * ((soc_matrix(ist1, ist - nsinglet_lz) / (AUTOCM**2)) / grad_diff))
                  write (*, *) "Hop probability: ", prob(ist1)
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
            write (*, *) "Nonadiabatic (S/T) HOP! (", trim(fmt_in), "->", trim(fmt_out), ")dGrad/a.u.", grad_diff, &
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
            !vel_rescale = sqrt(1-(Epot/Ekin))
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
      use mod_utils, only: abinerror, lowertoupper
      use mod_files, only: MAXUNITS

      real(DP), intent(inout) :: soc_matrix(:, :)
      character(len=*), intent(in) :: chpot

      integer :: itest, iost, ISTATUS, system, row, col
      character(len=20) :: chSOC = 'SOC.dat'
      character(len=100) :: chsystem
      logical :: file_exists

      !Perform SOC matrix computation
      chsystem = './'//trim(LowerToUpper(chpot))//'/r.'//trim(chpot)//'.soc' !r.pot.soc

      inquire (FILE=chsystem, EXIST=file_exists)
      if (.not. file_exists) then
         write (*, *) 'File ', chsystem
         write (*, *) 'does not exist! Exiting...'
         call abinerror('force_abin')
      end if

      write (chsystem, '(A30,I4.3)') chsystem, 1 !iw
      write (chsystem, '(A40,A14)') chsystem, ' < ./state.dat'
      !CALL TO SYSTEM
      ISTATUS = system(chsystem)
      if (ISTATUS /= 0 .and. ISTATUS /= 256) then
         write (*, *) 'ERROR: Something went wrong during the execution of the ab initio external program.'
         write (*, *) 'See the approprite output files in&
         & folder '//trim(LowerToUpper(chpot))//"/"
         write (*, *) 'CALL:', chsystem
         call abinerror('force_abin')
      end if

      !make sure that the file exist and flush the disc buffer
      itest = 0
      inquire (FILE=chSOC, EXIST=file_exists)
      do while (.not. file_exists .and. itest < 10)
         write (*, *) 'WARNING:File ', chSOC, ' does not exist. Waiting..'
         ISTATUS = system('sync') !mel by zajistit flush diskoveho bufferu
         !call system('sync')
         inquire (FILE=chSOC, EXIST=file_exists)
         itest = itest + 1
      end do

      open (unit=MAXUNITS + 100, file=chSOC, status='old', ACTION='READ', IOSTAT=iost)
      if (iost /= 0) then
         write (*, *) 'Fatal problem when trying to open the file ', chSOC
         call abinerror('force_abin')
      end if
      !Read the SOC values H^2 [cm^-1]
      !     t1   t2
      ! s1  ..   ..
      ! s2  ..   ..
      do row = 1, nsinglet_lz
         read (MAXUNITS + 100, *) (soc_matrix(row, col), col=1, ntriplet_lz)
      end do
      close (unit=MAXUNITS + 100)
   end subroutine lz_getsoc

   subroutine check_energy_lz(en_array_lz, istate_lz, energydifthr_lz)
      use mod_utils, only: abinerror
      use mod_const, only: AUTOEV
      use mod_general, only: it

      real(DP), intent(in) :: en_array_lz(:, :)
      integer, intent(in) :: istate_lz
      real(DP), intent(in) :: energydifthr_lz

      ! Did the PES changed abruptly? LZ might induce unphysical hops.
      if ((abs(en_array_lz(istate_lz, 1) - en_array_lz(istate_lz, 2)) * AUTOEV) >= energydifthr_lz .and. it > 1) then
         write (*, *) "!!!QM energy of current state changed more than", energydifthr_lz, "eV."
         write (*, *) "State ", istate_lz, ", change: ", &
                    & (abs(en_array_lz(istate_lz, 1) - en_array_lz(istate_lz, 2))) * AUTOEV, " eV"
         write (*, *) "LZ unstable. Terminating."
         call abinerror('check_energy_lz')
      end if

   end subroutine check_energy_lz

   subroutine lz_rewind(en_array_lz)
      !Routine for deleting interim calculations from LZ energy history (e.g. initialization of forces)
      real(DP), intent(inout) :: en_array_lz(:, :)

      en_array_lz(:, 1) = en_array_lz(:, 2)
      en_array_lz(:, 2) = en_array_lz(:, 3)
   end subroutine lz_rewind

   subroutine lz_getgraddiff(grad_diff, ist, ist1, x, y, z, vx, vy, vz, fxc, fyc, fzc)
      !Computes gradient difference * velocity to predict dEgap/dt
      use mod_general, only: natom, pot
      use mod_interfaces, only: force_clas

      real(DP), intent(inout) :: grad_diff
      integer, intent(in) :: ist, ist1
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: fxc(:, :), fyc(:, :), fzc(:, :)
      real(DP) :: eclas = 0.0D0
      integer :: i, iat

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

   subroutine lz_restout(fileunit)
      integer, intent(in) :: fileunit
      integer :: ist

      write (fileunit, *) istate_lz
      do ist = 1, nstate_lz
         write (fileunit, *) en_array_lz(ist, 1), en_array_lz(ist, 2), en_array_lz(ist, 3)
      end do

   end subroutine lz_restout

   subroutine lz_restin(fileunit, x, y, z, vx, vy, vz)
      integer, intent(in) :: fileunit
      real(DP), intent(out) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: vx(:, :), vy(:, :), vz(:, :)
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
         istate_init = istate_lz
         istate = istate_lz
      end if

      do ist = 1, nstate_lz
         read (fileunit, *) en_array_lz(ist, 1), en_array_lz(ist, 2), en_array_lz(ist, 3)
      end do

      !Store previous step values
      x_prev = x; y_prev = y; z_prev = z; vx_prev = vx; vy_prev = vy; vz_prev = vz; 
   end subroutine lz_restin

   subroutine lz_finalize()
      ! Deallocate arrays
      deallocate (en_array_lz, en_array_lz_backup, fx_old, fy_old, fz_old)
      deallocate (px_temp, py_temp, pz_temp, x_prev, y_prev, z_prev, vx_prev, vy_prev, vz_prev)
   end subroutine lz_finalize

end module mod_lz
