! Driver routines for Surface Hopping dynamics
! D. Hollas, O. Svoboda, P. Slavicek, M. Oncak
module mod_sh
   use mod_const, only: DP
   use mod_array_size, only: NSTMAX, NTRAJMAX
   use mod_utils, only: abinerror
   use mod_sh_integ
   implicit none
   private
   ! TODO: We should make some of these private
   ! We would need to read sh namelist inside this module
   ! and we should check input sanity here, not in input.F90
   public :: istate_init, substep, deltaE, inac, nohop, decoh_alpha, popthr, nac_accu1, nac_accu2
   public :: surfacehop, sh_init
   public :: istate, ntraj, tocalc, en_array
   public :: nacx, nacy, nacz
   public :: move_vars, get_nacm, write_nacmrest, read_nacmrest
   public :: energydifthr, energydriftthr, dE_S0S1_thr, adjmom, revmom
   public :: check_CIVector
   public :: ignore_state

   integer, parameter :: ntraj = NTRAJMAX ! Currently not in use, set to 1
   ! Initial electronic state
   integer :: istate_init = 1
   ! Number of substeps for integrating electronic Schrodinger eq.
   integer :: substep = 100
   ! Controls calculations of Non-adiabatic Couplings (NAC)
   ! 0 - Analytical NAC
   ! 1 - Numerical Hammers-Schffer-Tully model (currently not implemented)
   ! 2 - Do not compute couplings
   integer :: inac = 0
   ! 1 - Turn OFF hopping
   integer :: nohop = 0
   !
   integer :: adjmom = 0
   ! 1 - Reverse momentum direction after frustrated hop
   integer :: revmom = 0
   ! Numerical accuracy of MOLPRO NAC
   integer :: nac_accu1 = 7
   integer :: nac_accu2 = 5 !7 is MOLPRO default
   ! Decoherence correction parameter (a.u.)
   real(DP) :: decoh_alpha = 0.1D0
   real(DP) :: deltae = 5.0D0, popthr = 0.001D0
   real(DP) :: energydifthr = 1.0D0, energydriftthr = 1.0D0 !eV
   ! Special case for adiabatic TDDFT, terminate when close to S1-S0 crossing
   real(DP) :: dE_S0S1_thr = 0.0D0 !eV
   ! NA Couplings
   real(DP), allocatable :: nacx(:, :, :, :), nacy(:, :, :, :), nacz(:, :, :, :)
   ! *old variables holds data from the previous step
   real(DP), allocatable :: nacx_old(:, :, :, :), nacy_old(:, :, :, :), nacz_old(:, :, :, :)
   real(DP), allocatable :: en_array(:, :), en_array_old(:, :)
   ! Initial absolute electronic energy, needed for monitoring energy drift
   real(DP) :: entot0
   ! nstatexnstate matrix, off-diagonal elements determine NAC, diagonal elements = electronic gradients
   integer, allocatable :: tocalc(:, :)
   integer :: istate(NTRAJMAX)
   ! for ethylene 2-state SA3 dynamics
   integer :: ignore_state = 0
   save

contains

   ! Initialization routine (ugly, but works)
   subroutine sh_init(x, y, z, vx, vy, vz)
      use mod_const, only: AUtoEV
      use mod_general, only: irest, natom, pot, ipimd
      use mod_interfaces, only: force_clas
      use mod_kinetic, only: ekin_v
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: dum_fx(size(x, 1), size(x, 2))
      real(DP) :: dum_fy(size(y, 1), size(y, 2))
      real(DP) :: dum_fz(size(z, 1), size(z, 2))
      real(DP) :: dum_eclas
      integer :: itrj

      deltaE = deltaE / AUtoEV

      allocate (nacx(natom, NTRAJMAX, nstate, nstate))
      allocate (nacy(natom, NTRAJMAX, nstate, nstate))
      allocate (nacz(natom, NTRAJMAX, nstate, nstate))
      allocate (nacx_old(natom, NTRAJMAX, nstate, nstate))
      allocate (nacy_old(natom, NTRAJMAX, nstate, nstate))
      allocate (nacz_old(natom, NTRAJMAX, nstate, nstate))
      nacx = 0.0D0
      nacy = 0.0D0
      nacz = 0.0D0
      nacx_old = nacx
      nacy_old = nacx
      nacz_old = nacx

      allocate (en_array(nstate, ntraj))
      allocate (en_array_old(nstate, ntraj))
      en_array = 0.0D0
      en_array_old = en_array

      ! Determining the initial state
      if (irest /= 1) then
         do itrj = 1, ntraj
            istate(itrj) = istate_init
         end do
      end if

      ! computing only energies, used for subsequent
      ! determination of TOCALC according to deltaE threshold
      allocate (tocalc(nstate, nstate))
      tocalc = 0
      dum_eclas = 0.0D0
      dum_fx = 0.0D0; dum_fy = 0.0D0; dum_fz = 0.0D0

      call force_clas(dum_fx, dum_fy, dum_fz, x, y, z, dum_eclas, pot)

      do itrj = 1, ntraj

         call sh_set_initialwf(istate(itrj), en_array(1, itrj), itrj)

         ! WARNING: entot0 does not respect itrj, same for tocalc
         entot0 = en_array(istate(itrj), itrj) + ekin_v(vx, vy, vz)
         entot0 = entot0 * AUtoEV

         call set_tocalc(itrj)

      end do

      if (irest == 1 .and. ipimd /= 5) then
         call read_nacmrest()
      end if

   end subroutine sh_init

   subroutine get_nacm(itrj)
      use mod_general, only: pot
      integer, intent(in) :: itrj
      integer :: iost, ipom, ist1, ist2

      ! In TeraChem SH interface, we already got NACME
      if (pot == '_tera_') return

      ! Check whether we even want any NACME
      ipom = 0
      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate
            ipom = ipom + tocalc(ist1, ist2)
         end do
      end do

      if (inac == 0 .and. ipom > 0) then
         ! Calculate NACME using default accuracy

         call calc_nacm(itrj, nac_accu1)

         iost = read_nacm(itrj)

         if (iost /= 0 .and. nac_accu1 > nac_accu2) then
!        if NACME NOT COMPUTED: TRY TO DECREASE ACCURACY
            call calc_nacm(itrj, nac_accu2)
            iost = read_nacm(itrj)
         end if

         if (iost /= 0) then
            write (*, *) 'Some NACMEs not read. Exiting...'
            call abinerror('get_nacm')
         end if
         ! we always have to set tocalc because we change it in readnacm
         call set_tocalc(itrj)
      end if
   end subroutine get_nacm

   ! In this routine, we decide which gradients and which NACME we need to compute
   ! TOCALC is a symmetric matrix, upper-triangle defines NACME,
   ! the diagonal defines gradients
   subroutine set_tocalc(itrj)
      ! WARNING: tocalc array is currently the same for all trajectories
      integer, intent(in) :: itrj
      integer :: ist1, ist2
      real(DP) :: pop, pop2

      tocalc = 0

      if (inac /= 2) then ! for ADIABATIC dynamics, do not calculate NACME

         do ist1 = 1, nstate - 1
            do ist2 = ist1 + 1, nstate
               if (abs(en_array(ist1, itrj) - en_array(ist2, itrj)) < deltaE) then
                  tocalc(ist1, ist2) = 1
               end if
            end do
         end do

      end if

      if (popthr > 0) then
         ! COMPUTE NACME only if population of the states is gt.popthr
         do ist1 = 1, nstate - 1
            pop = el_pop(ist1, itrj)
            do ist2 = ist1 + 1, nstate
               pop2 = el_pop(ist2, itrj)
               if (pop < popthr .and. pop2 < popthr .and. ist1 /= istate(itrj) .and. ist2 /= istate(itrj)) then
                  tocalc(ist1, ist2) = 0
               end if
            end do
         end do
      end if

      ! The diagonal holds information about gradients that we need
      ! for SH, we just need the gradient of the current state
      tocalc(istate(itrj), istate(itrj)) = 1
   end subroutine set_tocalc

   subroutine write_nacmrest()
      use mod_general, only: narchive, it
      use mod_qmmm, only: natqm
      use mod_utils, only: archive_file
      integer :: ist1, ist2, iat, itrj
      integer :: iunit1
      logical :: file_exists
      character(len=*), parameter :: restart_file='restart_sh.bin'
      character(len=200) :: chsystem

      inquire (FILE=restart_file, EXIST=file_exists)
      chsystem = 'mv '//trim(restart_file)//'  '//trim(restart_file)//'.old'
      if (file_exists) call system(chsystem)

      open (newunit=iunit1, file=restart_file, action='write', status="new", access="sequential", form="unformatted")

      do itrj = 1, ntraj

         do ist1 = 1, nstate - 1

            do ist2 = ist1 + 1, nstate

               if (inac == 0) then

                  do iat = 1, natqm
                     write (iunit1) nacx(iat, itrj, ist1, ist2), &
                                  & nacy(iat, itrj, ist1, ist2), &
                                  & nacz(iat, itrj, ist1, ist2)
                  end do

               end if

            end do
         end do

         ! Printing total energy at t=0, so that we can safely restart
         ! and we do not break checking for energy drift
         ! For now, energy jump check is not handled well.
         write (iunit1) entot0

         if (phase == 1) then
            call sh_write_phase_bin(iunit1, itrj)
         end if

         ! ntraj enddo
      end do

      close (iunit1)

      if (modulo(it, narchive) == 0) then
         call archive_file(restart_file, it)
      end if

   end subroutine write_nacmrest

   subroutine read_nacmrest()
      use mod_general, only: it
      use mod_qmmm, only: natqm
      use mod_utils, only: archive_file
      character(len=*), parameter :: restart_file='restart_sh.bin'
      integer :: iost, ist1, ist2, iat, itrj
      integer :: iunit1
      logical :: file_exists
      character(len=200) :: chmsg

      print*,'Reading Surface Hopping restart data from file '//trim(restart_file)
      inquire (file=restart_file, exist=file_exists)
      if (.not. file_exists) then
         write (*, *) 'ERROR: Surface Hopping restart file does not exist! '//trim(restart_file)
         call abinerror('read_nacmrest')
      end if

      open (newunit=iunit1, file=restart_file, action="read", status="old", access="sequential", form="unformatted")

      do itrj = 1, ntraj

         do ist1 = 1, nstate - 1

            do ist2 = ist1 + 1, nstate

               if (inac == 0) then

                  do iat = 1, natqm
                     read (iunit1, iomsg=chmsg, IOSTAT=iost) nacx(iat, itrj, ist1, ist2), &
                                                           & nacy(iat, itrj, ist1, ist2), &
                                                           & nacz(iat, itrj, ist1, ist2)
                     if (iost /= 0) then
                        write (*, *) 'Error reading NACME from restart file '//trim(restart_file)
                        write (*, *) chmsg
                        call abinerror('read_nacmrest')
                     end if
                     nacx(iat, itrj, ist2, ist1) = -nacx(iat, itrj, ist1, ist2)
                     nacy(iat, itrj, ist2, ist1) = -nacy(iat, itrj, ist1, ist2)
                     nacz(iat, itrj, ist2, ist1) = -nacz(iat, itrj, ist1, ist2)
                  end do

               end if

            end do
         end do

         ! Reading total energy at t=0, so that we can safely restart
         ! and we do not break checking for energy drift
         ! For now, energy jump check is not handled well.
         read (iunit1) entot0

         if (phase == 1) then
            call sh_read_phase_bin(iunit1, itrj)
         end if

      end do

      close (iunit1)
      call archive_file(restart_file, it)

   end subroutine read_nacmrest

   integer function read_nacm(itrj)
      use mod_qmmm, only: natqm
      integer :: iost, ist1, ist2, iat, itrj, iunit

      iost = 0 ! needed if each tocalc=0
      open (newunit=iunit, file='nacm.dat')
      do ist1 = 1, nstate - 1
         do ist2 = ist1 + 1, nstate

            if (tocalc(ist1, ist2) == 1) then

               do iat = 1, natqm ! reading only for QM atoms
                  read (iunit, *, IOSTAT=iost) nacx(iat, itrj, ist1, ist2), &
                                             & nacy(iat, itrj, ist1, ist2), &
                                             & nacz(iat, itrj, ist1, ist2)
                  if (iost == 0) then
                     tocalc(ist1, ist2) = 0 ! marking as read, useful if we do decreased accuracy
                     nacx(iat, itrj, ist2, ist1) = -nacx(iat, itrj, ist1, ist2)
                     nacy(iat, itrj, ist2, ist1) = -nacy(iat, itrj, ist1, ist2)
                     nacz(iat, itrj, ist2, ist1) = -nacz(iat, itrj, ist1, ist2)
                  else
                     close (iunit, status='delete')
                     write (*, *) 'WARNING: NACME between states', ist1, ist2, 'not read.'
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

   subroutine calc_nacm(itrj, nac_accu)
      use mod_utils, only: toupper
      use mod_general, only: it, pot
      integer, intent(in) :: itrj, nac_accu
      integer :: ist1, ist2, u
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
      ! TODO: move the following line somwhere else
      ! write(*,*)'WARNING: Some NACMs not computed. Trying with decreased accuracy...'
      ! write(*,*)'Calling script r.'//pot//'with accuracy:',nac_accu
      write (chsystem, '(A30,I13,I4.3,I3,A12)') chsystem, it, itrj, nac_accu, ' < state.dat'

      call system(chsystem)

      ! TODO: catch errors here

   end subroutine calc_nacm

   ! move arrays from new step to old step
   subroutine move_vars(vx, vy, vz, vx_old, vy_old, vz_old, itrj)
      use mod_general, only: natom
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(out) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      integer, intent(in) :: itrj
      integer :: ist1, ist2, iat

      do ist1 = 1, nstate
         en_array_old(ist1, itrj) = en_array(ist1, itrj)

         if (inac == 0) then ! nedelame, pokud nacitame rovnou dotprodukt
            do ist2 = 1, nstate
               do iat = 1, natom
                  nacx_old(iat, itrj, ist1, ist2) = nacx(iat, itrj, ist1, ist2)
                  nacy_old(iat, itrj, ist1, ist2) = nacy(iat, itrj, ist1, ist2)
                  nacz_old(iat, itrj, ist1, ist2) = nacz(iat, itrj, ist1, ist2)
               end do
            end do
         end if

      end do

      do iat = 1, natom
         vx_old(iat, itrj) = vx(iat, itrj)
         vy_old(iat, itrj) = vy(iat, itrj)
         vz_old(iat, itrj) = vz(iat, itrj)
      end do
   end subroutine move_vars

   !******************************
   ! This is the main SH routine !
   !******************************
   subroutine surfacehop(x, y, z, vx, vy, vz, vx_old, vy_old, vz_old, dt, eclas)
      use mod_const, only: ANG, AUTOFS
      use mod_general, only: natom, nwrite, idebug, it, sim_time
      use mod_system, only: names
      use mod_files, only: UPOP, UPROB, UPES, UWFCOEF, UWFCOEF, UNACME, UBKL, UPHASE, UDOTPROD
      use mod_qmmm, only: natqm
      use mod_random, only: vranf
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP), intent(inout) :: eclas
      real(DP), intent(in) :: dt
      real(DP) :: vx_int(size(vx, 1), size(vx, 2))
      real(DP) :: vy_int(size(vy, 1), size(vy, 2))
      real(DP) :: vz_int(size(vz, 1), size(vz, 2))
      real(DP) :: vx_newint(size(vx, 1), size(vx, 2))
      real(DP) :: vy_newint(size(vy, 1), size(vy, 2))
      real(DP) :: vz_newint(size(vz, 1), size(vz, 2))
      real(DP) :: en_array_int(NSTMAX, NTRAJMAX), en_array_newint(NSTMAX, NTRAJMAX)
      real(DP) :: nacx_int(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: nacy_int(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: nacz_int(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: nacx_newint(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: nacy_newint(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: nacz_newint(size(vx, 1), NTRAJMAX, NSTMAX, NSTMAX)
      real(DP) :: dotproduct_int(NSTMAX, NSTMAX, NTRAJMAX), dotproduct_newint(NSTMAX, NSTMAX, NTRAJMAX) !rename dotproduct_int
      real(DP) :: t(NSTMAX, NSTMAX) ! switching probabilities
      real(DP) :: t_tot(NSTMAX, NSTMAX) ! cumulative switching probabilities
      real(DP) :: ran(10)
      real(DP) :: popsum !populations
      integer :: iat, ist1, ist2, itrj, itp ! iteration counters
      integer :: ist ! =istate(itrj)
      real(DP) :: vect_olap, fr, frd
      real(DP) :: Ekin, dtp
      integer :: ihop
      real(DP) :: pop0, prob(NSTMAX), hop_rdnum, stepfs
      character(len=500) :: formt
      character(len=20) :: chist, chihop, chit

      do itrj = 1, ntraj

         call check_energy(vx_old, vy_old, vz_old, vx, vy, vz, itrj)
         call check_energydrift(vx, vy, vz, itrj)

         t_tot = 1.0D0

         ! FIRST, CALCULATE NACME
         if (inac == 0) then

            do ist1 = 1, nstate - 1
               do ist2 = ist1 + 1, nstate
                  if (tocalc(ist1, ist2) == 0) then
                     write (*, *) 'Not computing NACME between states', ist1, ist2
                     ! We need to flush these to zero
                     ! Need to do this here, since tocalc is changed in GET_NACME routine
                     do iat = 1, natqm
                        nacx(iat, itrj, ist1, ist2) = 0.0D0
                        nacy(iat, itrj, ist1, ist2) = 0.0D0
                        nacz(iat, itrj, ist1, ist2) = 0.0D0
                        nacx(iat, itrj, ist2, ist1) = 0.0D0
                        nacy(iat, itrj, ist2, ist1) = 0.0D0
                        nacz(iat, itrj, ist2, ist1) = 0.0D0
                     end do
                  end if
               end do
            end do

            ! This computes and reads NACME
            call get_nacm(itrj)

            ! TODO: move this to a separate routine
            ! Calculating overlap between nacmes
            ! This is crucial, since NACME vectors can change orientation 180 degrees between time steps
            ! and we need to correct that
            do ist1 = 1, nstate
               do ist2 = 1, nstate
                  vect_olap = 0.0D0
                  do iat = 1, natom
                     vect_olap = vect_olap + &
                        & nacx_old(iat, itrj, ist1, ist2) * nacx(iat, itrj, ist1, ist2)
                     vect_olap = vect_olap + &
                        & nacy_old(iat, itrj, ist1, ist2) * nacy(iat, itrj, ist1, ist2)
                     vect_olap = vect_olap + &
                        & nacz_old(iat, itrj, ist1, ist2) * nacz(iat, itrj, ist1, ist2)
                  end do

                  if (vect_olap < 0) then
                     do iat = 1, natom
                        nacx(iat, itrj, ist1, ist2) = -nacx(iat, itrj, ist1, ist2)
                        nacy(iat, itrj, ist1, ist2) = -nacy(iat, itrj, ist1, ist2)
                        nacz(iat, itrj, ist1, ist2) = -nacz(iat, itrj, ist1, ist2)
                     end do
                  end if

               end do
            end do

            ! INAC=0  endif
         end if

         ! smaller time step
         dtp = dt / substep

         ! MAIN LOOP
         ! Smaller time step for electronic population transfer
         do itp = 1, substep

            ist = istate(itrj)

            ! pop0 is later used for Tully's fewest switches
            pop0 = el_pop(ist, itrj)

            ! INTERPOLATION
            fr = real(itp, DP) / real(substep, DP)
            frd = 1.0D0 - fr

            call interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_newint, vy_newint, vz_newint, &
                             nacx_newint, nacy_newint, nacz_newint, en_array_newint, &
                             dotproduct_newint, fr, frd, itrj)

            fr = real(itp - 1, DP) / real(substep, DP)
            frd = 1.0D0 - fr

            call interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                             nacx_int, nacy_int, nacz_int, en_array_int, &
                             dotproduct_int, fr, frd, itrj)

            ! Integrate electronic wavefunction for one dtp time step
            call sh_integrate_wf(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, itrj, dtp)

            ! Check whether total population is 1.0
            popsum = check_popsum(itrj)

            ! Calculate switching probabilities according to Tully's fewest switches
            ! Probabilities are stored in matrix t
            ! (although at this point we calculate only the relevant part of the matrix)
            call sh_TFS_transmat(dotproduct_int, dotproduct_newint, itrj, istate(itrj), pop0, t, dtp)

            if (idebug > 1) then
               ! WARNING: this will not work for adaptive time step
               ! stepfs = (it*substep+itp-substep) * dt * AUtoFS / substep
               stepfs = (sim_time + dtp * itp - substep * dtp) * AUtoFS
               call sh_debug_wf(ist, itrj, stepfs, t)
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
            ! TODO: Refactor this to a separate functions in sh_util
            ! and reuse for Landau-Zener.
            if (nohop /= 1) then

               ! Auxiliary calculations of probabilities on a number line
               prob = 0.0D0
               if (ist == 1) then
                  ! If we are in the ground state, we cannot jump into the ground state :-)
                  prob(1) = 0.0D0
               else
                  ! Probability of jumping from the current state to the ground state
                  prob(1) = t(ist, 1)
               end if

               do ist1 = 2, nstate
                  if (ist1 /= ist) then
                     prob(ist1) = prob(ist1 - 1) + t(ist, ist1)
                  end if
                  if (ist1 == ist) then
                     prob(ist1) = prob(ist1 - 1)
                  end if
               end do

               ihop = 0
               ! Get one random number between 0 and 1
               call vranf(ran, 1)
               hop_rdnum = ran(1)

               ! determine, whether we hopped or not
               do ist1 = 1, nstate
                  if (ist1 == ist) cycle
                  if (hop_rdnum < prob(ist1)) then
                     ihop = ist1
                     exit
                  end if
               end do

               ! Did HOP occur?
               if (ihop /= 0) then
                  if (adjmom == 0) then
                     call hop(vx, vy, vz, ist, ihop, itrj, eclas)
                  else if (adjmom == 1) then
                     call hop_dot(vx, vy, vz, ist, ihop, itrj, eclas)
                  end if

                  if (idebug > 0) then
                     write (formt, '(A8,I3,A7)') '(A1,I10,', nstate + 1, 'E20.10)'
                     write (UPOP, *) '# Substep RandomNum   Probabilities'
                     write (UPOP, fmt=formt) '#', itp, hop_rdnum, (t(ist, ist1), ist1=1, nstate)
                  end if

                  ! TODO: It seems that we're writing this geometry even for frustrated hop?
                  ! Is that desired?
                  ! write current geometry
                  write (chist, *) ist
                  write (chihop, *) ihop
                  write (chit, *) it
                  ! TODO: Rename this file to something more specific,
                  ! e.g. geom_hop.from.to.timestep.xyz
                  formt = 'geom.'//trim(adjustl(chist))//'.'//trim(adjustl(chihop))//'.'//adjustl(chit)
                  open (500, file=trim(formt))
                  write (500, *) natom
                  write (500, *) ''
                  do iat = 1, natom
                     write (500, *) names(iat), x(iat, itrj) / ANG, y(iat, itrj) / ANG, z(iat, itrj) / ANG
                  end do
                  close (500)
               end if

               !nohop endif
            end if

            ! Apply decoherence correction from Persico et al
            if (decoh_alpha > 0) then

               Ekin = ekin_v(vx_int, vy_int, vz_int)
               if (Ekin > 1.0D-4) then ! Decoherence diverges for zero velocities
                  call sh_decoherence_correction(en_array_int, decoh_alpha, Ekin, istate(itrj), itrj, dtp)
               end if

            end if

            !itp loop
         end do

         ! set tocalc array for the next step
         call set_tocalc(itrj)

         popsum = check_popsum(itrj)

         call move_vars(vx, vy, vz, vx_old, vy_old, vz_old, itrj)

         ! TODO: Move this to a separate function write_sh_output()
         if (modulo(it, nwrite) == 0) then
            stepfs = sim_time * AUtoFS
            write (formt, '(A10,I3,A13)') '(F15.2,I3,', nstate, 'F10.5,1F10.7)'
            write (UPOP, fmt=formt) stepfs, istate(itrj), (el_pop(ist1, itrj), ist1=1, nstate), popsum

            t_tot = 1 - t_tot ! up to know, t_tot was the probability of not hopping
            write (formt, '(A10,I3,A6)') '(F15.2,I3,', nstate, 'F10.5)'
            write (UPROB, fmt=formt) stepfs, istate(itrj), (t_tot(ist, ist1), ist1=1, nstate)
            write (formt, '(A7,I3,A7)') '(F15.2,', nstate, 'E20.10)'
            write (UPES, fmt=formt) stepfs, (en_array(ist1, itrj), ist1=1, nstate)
            if (inac == 0) write (UNACME, *) 'Time step:', it
            do ist1 = 1, nstate - 1
               do ist2 = ist1 + 1, nstate

                  if (ist1 == 1 .and. ist2 == 2) then
                     write (UDOTPROD, '(F15.2,E20.10)', advance='no') stepfs, dotproduct_int(ist1, ist2, itrj)
                  else
                     write (UDOTPROD, '(E20.10)', advance='no') dotproduct_int(ist1, ist2, itrj)
                  end if

                  if (inac == 0) then
                     write (UNACME, *) 'NACME between states:', ist1, ist2
                     do iat = 1, natom
                        write (UNACME, '(3E20.10)') nacx(iat, itrj, ist1, ist2), &
                                                  & nacy(iat, itrj, ist1, ist2), &
                                                  & nacz(iat, itrj, ist1, ist2)
                     end do
                  end if

               end do
            end do
            write (UDOTPROD, *) ''

         end if

         ! ntraj enddo
      end do

   end subroutine surfacehop

   subroutine hop(vx, vy, vz, instate, outstate, itrj, eclas)
      use mod_general, only: natom, pot
      use mod_system, only: am
      use mod_files, only: UPOP
      use mod_kinetic, only: ekin_v
      use mod_arrays, only: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: itrj, instate, outstate
      real(DP) :: a_temp, b_temp, c_temp, g_temp
      real(DP) :: ekin, ekin_new
      integer :: iat

      ! Checking for frustrated hop

      a_temp = 0.D0
      b_temp = 0.D0

      ekin = ekin_v(vx, vy, vz)

      do iat = 1, natom
         a_temp = a_temp + nacx(iat, itrj, instate, outstate)**2 / am(iat)
         a_temp = a_temp + nacy(iat, itrj, instate, outstate)**2 / am(iat)
         a_temp = a_temp + nacz(iat, itrj, instate, outstate)**2 / am(iat)
         b_temp = b_temp + nacx(iat, itrj, instate, outstate) * vx(iat, itrj)
         b_temp = b_temp + nacy(iat, itrj, instate, outstate) * vy(iat, itrj)
         b_temp = b_temp + nacz(iat, itrj, instate, outstate) * vz(iat, itrj)
      end do
      a_temp = 0.5D0 * a_temp
      c_temp = b_temp**2 + 4 * a_temp * (en_array(instate, itrj) - en_array(outstate, itrj))

      if (c_temp < 0) then
         write (*, *) '# Not enough momentum in the direction of NAC vector.'
         ! Try, whether there is enough total kinetic energy and scale velocities.
         call hop_dot(vx, vy, vz, instate, outstate, itrj, eclas)
         return
      end if

      istate(itrj) = outstate
      eclas = en_array(outstate, itrj)
      write (*, '(A24,I3,A10,I3)') '# Hop occured from state ', instate, ' to state ', outstate

      ! Rescaling the velocities

      if (b_temp < 0) then
         g_temp = (b_temp + dsqrt(c_temp)) / 2.0D0 / a_temp
      end if

      if (b_temp >= 0) then
         g_temp = (b_temp - dsqrt(c_temp)) / 2.0D0 / a_temp
      end if

      do iat = 1, natom
         vx(iat, itrj) = vx(iat, itrj) - g_temp * nacx(iat, itrj, instate, outstate) / am(iat)
         vy(iat, itrj) = vy(iat, itrj) - g_temp * nacy(iat, itrj, instate, outstate) / am(iat)
         vz(iat, itrj) = vz(iat, itrj) - g_temp * nacz(iat, itrj, instate, outstate) / am(iat)
      end do
      ekin_new = ekin_v(vx, vy, vz)

      write (*, '(A31,2E20.10)') '# deltaE_pot     E_kin-total', &
                               & en_array(outstate, itrj) - en_array(instate, itrj), ekin
      write (*, '(A,2E20.10)') '# Total_energy_old   Total_energy_new :', &
                             & ekin + en_array(instate, itrj), &
                             & ekin_new + en_array(outstate, itrj)

      call set_tocalc(itrj)
      write (*, *) 'Calculating forces for the new state.'
      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

   end subroutine hop

   subroutine interpolate(vx, vy, vz, vx_old, vy_old, vz_old, vx_int, vy_int, vz_int, &
                          nacx_int, nacy_int, nacz_int, en_array_int, &
                          dotproduct_int, fr, frd, itrj)
      use mod_general, only: natom
      real(DP), intent(out) :: dotproduct_int(:, :, :)
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      real(DP), intent(out) :: vx_int(:, :), vy_int(:, :), vz_int(:, :)
      real(DP), intent(out) :: nacx_int(:, :, :, :)
      real(DP), intent(out) :: nacy_int(:, :, :, :)
      real(DP), intent(out) :: nacz_int(:, :, :, :)
      real(DP), intent(out) :: en_array_int(:, :)
      real(DP) :: fr, frd
      integer :: iat, ist1, ist2, itrj !iteration counters

      do ist1 = 1, nstate

         en_array_int(ist1, itrj) = en_array(ist1, itrj) * fr + en_array_old(ist1, itrj) * frd
         do ist2 = 1, nstate
            dotproduct_int(ist1, ist2, itrj) = 0.0D0
            do iat = 1, natom
               nacx_int(iat, itrj, ist1, ist2) = nacx(iat, itrj, ist1, ist2) * fr + &
                                               & nacx_old(iat, itrj, ist1, ist2) * frd
               nacy_int(iat, itrj, ist1, ist2) = nacy(iat, itrj, ist1, ist2) * fr + &
                                               & nacy_old(iat, itrj, ist1, ist2) * frd
               nacz_int(iat, itrj, ist1, ist2) = nacz(iat, itrj, ist1, ist2) * fr + &
                                               & nacz_old(iat, itrj, ist1, ist2) * frd
               vx_int(iat, itrj) = vx(iat, itrj) * fr + vx_old(iat, itrj) * frd
               vy_int(iat, itrj) = vy(iat, itrj) * fr + vy_old(iat, itrj) * frd
               vz_int(iat, itrj) = vz(iat, itrj) * fr + vz_old(iat, itrj) * frd
               dotproduct_int(ist1, ist2, itrj) = dotproduct_int(ist1, ist2, itrj) + &
                                               & vx_int(iat, itrj) * nacx_int(iat, itrj, ist1, ist2)
               dotproduct_int(ist1, ist2, itrj) = dotproduct_int(ist1, ist2, itrj) + &
                                               & vy_int(iat, itrj) * nacy_int(iat, itrj, ist1, ist2)
               dotproduct_int(ist1, ist2, itrj) = dotproduct_int(ist1, ist2, itrj) + &
                                               & vz_int(iat, itrj) * nacz_int(iat, itrj, ist1, ist2)
            end do
         end do
      end do

   end subroutine interpolate

   subroutine hop_dot(vx, vy, vz, instate, outstate, itrj, eclas)
      use mod_general, only: natom, pot
      use mod_kinetic, only: ekin_v
      use mod_arrays, only: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: itrj, instate, outstate
      integer :: iat
      real(DP) :: de, ekin, alfa, ekin_new

      ekin = 0.0D0
      ekin_new = 0.0D0

      dE = en_array(outstate, itrj) - en_array(instate, itrj)
      ekin = ekin_v(vx, vy, vz)

      if (ekin >= de) then

         alfa = dsqrt(1 - de / ekin)

         do iat = 1, natom
            vx(iat, itrj) = alfa * vx(iat, itrj)
            vy(iat, itrj) = alfa * vy(iat, itrj)
            vz(iat, itrj) = alfa * vz(iat, itrj)
         end do
         istate(itrj) = outstate
         eclas = en_array(outstate, itrj)
         ekin_new = ekin_v(vx, vy, vz)

         write (*, '(A24,I3,A10,I3)') '# Hop occured from state ', instate, ' to state ', outstate
         write (*, *) '# Adjusting velocities by simple scaling.'
         write (*, '(A,2E20.10)') '#Total_energy_old   Total_energy_new :', &
                               & ekin + en_array(instate, itrj), ekin_new + en_array(outstate, itrj)

         call set_tocalc(itrj)
         write (*, *) 'Calculating forces for the new state.'
         call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)

      else

         write (*, '(A35,I3,A10,I3)') '# Frustrated Hop occured from state ', &
                                    & instate, ' to state ', outstate
         if (revmom == 1) then
            write (*, *) '# Reversing momentum direction.'
            vx = -vx
            vy = -vy
            vz = -vz
         end if

      end if

      write (*, '(A31,2E20.10)') '# deltaE_pot  E_kin-total', dE, ekin

   end subroutine hop_dot

   subroutine check_energy(vx_old, vy_old, vz_old, vx, vy, vz, itrj)
      use mod_interfaces, only: finish
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP), intent(in) :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
      integer, intent(in) :: itrj
      real(DP) :: ekin, ekin_old, entot, entot_old, dE_S0S1

      ! Special case for running MD with TDDFT:
      ! End the simulation when S1-S0 energy difference drops below certain
      ! small threshold.
      if (nstate >= 2) then
         dE_S0S1 = en_array(2, itrj) - en_array(1, itrj)
         dE_S0S1 = dE_S0S1 * AUtoEV
         if (dE_S0S1 < dE_S0S1_thr) then
            write (*, *) 'S1 - S0 gap dropped below threshold!'
            write (*, *) dE_S0S1, ' < ', dE_S0S1_thr
            call finish(0)
            stop 0
         end if
      end if

      ekin = ekin_v(vx, vy, vz)
      ekin_old = ekin_v(vx_old, vy_old, vz_old)

      entot = (ekin + en_array(istate(itrj), itrj)) * AUtoEV

      ! TODO: But what if we hopped to another state?
      ! en_array_old(istate(itrj), itrj) would not point to the correct energy,
      ! right?
      ! We need istate_old array ... but we should just refactor the whole thing...
      ! and make traj_data structure or something
      entot_old = (ekin_old + en_array_old(istate(itrj), itrj)) * AUtoEV

      if (abs(entot - entot_old) > energydifthr) then
         write (*, *) 'ERROR:Poor energy conservation. Exiting...'
         write (*, *) 'Total energy difference [eV] is:', entot - entot_old
         write (*, *) 'The threshold was:', energydifthr
         write (*, *) 'Ekin_old, Ekin, Epot_old, E_pot', &
            ekin_old, ekin, en_array_old(istate(itrj), itrj), en_array(istate(itrj), itrj)
         call abinerror('check_energy')
      end if

   end subroutine check_energy

   subroutine check_energydrift(vx, vy, vz, itrj)
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      integer, intent(in) :: itrj
      real(DP) :: ekin, entot

      ekin = ekin_v(vx, vy, vz)

      entot = (ekin + en_array(istate(itrj), itrj)) * AUtoEV

      if (abs(entot - entot0) > energydriftthr) then
         write (*, *) 'ERROR: Energy drift exceeded threshold value. Exiting...'
         write (*, *) 'Total energy difference [eV] is:', entot - entot0
         write (*, *) 'The threshold was:', energydriftthr
         call abinerror('check_energy')
      end if

   end subroutine check_energydrift

   integer function check_CIVector(CIvecs, CIvecs_old, ci_len, nstates)
      use mod_const, only: AUtoFS
      use mod_files, only: UDOTPRODCI
      use mod_general, only: sim_time
      real(DP), allocatable, intent(in) :: CIvecs(:, :), CIvecs_old(:, :)
      integer, intent(in) :: ci_len, nstates
      real(DP) :: cidotprod(NSTMAX)
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
