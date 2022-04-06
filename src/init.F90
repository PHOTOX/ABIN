! init() subroutine performs many things:
!
! 1. Reading input
! 2. Input sanity check
! 3. Allocation of arrays
! 4. Reading restart OR reading input geometry
! 5. Initialize velocities
! 6. Initialize everything else
!
! Coordinate and velocity transformations for path integral MD are NOT performed here.
! Surface hopping is initialized in sh_init().
!
! In general, each module should have their own init function,
! that can be called from here.

module mod_init
   implicit none
   private
   public :: init
   ! For unit test
   public :: initialize_masses
contains   

subroutine init(dt)
   use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
   use mod_const
   use mod_interfaces, only: print_compile_info, omp_set_num_threads
   use mod_cmdline, only: get_cmdline
   use mod_error, only: fatal_error
   use mod_mpi, only: initialize_mpi, get_mpi_size, get_mpi_rank, mpi_barrier_wrapper
   use mod_files
   use mod_arrays
   use mod_array_size
   use mod_general
   use mod_system
   use mod_nhc
   use mod_estimators
   use mod_potentials
   use mod_sh_integ, only: nstate, integ, phase, popsumthr, correct_decoherence
   use mod_sh
   use mod_lz
   use mod_qmmm, only: natqm, natmm
   use mod_force_mm
   use mod_gle
   use mod_sbc, only: sbc_init, rb_sbc, kb_sbc, isbc, rho
   use mod_random
   use mod_splined_grid, only: initialize_spline, potential_file
   use mod_utils, only: toupper, tolower, normalize_atom_name, file_exists_or_exit
   use mod_vinit
   use mod_analyze_geometry
   use mod_shake
   use mod_minimize, only: gamm, gammthr
   use mod_analysis, only: restin
   use mod_water, only: watpot, check_water
   use mod_plumed, only: iplumed, plumedfile, plumed_init
   use mod_en_restraint
   use mod_transform, only: initialize_pi_masses
   use mod_cp2k
   use mod_remd
   use mod_terampi
   use mod_terampi_sh
   implicit none
   real(DP), intent(out) :: dt
   ! Input parameters for analytical potentials
   real(DP) :: lambda_dw = -1.0D0, D0_dw = -1.0D0, k_dw = -1.0D0, r0_dw = -1.D0
   real(DP) :: r0_morse = -1, d0_morse = -1, k_morse = -1
   real(DP) :: k = -1, r0 = -1
   real(DP) :: kx = -1, ky = -1, kz = -1
   ! Initial temperature (read from namelist nhcopt)
   real(DP) :: temp0 = -1
   real(DP) :: masses(MAXTYPES)
   integer ::  iw, iat, natom_xyz, iost
   integer :: shiftdihed
   ! Number of OpenMP processes, read from ABIN input
   ! WARNING: We do NOT use OMP_NUM_THREADS environment variable!
   integer :: nproc
   integer :: getPID
!$ integer, external :: omp_get_max_threads
   character(len=2) :: massnames(MAXTYPES)
   character(len=2) :: atom
   character(len=200) :: chinput, chcoords, chveloc
   character(len=200) :: chiomsg, chout
   character(len=20) :: xyz_units
   character(len=60) :: chdivider
   character(len=60) :: mdtype
   character(len=1024) :: tc_server_name
   logical :: file_exists
   logical :: rem_comvel, rem_comrot
   logical :: testing_mode
   integer :: my_rank, mpi_world_size

   ! ABIN input parameters are read from the input file (default 'input.in')
   ! in the form of the standard Fortran namelist syntax.
   ! The input parameters are grouped in several different namelists:
   !
   ! - general: Most basic MD settings + misc
   ! - nhcopt:  parameters for thermostats
   ! - remd:    parameters for Replica Exchange MD
   ! - sh:      parameters for Surface Hopping
   ! - system:  system-specific parameters for model potentials, masses, SHAKE constraints...
   ! - lz:      parameters for Landau-Zener excited state dynamics.
   ! - qmmm:    parameters for internal QMMM (not really tested).
   !
   ! All namelists need to be in a single input file, and the code
   ! in this subroutine must ensure that the namelists can be in any order.
   namelist /general/ natom, pot, ipimd, mdtype, istage, inormalmodes, nwalk, nstep, icv, ihess, imini, nproc, iqmmm, &
      nwrite, nwritex, nwritev, nwritef, dt, irandom, nabin, irest, nrest, anal_ext, &
      isbc, rb_sbc, kb_sbc, gamm, gammthr, conatom, mpi_sleep, narchive, xyz_units, &
      dime, ncalc, idebug, enmini, rho, iknow, watpot, iremd, iplumed, plumedfile, &
      en_restraint, en_diff, en_kk, restrain_pot, &
      pot_ref, nstep_ref, nteraservers, max_wait_time, cp2k_mpi_beads, testing_mode

   namelist /remd/ nswap, nreplica, deltaT, Tmax, temp_list

   namelist /nhcopt/ inose, temp, temp0, nchain, tau0, tau0_langevin, imasst, nrespnose, nyosh, &
      scaleveloc, readNHC, nmolt, natmolt, nshakemol, rem_comrot, rem_comvel, gle_test

   namelist /system/ masses, massnames, ndist, dist1, dist2, &
      nang, ang1, ang2, ang3, ndih, dih1, dih2, dih3, dih4, shiftdihed, &
      k, r0, kx, ky, kz, r0_morse, d0_morse, k_morse, D0_dw, lambda_dw, k_dw, r0_dw, &
      Nshake, ishake1, ishake2, shake_tol, potential_file

   namelist /sh/ istate_init, nstate, substep, deltae, integ, inac, nohop, phase, decoh_alpha, popthr, ignore_state, &
      nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, adjmom, revmom, &
      dE_S0S1_thr, correct_decoherence

   namelist /lz/ initstate_lz, nstate_lz, nsinglet_lz, ntriplet_lz, deltaE_lz, energydifthr_lz

   namelist /qmmm/ natqm, natmm, q, rmin, eps, attypes

   chcoords = 'mini.dat'
   xyz_units = 'angstrom'
   chinput = 'input.in'
   chveloc = ''
   tc_server_name = ''
   mdtype = ''
   dt = -1
   nproc = 1
   iplumed = 0
   shiftdihed = 1

   chdivider = "######################################################"

   call get_cmdline(chinput, chcoords, chveloc, tc_server_name)

   ! Reading main input parameters from namelist &general
   open (150, file=chinput, status='OLD', delim='APOSTROPHE', action="READ")
   read (150, general)
   rewind (150)
   pot = tolower(pot)

   if (pot == "_cp2k_" .or. pot_ref == "_cp2k_") then
      call init_cp2k()
   end if

   call initialize_mpi(pot, pot_ref, nteraservers)

   ! Set OpenMP parallelization
   ! Currently only used in PIMD for trivial
   ! parallelization over PI beads.
   ! Note that scaling is actually not so great
   ! since SCF timings will vary for different beads,
   ! which decreases thread utilization.
   if (pot == "_tera_" .and. nteraservers > 1) then
      nproc = nteraservers
   end if
!$ call omp_set_num_threads(nproc)

   my_rank = get_mpi_rank()
   mpi_world_size = get_mpi_size()

   if (iremd == 0 .and. mpi_world_size > 1) then
      call fatal_error(__FILE__, __LINE__, &
         & 'MPI parallelization available only for REMD')
   end if

   ! Redirect stdout to /dev/null for all non-primary MPI replicas
   if (my_rank > 0) then
      call stdout_to_devnull()
   end if

   if (mpi_world_size > 1) then
      write (stdout, '(A,I0)') 'Number of MPI processes = ', mpi_world_size
   end if

   ! We need to connect to TeraChem as soon as possible,
   ! because we want to shut down TeraChem nicely in case something goes wrong.
   if (pot == '_tera_' .or. restrain_pot == '_tera_' .or. pot_ref == '_tera_') then
      call initialize_terachem_interface(trim(tc_server_name))
   end if

   if (mdtype /= '') then
      mdtype = tolower(mdtype)
      select case (mdtype)
      case ('md')
         ipimd = 0
      case ('pimd')
         ipimd = 1
      case ('sh')
         ipimd = 2
      case ('minimization')
         ipimd = 3
      case ('landau_zener')
         ipimd = 5
      case default
         call fatal_error(__FILE__, __LINE__, 'invalid mdtype in '//trim(chinput))
      end select
   end if

   if (iremd == 1) then
      write (chcoords, '(A,I2.2)') trim(chcoords)//'.', my_rank
      if (chveloc /= '') then
         write (chveloc, '(A,I2.2)') trim(chveloc)//'.', my_rank
      end if
   end if

   call file_exists_or_exit(chcoords)
   if (chveloc /= '') then
      call file_exists_or_exit(chveloc)
   end if

   call print_basic_info()

   ! Get number of atoms from XYZ coordinates NOW so that we can allocate arrays

   open (111, file=chcoords, status="old", action="read")
   read (111, '(I50)', iostat=iost) natom_xyz
   !TODO following line does not work
   if (iost /= 0) call print_read_error(chcoords, "Expected number of atoms on the first line.", iost)
   if (natom_xyz /= natom .and. natom > 0) then
      write (*, '(A,A)') 'WARNING: Number of atoms specified in ', trim(chinput)
      write (*, '(A,A)') 'does not match with the XYZ geometry in ', trim(chcoords)
      write (*, *) 'Going forward anyway...'
   end if

   natom = natom_xyz

   if (natom < 1) then
      write (*, '(A,A)') 'ERROR: Wrong number of atoms on the first line of the XYZ &
      & file ', trim(chcoords)
      write (*, *) natom
      call abinerror('init')
   end if

   ! This line is super important,
   ! cause we actually use natqm in many parts of the code
   if (iqmmm == 0 .and. pot /= 'mm') then
      natqm = natom
   end if

   if (irest == 1) then
      readnhc = .true.
      scaleveloc = 0 !do not scale velocities when restarting a job
   else
      readnhc = .false.
      scaleveloc = 1
   end if

   if (chveloc /= '') then
      scaleveloc = 0
   end if

! for future multiple time step integration in SH
   dt0 = dt

   ! We have to initialize here, because we read them from input
   allocate (names(natom))
   names = ''
   attypes = ''
   massnames = ''
   masses = -1.0D0

   allocate (ishake1(natom * 3 - 6))
   allocate (ishake2(natom * 3 - 6))
   ishake1 = 0
   ishake2 = 0

   ! By default, remove COM translation and rotation
   ! Unless restarting or taking velocities from a file
   if (irest == 1 .or. trim(chveloc) /= '') then
      rem_comrot = .false.
      rem_comvel = .false.
   else
      rem_comrot = .true.
      rem_comvel = .true.
   end if

   ! Allocate all basic arrays and set them to 0.0d0
   call allocate_arrays(natom, nwalk)

   if (iplumed == 1) then
      call plumed_init()
   end if

   ! Read initial geometry
   ! TODO: Move to a separate function
   read (111, *)
   do iat = 1, natom

      read (111, *, iostat=iost) names(iat), x(iat, 1), y(iat, 1), z(iat, 1)
      if (iost /= 0) call print_read_error(chcoords, 'Could not read atom names and coordinates', iost)
      names(iat) = normalize_atom_name(names(iat))
      if (tolower(trim(xyz_units)) == "angstrom") then
         x(iat, 1) = x(iat, 1) * ANG
         y(iat, 1) = y(iat, 1) * ANG
         z(iat, 1) = z(iat, 1) * ANG
      else if (tolower(trim(xyz_units)) == "bohr") then
         continue
      else
         write (*, *) 'ERROR: Wrong XYZ units: ', trim(xyz_units)
      end if

   end do
   close (111)

   ! Initialize all PIMD beads to the same initial geometry
   do iw = 1, nwalk
      do iat = 1, natom
         x(iat, iw) = x(iat, 1)
         y(iat, iw) = y(iat, 1)
         z(iat, iw) = z(iat, 1)
      end do
   end do

   ! the namelist system does not need to be present
   read (150, system, iostat=iost, iomsg=chiomsg)
   rewind (150)
   !check, whether we hit End-Of-File or other error
   if (IS_IOSTAT_END(iost)) then !fortran intrinsic for EOF
      continue
   else if (iost /= 0) then
      write (*, *) 'ERROR when reading namelist "system".'
      write (*, *) chiomsg
      call abinerror('init')
   end if

   do iat = 1, MAXTYPES
      massnames(iat) = normalize_atom_name(massnames(iat))
   end do

   ! Determine atomic masses from periodic table
   call initialize_masses(masses, massnames, natom)
   ! Transform masses for PIMD
   ! Note that amt array is used throughout the code
   ! even for non-PI simulations.
   call initialize_pi_masses(amg, amt)

   allocate (natmolt(natom))
   allocate (nshakemol(natom))
   natmolt = 0
   natmolt(1) = natom ! default for global NHC thermostat
   nshakemol = 0

   read (150, nhcopt)
   rewind (150)

   if (ipimd == 2) then
      read (150, sh)
      rewind (150)
      integ = tolower(integ)
   end if

   if (ipimd == 5) then
      read (150, lz)
      rewind (150)
      call lz_init() !Init arrays for possible restart
   end if

   if (iremd == 1) then
      read (150, remd)
      rewind (150)
      call remd_init(temp, temp0)
   end if

   if (iqmmm > 0 .or. pot == 'mm') then
      read (150, qmmm)
      rewind (150)
   end if

   close (150)
   ! END OF READING INPUT

   if (pot == '_tera_' .or. restrain_pot == '_tera_') then
      call initialize_tc_servers()
      if (ipimd == 2 .or. ipimd == 5) then
         call init_terash(x, y, z)
      end if
   end if

   ! Check whether input parameters make sense
   call check_inputsanity()

   ! resetting number of walkers to 1 in case of classical simulation
   if (ipimd == 0) then
      write (stdout, *) 'Using velocity Verlet integrator'
      md = 2
      nwalk = 1
      nabin = 1 ! TODO:safety for respashake code
      ! We should probably copy shake to velocity verlet
      ! algorithm as well
   end if

   ! Selecting proper integrator for a given MD type
   ! (controlled by variable 'md')
   if (ipimd == 2 .or. ipimd == 5) then
      nwalk = 1
      md = 2 ! velocity verlet
      nabin = 1
   else if (ipimd == 1 .and. inormalmodes /= 1) then
      write (stdout, *) 'Using RESPA integrator.'
      md = 1
   else if (ipimd == 1 .and. inormalmodes == 1) then
      write (stdout, *) 'Using velocity Verlet propagator with analytical PI normal mode propagation.'
      nabin = 1
      md = 2
   end if

   if (nshake /= 0) then
      md = 3
   end if

   if (pot_ref /= '_none_') then
      md = 4
      write (*, '(A)') 'Using Multiple Time-Step RESPA integrator'
      write (*, '(A)') "Reference (cheap) potential is "//trim(pot_ref)
      write (*, '(A, F6.2)') "with timestep [fs] ", dt / nstep_ref * AUtoFS
      write (*, '(A)') "Full potential is "//trim(pot)
      write (*, '(A, F6.2)') "with timestep [fs] ", dt * AUtoFS
      if (ipimd /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'ab initio MTS is implemented only for classical MD')
      end if
   end if

   if (temp0 > 0) then
      write (stdout, *) 'Initial temperature [K] =', temp0
   else
      write (stdout, *) 'Initial temperature [K] =', temp
   end if
   if (inose /= 0) write (stdout, *) 'Target temperature [K] =', temp

   ! Convert temperature from Kelvins to atomic units
   temp = temp / AUtoK
   temp0 = temp0 / AUtoK

   if (ihess == 1) then
      allocate (cvhess_cumul(nwalk))
      cvhess_cumul = 0.0D0
   end if

   ! Initialize SHAKE, determine the constrained bond lengths
   if (nshake /= 0) then
      call shake_init(x, y, z)
   end if

   if (pot == '_mmwater_' .or. pot_ref == '_mmwater_') then
      call check_water(natom, names)
   end if

   ! Initialize pseudo-random number generator
   ! TODO: move this up in the init
   call initialize_prng(irandom, my_rank, testing_mode)

   ! Initialize thermostat
   if (inose == 1) then
      call nhc_init()
   else if (inose == 2) then
      call gle_init(dt * 0.5 / nabin / nstep_ref) !nabin is set to 1 unless ipimd=1
   else if (inose == 3) then
      call pile_init(dt * 0.5, tau0_langevin)
   else if (inose == 4) then
      call gle_init(dt * 0.5 / nstep_ref)
   else if (inose == 0) then
      write (stdout, '(A)') 'No thermostat. NVE ensemble.'
   else
      call fatal_error(__FILE__, __LINE__, 'Invalid "inose" value!')
   end if

   ! performing RESTART from restart.xyz
   if (irest == 1) then
      call restin(x, y, z, vx, vy, vz, it)
   end if

   if (nchain > 1) then
      f = 0 !what about nchain=1?
      ! what about massive thermostat?
   end if

   ! Set initial velocities according to the Maxwell-Boltzmann distribution
   if (irest == 0 .and. chveloc == '') then
      if (temp0 >= 0) then
         call vinit(temp0, am, vx, vy, vz)
      else
         call vinit(temp, am, vx, vy, vz)
      end if
   end if

!     Reading velocities from file
   if (chveloc /= '' .and. irest == 0) then
      ! TODO: move the following to a separate function
      open (500, file=chveloc, status='OLD', action="READ")
      do iw = 1, nwalk
         read (500, *, IOSTAT=iost) natom_xyz
         if (iost /= 0) call print_read_error(chveloc, "Could not read velocities on line 1.", iost)
         if (natom_xyz /= natom) then
            write (*, '(A,A)') 'Nunmber of atoms in velocity input ', trim(chveloc)
            write (*, '(A,A)') 'does not match with XYZ coordinates in ', trim(chcoords)
            call abinerror('init')
         end if
         read (500, *, IOSTAT=iost)
         if (iost /= 0) call print_read_error(chveloc, "Could not read velocities on line 2.", iost)

         do iat = 1, natom
            read (500, *, IOSTAT=iost) atom, vx(iat, iw), vy(iat, iw), vz(iat, iw)
            if (iost /= 0) call print_read_error(chveloc, "Could not read velocities.", iost)
            atom = normalize_atom_name(atom)
            if (atom /= names(iat)) then
               write (*, *) 'Offending line:'
               write (*, *) atom, vx(iat, iw), vy(iat, iw), vz(iat, iw)
               call print_read_error(chveloc, "Inconsistent atom types in input velocities.", iost)
            end if
         end do
      end do

      close (500)

   end if

!     END OF READING VELOCITIES--------------------

   ! doing this here so that we can do it even when reading velocities from file
   if (rem_comvel) then
      call remove_comvel(vx, vy, vz, am, rem_comvel)
   end if
   if (rem_comrot) then
      call remove_rotations(x, y, z, vx, vy, vz, am, rem_comrot)
   end if

   if (conatom > 0) then
      call constrainP(vx, vy, vz, conatom)
   end if

   ! If scaleveloc=1, scale initial velocitites to match the temperature
   ! Otherwise, just print the temperature.
   call ScaleVelocities(vx, vy, vz)

   ! Initialize spherical boundary onditions
   if (isbc == 1) then
      call sbc_init(x, y, z)
   end if

   if (en_restraint >= 1) then
      call en_rest_init(natom)
   end if

   if (pot == '_doublewell_' .or. pot_ref == '_doublewell_') then
      call doublewell_init(natom, lambda_dw, d0_dw, k_dw, r0_dw, vy, vz)
   end if
   if (pot == '_harmonic_oscillator_' .or. pot_ref == '_harmonic_oscillator_') then
      call harmonic_oscillator_init(natom, kx, ky, kz, vx, vy, vz)
   end if
   if (pot == '_harmonic_rotor_' .or. pot_ref == '_harmonic_rotor_') then
      call harmonic_rotor_init(natom, k, r0)
   end if
   if (pot == '_morse_' .or. pot_ref == '_morse_') then
      call morse_init(natom, k_morse, r0_morse, d0_morse)
   end if
   if (pot == '_splined_grid_' .or. pot_ref == '_splined_grid_') then
      call initialize_spline(natom)
   end if

   if (pot == 'mm') then
      ! TODO: Move this to a single subroutine in force_mm.F90
      allocate (inames(natom))
      do iat = 1, MAXTYPES
         if (attypes(iat) == '') exit
         attypes(iat) = normalize_atom_name(attypes(iat))
      end do
      ! inames initialization for the MM part.
      ! We do this also because string comparison is very costly
      call inames_init()
      call ABr_init()
   end if

   if (my_rank == 0) then
      call print_simulation_parameters()
   end if

   if (mpi_world_size > 1) then
      call mpi_barrier_wrapper()
      write (*, '(A,I0,A,I0)') 'MPI rank: ', my_rank, ' PID: ', GetPID()
      call mpi_barrier_wrapper()
   else
      write (stdout, '(A,I0)') 'Process ID (PID): ', GetPID()
   end if

!$ write (stdout, '(A,I0)') 'Number of OpenMP threads: ', omp_get_max_threads()

   ! Open permanent files for writing
   call files_init(isbc, phase, ndist, nang, ndih)

   call flush (OUTPUT_UNIT)

contains

   subroutine check_inputsanity()
      use mod_chars, only: chknow
      use mod_utils, only: real_positive, real_nonnegative, &
                         & int_positive, int_nonnegative, int_switch
      integer :: error
!$    integer :: nthreads, omp_get_max_threads

      error = 0

      !  We should exclude all non-abinitio options, but whatever....
!$    nthreads = omp_get_max_threads()
!$    if (nthreads > 1 .and. (ipimd /= 1 .and. pot /= '_cp2k_')) then
!$       write (*, *) 'Number of threads is ', nthreads
!$       write (*, *) 'ERROR: Parallel execution is currently only supported with ab initio PIMD (ipimd=1)'
!$       call abinerror('init')
!$    end if

      if (nproc > 1) then
!$       if (.false.) then
            write (*, *) 'FATAL ERROR: This executable was not compiled with parallel support.'
            error = 1
!$       end if
      end if

      if (irest == 1 .and. chveloc /= '') then
         write (*, *) 'WARNING: Input velocities from file '//trim(chveloc)//' will be ignored!'
         write (*, *) 'Velocities will be taken from restart file because irest=1.'
      end if

      if (nstate > nstmax) then
         write (*, *) 'Maximum number of states is:'
         write (*, *) nstmax
         write (*, *) 'Adjust variable nstmax in modules.f90'
         error = 1
      end if

      if (pot /= '_cp2k_') then
         if (nproc > nwalk) then
            write (*, *) 'ERROR: Nproc greater than nwalk. That does not make sense.'
            write (*, *) 'Set nproc <= nwalk.'
            error = 1
         end if
         if (nproc <= 0) then
            write (*, *) 'ERROR: Nproc must be a positive integer.'
            error = 1
         end if
         if (modulo(nwalk, nproc) /= 0) then
            write (*, *) 'ERROR: Nwalk is not divisible by the number of OpenMP threads.'
            write (*, *) 'This is not a wise usage of your computer time.'
            error = 1
         end if
      end if
      if (pot == '_none_') then
         write (*, *) 'ERROR: Variable "pot" not specified.'
         error = 1
      end if
      if (ipimd == 1 .and. nwalk <= 1) then
         write (*, *) 'Number of walkers for PIMD (nwalk) <=1 !'
         write (*, *) 'Either set ipimd=0 for classical simulation or'
         write (*, *) 'set nwalk > 1'
         error = 1
      end if
      if (iqmmm < 0 .or. iqmmm > 1) then
         write (*, *) 'Error: iqmmm must be 0 or 1 for ONIOM.'
         error = 1
      end if
      if (integ /= 'euler' .and. integ /= 'rk4' .and. integ /= 'butcher') then
         write (*, *) 'integ must be "euler", "rk4" or "butcher".'
         error = 1
      end if
      if (integ /= 'butcher') then
         write (*, *) 'WARNING: variable integ is not "butcher", which is the default and most accurate.'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (deltae < 0) then
         write (*, *) 'Parameter deltae must be non-negative number.'
         error = 1
      end if
      if (popsumthr < 0) then
         write (*, *) 'Parameter popsumthr must be positive number.'
         error = 1
      end if
      if (energydifthr < 0) then
         write (*, *) 'Parameter energydifthr must be positive number in eV units.'
         error = 1
      end if
      if (energydriftthr < 0) then
         write (*, *) 'Parameter energydriftthr must be positive number in eV units.'
         error = 1
      end if
      if (shiftdihed /= 0 .and. shiftdihed /= 1) then
         write (*, *) 'Shiftdihed must be either 0 (for dihedrals -180:180) or 1 (for dihedrals 0:360)'
         error = 1
      end if

      if (shiftdihed == 0) shiftdih = 0.0D0
      if (shiftdihed == 1) shiftdih = 360.0D0

      if (icv == 1 .and. temp <= 0.0D0) then
         write (*, *) 'Cannot compute heat capacity for zero temperature.'
         error = 1
      end if
      if (icv == 1 .and. inose == 0) then
         write (*, *) 'Cannot compute heat capacity for NVE simulation.'
         error = 1
      end if
      if (ncalc > nwrite) then
         write (*, *) 'Ncalc greater than nwrite.Setting nwrite=ncalc'
         nwrite = ncalc
      end if

      if (ipimd == 1 .and. inose <= 0) then
         write (*, *) 'You have to use thermostat with PIMD! (inose>=0)'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (ipimd < 0 .or. ipimd > 5) then
         write (*, *) 'Invalid ipimd value'
         error = 1
      end if
      if (ipimd == 5 .and. pot == '_tera_' .and. ntriplet_lz > 0) then
         write (*, *) 'ERROR: Landau-Zener with Singlet-Triplet transitions not implemented with TeraChem over MPI.'
         error = 1
      end if
      if (ipimd == 5 .and. nwalk /= 1) then
         write (*, *) 'ERROR: LZ not implemented with multiple walkers.'
         error = 1
      end if
      if (inormalmodes < 0 .and. inormalmodes > 2) then
         write (*, *) 'ERROR: inormalmodes has to be 0, 1 or 2!'
         error = 1
      end if
      if (inac > 2 .or. inac < 0) then
         write (*, *) 'Parameter "inac" must be 0,1 or 2.' !be very carefull if you change this!
         error = 1
      end if
      if (adjmom > 1 .or. adjmom < 0) then
         write (*, *) 'Parameter "adjmom" must be 0 or 1.'
         error = 1
      end if
      if (adjmom == 0 .and. inac == 1) then
         write (*, *) 'Combination of adjmom=0 and inac=1 is not possible.'
         write (*, *) 'We dont have NAC vector if inac=1.'
         error = 1
      end if
      if (irest == 1 .and. scaleveloc == 1) then
         write (*, *) 'irest=1 AND scaleveloc=1.'
         write (*, *) 'You are trying to scale the velocities read from restart.xyz.'
         write (*, *) 'I assume this is an error in input. (set scaleveloc=0)'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (inose == 1 .and. (ipimd == 2) .and. en_restraint == 0) then
         write (*, *) 'Thermostating is not meaningful for surface hopping simulation.'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (inose == 1 .and. (ipimd == 5) .and. en_restraint == 0) then
         write (*, *) 'Thermostating is not meaningful for Landau-Zener MD.'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (istate_init > nstate) then
         write (*, *) 'Error:Initial state > number of computed states. Exiting...'
         error = 1
      end if
      if (ipimd == 5) then
         if (initstate_lz > nstate_lz) then
            write (*, *) initstate_lz, nstate_lz
            write (*, *) 'Error(LZ):Initial state > number of computed states. Exiting...'
            error = 1
         end if
         if (nstate_lz <= 0) then
            write (*, *) 'Error(LZ):No states to compute (nstate_lz<=0). Exiting...'
            error = 1
         end if
         if (nsinglet_lz == 0 .and. ntriplet_lz == 0 .and. nstate_lz > 0) then
            nsinglet_lz = nstate_lz !Assume singlet states
         end if
         if ((nsinglet_lz + ntriplet_lz) /= nstate_lz) then
            write (*, *) 'Error(LZ): Sum of singlet and triplet states must give total number of states. Exiting...'
            error = 1
         end if
      end if
      if (energydifthr_lz < 0) then
         write (*, *) 'Parameter energydifthr_lz must be positive number in eV units.'
         error = 1
      end if
      if (nac_accu1 <= 0 .or. nac_accu2 < 0) then
         write (*, *) 'Input error:NACME precision must be a positive integer.'
         write (*, *) 'The treshold is then 10^-(nac_accu).'
         error = 1
      end if
      if (nac_accu1 <= nac_accu2) then
         write (*, *) 'nac_accu1 < nac_accu2'
         write (*, *) 'I will compute NACME only with default accuracy:', nac_accu1
      end if
      if (istage == 1 .and. ipimd /= 1) then
         write (*, *) 'The staging transformation is only meaningful for PIMD'
         error = 1
      end if
      if (inormalmodes > 0 .and. ipimd /= 1) then
         write (*, *) 'The normal mode transformation is only meaningful for PIMD. Exiting...'
         error = 1
      end if
      if (istage == 0 .and. ipimd == 1 .and. inose /= 2 .and. inormalmodes == 0) then
         write (*, *) 'PIMD should be done with staging or normal mode transformation!'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (istage == 1 .and. inose == 2) then
         write (*, *) 'The staging transformation is not compatible with GLE thermostat.'
         error = 1
      end if
      if (nshake /= 0 .and. ipimd == 1) then
         write (*, *) 'PIMD with SHAKE cannot use massive thermostating!'
         error = 1
      end if
      if (nshake /= 0 .and. imasst == 1 .and. inose > 0) then
         write (*, *) 'SHAKE cannot use massive thermostating!'
         write (*, *) 'Set imasst=1 and nmolt, natmolt and nshakemol accordingly.'
         error = 1
      end if
      if ((natmm + natqm /= natom) .and. iqmmm > 0) then
         write (*, *) 'Natmm+natqm not equal to natom!'
         error = 1
      end if

      if (iremd == 1) then
         write (chout, '(A,I2.2)') 'movie.xyz.', my_rank
      else
         chout = 'movie.xyz'
      end if
      inquire (FILE=chout, EXIST=file_exists)
      if (file_exists) then
         if (irest == 0) then
            write (stdout, *) 'File '//trim(chout)//' exists. Please (re)move it or set irest=1.'
            error = 1
         else
            write (stdout, *) 'File "movie.xyz" exists and irest=1.Trajectory will be appended.'
         end if
      end if

      if (iremd == 1) then
         write (chout, '(A,I2.2)') 'restart.xyz.', my_rank
      else
         chout = 'restart.xyz'
      end if
      inquire (FILE=chout, EXIST=file_exists)
      if (file_exists) then
         if (irest == 0) then
            write (*, *) 'File ', trim(chout), ' exists. Please (re)move it or set irest=1.'
            error = 1
         end if
      else
         if (irest == 1) then
            write (*, *) 'File ', trim(chout), ' not found.'
            error = 1
         end if
      end if

      call real_nonnegative(temp, 'temp')
      call real_positive(dt, 'dt')

      call int_switch(irest, 'irest')
      call int_switch(icv, 'icv')
      call int_switch(istage, 'istage')

      call int_nonnegative(nstep, 'nstep')
      call int_nonnegative(imini, 'imini')
      call int_nonnegative(nwrite, 'nwrite')
      call int_nonnegative(nwritex, 'nwritex')
      call int_nonnegative(nwritev, 'nwritev')
      call int_nonnegative(nwritef, 'nwritef')
      call int_nonnegative(nrest, 'nrest')
      call int_nonnegative(narchive, 'narchive')

      call int_positive(nabin, 'nwalk')
      call int_positive(nwalk, 'nwalk')
      call int_positive(ncalc, 'ncalc')

      if (error == 1) then
         write (*, *) 'Input errors were found! Exiting now...'
         call abinerror('check_inputsanity')
      end if
   end subroutine check_inputsanity

   subroutine print_read_error(chfile, chmsg, iost)
      character(len=*), intent(in) :: chmsg, chfile
      integer, intent(in) :: iost
      write (*, *) trim(chmsg)
      write (*, '(A,A)') 'Error when reading file ', trim(chfile)
      write (*, *) 'Error code was', iost
      call abinerror('init')
   end subroutine print_read_error

   subroutine print_basic_info()
      write (stdout, *) 'Reading MD parameters from input file ', trim(chinput)
      write (stdout, *) 'Reading xyz coordinates from file ', trim(chcoords)
      write (stdout, *) 'XYZ Units = '//trim(xyz_units)
      if (chveloc /= '') then
         write (stdout, *) 'Reading initial velocities [a.u.] from file ', trim(chveloc)
      end if
      write (stdout, *) chdivider
      call print_logo()
      write (stdout, *) chdivider
      call print_compile_info()
      write (stdout, *) chdivider
      if (my_rank == 0) call print_runtime_info()
      write (stdout, *) chdivider
      write (stdout, *) '                                              '
      select case (ipimd)
      case (0)
         write (stdout, *) '              Classical MD                    '
      case (1)
         write (stdout, *) '            Path Integral MD                  '
      case (2)
         write (stdout, *) '           Surface Hopping MD                 '
      case (3)
         write (stdout, *) '              Minimization                    '
      case (5)
         write (stdout, *) '             Landau Zener MD                  '
      end select

      write (stdout, *) '    using potential: '//toupper(trim(pot))
      write (stdout, *)
      write (stdout, *) chdivider
   end subroutine print_basic_info

   subroutine print_simulation_parameters()
      write (stdout, *) chdivider
      write (stdout, *) ''
      write (stdout, *) '                SIMULATION PARAMETERS'
      write (stdout, nml=general, delim='APOSTROPHE')
      write (stdout, *)
      write (stdout, nml=system, delim='APOSTROPHE')
      write (stdout, *)
      if (inose >= 1) then
         write (stdout, nml=nhcopt, delim='APOSTROPHE')
         write (stdout, *)
      end if
      if (ipimd == 2) then
         write (stdout, nml=sh, delim='APOSTROPHE')
         write (stdout, *)
      end if
      if (ipimd == 5) then
         write (stdout, nml=lz, delim='APOSTROPHE')
         write (stdout, *)
      end if
      if (pot == 'mm') then
         write (stdout, nml=qmmm, delim='APOSTROPHE')
         write (stdout, *)
      end if
   end subroutine print_simulation_parameters

   subroutine print_logo()
      write (stdout, *) '                  _____     _    _      _ '
      write (stdout, *) '        /\       |  _  \   | |  |  \   | |'
      write (stdout, *) '       /  \      | |_|  |  | |  | | \  | |'
      write (stdout, *) '      / /\ \     |     /   | |  | |\ \ | |'
      write (stdout, *) '     / /__\ \    |=====|   | |  | | \ \| |'
      write (stdout, *) '    / /____\ \   |  _   \  | |  | |  \ | |'
      write (stdout, *) '   / /      \ \  | |_|  |  | |  | |   \  |'
      write (stdout, *) '  /_/        \_\ |_____/   |_|  |_|    \_|'
      write (stdout, *) ' '
      write (stdout, *)
      write (stdout, *) ' D. Hollas, J. Suchan, O. Svoboda, M. Oncak, P. Slavicek'
      write (stdout, *) ' '
   end subroutine print_logo

end subroutine init

   ! Subroutine initialize_masses() populates the global am() array,
   ! based on the atom names from names() array.
   ! User can also specify non-standard isotopes/elements.
   subroutine initialize_masses(masses, massnames, natom)
      use mod_const, only: DP, AMU
      use mod_files, only: stdout
      use mod_error, only: fatal_error
      use mod_system, only: am, names
      real(DP), intent(in) :: masses(:)
      character(len=2), intent(in) :: massnames(:)
      integer, intent(in) :: natom
      character(len=100) :: error_msg
      integer :: i, j
      ! Accurate values for H1 and H2 taken from:
      ! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
      ! Other atomic weights taken from Handbook of Chemistry and Physics, 2013
      ! Original citation: Wieser, M. E., et al., Pure Appl. Chem. 85, 1047, 2013
      allocate (am(natom))
      am = -1.0D0
      do i = 1, natom
         select case (names(i))
         case ('H')
            am(i) = 1.008D0
         case ('H1')
            am(i) = 1.00782503207D0
         case ('H2', 'D')
            am(i) = 2.01410177785D0
         case ('O')
            am(i) = 15.999D0
         case ('S')
            am(i) = 32.06D0
         case ('Se')
            am(i) = 78.971D0
         case ('Te')
            am(i) = 127.60D0
         case ('N')
            am(i) = 14.007D0
         case ('P')
            am(i) = 30.973761998D0
         case ('As')
            am(i) = 74.921595D0
         case ('Sb')
            am(i) = 121.760D0
         case ('Bi')
            am(i) = 208.98040D0
         case ('F')
            am(i) = 18.998403163D0
         case ('Cl')
            am(i) = 35.45D0
         case ('Br')
            am(i) = 79.904D0
         case ('I')
            am(i) = 126.90447D0
         case ('Li')
            am(i) = 6.94D0
         case ('Na')
            am(i) = 22.98976928D0
         case ('K')
            am(i) = 39.0983D0
         case ('Be')
            am(i) = 9.0121831D0
         case ('Mg')
            am(i) = 24.305D0
         case ('Ca')
            am(i) = 40.078D0
         case ('B')
            am(i) = 10.81D0
         case ('Al')
            am(i) = 26.9815385D0
         case ('C')
            am(i) = 12.0110D0
         case ('Si')
            am(i) = 28.085D0
         case ('Ge')
            am(i) = 72.630D0
         case ('Sn')
            am(i) = 118.710D0
         case ('Pb')
            am(i) = 207.2D0
         case ('He')
            am(i) = 4.002602D0
         case ('Ne')
            am(i) = 20.1797D0
         case ('Ar')
            am(i) = 39.948D0
         case ('Kr')
            am(i) = 83.798D0
         case ('Xe')
            am(i) = 131.293D0
         case ('Fe')
            am(i) = 55.845D0
         case ('Ti')
            am(i) = 47.867D0
         case ('V')
            am(i) = 50.9415D0
         case ('Cr')
            am(i) = 51.9961D0
         case ('Mn')
            am(i) = 54.938044D0
         case ('Co')
            am(i) = 58.933194D0
         case ('Ni')
            am(i) = 58.6934D0
         case ('Cu')
            am(i) = 63.546D0
         case ('Zn')
            am(i) = 65.38D0
         case ('Ag')
            am(i) = 107.8682D0
         case ('Au')
            am(i) = 196.966569D0
         case ('Pt')
            am(i) = 195.084D0
         case ('Cd')
            am(i) = 112.414D0
         case ('Hg')
            am(i) = 200.592D0
         case ('U')
            am(i) = 238.02891D0
         case ('Tl')
            am(i) = 204.38D0
         case ('Ba')
            am(i) = 137.327D0
         case ('Ce')
            am(i) = 140.116D0
         case ('Cs')
            am(i) = 132.90545196D0
         case ('Dy')
            am(i) = 162.500D0
         case ('Er')
            am(i) = 167.259D0
         case ('Eu')
            am(i) = 151.964D0
         case ('Gd')
            am(i) = 157.25D0
         case ('Ga')
            am(i) = 69.723D0
         case ('Hf')
            am(i) = 178.49D0
         case ('Ho')
            am(i) = 164.93033D0
         case ('In')
            am(i) = 114.818D0
         case ('Ir')
            am(i) = 192.217D0
         case ('La')
            am(i) = 138.90547D0
         case ('Lu')
            am(i) = 174.9668D0
         case ('Mo')
            am(i) = 95.95D0
         case ('Nd')
            am(i) = 144.242D0
         case ('Nb')
            am(i) = 92.90637D0
         case ('Os')
            am(i) = 190.23D0
         case ('Pd')
            am(i) = 106.42D0
         case ('Pr')
            am(i) = 140.90766D0
         case ('Pa')
            am(i) = 231.03588D0
         case ('Re')
            am(i) = 186.207D0
         case ('Rh')
            am(i) = 102.90550D0
         case ('Rb')
            am(i) = 85.4678D0
         case ('Ru')
            am(i) = 101.07D0
         case ('Sm')
            am(i) = 150.36D0
         case ('Sc')
            am(i) = 44.955908D0
         case ('Sr')
            am(i) = 87.62D0
         case ('Ta')
            am(i) = 180.94788D0
         case ('Tb')
            am(i) = 158.92535D0
         case ('Th')
            am(i) = 232.0377D0
         case ('Tm')
            am(i) = 168.93422D0
         case ('W')
            am(i) = 183.84D0
         case ('Yb')
            am(i) = 173.054D0
         case ('Y')
            am(i) = 88.90584D0
         case ('Zr')
            am(i) = 91.224D0
         case DEFAULT
            write (stdout, *) 'Atom name ', names(i), ' was not found in the library.'
            write (stdout, *) 'Using user-defined mass from input file'
         end select
      end do

      ! Check for duplicate user defined atom names
      do i = 1, size(massnames)
         do j = i + 1, size(massnames)
            if (massnames(i) == massnames(j) .and. trim(massnames(i)) /= '') then
               error_msg = 'Duplicate atom names in input array "massnames"'
               call fatal_error(__FILE__, __LINE__, error_msg)
               return
            end if
         end do
      end do

      ! Here we overwrite library values if user provided alternative value
      ! in the input file, or we define new atom names that are not in the library.
      do i = 1, natom
         do j = 1, size(massnames)
            if (names(i) == massnames(j)) then

               if (masses(j) <= 0.0D0) then
                  call fatal_error(__FILE__, __LINE__, &
                     & 'Mass cannot be negative. Fix arrays "masses" in the input file.')
                  return
               end if

               write (stdout, *) 'Defining new atom ', names(i), ' with mass=', masses(j)

               am(i) = masses(j)

            end if
         end do
      end do

      ! Catch any undefined masses
      do i = 1, natom
         if (am(i) <= 0) then
            write (error_msg, '(a,i0,a)') 'Atomic mass for atom '//names(i)//' was not specified'
            call fatal_error(__FILE__, __LINE__, error_msg)
            return
         end if
      end do

      if (natom <= 20) then
         write (stdout, *) ''
         write (stdout, *) '                        ATOMIC MASSES / a.u.'
         do i = 1, natom
            write (stdout, '(A2, A1)', advance='no') names(i), ' '
         end do
         write (stdout, *)
         write (stdout, *) 'The corresponding relative atomic masses are:'
         write (stdout, *) (am(i), i=1, natom)
         write (stdout, *) ''
      end if

      ! Finally, convert masses to atomic units
      am = am * AMU

   end subroutine initialize_masses

! NOTE: This subroutine is outside of init()
! due to a conflict of 'call system()` with namelist 'system',
! which some compilers do not like.
subroutine print_runtime_info()
   use mod_files, only: stdout
   character(len=1024) :: cmdline
   write (stdout, *) ''
   write (stdout, *) '          RUNTIME INFO'
   write (stdout, *) ' '
   write (stdout, *) "Running on node: "
   call system('uname -n')
   write (stdout, '(A)') 'Working directory: '
   call system('pwd')
   write (stdout, *)
   call get_command(cmdline)
   write (stdout, *) trim(cmdline)
   call flush (stdout)
   call get_command_argument(0, cmdline)
   write (stdout, *)
   call system('ldd '//cmdline)
   write (stdout, *) ''
end subroutine print_runtime_info

end module mod_init

! We cannot include finish in the module, since it would
! result in circular dependencies.
subroutine finish(error_code)
   use mod_arrays, only: deallocate_arrays
   use mod_general, only: pot, pot_ref, ipimd, inormalmodes
   use mod_files, only: close_files
   use mod_nhc, only: inose, finalize_nhc
   use mod_gle, only: finalize_gle, finalize_pile
   use mod_lz, only: lz_finalize
   use mod_transform, only: finalize_normalmodes
   use mod_cp2k, only: finalize_cp2k
   use mod_plumed, only: iplumed, finalize_plumed
   use mod_terampi, only: finalize_terachem
   use mod_terampi_sh, only: finalize_terash
   use mod_splined_grid, only: finalize_spline
   use mod_mpi, only: finalize_mpi
   implicit none
   integer, intent(in) :: error_code

   if (pot == '_tera_' .or. pot_ref == '_tera_') then
      if (ipimd == 2) then
         call finalize_terash()
      end if
      call finalize_terachem(error_code)
   end if

   call deallocate_arrays()
   call close_files()

   if (inormalmodes > 0) then
      call finalize_normalmodes()
   end if

   if (inose == 1) then
      call finalize_nhc()
   else if (inose == 2 .or. inose == 4) then
      call finalize_gle()
   else if (inose == 3) then
      call finalize_pile()
   end if

   if (iplumed == 1) then
      call finalize_plumed()
   end if

   ! Cleanup Landau-Zener
   if (ipimd == 5) then
      call lz_finalize()
   end if

   if (pot == '_splined_grid_') then
      call finalize_spline()
   end if
   ! MPI_FINALIZE is called in this routine as well
   if (pot == '_cp2k_') then
      call finalize_cp2k()
   end if

   ! This must be the last call
   call finalize_mpi(error_code)

end subroutine finish
