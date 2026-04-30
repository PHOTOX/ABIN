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
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
   use mod_utils, only: toupper, tolower, normalize_atom_name, &
                      & open_file_for_reading
   implicit none
   private
   public :: init
   ! For unit test
   public :: initialize_masses
   public :: read_xyz_file, read_atom_names
contains

   subroutine init(dt)
      use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
      use mod_const, only: ANG, AUtoK
      use mod_interfaces, only: print_compile_info, omp_set_num_threads
      use mod_cmdline, only: get_cmdline
      use mod_mpi, only: initialize_mpi, get_mpi_size, get_mpi_rank, mpi_barrier_wrapper
      use mod_files, only: files_init, stdout_to_devnull
      use mod_arrays
      use mod_general
      use mod_system, only: set_atom_names, conatom, am, f, dime
      use mod_nhc
      use mod_estimators
      use mod_potentials
      use mod_potentials_sh
      use mod_sh_integ, only: phase, nstate
      use mod_sh, only: read_sh_input, print_sh_input
      use mod_lz, only: read_lz_input, print_lz_input, lz_init
      use mod_qmmm, only: natqm, natmm
      use mod_force_mm, only: initialize_mm
      use mod_force_h2o, only: initialize_h2o_pot, h2opot
      use mod_gle
      use mod_sbc, only: sbc_init, rb_sbc, kb_sbc, isbc, rho
      use mod_prng_init, only: initialize_prng
      use mod_splined_grid, only: initialize_spline, potential_file
      use mod_utils, only: append_rank, file_exists_or_exit
      use mod_vinit
      use mod_analyze_geometry
      use mod_shake
      use mod_minimize, only: gamm, gammthr
      use mod_analysis, only: restin
      use mod_water, only: watpot, check_water
      use mod_plumed, only: iplumed, plumedfile, plumed_init
      use mod_en_restraint, only: en_rest_init, en_diff, en_kk, restrain_pot
      use mod_transform, only: initialize_pi_masses
      use mod_cp2k
      use mod_remd
      use mod_force_tcpb, only: initialize_tcpb
      use mod_terampi
      use mod_terampi_sh
      use mod_mace_mpi
      use mod_mdstep, only: initialize_integrator, nabin, nstep_ref
      real(DP), intent(out) :: dt
      ! Input parameters for analytical potentials
      ! TODO: Initialize these variable in the code not here
      real(DP), save :: lambda_dw = -1.0D0, D0_dw = -1.0D0, k_dw = -1.0D0, r0_dw = -1.D0
      real(DP), save :: r0_morse = -1, d0_morse = -1, k_morse = -1
      real(DP), save :: k = -1, r0 = -1
      real(DP), save :: kx = -1, ky = -1, kz = -1
      ! Lennard-Jones parameteres and Coulomb charges (pot=_mm_)
      ! All input parameters are expected to be in atomic units,
      ! except LJ_rmin which should be in angstroms.
      ! User-specified atomic types read from the input file
      character(len=2), allocatable :: mm_types(:)
      ! Coulomb charges correspoding to the atomic types defined above
      real(DP), allocatable :: q(:)
      ! L-J parameters
      real(DP), allocatable :: LJ_rmin(:), LJ_eps(:)
      ! Initial temperature (read from namelist nhcopt)
      real(DP), save :: temp0 = -1
      ! User-defined masses in relative atomic units
      real(DP), allocatable :: masses(:)
      integer :: iw, iat, natom_xyz, iost
      integer :: shiftdihed
      ! Random number seed
      ! Negative value means we get the seed from /dev/urandom
      integer, save :: irandom = -1
      ! Number of OpenMP processes, read from ABIN input
      ! WARNING: We do NOT use OMP_NUM_THREADS environment variable!
      integer :: nproc
      integer :: getPID
!$    integer, external :: omp_get_max_threads
      character(len=2), allocatable :: atnames(:)
      character(len=2), allocatable :: massnames(:)
      character(len=200) :: chinput, chcoords, chveloc
      character(len=200) :: chiomsg, chout
      character(len=20) :: xyz_units
      character(len=60) :: chdivider
      character(len=60) :: mdtype
      character(len=60) :: therm
      character(len=1024) :: tc_server_name, tcpb_input_file, tcpb_host
      integer :: tcpb_port
      logical :: file_exists
      logical :: rem_comvel, rem_comrot
      integer :: my_rank, mpi_world_size
      integer :: uinput, ucoord, uvel

      ! ABIN input parameters are read from the input file (default 'input.in')
      ! in the form of the standard Fortran namelist syntax.
      ! The input parameters are grouped in several different namelists:
      !
      ! - general:      Most basic MD settings + misc
      ! - thermostat:   parameters for thermostats (previously nhcopt)
      ! - remd:         parameters for Replica Exchange MD
      ! - sh:           parameters for Surface Hopping, moved to mod_sh module
      ! - system:       system-specific parameters for model potentials, masses, SHAKE constraints...
      ! - qmmm:         parameters for internal QMMM (not really tested).
      !
      ! All namelists need to be in a single input file, and the code
      ! in this subroutine must ensure that the namelists can be in any order.
      namelist /general/ pot, ipimd, mdtype, istage, inormalmodes, nwalk, nstep, icv, ihess, imini, nproc, iqmmm, &
         nwrite, nwritex, nwritev, nwritef, dt, irandom, nabin, irest, nrest, anal_ext, &
         isbc, rb_sbc, kb_sbc, gamm, gammthr, conatom, mpi_milisleep, narchive, xyz_units, &
         dime, idebug, enmini, rho, iknow, watpot, h2opot, iremd, iplumed, plumedfile, &
         en_restraint, en_diff, en_kk, restrain_pot, &
         pot_ref, nstep_ref, nteraservers, max_mpi_wait_time, cp2k_mpi_beads

      namelist /remd/ nswap, nreplica, deltaT, Tmax, temp_list

      ! new namelist section for thermostat, previously nhcopt
      namelist /thermostat/ inose, therm, temp, temp0, nchain, tau0, tau0_langevin, imasst, nrespnose, nyosh, &
         scaleveloc, readNHC, nmolt, natmolt, nshakemol, rem_comrot, rem_comvel, gle_test

      ! old namelist section for thermostat, kept for backwards compatibility
      namelist /nhcopt/ inose, therm, temp, temp0, nchain, tau0, tau0_langevin, imasst, nrespnose, nyosh, &
         scaleveloc, readNHC, nmolt, natmolt, nshakemol, rem_comrot, rem_comvel, gle_test

      namelist /system/ masses, massnames, ndist, dist1, dist2, &
         nang, ang1, ang2, ang3, ndih, dih1, dih2, dih3, dih4, shiftdihed, &
         k, r0, kx, ky, kz, r0_morse, d0_morse, k_morse, D0_dw, lambda_dw, k_dw, r0_dw, &
         Nshake, ishake1, ishake2, shake_tol, potential_file

      namelist /qmmm/ natqm, natmm, q, LJ_rmin, LJ_eps, mm_types

      namelist /mace/ mace_model, mace_device, mace_default_dtype, &
         mace_batch_size, mace_compute_stress, mace_return_contributions, &
         mace_info_prefix, mace_head, mace_max_mpi_wait_time, mace_mpi_milisleep

      chcoords = 'mini.dat'
      xyz_units = 'angstrom'
      chinput = 'input.in'
      chveloc = ''
      tc_server_name = ''
      tcpb_host = 'localhost'
      tcpb_input_file = ''
      tcpb_port = -1
      mdtype = ''
      therm = ''
      dt = -1
      nproc = 1
      iplumed = 0
      shiftdihed = 1

      chdivider = "######################################################"

      call get_cmdline(chinput, chcoords, chveloc, tc_server_name, tcpb_input_file, tcpb_host, tcpb_port)

      ! Reading main input parameters from namelist &general
      open (newunit=uinput, file=chinput, status='OLD', delim='APOSTROPHE', action="READ")
      read (uinput, general)
      rewind (uinput)

      pot = tolower(pot)
      pot_ref = tolower(pot_ref)

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
!$    call omp_set_num_threads(nproc)

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

      ! Connect to MACE server as early as possible
      if (pot == '_mace_') then
         call initialize_mace_interface()
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
         chcoords = append_rank(chcoords)
         if (chveloc /= '') then
            chveloc = append_rank(chveloc)
         end if
      end if

      call file_exists_or_exit(chcoords)
      if (chveloc /= '') then
         call file_exists_or_exit(chveloc)
      end if

      call print_basic_info()

      ! Initialize pseudo-random number generator
      call initialize_prng(seed=irandom, mpi_rank=my_rank)

      ! Get number of atoms and atom names from XYZ coordinates NOW
      ! so that we can allocate arrays
      call read_atom_names(chcoords, natom_xyz, atnames)

      ! Set global variable "natom"
      call set_natom(natom_xyz)
      ! Set global variable "names"
      call set_atom_names(atnames, natom_xyz)

      ! This line is super important,
      ! cause we actually use natqm in many parts of the code
      if (iqmmm == 0 .and. pot /= '_mm_') then
         natqm = natom_xyz
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

      ! We have these allocate these arrays here
      ! because we read them from input file.
      allocate (massnames(natom))
      allocate (masses(natom), source=-1.0_DP)

      massnames = ''

      ! Lennard-Jones / Coulomb parameters
      allocate (mm_types(natom))
      allocate (LJ_rmin(natom), source=-1.0D0)
      allocate (q(natom), source=0.0D0)
      allocate (LJ_eps(natom), source=-1.0D0)
      mm_types = ''

      allocate (ishake1(natom * 3 - 6), source=0)
      allocate (ishake2(natom * 3 - 6), source=0)

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

      ! Read initial geometry
      ucoord = open_file_for_reading(chcoords)
      call read_xyz_file(ucoord, chcoords, atnames, natom_xyz, 1, x, y, z)
      close (ucoord)

      if (tolower(trim(xyz_units)) == "angstrom") then
         x(:, 1) = x(:, 1) * ANG
         y(:, 1) = y(:, 1) * ANG
         z(:, 1) = z(:, 1) * ANG
      else if (tolower(trim(xyz_units)) == "bohr") then
         continue
      else
         write (*, *) 'ERROR: Wrong XYZ units: ', trim(xyz_units)
      end if

      ! Initialize all PIMD beads to the same initial geometry
      do iw = 1, nwalk
         x(:, iw) = x(:, 1)
         y(:, iw) = y(:, 1)
         z(:, iw) = z(:, 1)
      end do

      ! the namelist system does not need to be present
      read (uinput, system, iostat=iost, iomsg=chiomsg)
      rewind (uinput)
      !check, whether we hit End-Of-File or other error
      if (IS_IOSTAT_END(iost)) then !fortran intrinsic for EOF
         continue
      else if (iost /= 0) then
         write (stderr, *) chiomsg
         call fatal_error(__FILE__, __LINE__, 'Could not read namelist "system".')
      end if

      do iat = 1, size(massnames)
         massnames(iat) = normalize_atom_name(massnames(iat))
      end do

      ! Determine atomic masses from periodic table
      ! names - array of atomic names read from the XYZ structure
      ! masses - user-defined atomic masses
      ! massnames - atomic symbols of user-defined masses
      ! am - the output array of atomic masses in atomic units
      call initialize_masses(atnames, masses, massnames, natom, am)

      ! Transform masses for PIMD
      ! Note that amt array is used throughout the code
      ! even for non-PI simulations.
      call initialize_pi_masses(am, amg, amt)

      allocate (natmolt(natom))
      allocate (nshakemol(natom))
      natmolt = 0
      natmolt(1) = natom ! default for global NHC thermostat
      nshakemol = 0

      ! read &thermostat section
      read (uinput, thermostat, iostat=iost)
      rewind (uinput)
      ! if &thermostat is not found, try the older &nhcopt
      if (iost /= 0) then
         read (uinput, nhcopt)
         rewind (uinput)
      end if

      if (ipimd == 2) then
         call read_sh_input(uinput)
      end if

      if (ipimd == 5) then
         call read_lz_input(uinput)
         call lz_init(pot)
      end if

      if (iremd == 1) then
         read (uinput, remd)
         rewind (uinput)
         call remd_init(temp, temp0)
      end if

      if (iqmmm > 0 .or. pot == '_mm_') then
         read (uinput, qmmm)
         rewind (uinput)
      end if

      ! Read &mace namelist when using MACE potential
      if (pot == '_mace_') then
         read (uinput, mace, iostat=iost, iomsg=chiomsg)
         rewind (uinput)
         if (iost /= 0) then
            write (stderr, *) chiomsg
            call fatal_error(__FILE__, __LINE__, &
               & 'Could not read namelist "mace". Required when pot="_mace_".')
         end if
      end if

      close (uinput)
      ! END OF READING INPUT

      if (pot == '_tera_' .or. restrain_pot == '_tera_') then
         call initialize_tc_servers()
         if (ipimd == 2 .or. ipimd == 5) then
            call init_terash(x, y, z)
         end if
      end if

      if (pot == '_mace_') then
         call initialize_mace_server()
      end if

      if (pot == '_tcpb_' .or. restrain_pot == '_tcpb_' .or. pot_ref == '_tcpb_') then
         call initialize_tcpb(natqm, atnames, tcpb_port, tcpb_host, tcpb_input_file)
      end if

      ! first check if 'therm' keyword was used and then transfer to inose
      if (therm /= '') then
         therm = tolower(therm)
         select case (therm)
         case ('none')
            inose = 0
         case ('nhc')
            inose = 1
         case ('gle')
            inose = 2
         case ('langevine')
            inose = 3
         case default
            call fatal_error(__FILE__, __LINE__, 'invalid thermostat (therm) in '//trim(chinput))
         end select
      end if

      ! Check whether input parameters make sense
      call check_inputsanity()

      call initialize_integrator(dt, ipimd, inormalmodes, nshake, pot_ref, pot_ref)

      if (iplumed == 1) then
         call plumed_init(natom, irest, dt0, nrest)
      end if

      if (temp0 < 0) then
         temp0 = temp
      end if
      write (stdout, *) 'Initial temperature [K] =', temp0
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
         call check_water(natom_xyz, atnames)
      end if

      ! Initialize thermostat
      if (inose == 1) then
         call nhc_init()
      else if (inose == 2) then
         call gle_init(dt * 0.5D0 / nabin / nstep_ref) !nabin is set to 1 unless ipimd=1
      else if (inose == 3) then
         call pile_init(dt * 0.5D0, tau0_langevin)
      else if (inose == 4) then
         call gle_init(dt * 0.5D0 / nstep_ref)
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
         call vinit(temp0, am, vx, vy, vz)
      end if

      ! Read velocities from file (optional)
      if (chveloc /= '' .and. irest == 0) then
         uvel = open_file_for_reading(chveloc)
         do iw = 1, nwalk
            call read_xyz_file(uvel, chveloc, atnames, natom_xyz, iw, vx, vy, vz)
         end do
         close (uvel)
      end if

      ! Doing this here so that we can do it even when reading velocities from file
      if (rem_comvel) then
         call remove_comvel(vx, vy, vz, am)
      end if
      if (rem_comrot) then
         call remove_rotations(x, y, z, vx, vy, vz, am)
      end if

      if (conatom > 0) then
         call constrainP(vx, vy, vz, conatom)
      end if

      ! If scaleveloc=1, scale initial velocitites to match the temperature
      ! Otherwise, just print the temperature.
      call scale_velocities(vx, vy, vz)

      ! Initialize spherical boundary onditions
      if (isbc == 1) then
         call sbc_init(x, y, z)
      end if

      if (en_restraint >= 1) then
         call en_rest_init(natom)
      end if

      ! Initialize in-built analytical potentials
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
      if (pot == '_h2o_' .or. pot_ref == '_h2o_') then
         call initialize_h2o_pot(natom, atnames)
      end if
      if (pot == '_mm_' .or. pot_ref == '_mm_') then
         call initialize_mm(natom, atnames, mm_types, q, LJ_rmin, LJ_eps)
      end if
      if (pot == '_nai_' .or. pot_ref == '_nai_') then
         call nai_init(natom, nwalk, ipimd, nstate)
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

!$    write (stdout, '(A,I0)') 'Number of OpenMP threads: ', omp_get_max_threads()

      ! Open permanent files for writing
      call files_init(isbc, phase, ndist, nang, ndih)

      call flush (OUTPUT_UNIT)

   contains

      subroutine check_inputsanity()
         use mod_chars, only: chknow
         use mod_utils, only: real_positive, real_nonnegative, &
                            & int_positive, int_nonnegative, int_switch
         integer :: error
!$       integer :: nthreads, omp_get_max_threads

         error = 0

         !  We should exclude all non-abinitio options, but whatever....
!$       nthreads = omp_get_max_threads()
!$       if (nthreads > 1 .and. (ipimd /= 1 .and. pot /= '_cp2k_')) then
!$          write (*, *) 'Number of threads is ', nthreads
!$          call fatal_error(__FILE__, __LINE__, &
!$             & 'Parallel execution is currently only supported with ab initio PIMD (ipimd=1)')
!$       end if

         if (nproc > 1) then
!$          if (.false.) then
               write (*, *) 'FATAL ERROR: This executable was not compiled with parallel support.'
               error = 1
!$          end if
         end if

         if (irest == 1 .and. chveloc /= '') then
            write (stderr, *) 'WARNING: Input velocities from file '//trim(chveloc)//' will be ignored!'
            write (stderr, *) 'Velocities will be taken from restart file because irest=1.'
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
            write (*, *) 'Number of walkers for PIMD (nwalk) mus be >= 1!'
            write (*, *) 'Either set ipimd=0 for classical simulation'
            write (*, *) 'or set nwalk > 1'
            error = 1
         end if
         if (iqmmm < 0 .or. iqmmm > 1) then
            write (*, *) 'Error: iqmmm must be 0 or 1 for ONIOM.'
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

         if (ipimd == 1 .and. inose <= 0) then
            write (*, *) 'You have to use thermostat with PIMD! (inose>=0)'
            write (*, *) chknow
            if (iknow /= 1) error = 1
         end if
         if (ipimd < 0 .or. ipimd > 5) then
            write (*, *) 'Invalid ipimd value'
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
         if (nshake /= 0 .and. (inose == 2 .or. inose == 4)) then
            write (*, *) 'SHAKE is not compatible with GLE thermostat!'
            error = 1
         end if
         if (nshake /= 0 .and. inose == 3) then
            write (*, *) 'SHAKE is currently not compatible with Langeving thermostat!'
            error = 1
         end if
         if ((natmm + natqm /= natom)) then
            write (*, *) 'Natmm+natqm /= natom!'
            error = 1
         end if

         chout = append_rank('movie.xyz')
         inquire (file=chout, exist=file_exists)
         if (file_exists) then
            if (irest == 0) then
               write (stderr, *) 'File '//trim(chout)//' exists. Please (re)move it or set irest=1.'
               error = 1
            else
               write (stdout, *) 'File '//trim(chout)//' exists and irest=1. Trajectory will be appended.'
            end if
         end if

         chout = append_rank('restart.xyz')
         inquire (file=chout, exist=file_exists)
         if (file_exists .and. irest == 0) then
            write (stderr, *) 'File ', trim(chout), ' exists. Please (re)move it or set irest=1.'
            error = 1
         end if
         if (irest == 1 .and. .not. file_exists) then
            write (stderr, *) 'File ', trim(chout), ' not found.'
            error = 1
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

         call int_positive(nwalk, 'nwalk')

         if (error == 1) then
            call fatal_error(__FILE__, __LINE__, 'Invalid input parameters')
         end if
      end subroutine check_inputsanity

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
            call print_sh_input()
            write (stdout, *)
         end if
         if (ipimd == 5) then
            call print_lz_input()
            write (stdout, *)
         end if
         if (iqmmm > 0 .or. pot == '_mm_') then
            write (stdout, nml=qmmm, delim='APOSTROPHE')
            write (stdout, *)
         end if
         if (pot == '_mace_') then
            write (stdout, nml=mace, delim='APOSTROPHE')
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
         write (stdout, *) ' D. Hollas, J. Suchan, J. Janos, M. Oncak and P. Slavicek'
         write (stdout, *) ' '
         write (stdout, *) ' with contributions by O. Svoboda, S. Srsen, Jan Postulka,'
         write (stdout, *) '      V. Juraskova, M. Barnfield, J. Chalabala and others'
         write (stdout, *) ' '
      end subroutine print_logo

   end subroutine init

   ! Read atom names from XYZ file, return number of atoms
   subroutine read_atom_names(coordfile, natom_xyz, atnames)
      character(len=*), intent(in) :: coordfile
      integer, intent(out) :: natom_xyz
      character(len=2), allocatable, intent(out) :: atnames(:)
      integer :: u, iost, iat
      real(DP) :: xtmp, ytmp, ztmp

      u = open_file_for_reading(coordfile)

      read (u, '(I50)', iostat=iost) natom_xyz
      if (iost /= 0) then
         close (u)
         call fatal_error(__FILE__, __LINE__, &
            &'Could not read number of atoms from the first line of file '//trim(coordfile))
         return
      end if

      if (natom_xyz < 1) then
         close (u)
         call fatal_error(__FILE__, __LINE__, &
            & 'Invalid number of atoms on the first line of the XYZ file '//trim(coordfile))
         return
      end if

      allocate (atnames(natom_xyz))
      atnames = ''

      ! Ignore comment line
      read (u, *)

      do iat = 1, natom_xyz
         ! Ignore the positions for now, only read atom names
         read (u, *, iostat=iost) atnames(iat), xtmp, ytmp, ztmp
         if (iost /= 0) then
            close (u)
            call fatal_error(__FILE__, __LINE__, &
               & 'Invalid line in file '//trim(coordfile))
            return
         end if
         atnames(iat) = normalize_atom_name(atnames(iat))
      end do
      close (u)
   end subroutine read_atom_names

   subroutine read_xyz_file(u, fname, atnames, num_atom, iw, x, y, z)
      integer, intent(in) :: u
      character(len=*), intent(in) :: fname
      character(len=2), intent(in) :: atnames(:)
      integer, intent(in) :: num_atom
      ! Bead index
      integer, intent(in) :: iw
      real(DP), dimension(:, :), intent(inout) :: x, y, z
      character(len=2) :: atom
      integer :: natom_xyz
      integer :: iat, iost

      ! Due to a bug in GFortran 7.0, we cannot use this
      ! so we are passing fname explicitly
      ! inquire (unit=u, name=fname)

      ! Verify that number of atoms matches what we expect
      read (u, *, IOSTAT=iost) natom_xyz
      if (iost /= 0) then
         close (u)
         call fatal_error(__FILE__, __LINE__, &
            &'Could not read number of atoms from the first line of file '//trim(fname))
         return
      end if

      if (natom_xyz /= num_atom) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Inconsistent number of atoms on the first line of file '//trim(fname))
         return
      end if

      ! Skip comment line
      read (u, *)

      do iat = 1, num_atom
         read (u, *, iostat=iost) atom, x(iat, iw), y(iat, iw), z(iat, iw)
         if (iost /= 0) then
            call fatal_error(__FILE__, __LINE__, 'Invalid line in file '//trim(fname))
            return
         end if
         if (normalize_atom_name(atom) /= atnames(iat)) then
            write (stderr, *) 'Offending line:'
            write (stderr, *) atom, x(iat, iw), y(iat, iw), z(iat, iw)
            call fatal_error(__FILE__, __LINE__, 'Inconsistent atom type in file '//trim(fname))
            return
         end if
      end do
   end subroutine read_xyz_file

   ! Subroutine initialize_masses() populates the global am() array,
   ! based on the atom names from names() array.
   ! User can also specify non-standard isotopes/elements.
   subroutine initialize_masses(names, masses, massnames, natom, am)
      use mod_const, only: DP, AMU
      use mod_files, only: stdout
      ! Atomic names of simulated structure
      character(len=2), intent(in) :: names(natom)
      real(DP), intent(in) :: masses(:)
      character(len=2), intent(in) :: massnames(:)
      integer, intent(in) :: natom
      real(DP), allocatable, intent(out) :: am(:)
      character(len=100) :: error_msg
      integer :: i, j

      allocate (am(natom))
      am = -1.0D0

      ! Accurate values for H1 and H2 taken from:
      ! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
      ! Other atomic weights taken from Handbook of Chemistry and Physics, 2013
      ! Original citation: Wieser, M. E., et al., Pure Appl. Chem. 85, 1047, 2013
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
            ! The following values were taken from QCelemental Python library
         case ('Tc')
            am(i) = 97.9072124D0
         case ('Pm')
            am(i) = 144.9127559D0
         case ('Po')
            am(i) = 208.9824308D0
         case ('At')
            am(i) = 209.9871479D0
         case ('Rn')
            am(i) = 222.0175782D0
         case ('Fr')
            am(i) = 223.019736D0
         case ('Ra')
            am(i) = 226.0254103D0
         case ('Ac')
            am(i) = 227.0277523D0
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

   subroutine print_runtime_info()
      use mod_files, only: stdout
      character(len=1024) :: cmdline
      write (stdout, *) ''
      write (stdout, *) '          RUNTIME INFO'
      write (stdout, *) ' '
      write (stdout, *) "Running on node: "
      call execute_command_line('uname -n')
      write (stdout, '(A)') 'Working directory: '
      call execute_command_line('pwd')
      write (stdout, *)
      call get_command(cmdline)
      write (stdout, *) trim(cmdline)
      call flush (stdout)
      call get_command_argument(0, cmdline)
      write (stdout, *)
      call execute_command_line('ldd '//cmdline)
      write (stdout, *) ''
   end subroutine print_runtime_info

end module mod_init

! We cannot include finish in the module, since it would
! result in circular dependencies.
! allow(procedure-not-in-module) ! fortitude linter
subroutine finish(error_code)
   use mod_arrays, only: deallocate_arrays
   use mod_general, only: pot, pot_ref, ipimd, inormalmodes, en_restraint
   use mod_files, only: close_files
   use mod_nhc, only: inose, finalize_nhc
   use mod_gle, only: finalize_gle, finalize_pile
   use mod_lz, only: lz_finalize
   use mod_en_restraint, only: en_rest_finalize
   use mod_transform, only: finalize_normalmodes
   use mod_cp2k, only: finalize_cp2k
   use mod_plumed, only: iplumed, finalize_plumed
   use mod_terampi, only: finalize_terachem
   use mod_terampi_sh, only: finalize_terash
   use mod_mace_mpi, only: finalize_mace
   use mod_splined_grid, only: finalize_spline
   use mod_force_mm, only: finalize_mm
   use mod_force_tcpb, only: finalize_tcpb
   use mod_mpi, only: finalize_mpi
   implicit none
   integer, intent(in) :: error_code

   if (pot == '_tera_' .or. pot_ref == '_tera_') then
      if (ipimd == 2) then
         call finalize_terash()
      end if
      call finalize_terachem(error_code)
   end if

   if (pot == '_mace_') then
      call finalize_mace(error_code)
   end if

   if (pot == '_tcpb_') then
      call finalize_tcpb()
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

   if (en_restraint >= 1) then
      call en_rest_finalize()
   end if

   if (pot == '_splined_grid_') then
      call finalize_spline()
   else if (pot == '_mm_') then
      call finalize_mm()
   end if

   ! MPI_FINALIZE is called in this routine as well
   if (pot == '_cp2k_') then
      call finalize_cp2k()
   end if

   ! This must be the last call
   call finalize_mpi(error_code)

end subroutine finish
