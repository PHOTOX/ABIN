! This messy function performs many things, among others:
! 1. Reading input
! 2. Input sanity check
! 3. Allocation of arrays
! 4. Reading restart OR reading input geometry
! 5. Initialize velocities
! 6. Initialize everything else
! At this moment, coordinate and velocity transformations are NOT performed here
! Surface hopping is NOT initialized here.

subroutine init(dt)
   use mod_const
   use mod_interfaces, only: print_compile_info
   use mod_cmdline, only: get_cmdline
   use mod_files
   use mod_arrays
   use mod_array_size
   use mod_general
   use mod_system
   use mod_nhc
   use mod_estimators
   use mod_harmon
   use mod_sh_integ, only: nstate, integ, phase, popsumthr, correct_decoherence
   use mod_sh
   use mod_lz
   use mod_qmmm, only: natqm, natmm
   use mod_force_mm
   use mod_gle
   use mod_sbc, only: sbc_init, rb_sbc, kb_sbc, isbc, rho
   use mod_random
   use mod_splined_grid, only: initialize_spline
   use mod_utils, only: lowertoupper, uppertolower, file_exists_or_exit
   use mod_vinit
   use mod_analyze_geometry
   use mod_shake
   use mod_minimize, only: gamm, gammthr
   use mod_analysis, only: restin
   use mod_water, only: watpot, check_water
   use mod_plumed, only: iplumed, plumedfile, plumed_init
   use mod_en_restraint
   use mod_transform, only: init_mass
   use mod_cp2k
   use mod_remd
   use mod_terampi
   use mod_terampi_sh
#ifdef USE_MPI
   use mpi, only: MPI_COMM_WORLD, MPI_Init, MPI_Comm_Rank, MPI_Comm_Size, MPI_Barrier
#endif
   implicit none
   real(DP), intent(out) :: dt
   real(DP) :: masses(MAXTYPES)
   real(DP) :: rans(10)
   integer :: iw, iat, natom_xyz, imol, shiftdihed = 1, iost
   integer :: error, getpid, nproc = 1, ipom
   character(len=2) :: massnames(MAXTYPES), atom
   character(len=200) :: chinput, chcoords, chveloc
   character(len=200) :: chiomsg, chout
   character(len=20) :: xyz_units = 'angstrom'
   character(len=60) :: chdivider
   character(len=60) :: mdtype
   logical :: file_exists
   logical :: rem_comvel, rem_comrot
!  real(DP) :: wnw=5.0d-5
   ! Used for MPI calls
   integer :: ierr
   integer :: irand
!$ integer :: nthreads, omp_get_max_threads
!  wnw "optimal" frequency for langevin (inose=3)

   namelist /general/ natom, pot, ipimd, mdtype, istage, inormalmodes, nwalk, nstep, icv, ihess, imini, nproc, iqmmm, &
      nwrite, nwritex, nwritev, nwritef, dt, irandom, nabin, irest, nrest, anal_ext, &
      isbc, rb_sbc, kb_sbc, gamm, gammthr, conatom, mpi_sleep, narchive, xyz_units, &
      dime, ncalc, idebug, enmini, rho, iknow, watpot, iremd, iplumed, plumedfile, &
      en_restraint, en_diff, en_kk, restrain_pot, &
      pot_ref, nstep_ref, nteraservers, cp2k_mpi_beads

#ifdef USE_MPI
   namelist /remd/ nswap, nreplica, deltaT, Tmax, temp_list
#endif

   namelist /nhcopt/ inose, temp, temp0, nchain, ams, tau0, tau0_langevin, imasst, nrespnose, nyosh, &
      scaleveloc, readNHC, readQT, initNHC, nmolt, natmolt, nshakemol, rem_comrot, rem_comvel

   namelist /system/ masses, massnames, ndist, dist1, dist2, &
      nang, ang1, ang2, ang3, ndih, dih1, dih2, dih3, dih4, shiftdihed, &
      k, r0, k1, k2, k3, De, a, D0_dw, lambda_dw, k_dw, r0_dw, &
      Nshake, ishake1, ishake2, shake_tol

   namelist /sh/ istate_init, nstate, substep, deltae, integ, inac, nohop, phase, decoh_alpha, popthr, ignore_state, &
      nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, adjmom, revmom, natmm_tera, &
      dE_S0S1_thr, correct_decoherence

   namelist /lz/ initstate_lz, nstate_lz, nsinglet_lz, ntriplet_lz, deltaE_lz, energydifthr_lz

   namelist /qmmm/ natqm, natmm, q, rmin, eps, attypes

   chcoords = 'mini.dat'
   chinput = 'input.in'
   chveloc = ''
   mdtype = ''
   dt = -1
   error = 0
   iplumed = 0

   chdivider = "######################################################"

   call get_cmdline(chinput, chcoords, chveloc, teraport)

   ! READING MAIN INPUT
   open (150, file=chinput, status='OLD', delim='APOSTROPHE', action="READ")
   read (150, general)
   rewind (150)
   pot = UpperToLower(pot)

   if (pot == 'splined_grid') then
      natom = 1
      dime = 1
      f = 0
      call initialize_spline()
   end if

   if (pot == "_cp2k_" .or. pot_ref == "_cp2k_") then
      call init_cp2k()
#ifdef USE_MPI
   else
      call MPI_INIT(ierr)
      if (ierr /= 0) then
         write (*, *) 'Bad signal from MPI_INIT:', ierr
         stop 1
      end if
#endif
   end if

   ! We need to connect to TeraChem as soon as possible
   ! because we want to shut down TeraChem nicely in case something goes wrong.
#ifdef USE_MPI
   ! TODO: Add explicit checks for ierr for all MPI calls!
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)
   if (my_rank == 0 .and. mpi_world_size > 1) then
      write (*, '(A,I3)') 'Number of MPI processes = ', mpi_world_size
   end if
   if (pot == '_tera_' .or. restrain_pot == '_tera_') then
      if (nwalk > 1) then
         write (*, *) 'WARNING: You are using PIMD with direct TeraChem interface.'
         write (*, *) 'You should have "integrator regular" in &
& the TeraChem input file'
      end if
      write (*, *) 'Number of TeraChem servers = ', nteraservers
      do ipom = 1, nteraservers
         call connect_terachem(ipom)
      end do

      if (nproc /= nteraservers) then
         write (*, *) 'WARNING: parameter "nproc" must equal "nteraservers"'
         write (*, *) 'Setting nproc = ', nteraservers
         nproc = nteraservers
      end if
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

#else
   if (pot == '_tera_') then
      write (*, *) 'FATAL ERROR: This version was not compiled with MPI support.'
      write (*, *) 'You cannot use the direct MPI interface to TeraChem.'
      call abinerror('init')
   end if
#endif

   if (en_restraint >= 1) then
      call en_rest_init()

      if (en_restraint == 1) then
         write (*, *) 'Energy restraint is ON(1): Using method of Lagrange multipliers.'
      else if (en_restraint == 2 .and. en_kk >= 0) then
         write (*, *) 'Energy restraint is ON(2): Using quadratic potential restraint.'
      else
         write (*, *) 'FATAL ERROR: en_restraint must be either 0 or 1(Lagrange multipliers),2(umbrella, define en_kk)'
         call abinerror('init')
      end if
   end if

   if (mdtype /= '') then
      mdtype = UpperToLower(mdtype)
      select case (mdtype)
      case ('md')
         ipimd = 0
      case ('pimd')
         ipimd = 1
      case ('sh')
         ipimd = 2
      case ('minimization')
         ipimd = 3
      case ('ehrenfest')
         ipimd = 4
      case ('landau_zener')
         ipimd = 5
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

   if (my_rank == 0) then
      write (*, *) 'Reading MD parameters from input file ', trim(chinput)
      write (*, *) 'Reading xyz coordinates from file ', trim(chcoords)
      write (*, *) 'XYZ Units = '//trim(xyz_units)
      if (chveloc /= '') then
         write (*, *) 'Reading initial velocities [a.u.] from file ', trim(chveloc)
      end if
      print '(a)', chdivider
      call print_logo()
      print '(a)', chdivider
      call print_compile_info()
      print '(a)', chdivider
      call print_runtime_info()
      print '(a)', chdivider
      print '(a)', '                                              '
      select case (ipimd)
      case (0)
         print '(a)', '              Classical MD                    '
      case (1)
         print '(a)', '            Path Integral MD                  '
      case (2)
         print '(a)', '           Surface Hopping MD                 '
      case (3)
         print '(a)', '              Minimization                    '
      case (4)
         print '(a)', '              Ehrenfest MD                    '
      case (5)
         print '(a)', '             Landau Zener MD                  '
      end select

      write (*, *) '    using potential: ', LowerToUpper(pot)
      print '(a)', '                                              '
      print '(a)', chdivider
   end if

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
   if (iqmmm == 0 .and. pot /= 'mm') natqm = natom

   if (irest == 1) then
      readnhc = 1 !readnhc has precedence before initNHC
      readQT = 1
      initNHC = 0 !i.e. if(readnhc.eq.1.and.initNHC.eq.1)
      scaleveloc = 0 !do not scale velocities when restarting a job
   else !then nhc momenta from restart.xyz will be used
      readnhc = 0
      readQT = 0
      initNHC = 1
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

#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4
   allocate (ishake1(natom * 3 - 6))
   allocate (ishake2(natom * 3 - 6))
#endif
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

   ! allocate all basic arrays and set them to 0.0d0
   call allocate_arrays(natom, nwalk + 1)
   ! Ehrenfest require larger array since gradients for all of the states are need
   ! TODO: We should really make this differently..
   if (ipimd == 4) call allocate_ehrenfest(natom, nstate)

   if (iplumed == 1) then
      call plumed_init()
   end if

   ! READING GEOMETRY
   read (111, *)
   do iat = 1, natom

      read (111, *, iostat=iost) names(iat), x(iat, 1), y(iat, 1), z(iat, 1)
      if (iost /= 0) call print_read_error(chcoords, 'Could not read atom names and coordinates', iost)
      names(iat) = LowerToUpper(names(iat))
      if (UpperToLower(trim(xyz_units)) == "angstrom") then
         x(iat, 1) = x(iat, 1) * ANG
         y(iat, 1) = y(iat, 1) * ANG
         z(iat, 1) = z(iat, 1) * ANG
      else if (UpperToLower(trim(xyz_units)) == "bohr") then
         continue
      else
         write (*, *) 'ERROR: Wrong XYZ units: ', trim(xyz_units)
      end if

   end do
   close (111)

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
   ! check, whether we hit End-Of-File or other error
   if (IS_IOSTAT_END(iost)) then !fortran intrinsic for EOF
      continue
   else if (iost /= 0) then
      write (*, *) 'ERROR when reading namelist "system".'
      write (*, *) chiomsg
      call abinerror('init')
   end if

   do iat = 1, MAXTYPES
      massnames(iat) = LowerToUpper(massnames(iat))
   end do

   ! Determine atomic masses from periodic table
   call mass_init(masses, massnames)
   ! Transform masses for PIMD
   ! TODO: rename this function
   call init_mass(amg, amt)
   ! Lower the second character of atom name.
   ! This is because of TeraChem.
   do iat = 1, natom
      names(iat) (2:2) = UpperToLower(names(iat) (2:2))
   end do

#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4
   allocate (natmolt(natom))
   allocate (nshakemol(natom))
#endif
   natmolt = 0
   natmolt(1) = natom ! default for global NHC thermostat
   nshakemol = 0

   read (150, nhcopt)
   rewind (150)

   if (ipimd == 2 .or. ipimd == 4) then
      read (150, sh)
      rewind (150)
      integ = UpperToLower(integ)
   end if

   if (ipimd == 5) then
      read (150, lz)
      rewind (150)
      call lz_init() !Init arrays for possible restart
   end if

   if (iremd == 1) then
#ifdef USE_MPI
      read (150, remd)
      rewind (150)
      call remd_init(temp, temp0)
#else
      write (*, *) 'FATAL ERROR: This version was not compiled with MPI support.'
      write (*, *) 'You cannot do REMD.'
      call abinerror('init')
#endif
   end if

   if (iqmmm > 0 .or. pot == 'mm') then
      read (150, qmmm)
      rewind (150)
   end if

   close (150)
   ! END OF READING INPUT

!$ call OMP_set_num_threads(nproc)
!$ nthreads = omp_get_max_threads()

#ifdef USE_MPI
   if (pot == '_tera_' .or. restrain_pot == '_tera_') then
      call initialize_terachem()
      if (ipimd == 2 .or. ipimd == 4 .or. ipimd == 5) call init_terash(x, y, z)
   end if
#endif

   ! HERE WE CHECK FOR ERRORS IN INPUT
   call check_inputsanity()

   ! resetting number of walkers to 1 in case of classical simulation
   if (ipimd == 0) then
      if (my_rank == 0) then
         write (*, *) 'Using velocity Verlet integrator'
      end if
      md = 2
      nwalk = 1
      nabin = 1 ! TODO:safety for respashake code
      ! We should probably copy shake to velocity verlet
      ! algorithm as well
   end if

   ! for surface hopping and ehrenfest
   if (ipimd == 2 .or. ipimd == 4 .or. ipimd == 5) then
      nwalk = ntraj !currently 1
      md = 2 ! velocity verlet
      nabin = 1
   else if (ipimd == 1 .and. inormalmodes /= 1) then
      if (my_rank == 0) write (*, *) 'Using RESPA integrator.'
      md = 1
   else if (ipimd == 1 .and. inormalmodes == 1) then
      md = 2
   end if

   ! we should include shake into the verlet routine
   if (nshake /= 0) then
      md = 3
   end if

   if (pot_ref /= 'none') then
      md = 4
      write (*, '(A)') 'Using Multiple Time-Step RESPA integrator!'
      write (*, '(A)') "Reference (cheap) potential is "//trim(pot_ref)
      write (*, '(A, F6.2)') "with timestep [fs] ", dt / nstep_ref * AUtoFS
      write (*, '(A)') "Full potential is "//trim(pot)
      write (*, '(A, F6.2)') "with timestep [fs] ", dt * AUtoFS
   end if

   if (my_rank == 0) then
      if (temp0 > 0) then
         write (*, *) 'Initial temperature [K] =', temp0
      else
         write (*, *) 'Initial temperature [K] =', temp
      end if
      if (inose /= 0) write (*, *) 'Target temperature [K] =', temp
   end if

   ! conversion of temperature from K to au
   temp = temp / AUtoK
   temp0 = temp0 / AUtoK

   if (ihess == 1) then
      allocate (hess(natom * 3, natom * 3, nwalk))
      allocate (cvhess_cumul(nwalk))
      cvhess_cumul = 0.0D0
   end if

   ! SHAKE initialization, determining the constrained bond lenghts
   if (nshake /= 0) then
      if (my_rank == 0) then
         write (*, *) 'Setting distances for SHAKE from XYZ coordinates'
      end if
      call shake_init(x, y, z)
   end if

   if (pot == 'mmwater' .or. pot_ref == 'mmwater') call check_water(natom, names)

   ! MUST BE BEFORE RESTART DUE TO ARRAY ALOCATION
   if (my_rank /= 0) then
      call srand(irandom)
      do ipom = 0, my_rank
         irandom = irand()
      end do
   end if

   ! initialize prng, maybe rewritten during restart
   ! call vranf(rans,0,irandom)
   call gautrg(rans, 0, irandom)

   ! THERMOSTAT INITIALIZATION
   if (inose == 1) then
      call nhc_init()
   else if (inose == 2) then
      call gle_init(dt * 0.5 / nabin / nstep_ref) !nabin is set to 1 unless ipimd=1
   else if (inose == 3) then
      call pile_init(dt * 0.5, tau0_langevin)
   else if (inose == 0) then
      write (*, '(A)') 'No thermostat. NVE ensemble.'
   else
      write (*, '(A)') 'ERROR: Invalid "inose" value!'
      call abinerror('init')
   end if

   ! performing RESTART from restart.xyz
   if (irest == 1) call restin(x, y, z, vx, vy, vz, it)

   if (pot == '2dho') then
      f = 0 !temporary hack
   end if
   if (nchain > 1) then
      f = 0 ! what about nchain=1?
      ! what about massive thermostat?
   end if

   ! SETTING initial velocities according to the Maxwell-Boltzmann distribution
   if (irest == 0 .and. chveloc == '') then
      ! TODO: GLE thermostat, initialize momenta in gle_init
      if (temp0 >= 0) then
         call vinit(temp0, am, vx, vy, vz)
      else
         call vinit(temp, am, vx, vy, vz)
      end if
   end if

   ! Reading velocities from file
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
            atom = LowerToUpper(atom)
            if (atom /= LowerToUpper(names(iat))) then
               write (*, *) 'Offending line:'
               write (*, *) atom, vx(iat, iw), vy(iat, iw), vz(iat, iw)
               call print_read_error(chveloc, "Inconsistent atom types in input velocities.", iost)
            end if
         end do
      end do

      close (500)

   end if

   ! END OF READING VELOCITIES

   ! doing this here so that we can do it even when reading velocities from file
   if (rem_comvel) call remove_comvel(vx, vy, vz, am, rem_comvel)
   if (rem_comrot) call remove_rotations(x, y, z, vx, vy, vz, am, rem_comrot)

   if (conatom > 0) call constrainP(vx, vy, vz, conatom)

   ! If scaleveloc=1, scale initial velocitites to match the temperature
   ! Otherwise, just print the temperature.
   call ScaleVelocities(vx, vy, vz)

   ! Initialize spherical boundary onditions
   if (isbc == 1) call sbc_init(x, y, z)

   ! inames initialization for the MM part.
   ! We do this also because string comparison is very costly
   if (iqmmm == 3 .or. pot == 'mm') allocate (inames(natom))

   if (iqmmm == 3 .or. pot == 'mm') then
      do iat = 1, MAXTYPES
         if (attypes(iat) == '') exit
         attypes(iat) = LowerToUpper(attypes(iat))
      end do
      call inames_init()
      call ABr_init()
   end if

   if (my_rank == 0) then
      write (*, *)
      write (*, *) '--------------SIMULATION PARAMETERS--------------'
      write (*, nml=general)
      write (*, *)
      write (*, nml=system)
      write (*, *)
      if (inose >= 1) write (*, nml=nhcopt)
      write (*, *)
      if (ipimd == 2 .or. ipimd == 4) write (*, nml=sh)
      write (*, *)
      if (ipimd == 5) write (*, nml=lz)
      write (*, *)
      if (iqmmm == 3 .or. pot == 'mm') write (*, nml=qmmm)
      write (*, *)
   end if
   call flush (6)
#ifdef USE_MPI
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
   pid = GetPID()
   if (my_rank == 0) write (*, *) 'Pid of the current proccess is:', pid

   ! Open files for writing
   ! TODO: It's strange that we're passing these random params here...
   call files_init(isbc, phase, ndist, nang, ndist)

   call flush (6)

contains

   subroutine check_inputsanity()
      use mod_chars, only: chknow

      !  We should exclude all non-abinitio options, but whatever....
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
         ! write(*,*)'ERROR: Input velocities are not compatible with irest=1.'
         write (*, *) 'WARNING: Input velocities from file'//trim(chveloc)//' will be ignored!'
         write (*, *) 'Velocities will be taken from restart file because irest=1.'
         ! write(*,*)chknow
         ! if(iknow.ne.1) error=1
      end if

      ! Check, whether input variables don't exceeds array limits
      if (ntraj > ntrajmax) then
         write (*, *) 'Maximum number of trajectories is:'
         write (*, *) ntrajmax
         write (*, *) 'Adjust variable ntrajmax in modules.f90'
         error = 1
      end if
      if (nstate > nstmax) then
         write (*, *) 'Maximum number of states is:'
         write (*, *) nstmax
         write (*, *) 'Adjust variable nstmax in modules.f90'
         error = 1
      end if
      if (nchain > maxchain) then
         write (*, *) 'Maximum number of Nose-Hoover chains is:'
         write (*, *) maxchain
         write (*, *) 'Adjust variable maxchain in modules.f90'
         error = 1
      end if
      if (ndist >= ndistmax) then
         write (*, *) 'Maximum number of bonds for printing is:'
         write (*, *) ndistmax
         write (*, *) 'Adjust variable ndistmax in modules.f90'
         error = 1
      end if
      if (nang >= ndistmax) then
         write (*, *) 'Maximum number of angles (nang) for printing is:'
         write (*, *) ndistmax
         write (*, *) 'Adjust variable ndistmax in modules.f90'
         error = 1
      end if
      if (ndih >= ndistmax) then
         write (*, *) 'Maximum number of dihedral angles (ndih) for printing is:'
         write (*, *) ndistmax
         write (*, *) 'Adjust variable ndistmax in modules.f90'
         error = 1
      end if
      ! HERE we check for errors in input.
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
            write (*, *) 'ERROR: Nwalk is not divisible by nproc. This is not a wise usage of your computer time.'
            error = 1
         end if
      end if
      if (pot == 'none') then
         write (*, *) 'FATAL: Variable "pot" not specified.Exiting now...'
         error = 1
      end if
      if (ipimd == 1 .and. nwalk <= 1) then
         write (*, *) 'Number of walkers for PIMD (nwalk) <=1 !'
         write (*, *) 'Either set ipimd=0 for classical simulation or'
         write (*, *) 'set nwalk > 1'
         error = 1
      end if
      if (iqmmm < 0 .or. iqmmm > 3) then
         write (*, *) 'Error: iqmmm must be 0, 1, 2 or 3.'
         error = 1
      end if
      if (integ /= 'euler' .and. integ /= 'rk4' .and. integ /= 'butcher') then
         write (*, *) 'integ must be "euler","rk4" or "butcher".'
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

      if (imini < 0) then
         write (*, *) 'Input error: imini must be positiv or zero.'
         error = 1
      end if
      if (nstep < 0) then
         write (*, *) 'Input error: nstep must be positive.'
         error = 1
      end if
      if (nwrite <= 0) then
         write (*, *) 'Input error: nwrite must be positive.'
         error = 1
      end if
      if (nwritex <= 0) then
         write (*, *) 'Input error: nwritex must be positive.'
         error = 1
      end if
      if (nrest <= 0) then
         write (*, *) 'Input error: nrest must be positive.'
         error = 1
      end if
      if (nabin <= 0) then
         write (*, *) 'Input error: nabin must be positive.'
         error = 1
      end if
      if (icv /= 0 .and. icv /= 1) then
         write (*, *) 'Input error: icv must be 1 or zero.'
         error = 1
      end if
      if (temp < 0) then
         write (*, *) 'Input error: temp must be positive.'
         error = 1
      end if
      if (dt <= 0) then
         write (*, *) 'Time step negative or undefined!'
         write (*, *) 'Modify variable "dt" in input the general input section.'
         error = 1
      end if
      if (ncalc <= 0) then
         write (*, *) 'Ncalc must be positive integer number!'
         error = 1
      end if
      if (ncalc > nwrite) then
         write (*, *) 'Ncalc greater than nwrite.Setting nwrite=ncalc'
         nwrite = ncalc
      end if

      if (ipimd == 1 .and. inose <= 0) then
         write (*, *) 'You have to use thermostat with PIMD!(inose>=0)'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (ipimd < 0 .or. ipimd > 3) then
         if (ipimd /= 5) then
            write (*, *) 'ipimd has to be 0,1,2 or 3.'
            error = 1
         end if
      end if
      if (ipimd == 5 .and. pot == '_tera_' .and. ntriplet_lz > 0) then
         write (*, *) 'ERROR: Landau-Zener with Singlet-Triplet transitions not implemented with TeraChem over MPI.'
         error = 1
      end if
      if (ipimd == 5 .and. nwalk /= 1) then
         write (*, *) 'ERROR: LZ not implemented with multiple walkers.'
         error = 1
      end if
      if (istage /= 1 .and. istage /= 0) then
         write (*, *) 'ERROR: istage has to be 0 or 1'
         error = 1
      end if
      if (inormalmodes < 0 .and. inormalmodes > 2) then
         write (*, *) 'ERROR: inormalmodes has to be 0, 1 or 2!'
         error = 1
      end if
      if (readnhc == 1 .and. initNHC == 1 .and. irest == 1) then
         write (*, *) 'Warning: Conflicting keywords readnhc and initNHC set to 1.'
         write (*, *) 'Momenta from restart.xyz will be used.'
      end if
      if (readnhc == 1 .and. irest == 0) then
         write (*, *) 'Ignoring readnhc=1 since irest=0.'
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
      if (inose == 1 .and. (ipimd == 2 .or. ipimd == 4) .and. en_restraint == 0) then
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
      if (imasst /= 0 .and. imasst /= 1) then
         write (*, *) 'Input error: imasst must be 1 or zero.'
         error = 1
      end if
      if (imasst == 0 .and. ipimd == 1) then
         write (*, *) 'PIMD simulations must use massive thermostat ( imasst=1)! '
         error = 1
      end if
      if (imasst == 0 .and. nmolt <= 0) then
         write (*, *) 'Number of molecules coupled to separate NH chains not specified!Set nmolt > 0.'
         error = 1
      end if
      if (nmolt > natom) then
         write (*, *) 'Input error: nmolt > natom, which is not possible. Consult the manual.'
         error = 1
      end if
      if (imasst == 0) then
         do imol = 1, nmolt
            if (natmolt(imol) <= 0) then
               write (*, *) 'Number of atoms in molecules not specified! Set array natmolt properly.'
               error = 1
            end if
         end do
      end if
      if (inose < 0 .and. inose > 3) then
         write (*, *) 'inose has to be 0,1,2 or 3.'
         error = 1
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
         write (*, *) 'PIMD should be done with staging or normal mode transformation! Exiting...'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (istage == 1 .and. inose == 2) then
         write (*, *) 'The staging transformation is not compatible with GLE thermostat.'
         error = 1
      end if
      if (nyosh /= 1 .and. nyosh /= 3 .and. nyosh /= 7) then
         write (*, *) 'Variable nyosh(order of Suzuki-Yoshiga scheme) must be 1,3 or 7'
         error = 1
      end if
      if (nyosh <= 1 .and. inose == 1) then
         write (*, *) 'It is strongly reccommended to use Suzuki-Yoshida scheme'//&
                     &' when using Nose-Hoover thermostat (nyosh=3 or 7).'
         write (*, *) iknow, error, chknow
         if (iknow /= 1) error = 1
      end if
      if (nrespnose < 3 .and. inose == 1) then
         write (*, *) 'Variable nrespnose < 3! Assuming this is an error in input and exiting.'
         write (*, *) 'Such low value would probably not produce stable results.'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if (nrespnose <= 0) then
         write (*, *) 'Variable nrespnose must be positive integer'
         error = 1
      end if
      if (irest /= 1 .and. irest /= 0) then
         write (*, *) 'ERROR:irest has to be 1 or 0'
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

      if (pot == '2dho' .and. natom > 1) then
         write (*, *) 'Only 1 particle is allowed for 2D harmonic oscillator!'
         error = 1
      end if
      if (pot == 'mm' .and. iqmmm > 0) then
         write (*, *) 'Pot="mm"is not compatible with iqmmm>0!'
         error = 1
      end if
      if (iqmmm > 1) then
         write (*, *) 'WARNING: QMMM is higly experimental at this point. Use with care!'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if
      if ((natmm + natqm /= natom) .and. iqmmm > 0) then
         write (*, *) 'Natmm+natqm not equal to natom!'
         error = 1
      end if

      if (inose == 1 .and. imasst == 0) then
         ipom = 0
         do iat = 1, nmolt
            ipom = ipom + natmolt(iat)
         end do
         if (ipom /= natom) then
            write (*, *) 'Number of atoms in thermostated molecules(natmolt) doesnt match with natom.'
            write (*, *) 'This is probably mistake in input.Exiting...'
            write (*, *) chknow
            if (iknow /= 1) error = 1
         end if
      end if

      if (temp < 1 .and. inose >= 1) then
         write (*, *) 'WARNING!:Temperature below 1K. Are you sure?'
         write (*, *) 'This is probably mistake in input.Exiting...'
         write (*, *) chknow
         if (iknow /= 1) error = 1
      end if

      if (iremd == 1) then
         write (chout, '(A,I2.2)') 'movie.xyz.', my_rank
      else
         chout = 'movie.xyz'
      end if
      inquire (FILE=chout, EXIST=file_exists)
      if (file_exists) then
         if (irest == 0) then
            if (my_rank == 0) write (*, *) 'File '//trim(chout)//' exists. Please (re)move it or set irest=1.'
            error = 1
         else
            if (my_rank == 0) write (*, *) 'File "movie.xyz" exists and irest=1.Trajectory will be appended.'
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

   subroutine print_logo()
      print '(a)', '                  _____     _    _      _ '
      print '(a)', '        /\       |  _  \   | |  |  \   | |'
      print '(a)', '       /  \      | |_|  |  | |  | | \  | |'
      print '(a)', '      / /\ \     |     /   | |  | |\ \ | |'
      print '(a)', '     / /__\ \    |=====|   | |  | | \ \| |'
      print '(a)', '    / /____\ \   |  _   \  | |  | |  \ | |'
      print '(a)', '   / /      \ \  | |_|  |  | |  | |   \  |'
      print '(a)', '  /_/        \_\ |_____/   |_|  |_|    \_|'
      print '(a)', ' '

      print '(a)', ' version 1.1'
      print '(a)', ' D. Hollas, J. Suchan, O. Svoboda, M. Oncak, P. Slavicek'
      print '(a)', ' '
   end subroutine print_logo

   subroutine print_runtime_info()
      character(len=1024) :: cmdline
      print '(a)', ''
      print '(a)', '          RUNTIME INFO'
      print '(a)', ' '
      write (*, '(A17)') "Running on node: "
      call system('uname -n')
      write (*, '(A19)') 'Working directory: '
      call system('pwd')
      write (*, *)
      call get_command(cmdline)
      write (*, *) trim(cmdline)
      call flush (6)
      call get_command_argument(0, cmdline)
      write (*, *)
      call system('ldd '//cmdline)
      print '(a)', ' '
   end subroutine print_runtime_info

end subroutine init

! TODO: Move to a separate file.
subroutine finish(error_code)
   use mod_arrays, only: deallocate_arrays
   use mod_general
   use mod_files, only: MAXUNITS
   use mod_nhc, only: inose, finalize_nhc
   use mod_gle, only: finalize_gle
   use mod_estimators, only: h
   use mod_harmon, only: hess
   use mod_lz, only: lz_finalize
   use mod_transform, only: finalize_normalmodes
   use mod_cp2k, only: finalize_cp2k

   use mod_plumed, only: iplumed, finalize_plumed

#ifdef USE_MPI
   use mod_terampi, only: finalize_terachem
   use mod_terampi_sh, only: finalize_terash
   use mpi, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_Finalize, MPI_Abort
#endif
   implicit none
   integer, intent(in) :: error_code
   integer :: i, ierr
   logical :: lopen

#ifdef USE_MPI
   if (pot == '_tera_') then
      if (ipimd == 2) then
         call finalize_terash()
      end if
      call finalize_terachem(error_code)
   end if
#endif

   if (my_rank == 0) then
      write (*, *) ''
      if (error_code == 0) then
         write (*, '(A)') 'Job finished!'
      end if
   end if

   call deallocate_arrays()

   ! TODO: Move this to a subroutine in mod_files
   do i = 2, MAXUNITS
      inquire (unit=i, opened=lopen)
      ! TODO: This is not portable, do not hardcode 5 and 6!
      if (lopen .and. i /= 5 .and. i /= 6) then
         close (i)
      end if
   end do

   if (allocated(hess)) then
      deallocate (hess)
   end if
   if (allocated(h)) then
      deallocate (h)
   end if

   if (inormalmodes > 0) then
      call finalize_normalmodes()
   end if

   if (inose == 1) then
      call finalize_nhc()
   end if
   if (inose > 1 .and. inose < 5) then
      call finalize_gle()
   end if

   if (iplumed == 1) then
      call finalize_plumed()
   end if

   ! Cleanup Landau-Zener
   if (ipimd == 5) then
      call lz_finalize()
   end if

   ! MPI_FINALIZE is called in this routine as well
   if (pot == '_cp2k_') then
      call finalize_cp2k()
   end if

   ! At last, we terminate MPI processes
   ! TODO: We should have an MPI module handling this
   ! TODO: We should check whether MPI was initialized with MPI_Init
   ! before we attempt to call MPI_Finalize().
#ifdef USE_MPI
   if (iremd == 1 .or. pot == '_tera_' .or. pot == '_cp2k_') then
      if (error_code == 0 .and. pot /= "_cp2k_") then
         call MPI_Finalize(ierr)
         if (ierr /= MPI_SUCCESS) then
            write (*, '(A)') 'Bad signal from MPI_FINALIZE: ', ierr
            ! Let's try to continue
         end if
      else if (error_code > 0) then
         call MPI_Abort(MPI_COMM_WORLD, error_code, ierr)
      end if
   end if
#endif

end subroutine finish
