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
   use mod_sbc,      only: sbc_init, rb_sbc, kb_sbc, isbc, rho
   use mod_random
   use mod_splined_grid, only: initialize_spline
   use mod_utils, only: lowertoupper, uppertolower, file_exists_or_exit
   use mod_vinit
   use mod_analyze_geometry
   use mod_shake
   use mod_minimize, only: gamm, gammthr
   use mod_analysis, only: restin
   use mod_water,    only: watpot, check_water
   use mod_plumed,   only: iplumed, plumedfile, plumed_init
   use mod_en_restraint
   use mod_transform, only:init_mass
   use mod_cp2k
   use mod_remd
   use mod_terampi
   use mod_terampi_sh
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   real(DP),intent(out) :: dt
   real(DP) :: masses(MAXTYPES)
   real(DP)  :: rans(10)
   integer :: iw, iat, natom_xyz, imol, shiftdihed = 1, iost
   integer :: error, getpid, nproc=1, ipom, i
   character(len=2)    :: massnames(MAXTYPES), atom
   character(len=200)  :: chinput, chcoords, chveloc
   character(len=200)  :: chiomsg, chout
   character(len=20)   :: xyz_units='angstrom'
   character(len=60)   :: chdivider
   character(len=60)   :: mdtype
   LOGICAL :: file_exists
   logical :: rem_comvel, rem_comrot
   integer :: ierr
   integer :: irand

   namelist /general/ natom, pot, ipimd, mdtype, istage, inormalmodes, nwalk, nstep, icv, ihess, imini, nproc, iqmmm, &
            nwrite,nwritex,nwritev, nwritef, dt,irandom,nabin,irest,nrest,anal_ext,  &
            isbc,rb_sbc,kb_sbc,gamm,gammthr,conatom,mpi_sleep,narchive,xyz_units, &
            dime,ncalc,idebug, enmini, rho, iknow, watpot, iremd, iplumed, plumedfile, &
            en_restraint, en_diff, en_kk, restrain_pot, &
            pot_ref, nstep_ref, nteraservers, max_wait_time, cp2k_mpi_beads

#ifdef USE_MPI
   namelist /remd/   nswap, nreplica, deltaT, Tmax, temp_list
#endif

   namelist /nhcopt/ inose,temp,temp0,nchain,ams,tau0, tau0_langevin, imasst,nrespnose,nyosh,      &
                     scaleveloc,readNHC,readQT,initNHC,nmolt,natmolt,nshakemol,rem_comrot,rem_comvel

   namelist /system/ masses,massnames,ndist,dist1,dist2, &
                     nang,ang1,ang2,ang3, ndih,dih1,dih2,dih3,dih4,shiftdihed, &
                     k,r0,k1,k2,k3,De,a,D0_dw,lambda_dw,k_dw, r0_dw, &
                     Nshake,ishake1,ishake2,shake_tol

   namelist /sh/     istate_init,nstate,substep,deltae,integ,inac,nohop,phase,decoh_alpha,popthr,ignore_state, &
                     nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, adjmom, revmom, &
                     dE_S0S1_thr, correct_decoherence
             
   namelist /lz/     initstate_lz, nstate_lz, nsinglet_lz, ntriplet_lz, deltaE_lz, energydifthr_lz 

   namelist /qmmm/   natqm,natmm,q,rmin,eps,attypes


   chcoords = 'mini.dat'
   chinput = 'input.in'
   chveloc = ''
   mdtype = ''
   dt = -1
   error = 0
   iplumed = 0

   chdivider = "######################################################"

   call get_cmdline(chinput, chcoords, chveloc, tc_server_name)

   ! READING MAIN INPUT
   open(150,file=chinput, status='OLD', delim='APOSTROPHE', action = "READ")
   read(150, general)
   rewind(150)
   pot = UpperToLower(pot)

   if(pot.eq.'splined_grid')then
      natom = 1
      dime = 1
      f = 0
      call initialize_spline()
   end if

   if(pot.eq."_cp2k_".or.pot_ref.eq."_cp2k_")then
      call init_cp2k()
#ifdef USE_MPI
   else
      if (pot == "_tera_" .and. nteraservers > 1) then
         ! We will be calling TS servers concurently
         ! via OpenMP parallelization, hence we need MPI_Init_thread().
         ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm
         call MPI_Init_thread(MPI_THREAD_MULTIPLE, i, ierr)
         if (i /= MPI_THREAD_MULTIPLE) then
            write (*, *) 'Provided safety level is not MPI_THREAD_MULTIPLE'
            write (*, '(A,I1,A,I1)') 'Requested ', MPI_THREAD_MULTIPLE, 'got:', i
            call abinerror('init')
         end if
         ! nproc is used to initialize OpenMP threads below.
         if (nproc /= nteraservers) then
            nproc = nteraservers
         end if
      else
         ! TODO: Check whether MPI is already initialized
         ! This can happen when using internal CP2K interface.
         call MPI_Init(ierr)
      end if
      if (ierr.ne.0)then
         write(*,*)'Bad signal from MPI_INIT:', ierr
         stop 1
      end if
#endif
   end if

   ! Set OpenMP parallelization
   ! Currently only used in PIMD for trivial
   ! parallelization over PI beads.
   ! Note that scaling is actually not so great
   ! since SCF timings will vary for different beads,
   ! which decreases thread utilization.
!$ call OMP_set_num_threads(nproc)

! We need to connect to TeraChem as soon as possible,
! because we want to shut down TeraChem nicely in case something goes wrong.
#ifdef USE_MPI
   ! TODO: Move this to an mpi_wrapper module
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)
   ! TODO: allow mpi_world_size > 1 only for REMD
   if (my_rank.eq.0.and.mpi_world_size.gt.1)then
      write(*,'(A,I3)')'Number of MPI processes = ', mpi_world_size
   end if
   if(pot.eq.'_tera_'.or.restrain_pot.eq.'_tera_')then
      call initialize_terachem_interface()
   end if

   ! TODO: Do we need a barrier here?
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

#else
   if(pot.eq.'_tera_')then
      write(*,*)'FATAL ERROR: This version was not compiled with MPI support.'
      write(*,*)'You cannot use the direct MPI interface to TeraChem.'
      call abinerror('init')
   end if
#endif


   if (en_restraint.ge.1) then
      call en_rest_init()
 
      if (en_restraint.eq.1)then
         write(*,*) 'Energy restraint is ON(1): Using method of Lagrange multipliers.'
      else if (en_restraint.eq.2.and.en_kk.ge.0)then
         write(*,*) 'Energy restraint is ON(2): Using quadratic potential restraint.'
      else
         write(*,*) 'FATAL ERROR: en_restraint must be either 0 or 1(Lagrange multipliers),2(umbrella, define en_kk)'
         call abinerror('init')
      end if
   end if

   if (mdtype.ne.'')then
      mdtype = UpperToLower(mdtype)
      SELECT CASE (mdtype)
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
      END SELECT
   end if

   if(iremd.eq.1)then
      write(chcoords,'(A,I2.2)')trim(chcoords)//'.',my_rank
      if(chveloc.ne.'')then
         write(chveloc,'(A,I2.2)')trim(chveloc)//'.',my_rank
      end if
   end if

   call file_exists_or_exit(chcoords)
   if (chveloc.ne.'')then
      call file_exists_or_exit(chveloc)
   end if
   
   if (my_rank.eq.0)then
      write(*,*)'Reading MD parameters from input file ', trim(chinput)
      write(*,*)'Reading xyz coordinates from file ',trim(chcoords)
      write(*,*)'XYZ Units = '//trim(xyz_units)
      if(chveloc.ne.'')then
         write(*,*)'Reading initial velocities [a.u.] from file ', trim(chveloc)
      end if
      print '(a)', chdivider
      call print_logo()
      print '(a)', chdivider
      call print_compile_info()
      print '(a)', chdivider
      call print_runtime_info()
      print '(a)', chdivider
      print '(a)','                                              '
      SELECT CASE (ipimd)
         case (0)
            print '(a)','              Classical MD                    '
         case (1)
            print '(a)','            Path Integral MD                  '
         case (2)
            print '(a)','           Surface Hopping MD                 '
         case (3)
            print '(a)','              Minimization                    '
         case (4)
            print '(a)','              Ehrenfest MD                    '
         case (5)
            print '(a)','             Landau Zener MD                  '
      END SELECT

            write(*,*)'    using potential: ', LowerToUpper(pot)
               print '(a)','                                              '
            print '(a)', chdivider
   end if

   ! Get number of atoms from XYZ coordinates NOW so that we can allocate arrays


   open(111, file = chcoords, status = "old", action = "read")
   read(111, '(I50)', iostat = iost)natom_xyz
   !TODO following line does not work
   if (iost.ne.0) call print_read_error(chcoords,"Expected number of atoms on the first line.", iost)
   if(natom_xyz.ne.natom.and.natom.gt.0)then
     write(*,'(A,A)')'WARNING: Number of atoms specified in ', trim(chinput)
     write(*,'(A,A)')'does not match with the XYZ geometry in ', trim(chcoords)
     write(*,*)'Going forward anyway...'
   endif

   natom = natom_xyz

   if (natom.lt.1)then 
      write(*,'(A,A)')'ERROR: Wrong number of atoms on the first line of the XYZ &
      & file ',trim(chcoords)
      write(*,*)natom
      call abinerror('init')
   end if

   ! This line is super important,
   ! cause we actually use natqm in many parts of the code
   if(iqmmm.eq.0.and.pot.ne.'mm') natqm = natom

   if(irest.eq.1)then
    readnhc=1   !readnhc has precedence before initNHC
    readQT=1
    initNHC=0   !i.e. if(readnhc.eq.1.and.initNHC.eq.1)
    scaleveloc=0  !do not scale velocities when restarting a job
   else         !then nhc momenta from restart.xyz will be used
    readnhc=0
    readQT=0
    initNHC=1
    scaleveloc=1
   endif
   
   if (chveloc.ne.'')then
      scaleveloc=0
   end if

! for future multiple time step integration in SH
   dt0 = dt

   ! We have to initialize here, because we read them from input
   allocate(names(natom))
   names     = ''
   attypes   = ''
   massnames = ''
   masses    = -1.0d0

#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
   allocate( ishake1(natom*3-6) )
   allocate( ishake2(natom*3-6) )
#endif
   ishake1 = 0
   ishake2 = 0

   ! By default, remove COM translation and rotation
   ! Unless restarting or taking velocities from a file
   if(irest.eq.1.or.trim(chveloc).ne.'')then
      rem_comrot=.false.
      rem_comvel=.false.
   else
      rem_comrot=.true.
      rem_comvel=.true.
   end if


!  allocate all basic arrays and set them to 0.0d0
   call allocate_arrays(natom, nwalk+1)
!  Ehrenfest require larger array since gradients for all of the states are need   
!  TODO: We should really make this differently..
   if(ipimd.eq.4) call allocate_ehrenfest(natom, nstate)

   if(iplumed.eq.1) then
      call plumed_init()
   endif


!  READING GEOMETRY
   read(111, *)
   do iat = 1, natom

     read(111,*, iostat=iost)names(iat),x(iat,1),y(iat,1),z(iat,1)
     if(iost.ne.0) call print_read_error(chcoords,'Could not read atom names and coordinates', iost)
     names(iat) = LowerToUpper(names(iat))
     if (UpperToLower(trim(xyz_units)).eq."angstrom")then
         x(iat,1) = x(iat,1) * ANG
         y(iat,1) = y(iat,1) * ANG
         z(iat,1) = z(iat,1) * ANG
     else if (UpperToLower(trim(xyz_units)).eq."bohr")then
         continue
     else
         write(*,*)'ERROR: Wrong XYZ units: ', trim(xyz_units)
     end if

   enddo 
   close(111)

   do iw=1,nwalk
      do iat=1,natom
         x(iat,iw) = x(iat,1)
         y(iat,iw) = y(iat,1)
         z(iat,iw) = z(iat,1)
      enddo
   enddo


      ! the namelist system does not need to be present
      read(150,system,iostat=iost,iomsg=chiomsg)
      rewind(150)
      !check, whether we hit End-Of-File or other error
      if(IS_IOSTAT_END(iost))then  !fortran intrinsic for EOF
         continue
      else if (iost.ne.0)then
         write(*,*)'ERROR when reading namelist "system".'
         write(*,*)chiomsg
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
      do iat=1,natom
         names(iat)(2:2)=UpperToLower(names(iat)(2:2))
      end do

#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
      allocate ( natmolt(natom)  )
      allocate ( nshakemol(natom)  )
#endif
      natmolt   = 0
      natmolt(1) = natom ! default for global NHC thermostat
      nshakemol = 0

      read(150,nhcopt)
      rewind(150)

      if(ipimd.eq.2.or.ipimd.eq.4)then
         read(150, sh)
         rewind(150)
         integ = UpperToLower(integ)
      end if

      if(ipimd.eq.5)then
         read(150, lz)
         rewind(150)
         call lz_init() !Init arrays for possible restart
      end if

   if(iremd.eq.1)then
#ifdef USE_MPI
      read(150, remd)
      rewind(150)
      call remd_init(temp, temp0)
#else
      write(*,*)'FATAL ERROR: This version was not compiled with MPI support.'
      write(*,*)'You cannot do REMD.'
      call abinerror('init')
#endif
      end if

      if(iqmmm.gt.0.or.pot.eq.'mm')then
         read(150, qmmm)
         rewind(150)
      end if

      close(150)
!--END OF READING INPUT---------------


#ifdef USE_MPI
   if(pot.eq.'_tera_'.or.restrain_pot.eq.'_tera_')then
      call initialize_tc_servers()
      if (ipimd.eq.2.or.ipimd.eq.4.or.ipimd.eq.5)then
         call init_terash(x, y, z)
      end if
   end if
#endif

!-----HERE WE CHECK FOR ERRORS IN INPUT-----------
      call check_inputsanity()

!     resetting number of walkers to 1 in case of classical simulation      
      if(ipimd.eq.0)then
         if(my_rank.eq.0)then
            write(*,*)'Using velocity Verlet integrator'
         end if
         md = 2
         nwalk = 1
         nabin = 1 ! TODO:safety for respashake code
                   ! We should probably copy shake to velocity verlet
                   ! algorithm as well
      endif

!for surface hopping and ehrenfest     
      if(ipimd.eq.2.or.ipimd.eq.4.or.ipimd.eq.5)then
         nwalk = ntraj !currently 1
         md = 2 ! velocity verlet
         nabin = 1
      else if(ipimd.eq.1.and.inormalmodes.ne.1)then
         if (my_rank.eq.0) write(*,*)'Using RESPA integrator.'
         md = 1
      else if(ipimd.eq.1.and.inormalmodes.eq.1)then
         md = 2
      end if

      ! we should include shake into the verlet routine
      if(nshake.ne.0)then
         md = 3
      end if

      if(pot_ref.ne.'_none_')then
         md = 4
         write(*, '(A)')'Using Multiple Time-Step RESPA integrator!'
         write(*, '(A)')"Reference (cheap) potential is "//trim(pot_ref)
         write(*, '(A, F6.2)')"with timestep [fs] ", dt / nstep_ref * AUtoFS
         write(*, '(A)')"Full potential is "//trim(pot)
         write(*, '(A, F6.2)')"with timestep [fs] ", dt * AUtoFS
      end if

      if (my_rank.eq.0)then
         if (temp0.gt.0)then
            write(*,*)'Initial temperature [K] =', temp0
         else
            write(*,*)'Initial temperature [K] =', temp
         end if
         if (inose.ne.0) write(*,*)'Target temperature [K] =', temp
      end if

      ! conversion of temperature from K to au
      temp = temp / AUtoK
      temp0 = temp0 / AUtoK

      if (ihess.eq.1)then 
         allocate ( hess(natom*3,natom*3,nwalk) )
         allocate ( cvhess_cumul(nwalk) )
         cvhess_cumul=0.0d0
      endif

!     SHAKE initialization, determining the constrained bond lenghts
      if(nshake.ne.0)then
         if (my_rank.eq.0)then 
            write(*,*)'Setting distances for SHAKE from XYZ coordinates'
         end if
         call shake_init(x,y,z)
      endif

      if (pot.eq.'mmwater'.or.pot_ref.eq.'mmwater') call check_water(natom, names)

!    MUST BE BEFORE RESTART DUE TO ARRAY ALOCATION
     if (my_rank .ne. 0) then
        call srand(irandom)
        do ipom = 0, my_rank
           irandom = irand()
        end do
     end if

!    call vranf(rans,0,IRandom)  !initialize prng,maybe rewritten during restart
     call gautrg(rans, 0, IRandom)  !initialize prng, maybe rewritten during restart

!    THERMOSTAT INITIALIZATION
     if (inose.eq.1)then
        call nhc_init()
     else if (inose.eq.2)then
        call gle_init(dt*0.5/nabin/nstep_ref) !nabin is set to 1 unless ipimd=1
     else if (inose.eq.3)then
        call pile_init(dt * 0.5, tau0_langevin)
     else if (inose.eq.0)then
        write(*, '(A)')'No thermostat. NVE ensemble.'
     else
        write(*,'(A)')'ERROR: Invalid "inose" value!'
        call abinerror('init')
     end if


!    performing RESTART from restart.xyz
     if(irest.eq.1) call restin(x, y, z, vx, vy, vz, it)


      if(pot.eq.'2dho')then
         f=0 !temporary hack
      endif
      if(nchain.gt.1)then
         f=0 !what about nchain=1?
         ! what about massive therm?
      endif

!     SETTING initial velocities according to the Maxwell-Boltzmann distribution
      if(irest.eq.0.and.chveloc.eq.'')then
         ! TODO: GLE thermostat, initialize momenta in gle_init
         if (temp0.ge.0)then
            call vinit(temp0, am, vx, vy, vz)
         else
            call vinit(temp, am, vx, vy, vz)
         end if
      end if

!     Reading velocities from file
      if (chveloc.ne.''.and.irest.eq.0)then
         ! TODO: move the following to a separate function
         open(500,file=chveloc, status='OLD', action = "READ")
         do iw=1,nwalk
            read(500,*, IOSTAT=iost)natom_xyz
            if (iost.ne.0) call print_read_error(chveloc,"Could not read velocities on line 1.", iost)
            if(natom_xyz.ne.natom)then
               write(*,'(A,A)')'Nunmber of atoms in velocity input ', trim(chveloc)
               write(*,'(A,A)')'does not match with XYZ coordinates in ', trim(chcoords)
               call abinerror('init')
            endif
            read(500,*, IOSTAT=iost)
            if (iost.ne.0) call print_read_error(chveloc,"Could not read velocities on line 2.", iost)
          
            do iat=1,natom
               read(500,*, IOSTAT=iost)atom, vx(iat,iw), vy(iat,iw), vz(iat, iw)
               if (iost.ne.0) call print_read_error(chveloc,"Could not read velocities.", iost)
               atom = LowerToUpper(atom)
               if (atom.ne.LowerToUpper(names(iat)))then
                  write(*,*)'Offending line:'
                  write(*,*)atom, vx(iat,iw), vy(iat,iw), vz(iat, iw)
                  call print_read_error(chveloc,"Inconsistent atom types in input velocities.", iost)
               end if
            end do
         end do

         close(500)

      end if

!     END OF READING VELOCITIES--------------------

      ! doing this here so that we can do it even when reading velocities from file
      if(rem_comvel) call remove_comvel(vx, vy, vz, am, rem_comvel)
      if(rem_comrot) call remove_rotations(x, y, z, vx, vy, vz, am, rem_comrot)

      if(conatom.gt.0) call constrainP(vx, vy, vz, conatom)

      ! If scaleveloc=1, scale initial velocitites to match the temperature
      ! Otherwise, just print the temperature.
      call ScaleVelocities(vx, vy, vz)

      ! Initialize spherical boundary onditions
      if(isbc.eq.1) call sbc_init(x,y,z)

      ! inames initialization for the MM part. 
      ! We do this also because string comparison is very costly
      if(iqmmm.eq.3.or.pot.eq.'mm') allocate( inames(natom) )

      if(iqmmm.eq.3.or.pot.eq.'mm')then 
         do iat = 1, MAXTYPES
            if(attypes(iat).eq.'') exit
            attypes(iat) = LowerToUpper(attypes(iat))
         end do
         call inames_init()
         call ABr_init()
      endif

   if(my_rank.eq.0)then
      write(*,*)
      write(*,*)'--------------SIMULATION PARAMETERS--------------'
      write(*,nml=general)
      write(*,*)
      write(*,nml=system)
      write(*,*)
      if ( inose.ge.1 ) write(*,nml=nhcopt)
      write(*,*)
      if ( ipimd.eq.2.or.ipimd.eq.4 ) write(*,nml=sh)
      write(*,*)
      if ( ipimd.eq.5 ) write(*,nml=lz)
      write(*,*)
      if (iqmmm.eq.3.or.pot.eq.'mm') write(*,nml=qmmm)
      write(*,*)
   end if
   call flush(6)
#ifdef USE_MPI
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
   pid = GetPID()
   ! TODO: Print pid together with my_rank, need to be part of a single write statement
   write (*, '(A,I0)') 'Pid of the current proccess is: ', pid


   ! Open files for writing
   ! TODO: It's strange that we're passing these random params here...
   call files_init(isbc, phase, ndist, nang, ndist)

   call flush(6)

   CONTAINS

   subroutine check_inputsanity()
   use mod_chars, only: chknow
!$ integer :: nthreads, omp_get_max_threads

   !  We should exclude all non-abinitio options, but whatever....
!$    nthreads = omp_get_max_threads()
!$    if(nthreads.gt.1.and.(ipimd.ne.1.and.pot.ne.'_cp2k_'))then
!$     write(*,*)'Number of threads is ', nthreads
!$     write(*,*)'ERROR: Parallel execution is currently only supported with ab initio PIMD (ipimd=1)'
!$     call abinerror('init')
!$    endif

      if(nproc.gt.1)then
!$       if(.false.)then
         write(*,*)'FATAL ERROR: This executable was not compiled with parallel support.'
         error=1
!$       end if
      end if

      if(irest.eq.1.and.chveloc.ne.'')then
      !   write(*,*)'ERROR: Input velocities are not compatible with irest=1.'
         write(*,*)'WARNING: Input velocities from file'//trim(chveloc) //' will be ignored!'
         write(*,*)'Velocities will be taken from restart file because irest=1.'
      !   write(*,*)chknow
      !   if(iknow.ne.1) error=1
      end if

      !-----Check, whether input variables don't exceeds array limits
      if(ntraj.gt.ntrajmax)then
         write(*,*)'Maximum number of trajectories is:'
         write(*,*)ntrajmax
         write(*,*)'Adjust variable ntrajmax in modules.f90'
         error=1
      endif
      if(nstate.gt.nstmax)then
         write(*,*)'Maximum number of states is:'
         write(*,*)nstmax
         write(*,*)'Adjust variable nstmax in modules.f90'
         error=1
      endif
      if(nchain.gt.maxchain)then
       write(*,*)'Maximum number of Nose-Hoover chains is:'
       write(*,*)maxchain
       write(*,*)'Adjust variable maxchain in modules.f90'
       error=1
      endif
      if(ndist.ge.ndistmax)then
       write(*,*)'Maximum number of bonds for printing is:'
       write(*,*)ndistmax
       write(*,*)'Adjust variable ndistmax in modules.f90'
       error=1
      endif
      if(nang.ge.ndistmax)then
       write(*,*)'Maximum number of angles (nang) for printing is:'
       write(*,*)ndistmax
       write(*,*)'Adjust variable ndistmax in modules.f90'
       error=1
      endif
      if(ndih.ge.ndistmax)then
       write(*,*)'Maximum number of dihedral angles (ndih) for printing is:'
       write(*,*)ndistmax
       write(*,*)'Adjust variable ndistmax in modules.f90'
       error=1
      endif
!----------HERE we check for errors in input.      
      if (pot.ne.'_cp2k_')then
         if (nproc.gt.nwalk)then
            write(*,*)'ERROR: Nproc greater than nwalk. That does not make sense.'
            write(*,*)'Set nproc <= nwalk.'
            error=1
         end if
         if (nproc.le.0)then
            write(*,*)'ERROR: Nproc must be a positive integer.'
            error=1
         end if
         if (modulo(nwalk,nproc).ne.0)then
            write(*,*)'ERROR: Nwalk is not divisible by the number of OpenMP threads.'
            write(*,*)'This is not a wise usage of your computer time.'
            error=1
         end if
      end if
      if(pot.eq.'_none_')then
       write(*,*)'FATAL: Variable "pot" not specified.Exiting now...'
       error=1
      endif
      if(ipimd.eq.1.and.nwalk.le.1)then
       write(*,*)'Number of walkers for PIMD (nwalk) <=1 !'
       write(*,*)'Either set ipimd=0 for classical simulation or'
       write(*,*)'set nwalk > 1'
       error=1
      endif
      if(iqmmm.lt.0.or.iqmmm.gt.3)then
       write(*,*)'Error: iqmmm must be 0, 1, 2 or 3.'
       error=1
      endif
      if(integ.ne.'euler'.and.integ.ne.'rk4'.and.integ.ne.'butcher')then
       write(*,*)'integ must be "euler","rk4" or "butcher".'
       error=1
      endif
      if(integ.ne.'butcher')then
         write(*,*)'WARNING: variable integ is not "butcher", which is the default and most accurate.'
         write(*,*)chknow
         if(iknow.ne.1) error=1
      end if
      if(deltae.lt.0)then
       write(*,*)'Parameter deltae must be non-negative number.'
       error=1
      endif
      if(popsumthr.lt.0)then
       write(*,*)'Parameter popsumthr must be positive number.'
       error=1
      endif
      if(energydifthr.lt.0)then
       write(*,*)'Parameter energydifthr must be positive number in eV units.'
       error=1
      endif
      if(energydriftthr.lt.0)then
       write(*,*)'Parameter energydriftthr must be positive number in eV units.'
       error=1
      endif
      if(shiftdihed.ne.0.and.shiftdihed.ne.1)then
       write(*,*)'Shiftdihed must be either 0 (for dihedrals -180:180) or 1 (for dihedrals 0:360)'
       error=1
      endif

      if(shiftdihed.eq.0) shiftdih=0.0d0
      if(shiftdihed.eq.1) shiftdih=360.0d0

      if(imini.lt.0)then
              write(*,*)'Input error: imini must be positiv or zero.'
              error=1
      endif
      if(nstep.lt.0)then
              write(*,*)'Input error: nstep must be positive.'
              error=1
      endif
      if(nwrite.le.0)then
              write(*,*)'Input error: nwrite must be positive.'
              error=1
      endif
      if(nwritex.le.0)then
              write(*,*)'Input error: nwritex must be positive.'
              error=1
      endif
      if(nrest.le.0)then
              write(*,*)'Input error: nrest must be positive.'
              error=1
      endif
      if(nabin.le.0)then
              write(*,*)'Input error: nabin must be positive.'
              error=1
      endif
      if(icv.ne.0.and.icv.ne.1)then
              write(*,*)'Input error: icv must be 1 or zero.'
              error=1
      endif
      if(temp.lt.0)then
              write(*,*)'Input error: temp must be positive.'
              error=1
      endif
      if(dt.le.0)then
              write(*,*)'Time step negative or undefined!'
              write(*,*)'Modify variable "dt" in input the general input section.'
              error=1
      endif
      if(ncalc.le.0)then
              write(*,*)'Ncalc must be positive integer number!'
              error=1
      endif
      if(ncalc.gt.nwrite)then
              write(*,*)'Ncalc greater than nwrite.Setting nwrite=ncalc'
              nwrite=ncalc
      endif

      if(ipimd.eq.1.and.inose.le.0)then
       write(*,*)'You have to use thermostat with PIMD!(inose>=0)'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif
      if(ipimd.lt.0.or.ipimd.gt.3)then
       if(ipimd.ne.5)then
          write(*,*)'ipimd has to be 0,1,2 or 3.'
          error=1
       endif
      endif
      if(ipimd.eq.5.and.pot.eq.'_tera_'.and.ntriplet_lz.gt.0)then
         write(*,*)'ERROR: Landau-Zener with Singlet-Triplet transitions not implemented with TeraChem over MPI.'
         error=1
      endif
      if(ipimd.eq.5.and.nwalk.ne.1)then
         write(*,*)'ERROR: LZ not implemented with multiple walkers.'
         error=1
      endif
      if(istage.ne.1.and.istage.ne.0)then
         write(*,*)'ERROR: istage has to be 0 or 1'
         error=1 
      endif
      if(inormalmodes.lt.0.and.inormalmodes.gt.2)then
         write(*,*)'ERROR: inormalmodes has to be 0, 1 or 2!'
         error=1 
      endif
      if(readnhc.eq.1.and.initNHC.eq.1.and.irest.eq.1)then
       write(*,*)'Warning: Conflicting keywords readnhc and initNHC set to 1.'
       write(*,*)'Momenta from restart.xyz will be used.'
      endif
      if(readnhc.eq.1.and.irest.eq.0)then
       write(*,*)'Ignoring readnhc=1 since irest=0.'
      endif
      if(inac.gt.2.or.inac.lt.0)then
       write(*,*)'Parameter "inac" must be 0,1 or 2.'   !be very carefull if you change this!
       error=1
      endif
      if(adjmom.gt.1.or.adjmom.lt.0)then
       write(*,*)'Parameter "adjmom" must be 0 or 1.' 
       error=1
      endif
      if(adjmom.eq.0.and.inac.eq.1)then
       write(*,*)'Combination of adjmom=0 and inac=1 is not possible.' 
       write(*,*)'We dont have NAC vector if inac=1.' 
       error=1
      endif
      if(irest.eq.1.and.scaleveloc.eq.1)then
       write(*,*)'irest=1 AND scaleveloc=1.'
       write(*,*)'You are trying to scale the velocities read from restart.xyz.'
       write(*,*)'I assume this is an error in input. (set scaleveloc=0)'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif
      if(inose.eq.1.and.(ipimd.eq.2.or.ipimd.eq.4).and.en_restraint.eq.0)then
       write(*,*)'Thermostating is not meaningful for surface hopping simulation.'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif
      if(inose.eq.1.and.(ipimd.eq.5).and.en_restraint.eq.0)then
       write(*,*)'Thermostating is not meaningful for Landau-Zener MD.'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif
      if(istate_init.gt.nstate)then
       write(*,*)'Error:Initial state > number of computed states. Exiting...'
       error=1
      endif
      if(ipimd.eq.5)then
          if(initstate_lz.gt.nstate_lz)then
             write(*,*) initstate_lz, nstate_lz
             write(*,*)'Error(LZ):Initial state > number of computed states. Exiting...'
             error=1
          endif
          if(nstate_lz.le.0)then
             write(*,*)'Error(LZ):No states to compute (nstate_lz<=0). Exiting...'
             error=1
          endif
          if(nsinglet_lz.eq.0.and.ntriplet_lz.eq.0.and.nstate_lz.gt.0)then
             nsinglet_lz = nstate_lz                !Assume singlet states
          endif
          if((nsinglet_lz+ntriplet_lz).ne.nstate_lz)then
             write(*,*)'Error(LZ): Sum of singlet and triplet states must give total number of states. Exiting...'
             error=1
          endif
      endif
      if(energydifthr_lz.lt.0)then
       write(*,*)'Parameter energydifthr_lz must be positive number in eV units.'
       error=1
      endif
      if(nac_accu1.le.0.or.nac_accu2.lt.0)then
       write(*,*)'Input error:NACME precision must be a positive integer.'
       write(*,*)'The treshold is then 10^-(nac_accu).'
       error=1
      endif
      if(nac_accu1.le.nac_accu2)then
       write(*,*)'nac_accu1 < nac_accu2'
       write(*,*)'I will compute NACME only with default accuracy:',nac_accu1
      endif
      if(imasst.ne.0.and.imasst.ne.1)then
              write(*,*)'Input error: imasst must be 1 or zero.'
              error=1
      endif
      if(imasst.eq.0.and.ipimd.eq.1)then
              write(*,*)'PIMD simulations must use massive thermostat ( imasst=1)! '
              error=1
      endif
      if(imasst.eq.0.and.nmolt.le.0)then
              write(*,*)'Number of molecules coupled to separate NH chains not specified!Set nmolt > 0.'
              error=1
      endif
      if(nmolt.gt.natom)then
              write(*,*)'Input error: nmolt > natom, which is not possible. Consult the manual.'
              error=1
      endif
      if(imasst.eq.0)then
       do imol=1,nmolt
        if(natmolt(imol).le.0)then
         write(*,*)'Number of atoms in molecules not specified! Set array natmolt properly.'
         error=1
        endif
       enddo
      endif
      if(inose.lt.0.and.inose.gt.3)then
       write(*,*)'inose has to be 0,1,2 or 3.'
       error=1
      endif
      if(istage.eq.1.and.ipimd.ne.1)then
         write(*,*)'The staging transformation is only meaningful for PIMD'
         error=1
      endif
      if(inormalmodes.gt.0.and.ipimd.ne.1)then
         write(*,*)'The normal mode transformation is only meaningful for PIMD. Exiting...'
         error=1
      endif
      if(istage.eq.0.and.ipimd.eq.1.and.inose.ne.2.and.inormalmodes.eq.0)then
       write(*,*)'PIMD should be done with staging or normal mode transformation! Exiting...'
       write(*,*)chknow
       if (iknow.ne.1) error=1
      endif
      if(istage.eq.1.and.inose.eq.2)then
         write(*,*)'The staging transformation is not compatible with GLE thermostat.'
         error=1
      endif
      if(nyosh.ne.1.and.nyosh.ne.3.and.nyosh.ne.7)then
       write(*,*)'Variable nyosh(order of Suzuki-Yoshiga scheme) must be 1,3 or 7'
       error=1
      endif
      if(nyosh.le.1.and.inose.eq.1)then
       write(*,*)'It is strongly reccommended to use Suzuki-Yoshida scheme when using Nose-Hoover thermostat (nyosh=3 or 7).'
       write(*,*)iknow, error, chknow
       if (iknow.ne.1) error=1
      endif
      if(nrespnose.lt.3.and.inose.eq.1)then
       write(*,*)'Variable nrespnose < 3! Assuming this is an error in input and exiting.'
       write(*,*)'Such low value would probably not produce stable results.'
       write(*,*)chknow
       if (iknow.ne.1)error=1
      endif
      if(nrespnose.le.0)then
       write(*,*)'Variable nrespnose must be positive integer'
       error=1
      endif
      if(irest.ne.1.and.irest.ne.0)then
      write(*,*)'ERROR:irest has to be 1 or 0'
       error=1
      endif
      if(nshake.ne.0.and.ipimd.eq.1)then
       write(*,*)'PIMD with SHAKE cannot use massive thermostating!'
       error=1
      endif
      if(nshake.ne.0.and.imasst.eq.1.and.inose.gt.0)then
       write(*,*)'SHAKE cannot use massive thermostating!'
       write(*,*)'Set imasst=1 and nmolt, natmolt and nshakemol accordingly.'
       error=1
      endif

      if(pot.eq.'2dho'.and.natom.gt.1)then
       write(*,*)'Only 1 particle is allowed for 2D harmonic oscillator!'
       error=1
      endif
      if(pot.eq.'mm'.and.iqmmm.gt.0)then
       write(*,*)'Pot="mm"is not compatible with iqmmm>0!'
       error=1
      endif
      if(iqmmm.gt.1)then
       write(*,*)'WARNING: QMMM is higly experimental at this point. Use with care!'
       write(*,*)chknow
       if (iknow.ne.1) error=1
      endif
      if((natmm+natqm.ne.natom).and.iqmmm.gt.0)then
       write(*,*)'Natmm+natqm not equal to natom!'
       error=1
      endif

      if(inose.eq.1.and.imasst.eq.0)then
       ipom=0
       do iat=1,nmolt
        ipom=ipom+natmolt(iat)
       enddo
       if(ipom.ne.natom)then
        write(*,*)'Number of atoms in thermostated molecules(natmolt) doesnt match with natom.'
        write(*,*)'This is probably mistake in input.Exiting...'
       write(*,*)chknow
        if(iknow.ne.1) error=1
       endif
      endif

      if(temp.lt.1.and.inose.ge.1)then
       write(*,*)'WARNING!:Temperature below 1K. Are you sure?'
       write(*,*)'This is probably mistake in input.Exiting...'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif

      if(iremd.eq.1)then
         write(chout, '(A,I2.2)')'movie.xyz.',my_rank
      else
         chout='movie.xyz'
      end if
      INQUIRE(FILE=chout, EXIST=file_exists)
      if(file_exists)then
       if(irest.eq.0)then
        if (my_rank.eq.0) write(*,*)'File '//trim(chout)//' exists. Please (re)move it or set irest=1.'
        error=1
       else
         if (my_rank.eq.0) write(*,*)'File "movie.xyz" exists and irest=1.Trajectory will be appended.'
       endif
      endif

      if(iremd.eq.1)then
         write(chout, '(A,I2.2)')'restart.xyz.',my_rank
      else
         chout='restart.xyz'
      end if
      INQUIRE(FILE=chout, EXIST=file_exists)
      if(file_exists)then
       if(irest.eq.0)then
        write(*,*)'File ',trim(chout),' exists. Please (re)move it or set irest=1.'
        error=1
       endif
      else
       if(irest.eq.1)then
         write(*,*)'File ', trim(chout), ' not found.' 
         error=1
       endif 
      endif

      if(error.eq.1)then
         write(*,*)'Input errors were found! Exiting now...'
         call abinerror('check_inputsanity')
      endif
   end subroutine check_inputsanity


   subroutine print_read_error(chfile, chmsg, iost)
      character(len=*), intent(in)  :: chmsg, chfile
      integer, intent(in)  :: iost
      write(*,*) trim(chmsg)
      write(*,'(A,A)')'Error when reading file ', trim(chfile)
      write(*,*)'Error code was', iost
      call abinerror('init')
   end subroutine print_read_error


   subroutine print_logo()
print '(a)','                  _____     _    _      _ '
print '(a)','        /\       |  _  \   | |  |  \   | |'
print '(a)','       /  \      | |_|  |  | |  | | \  | |'
print '(a)','      / /\ \     |     /   | |  | |\ \ | |'
print '(a)','     / /__\ \    |=====|   | |  | | \ \| |'
print '(a)','    / /____\ \   |  _   \  | |  | |  \ | |'
print '(a)','   / /      \ \  | |_|  |  | |  | |   \  |'
print '(a)','  /_/        \_\ |_____/   |_|  |_|    \_|'
print '(a)',' '

! TODO: Pass version as compiler parameter
print '(a)',' version 1.1'
print '(a)',' D. Hollas, J. Suchan, O. Svoboda, M. Oncak, P. Slavicek'
print '(a)',' '
   end subroutine print_logo


   subroutine print_runtime_info()
   character(len=1024) :: cmdline
   print '(a)',''
   print '(a)','          RUNTIME INFO'
   print '(a)',' '
   write(*,'(A17)')"Running on node: "
   call system('uname -n')
   write(*,'(A19)')'Working directory: '
   call system('pwd')
   write(*,*)
   call get_command(cmdline)
   write(*,*)trim(cmdline)
   call flush(6)
   call get_command_argument(0, cmdline)
   write(*,*)
   call system('ldd ' // cmdline)
   print '(a)',' '
   end  subroutine print_runtime_info

end subroutine init


subroutine finish(error_code)
   use mod_arrays, only: deallocate_arrays
   use mod_general
   use mod_files,  only: MAXUNITS
   use mod_nhc!,   only: finalize_nhc
   use mod_gle,    only: finalize_gle
   use mod_estimators, only: h
   use mod_harmon, only: hess
   use mod_lz,     only: lz_finalize
   use mod_transform, only: finalize_normalmodes
   use mod_cp2k,   only: finalize_cp2k

   use mod_plumed, only: iplumed, finalize_plumed

#ifdef USE_MPI
   use mod_terampi, only: finalize_terachem
   use mod_terampi_sh, only: finalize_terash
   use mpi, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_Finalize, MPI_Abort
#endif
   implicit none
   integer, intent(in) :: error_code
   integer  :: i, ierr
   logical  :: lopen

#ifdef USE_MPI
   if (pot.eq.'_tera_')then
      if (ipimd.eq.2) then
         call finalize_terash()
      end if
      call finalize_terachem(error_code)
   end if
#endif

   if (my_rank.eq.0)then
      write(*,*)''
      if (error_code.eq.0)then
         write(*,'(A)')'Job finished!'
      end if
   end if


   call deallocate_arrays( )

   ! TODO: Move this to a subroutine in mod_files
   do i=2,MAXUNITS
      inquire(unit=i,opened=lopen)
      ! TODO: This is not portable, do not hardcode 5 and 6!
      if (lopen.and.i.ne.5.and.i.ne.6)then
         close(i)
      end if
   end do

   if (allocated(hess))then
      deallocate ( hess )
   end if
   if (allocated(h))then
      deallocate ( h )
   end if

   if (inormalmodes.gt.0)then
      call finalize_normalmodes()
   end if

   if(inose.eq.1)then
      call finalize_nhc()
   end if
   if(inose.gt.1.and.inose.lt.5)then
      call finalize_gle()
   end if

   if (iplumed.eq.1) then
      call finalize_plumed()
   end if

   ! Cleanup Landau-Zener
   if(ipimd.eq.5)then
      call lz_finalize()
   end if

   ! MPI_FINALIZE is called in this routine as well
   if(pot.eq.'_cp2k_')then
      call finalize_cp2k()
   end if

! At last, we terminate MPI processes
! TODO: We should have an MPI module handling this
! TODO: We should check whether MPI was initialized with MPI_Init
! before we attempt to call MPI_Finalize().
#ifdef USE_MPI
if(iremd.eq.1.or.pot.eq.'_tera_'.or.pot.eq.'_cp2k_')then
   if (error_code.eq.0.and.pot.ne."_cp2k_".or.pot.eq.'_tera_')then
      call MPI_Finalize(ierr)
      if (ierr.ne.MPI_SUCCESS)then
         write(*,'(A)')'Bad signal from MPI_FINALIZE: ', ierr
         ! Let's try to continue
      end if
   else if (error_code.gt.0)then
      call MPI_Abort(MPI_COMM_WORLD, error_code, ierr)
   end if
end if
#endif

end subroutine finish
