! This messy function performs many things, among others:
! 1. Reading input
! 2. Input sanity check
! 3. Allocation of arrays
! 4. Reading restart OR reading input geometry
! 5. Initialize velocities
! 6. Initialize everything else
! At this moment, coordinate and velocity transformations are NOT performed here
! Surface hopping is NOT initialized here.

subroutine init(dt, values1)
   use mod_const
   use mod_files
   use mod_arrays
   use mod_array_size
   use mod_general
   use mod_system
   use mod_nhc
   use mod_estimators
   use mod_harmon
   use mod_sh
   use mod_qmmm
   use mod_gle
   use mod_nab
   use mod_sbc,      only: sbc_init, rb_sbc, kb_sbc, isbc, rho
   use mod_random
   use mod_guillot,  only: inames_guillot
   use mod_utils
   use mod_vinit
   use mod_density
   use mod_shake
   use mod_minimize, only: gamm, gammthr
   use mod_analysis, only: restin
   use mod_water,    only: watpot, check_water
   use mod_plumed,   only: iplumed, plumedfile, plumed_init
   use mod_transform, only:init_mass
#ifdef USEFFTW
   use mod_fftw3,    only: fftw_init
#endif
! #ifdef CP2K
   use mod_cp2k
! #endif
#ifdef MPI
   use mod_remd
#endif
   use mod_terampi
   use mod_terampi_sh
   implicit none
   real(DP),intent(out) :: dt
   integer,dimension(8) :: values1
   real(DP) :: masses(MAXTYPES)
   real(DP)  :: rans(10)
   integer :: iw,iat,natom1,imol,shiftdihed=1, iost
   integer :: error, getpid, nproc=1, ipom, ipom2=0
   character(len=2)    :: massnames(MAXTYPES), atom
   character(len=200)  :: chinput, chcoords, chveloc
   character(len=200)  :: chiomsg, chout
   LOGICAL :: file_exists
!  real(DP) :: wnw=5.0d-5
   integer :: ierr
!$ integer :: nthreads, omp_get_max_threads
!  wnw "optimal" frequency for langevin (inose=3) 
   REAL, POINTER, DIMENSION(:) :: VECPTR => NULL ()  !null pointer
!  REAL, POINTER, DIMENSION(:,:) :: VECPTR2 => NULL ()
   REAL, POINTER  :: REALPTR => NULL ()
   logical        :: rem_comvel, rem_comrot

   namelist /general/natom, pot, ipimd, istage, inormalmodes, nwalk, nstep, icv, ihess,imini,nproc,iqmmm, &
            nwrite,nwritex,nwritev, nwritef, dt,irandom,nabin,irest,nrest,anal_ext,  &
            isbc,rb_sbc,kb_sbc,gamm,gammthr,conatom,mpi_sleep,narchive, &
            parrespa,dime,ncalc,idebug, enmini, rho, iknow, watpot, iremd, iplumed, plumedfile, &
            pot_ref, nstep_ref, teraport, nteraservers, cp2k_mpi_beads

#ifdef MPI
   namelist /remd/   nswap, nreplica, deltaT, Tmax, temps
#endif

   namelist /nhcopt/ inose,temp,temp0,nchain,ams,tau0,imasst,nrespnose,nyosh,      &
                     scaleveloc,readNHC,readQT,initNHC,nmolt,natmolt,nshakemol,rem_comrot,rem_comvel

   namelist /system/ masses,massnames,nbin,nbin_ang,ndist,dist1,dist2,xmin,xmax,disterror, &
                     nang,ang1,ang2,ang3,ndih,dih1,dih2,dih3,dih4,shiftdihed, &
                     k,r0,k1,k2,k3,De,a,D0_dw,lambda_dw,k_dw, r0_dw, &
                     Nshake,ishake1,ishake2,shake_tol

   namelist /sh/     istate_init,nstate,substep,deltae,integ,inac,nohop,phase,alpha,popthr,ignore_state, &
                     nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, adjmom, revmom

   namelist /qmmm/   natqm,natmm,q,rmin,eps,attypes

   namelist /nab/    ipbc,alpha_pme,kappa_pme,cutoff,nsnb,ips,epsinf,natmol,nmol


   chcoords='mini.dat'
   chinput='input.in'
   chveloc=''
   dt=-1  
   error=0
   iplumed=0


   call Get_cmdline(chinput, chcoords, chveloc)

   !-READING INPUT----------------------------------------- 

   open(150,file=chinput, status='OLD', delim='APOSTROPHE', action = "READ") !here ifort has some troubles
   read(150,general)
   rewind(150)
   pot=UpperToLower(pot)


   if(pot.eq."_cp2k_".or.pot_ref.eq."_cp2k_")then
#ifdef CP2K
      call init_cp2k()
#else
      write(*,*)'FATAL ERROR: ABIN was not compiled with CP2K interface.'
      write(*,*)''
      call abinerror('init')
#endif
#ifdef MPI
   else
      call MPI_INIT ( ierr )
      if (ierr.ne.0)then
         write(*,*)'Bad signal from MPI_INIT:', ierr
         call abinerror('MPI_INIT')
      end if
#endif
   end if

   if(iqmmm.eq.0.and.pot.ne.'mm')then
      natqm=natom
   endif

! We need to connect to TeraChem as soon as possible,
! because we want to shut down TeraChem nicely in case something goes wrong.
#ifdef MPI
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, mpi_world_size, ierr)
   if(pot.eq.'_tera_')then
      if (nwalk.gt.1)then
         write(*,*)'WARNING: You are using PIMD with direct TeraChem interface.'
         write(*,*)'You should have "integrator regular" in &
& the TeraChem input file'
      end if
      write(*,*)'Number of TeraChem servers = ', nteraservers
      do ipom=1, nteraservers
         call connect_terachem(ipom)
      end do
      
      if(nproc.ne.nteraservers)then
         write(*,*)'WARNING: parameter "nproc" must equal "nteraservers"'
         write(*,*)'Setting nproc = ', nteraservers
         nproc = nteraservers
      end if
   end if
#else
   if(pot.eq.'_tera_')then
      write(*,*)'FATAL ERROR: This version was not compiled with MPI support.'
      write(*,*)'You cannot use the direct MPI interface to TeraChem.'
      call abinerror('init')
   end if
#endif

   if(iplumed.eq.1) then

#ifdef PLUM
      call plumed_init(dt)
      write(*,*) 'PLUMED is on'
      write(*,*) 'PLUMEDfile is ',plumedfile   
#else
      write(*,*)'FATAL ERROR: ABIN was not compiled with PLUMED.'
      stop 1
#endif

   endif
   
   if (my_rank.eq.0)then
      write(*,*)'Reading parameters from input file ',chinput
      write(*,*)'Reading xyz coordinates from file ',chcoords
      call PrintLogo(values1)

print '(a)','**********************************************'
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
   END SELECT

  write(*,*)'    using potential: ', LowerToUpper(pot)
print '(a)','                                              '
print '(a)','**********************************************'

   end if

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
   allocate( names( natom )     )
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
   if(irest.eq.1)then
      rem_comrot=.false.
      rem_comvel=.false.
   else
      rem_comrot=.true.
      rem_comvel=.true.
   end if


!  allocate all basic arrays and set them to 0.0d0
   call allocate_arrays( natom, nwalk+1 )

!-----READING GEOMETRY
   if(iremd.eq.1) write(chcoords,'(A,I2.2)')trim(chcoords)//'.',my_rank

   open(111,file=chcoords,status = "old", action = "read", iostat=iost) 
   read(111,*, iostat=iost)natom1
   !TODO following line does not work
   if (iost.ne.0) call err_read(chcoords,"Expected number of atoms on the first line.", iost)
   if(natom1.ne.natom)then
     write(*,'(A,A)')'No. of atoms in ',chinput
     write(*,'(A,A)')'and in ',chcoords
     write(*,*)'do not match.'
     call abinerror('init')
   endif
   read(111,*)
   do iat=1,natom

     read(111,*, iostat=iost)names(iat),x(iat,1),y(iat,1),z(iat,1)
     if(iost.ne.0) call err_read(chcoords,'Could not read atom names and coordinates', iost)
     names(iat) = LowerToUpper(names(iat))
     x(iat,1) = x(iat,1) * ANG
     y(iat,1) = y(iat,1) * ANG
     z(iat,1) = z(iat,1) * ANG

   enddo 
   close(111)

   do iw=1,nwalk
      do iat=1,natom
         x(iat,iw) = x(iat,1)
         y(iat,iw) = y(iat,1)
         z(iat,iw) = z(iat,1)
      enddo
   enddo
!-----END OF READING GEOMETRY      


      ! the namelist system does not need to be present
      read(150,system,iostat=iost,iomsg=chiomsg)
      rewind(150)
      !check, whether we hit End-Of-File or other error
      if(IS_IOSTAT_END(iost))then  !fortran intrinsic for EOF
         if (my_rank.eq.0) write(*,*)'Namelist "system" not found. Ignoring...'
      else if (iost.ne.0)then
         write(*,*)'ERROR when reading namelist "system".'
         write(*,*)chiomsg
         call abinerror('init')
      else
         do iat=1,MAXTYPE
            massnames(iat)=LowerToUpper(massnames(iat))
         end do
      end if

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
      nshakemol = 0

      read(150,nhcopt)
      rewind(150)

      if(ipimd.eq.2)then
         read(150, sh)
         rewind(150)
         integ = UpperToLower(integ)
      end if

   if(iremd.eq.1)then
#ifdef MPI
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

      if(iqmmm.eq.2.or.pot.eq.'nab'.or.pot_ref.eq.'nab')then
#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
         allocate ( natmol(natom) )
#endif
         natmol = 0
         read(150, nab)
         rewind(150)
      end if

      close(150)
!--END OF READING INPUT---------------

!$ call OMP_set_num_threads(nproc)
!$ nthreads = omp_get_max_threads()

#ifdef MPI
   if(pot.eq.'_tera_')then
      call initialize_terachem()
      if (ipimd.eq.2) call init_terash(x, y, z)
   end if
#endif

!-----HERE WE CHECK FOR ERRORS IN INPUT-----------
      call check_inputsanity()

!     resetting number of walkers to 1 in case of classical simulation      
      if(ipimd.eq.0)then
         if(my_rank.eq.0)then
            write(*,*)'ipimd=0,Resetting number of walkers to 1.'
            write(*,*)'ipimd=0,Using velocity Verlet integrator'
         end if
         md = 2
         nwalk = 1
         nabin = 1   ! TODO:safety for respashake code
                   ! We should probably copy shake to velocity verlet
                   ! algorithm as well
      endif

!for surface hopping      
      if(ipimd.eq.2)then
         nwalk = ntraj
         md = 2 ! velocity verlet
         nabin = 1
      else if(ipimd.eq.1.and.inormalmodes.ne.1)then
         if (my_rank.eq.0) write(*,*)'Using RESPA integrator.'
         md = 1
      else if(ipimd.eq.1.and.inormalmodes.eq.1)then
         md = 2
      end if

      ! we should include shake into the verlet routine
      if(nshake.gt.0)then
         md = 3
      end if

      if(pot_ref.ne.'none')then
         md = 4
      end if

#ifdef USEFFTW
      if(inormalmodes.gt.0) call fftw_init(nwalk)
#endif


      if (my_rank.eq.0)then
         if (temp0.gt.0)then
            write(*,*)'Initial temperature in Kelvins =', temp0
         else
            write(*,*)'Initial temperature in Kelvins =', temp
         end if
         if (inose.ne.0) write(*,*)'Target temperature in Kelvins =', temp
      end if

!-----conversion of temperature from K to au
      temp = temp / AUtoK
      temp0 = temp0 / AUtoK

      if (ihess.eq.1)then 
         allocate ( hess(natom*3,natom*3,nwalk) )
         allocate ( cvhess_cumul(nwalk) )
         cvhess_cumul=0.0d0
         if (pot.eq.'nab'.or.pot_ref.eq.'nab')then
!!$OMP PARALLEL
            allocate ( h(natom*natom*9) )
!!$OMP END PARALLEL
         endif
      endif

!----SHAKE initialization,determining the constrained bond lenghts
      if(nshake.ge.1)then
         if (my_rank.eq.0) write(*,*)'Setting distances for SHAKE from mini.dat'
         call shake_init(x,y,z)
      endif

      call dist_init() !zeroing distribution arrays

      if (pot.eq.'mmwater'.or.pot_ref.eq.'mmwater') call check_water(natom, names)

!----THERMOSTAT INITIALIZATION------------------ 
!----MUST BE BEFORE RESTART DUE TO ARRAY ALOCATION
     if (my_rank .ne. 0) then
        call srand(irandom)
        do ipom=0,my_rank
        irandom = irand()
        end do
     end if
!    call vranf(rans,0,IRandom,6)  !initialize prng,maybe rewritten during restart
     call gautrg(rans,0,IRandom,6)  !initialize prng,maybe rewritten during restart
     if (inose.eq.1) call nhc_init()
     if (inose.eq.2) call gle_init(dt*0.5/nabin/nstep_ref) !nabin is set to 1 unless ipimd=1
     if (inose.eq.3) call pile_init(dt*0.5,tau0)


!----performing RESTART from restart.xyz
     if(irest.eq.1)then
        call restin(x,y,z,vx,vy,vz,it)
     end if
!----END OF INPUT SECTION---------------------



!-----INITIALIZATION-----------------------


      if(pot.eq.'2dho')then
         f=0 !temporary hack
      endif
      if(nchain.gt.1)then
         f=0 !what about nchain=1?
         ! what about massive therm?
      endif

!-----SETTING initial velocities according to the Maxwell-Boltzmann distribution
      if(irest.eq.0.and.chveloc.eq.'')then
         if (temp0.ge.0)then
            call vinit(TEMP0, am, vx, vy, vz, rem_comvel, rem_comrot)
         else
            call vinit(TEMP, am, vx, vy, vz, rem_comvel, rem_comrot)
         end if
      end if

!     Reading velocities from file
      if (chveloc.ne.'')then
         if(iremd.eq.1) write(chveloc,'(A,I2.2)')trim(chveloc)//'.',my_rank
         write(*,*)'Reading initial velocities in a.u. from external file:'
         write(*,*)chveloc 
         open(500,file=chveloc, status='OLD', action = "READ")
         do iw=1,nwalk
            read(500,*, IOSTAT=iost)natom1
            if (iost.ne.0) call err_read(chveloc,"Could not read velocities on line 1.", iost)
            if(natom1.ne.natom)then
               write(*,'(A,A)')'No. of atoms in ',chinput
               write(*,'(A,A)')'and in ',chcoords
               write(*,*)'do not match.'
               call abinerror('init')
            endif
            read(500,*, IOSTAT=iost)
            if (iost.ne.0) call err_read(chveloc,"Could not read velocities on line 2.", iost)
          
            do iat=1,natom
               read(500,*, IOSTAT=iost)atom, vx(iat,iw), vy(iat,iw), vz(iat, iw)
               if (iost.ne.0) call err_read(chveloc,"Could not read velocities.", iost)
               if (atom.ne.names(iat)) call err_read(chveloc,"Inconsistent atom types in input velocities.", iost)
            end do
         end do

         close(500)

      end if

!     END OF READING VELOCITIES--------------------

      ! doing this here so that we can do it even when reading velocities from file
      if(rem_comvel) call remove_comvel(vx, vy, vz, am, rem_comvel)
      if(rem_comrot) call remove_rotations(x, y, z, vx, vy, vz, am, rem_comrot)

      if(conatom.gt.0) call constrainP(vx,vy,vz)

      ! If scaleveloc=1, scale initial velocitites to match the temperature
      ! Otherwise, just print the temperature.
      call ScaleVelocities(vx, vy, vz)

       
!-----some stuff for spherical boundary onditions
      if(isbc.eq.1) call sbc_init(x,y,z)

!-----inames initialization for guillot rm MM part. 
!-----We do this also because string comparison is very costly
      if(iqmmm.eq.3.or.pot.eq.'mm'.or.pot=='guillot') allocate( inames(natom) )
      if(pot.eq.'guillot') call inames_guillot()

      if(iqmmm.eq.3.or.pot.eq.'mm')then 
         do iat=1, natom
            attypes(iat)=LowerToUpper(attypes(iat))
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
      if ( ipimd.eq.2 ) write(*,nml=sh)
      write(*,*)
      if (iqmmm.eq.3.or.pot.eq.'mm') write(*,nml=qmmm)
      write(*,*)
   end if
   call flush(6)
#ifdef MPI
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
   pid=GetPID()
   if(my_rank.eq.0) write(*,*)'Pid of the current proccess is:',pid


#ifdef NAB
      if (pot.eq."nab".or.pot_ref.eq."nab".or.iqmmm.eq.2)then
         if (alpha_pme.lt.0) alpha_pme = pi / cutoff
         if (kappa_pme.lt.0) kappa_pme = alpha_pme
         call nab_init(alpha_pme, cutoff, nsnb, ipbc, ips, iqmmm) !C function...see nabinit.c

         if(ipbc.eq.1)then
            allocate ( charges(natom) )
         if (nchain.eq.1) f=3
            call nab_getbox(boxx,boxy,boxz) !see nabinit.c
            call nab_getcharges(charges) !see nabinit.c

            do iat=1,natom
               charges(iat)=charges(iat)*ambtoau*sqrt(167100.75d0)  !for macsimus units
            enddo
            write(*,*)'Box sizes[Ang]: ',boxx,boxy,boxz
            write(*,*)'Half of box size: ',boxx2,boxy2,boxz2
!           write(*,*)x(iat,1)/ang,y(iat,1)/ang,z(iat,1)/ang,charges(iat)
            call ewald(VECPTR,VECPTR,charges,REALPTR,boxx,boxy,boxz,cutoff,alpha_pme,kappa_pme,epsinf,natom,ipom2)
            boxx = boxx * ANG
            boxy = boxy * ANG
            boxz = boxz * ANG
            boxx2 = 0.5d0 * boxx
            boxy2 = 0.5d0 * boxy
            boxz2 = 0.5d0 * boxz
         endif
      endif
#endif

!  Open files for writing
   call files_init(isbc, phase)

!--------END OF INITIALIZATION-------------------
   call flush(6)

   CONTAINS

   subroutine check_inputsanity()
   use mod_chars, only: chknow

   !  We should exclude all non-abinitio options, but whatever....
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

#ifndef NAB
      if(pot.eq.'nab')then
         write(*,*)'FATAL ERROR: The program was not compiled with NAB libraries.'
         call abinerror('init')
      end if
#endif

#ifndef USEFFTW
      if(inormalmodes.gt.0)then
         write(*,*)'FATAL ERROR: The program was not compiled with FFTW libraries.'
         write(*,*)'Normal mode transformations cannot be performed.'
         call abinerror('init')
      end if
#endif
      if(irest.eq.1.and.chveloc.ne.'')then
         write(*,*)'ERROR: Input velocities are not compatible with irest=1.'
         write(*,*)chknow
         if(iknow.ne.1) error=1
      end if

      !-----Check,whether input variables don't exceeds array limits
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
       write(*,*)'Maximum number of bonds for binning is:'
       write(*,*)ndistmax
       write(*,*)'Adjust variable ndistmax in modules.f90'
       error=1
      endif
      if(nbin.gt.nbinmax)then
       write(*,*)'Maximum number of bins for densities is:'
       write(*,*)nbinmax
       write(*,*)'Adjust variable nbinmax in modules.f90'
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
            write(*,*)'ERROR: Nwalk is not divisible by nproc. This is not a wise usage of your computer time.'
            error=1
         end if
      end if
      if(pot.eq.'none')then
       write(*,*)'FATAL: Variable "pot" not specified.Exiting now...'
       error=1
      endif
      if(ipimd.eq.1.and.nwalk.le.1)then
       write(*,*)'Number of walkers for PIMD (nwalk) <=1 !'
       write(*,*)'Either set ipimd=0 for classical simulation or'
       write(*,*)'set nwalk > 1'
       error=1
      endif
      if(ipbc.eq.1.and.nmol.le.1)then
       write(*,*)'You have to specify number of molecules(nmol=x) for PBC!'
       write(*,*)'Also dont forget to specify number of atoms in molecules(array natmol)'
       write(*,*)'These are used to wrap molecules back to the box'
       error=1
      endif
      if(iqmmm.lt.0.or.iqmmm.gt.3)then
       write(*,*)'Error: iqmmm must be 0, 1, 2 or 3.'
       error=1
      endif
      if(iqmmm.eq.2.and.ipbc.eq.1)then
       write(*,*)'QM/MM with PBC not supported !'
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
       write(*,*)'ipimd has to be 0,1,2 or 3.'
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
      if(inose.eq.1.and.ipimd.eq.2)then
       write(*,*)'Thermostating is not meaningful for surface hopping simulation.'
       write(*,*)chknow
       if(iknow.ne.1) error=1
      endif
      if(istate_init.gt.nstate)then
       write(*,*)'Error:Initial state > number of computed states. Exiting...'
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
      if(nshake.gt.0.and.ipimd.eq.1)then
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
       write(*,*)'PIMD with SHAKE cannot use massive thermostating!Exiting... !'
       error=1
      endif
      if(nshake.gt.0.and.imasst.eq.1)then
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
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

      if(isbc.eq.1.and.ipbc.eq.1)then
       write(*,*)'Spherical boundary conditions not compatible with periodic boundary conditions!'
       error=1
      endif

      if(ipbc.eq.1)then
       ipom=0
       do iat=1,nmol
        ipom=ipom+natmol(iat)
       enddo
       if(ipom.ne.natom)then
        write(*,*)'Number of atoms in molecules(natmol) doesnt match with natom.'
        error=1
       endif
      endif

      if(inose.eq.1.and.imasst.eq.0)then
       ipom=0
       do iat=1,nmolt
        ipom=ipom+natmolt(iat)
       enddo
       if(ipom.ne.natom)then
        write(*,*)'Number of atoms in thermostated molecules(natmol) doesnt match with natom.'
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
         call abinerror('init')
      endif



   end subroutine check_inputsanity

   subroutine err_read(chfile, chmsg, iost)
      character(len=*), intent(in)  :: chmsg, chfile
      integer, intent(in)  :: iost
      write(*,*) trim(chmsg)
      write(*,'(A,A)')'Error when reading file ', trim(chfile)
      write(*,*)'Error code was', iost
      call abinerror('init')
   end subroutine err_read

   subroutine PrintLogo(values1)
   integer,dimension(8),intent(out) :: values1
   call date_and_time(VALUES=values1)

print '(a)','                    _____      _     _      _ '
print '(a)','        /\         |  _  \    | |   |  \   | |'
print '(a)','       /  \        | |_|  |   | |   | | \  | |'
print '(a)','      / /\ \       |     /    | |   | |\ \ | |'
print '(a)','     / /__\ \      |=====|    | |   | | \ \| |'
print '(a)','    / /____\ \     |  _   \   | |   | |  \ | |'
print '(a)','   / /      \ \    | |_|  |   | |   | |   \  |'
print '(a)','  /_/        \_\   |_____/    |_|   |_|    \_|'
print '(a)',' '

print '(a)','     version 1.0'
print '(a)',' D. Hollas, O.Svoboda, M. Oncak, P. Slavicek 2011-2015'
print '(a)',' '

call print_compile_info()

write(*,*)'Job started at:'
write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")values1(5),':', &
        values1(6),':',values1(7),'  ',values1(3),'.',values1(2),'.',&
        values1(1)

   end subroutine PrintLogo

   subroutine PrintHelp()
   implicit none
    print '(a)', ''
    print '(a)', 'ABIN: Multipurpose ab initio MD program.'
    print '(a)', ''
    call print_compile_info()
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help               print help and exit'
    print '(a)', '  -i <input_parameters>    default: input.in'
    print '(a)', '  -x <input_coordinates>   default: mini.dat'
    print '(a)', '  -v <input_velocities>    no default'
    print '(a)', ''
   end subroutine PrintHelp

   subroutine Get_cmdline(chinput, chcoords, chveloc )
   character(len=*),intent(inout)   :: chinput, chcoords, chveloc
   character(len=len(chinput))   :: arg
   integer            :: i
   logical            :: lexist
   
   i=0
   do while (i < command_argument_count())

     i=i+1
     call get_command_argument(i, arg)
   
      select case (arg)
      case ('-h', '--help')
         call PrintHelp()
         stop 0

      case ('-i')
         i=i+1
         call get_command_argument(i, arg)
         !-format specifier is needed here in case of slashes
         read(arg,'(A)')chinput
         chinput=trim(chinput)

      case ('-x')
         i=i+1
         call get_command_argument(i, arg)
         read(arg,'(A)')chcoords
         chcoords=trim(chcoords)
      case ('-v')
         i=i+1
         call get_command_argument(i, arg)
         read(arg,'(A)')chveloc
         chveloc=trim(chveloc)
      case default
         write(*,*)'Invalid command line argument!'
         call abinerror('Get_cmdline')
      end select

   end do
   !check for existence of input files


   inquire(file=chinput,exist=lexist)
   if (.not.lexist)then
      write(*,*)'FATAL: The following input file does not exists!'
      write(*,*)chinput
      call abinerror('Get_Cmdline')
   end if

#ifndef MPI
   inquire(file=chcoords,exist=lexist)
   if (.not.lexist)then
      write(*,*)'FATAL: Input file does not exists!'
      write(*,*)chcoords
      call abinerror('Get_Cmdline')
   end if
   if (chveloc.ne.'')then
      inquire(file=chveloc,exist=lexist)
      if (.not.lexist)then
         write(*,*)'FATAL: The following input file does not exists!'
         write(*,*)chveloc
         call abinerror('Get_Cmdline')
      end if
   end if
#endif

   end subroutine Get_cmdline

   subroutine print_compile_info()
   character(len=1024)   :: cmdline
   !   include 'date.inc'

   print *,'Compiled at  ', DATE
   print *,COMMIT
!$ print *,'Compiled with parallel OpenMP support for PIMD.'
#ifdef USEFFTW
   write(*,*)'Compiled with FFTW support.'
#endif
#ifdef CP2K
   write(*,*)'Compiled with in-built CP2K interface.'
#endif
#ifdef PLUM
   write(*,*)'Compiled with PLUMED (static lib).'
#endif
#ifdef MPI
   write(*,*)'Compiled with MPI support.'
   write(*,*)'(used for REMD and direct CP2K and TeraChem interfaces.)'
#endif
   print *,' '

#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 6
   print *, 'This program was compiled by ',  &
             compiler_version(), ' using the options: '
   print *,     compiler_options()
#endif

print '(a)','#################### RUNTIME INFO ####################'
call get_command(cmdline)
write(*,*)trim(cmdline)
call flush(6)
call get_command_argument(0, cmdline)
call system('ldd '//cmdline)
print '(a)','######################################################'

   end subroutine print_compile_info

end subroutine init


