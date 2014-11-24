!-----This messy function performs many thing, among others:
! 1. Read input
! 2. Input sanity check
! 3. Allocation of arrays
! 4. Reading restart OR reading input geometry
! 5. Initialize velocities
! 6. Initialize everything else
! At this moment, coordinate and velocity transformations are NOT performed here
! Surface hopping is NOT initialized here.
! Daniel Hollas,9.2.2012 

subroutine init(dt)
      use mod_const
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
      use mod_fftw3,    only: fftw_init
      use mod_random
      use mod_guillot,  only: inames_guillot
      use mod_utils
      use mod_vinit,    only: vinit, scalevelocities
      use mod_density
      use mod_shake
      use mod_minimize, only: gamm, gammthr
      use mod_analysis, only: restin
      use mod_water,    only: watpot, check_water
      implicit none
      real(DP),intent(out) :: dt
      real(DP) :: masses(MAXTYPES)
      real(DP)  :: rans(10)
      integer :: iw,iat,natom1,imol,shiftdihed=1, iost
      integer :: error, getpid, nproc=1, iknow=0, ipom, ipom2=0
      character(len=2) :: massnames(MAXTYPES)
      character(len=10)   :: chaccess
      character(len=200)  :: chinput, chcoords
      character(len=200)  :: chiomsg
      LOGICAL :: file_exists
      real(DP)  :: wnw=5.0d-5
!$    integer :: nthreads,omp_get_max_threads
! wnw "optimal" frequency for langevin (inose=3) 
      REAL, POINTER, DIMENSION(:) :: VECPTR => NULL ()  !null pointer
!      REAL, POINTER, DIMENSION(:,:) :: VECPTR2 => NULL ()
      REAL, POINTER  :: REALPTR => NULL ()

      namelist /general/natom, pot,ipimd,istage,nwalk,nstep,icv,ihess,imini,nproc,iqmmm, &
               nwrite,nwritex,nwritev,dt,irandom,nabin,irest,nrest,anal_ext,isbc,rb_sbc,kb_sbc,gamm,gammthr,conatom, &
               parrespa,dime,ncalc,idebug,enmini,rho,iknow, watpot

      namelist /nhcopt/ inose,temp,temp0,nchain,ams,tau0,imasst,wnw,nrespnose,nyosh,      &
                        scaleveloc,readNHC,readQT,initNHC,nmolt,natmolt,nshakemol
      namelist /system/ masses,massnames,nbin,nbin_ang,ndist,dist1,dist2,xmin,xmax,disterror, &
                        nang,ang1,ang2,ang3,ndih,dih1,dih2,dih3,dih4,shiftdihed, &
                        k,r0,k1,k2,k3,De,a, &
                        Nshake,ishake1,ishake2,shake_tol
      namelist /sh/     istate_init,nstate,substep,deltae,integ,inac,nohop,phase,alpha,popthr, &
                        nac_accu1, nac_accu2, popsumthr, energydifthr, energydriftthr, adjmom, revmom
      namelist /qmmm/   natqm,natmm,q,rmin,eps,attypes
      namelist /nab/    ipbc,alpha_pme,kappa_pme,cutoff,nsnb,ips,epsinf,natmol, nmol

      chcoords='mini.dat'
      chinput='input.in'
      call Get_cmdline(chinput, chcoords)
      write(*,*)'I will read parameters from input file ',chinput
      write(*,*)'I will read xyz coordinates from file ',chcoords

      dt=-1  
      error=0
      !-READING INPUT----------------------------------------- 

      open(150,file=chinput, status='OLD', delim='APOSTROPHE', action = "READ") !here ifort has some troubles
      read(150,general)
      rewind(150)

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


      ! the namelist system does not need to be present
      read(150,system,iostat=iost,iomsg=chiomsg)
      rewind(150)
      !check, whether we hit End-Of-File or other error
      if(IS_IOSTAT_END(iost))then  !fortran intrinsic for EOF
         write(*,*)'Namelist "system" not found. Ignoring...'
      else if (iost.ne.0)then
         write(*,*)'ERROR when reading namelist "system".'
         write(*,*)chiomsg
         call abinerror('init')
      else
         massnames=LowerToUpper(massnames)
      end if

#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
      allocate ( natmolt(natom)  )
      allocate ( nshakemol(natom)  )
#endif
      natmolt   = 0
      nshakemol = 0

      read(150,nhcopt)
      rewind(150)

      pot=UpperToLower(pot)
      if(ipimd.eq.2)then
         read(150,sh)
         rewind(150)
         integ=UpperToLower(integ)
      end if

      if(iqmmm.gt.0.or.pot.eq.'mm')then
         read(150,qmmm)
         rewind(150)
      end if

      if(iqmmm.eq.2.or.pot.eq.'nab')then
#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
         allocate ( natmol(natom) )
#endif
         natmol = 0
         read(150,nab)
         rewind(150)
      end if

      close(150)
!-----END OF READING INPUT----------------------------------------- 

!$    call OMP_set_num_threads(nproc)
!$    nthreads=omp_get_max_threads()

!-----HERE WE CHECK FOR ERRORS IN INPUT-----------------------------------------------
      call check_inputsanity()

!     resetting number of walkers to 1 in case of classical simulation      
      if(ipimd.eq.0)then
              write(*,*)'ipimd=0,Resetting number of walkers to 1.'
              write(*,*)'ipimd=0,Using velocity Verlet integrator'
              md=2
              nwalk=1
              nabin=1   !TODO:safety for respashake code
                        !we should probably copy shake to velocity verlet
                        !algorithm as well
      endif
      if(iqmmm.eq.0.and.pot.ne.'mm')then
              natqm=natom
      endif

!for surface hopping      
      if(ipimd.eq.2)then
              nwalk=ntraj
              md=2
              nabin=1
      endif

      if(ipimd.eq.1)then
              write(*,*)'ipimd=1,using RESPA integrator'
              md=1
      endif

!-------------------------INITIALIZATION OF FFTW ROUTINES-------
      if(istage.eq.2) call fftw_init(nwalk)


!     allocate all basic arrays and set them to 0.0d0
      call allocate_arrays( natom, nwalk+1 )

!-----READING GEOMETRY
      open(111,file=chcoords,status = "old", action = "read") 
      read(111,*)natom1
      if(natom1.ne.natom)then
        write(*,'(A,A)')'No. of atoms in ',chinput
        write(*,'(A,A)')'and in ',chcoords
        write(*,*)'do not match.'
        call abinerror('init')
      endif
      read(111,*)
      do iat=1,natom

        read(111,*)names(iat),x(iat,1),y(iat,1),z(iat,1)
        names(iat)=LowerToUpper(names(iat))
        x(iat,1)=x(iat,1)*ang
        y(iat,1)=y(iat,1)*ang
        z(iat,1)=z(iat,1)*ang

      enddo 
      close(111)

      do iw=1,nwalk
       do iat=1,natom
       x(iat,iw)=x(iat,1)
       y(iat,iw)=y(iat,1)
       z(iat,iw)=z(iat,1)
       enddo
      enddo
!-----END OF READING GEOMETRY      

!-----conversion of temperature from K to au
      if (inose.ne.0) write(*,*)'Target temperature in Kelvins =',temp
      temp=temp/autok

      if (ihess.eq.1)then 
       allocate ( hess(natom*3,natom*3,nwalk) )
       allocate ( cvhess_cumul(nwalk) )
       cvhess_cumul=0.0d0
       if (pot.eq.'nab')then
!!$OMP PARALLEL
        allocate ( h(natom*natom*9) )
!!$OMP END PARALLEL
       endif
      endif

!----SHAKE initialization,determining the constrained bond lenghts
      if(nshake.ge.1)then
       write(*,*)'Setting distances for SHAKE from mini.dat'
       call shake_init(x,y,z)
      else 
       do iat=1,natom
        nshakemol(iat)=0
       enddo
      endif

     call dist_init() !zeroing distribution arrays

!----In case of big systems, we don't want to manually set am and names arrays.
!----Currently supported for most of the elements
     call mass_init(masses, massnames)
     !lower the second character of atom name
     !this is because of terachem
      do iat=1,natom
         names(iat)(2:2)=UpperToLower(names(iat)(2:2))
      end do

      if (pot.eq.'mmwater') call check_water(natom, names)


!-----THERMOSTAT INITIALIZATION------------------ 
!----MUST BE BEFORE RESTART DUE TO ARRAY ALOCATION
!    call vranf(rans,0,IRandom,6)  !initialize prng,maybe rewritten during restart
     call gautrg(rans,0,IRandom,6)  !initialize prng,maybe rewritten during restart
     if (inose.eq.1) call nhc_init()
     if (inose.eq.2) call gle_init(dt*0.5/nabin) !nabin is set to 1 unless ipimd=1
     if (inose.eq.3) call wn_init(dt*0.5,wnw)


!----performing RESTART from restart.xyz
     if(irest.eq.1) call restin(x,y,z,vx,vy,vz,it)

     write(*,*)
     write(*,*)'-------SIMULATION PARAMETERS---------------'
     write(*,nml=general)
     write(*,*)
     write(*,nml=system)
     if ( inose.ge.1 ) write(*,nml=nhcopt)
     if ( ipimd.eq.2 ) write(*,nml=sh)
     write(*,*)

!------------END OF INPUT SECTION---------------------

!----------------INITIALIZATION-----------------------

!--Here we ensure, that previous files are deleted----
if(irest.eq.0)then
  open(10,file='movie_mini.xyz')
  close(10,status='delete')
  chaccess='SEQUENTIAL'
!  if (nwritev.gt.0)then !I think its not needed
!    open(13,file='vel.dat')
!    close(13,status='delete')
!  endif
else
  chaccess='APPEND'
endif

if (nwritev.gt.0)then
   open(13,file='vel.dat',access=chaccess,action='write')
endif

if (ipimd.ne.1)then
   open(1,file='energies.dat',access=chaccess,action='write')
   write(1,*)'#        Time[fs] E-potential           E-kinetic     E-Total    E-Total-Avg'
end if
if(ipimd.eq.1.or.icv.eq.1)then
  open(7,file='est_energy.dat',access=chaccess,action='write')
  write(7,*)'#     Time[fs] E-potential  E-primitive   E-virial  CumulAvg_prim  CumulAvg_vir'
endif

open(2,file='temper.dat',access=chaccess,action='write')
write(2,*)'#      Time[fs] Temperature T-Average Conserved_quantity_of_thermostat'

if(ipimd.eq.2)then
  open(3,file='pop.dat',access=chaccess,action='write')
  write(3,*)'#    Time[fs] CurrentState   Populations Sum-of-Populations'
  open(4,file='prob.dat',access=chaccess,action='write')
  write(4,*)'#    Time[fs] CurrentState   Probabilities'
  open(8,file='PES.dat',access=chaccess,action='write')
  write(8,*)'#    Time[fs] Potential energies'
  open(14,file='nacm_all.dat',access=chaccess,action='write')
  open(15,file='dotprod.dat',access=chaccess,action='write')
  write(15,*)'#    Time[fs] dotproduct(i,j) [i=1,nstate-1;j=i+1,nstate]'
  if(idebug.gt.1)then
   open(16,file='bkl.dat',access=chaccess,action='write')
   write(16,*)'# Hopping probabilities - bkl(i) [i=1,nstate]'
   open(17,file='coef.dat',access=chaccess,action='write')
   write(17,*)'# WF coefficients c_real(i),i=1,nstate c_imag(i),i=1,nstate'
   if(phase.eq.1)then
      open(18,file='phase.dat',access=chaccess,action='write')
      write(18,*)'# Lower triangular matrix of gamma (phase)  gamma(i,j) [i=1,nstate ;j=1,i-1]'
   end if
  endif
endif

if(isbc.eq.1)then
 open(11,file='r.dat',access=chaccess,action='write')
 write(11,*)'#Radius    density'
endif

if(icv.eq.1)then
 open(122,file='cv.dat',access=chaccess,action='write')
 write(122,*)'#         Time[fs]  Cv-prim   Cv-vir  Cv_cumul_prim  Cv_cumul_vir'
 close(122)
 if(ihess.eq.1)then
  open(123,file='cv_dcv.dat',access=chaccess,action='write')
  write(123,*)'#         Time[fs]  Cv-DCV   Cv-DCV'
  close(123)
 endif
endif

!------------------------------------------------------

      pid=GetPID()
      write(*,*)'Pid of the current proccess is:',pid



      if(pot.eq.'2dho')then
       f=0 !temporary hack
      endif
      if(nchain.gt.1)then
       f=0 !what about nchain=1?
       !what about massive therm?
      endif

!-----SETTING initial velocities according to the Maxwell-Boltzmann distribution
!     TODO: odstraneni rotace
      if(irest.eq.0)then
         if (temp0.ge.0)then
            call vinit(TEMP0,am,vx,vy,vz)
         else
            call vinit(TEMP,am,vx,vy,vz)
         end if
      end if

      if(conatom.gt.0)then
         write(*,*)'Removing initial velocity of constrained atoms.'
         call constrainP(vx,vy,vz)
      endif

      ! If scaleveloc=1, scale initial velocitites
      ! Otherwise, just print the temperature.
      call ScaleVelocities(vx, vy, vz)

       
!----some stuff for spherical boundary onditions
      if(isbc.eq.1) call sbc_init(x,y,z)

!-----inames initialization for guillot rm MM part. 
!-----We do this also because string comparison is very costly
      if(iqmmm.eq.3.or.pot.eq.'mm'.or.pot=='guillot') allocate( inames(natom) )
      if(pot.eq.'guillot') call inames_guillot()

      if(iqmmm.eq.3.or.pot.eq.'mm')then 
         attypes=LowerToUpper(attypes)
         call inames_init()
         call ABr_init()
      endif

      if (iqmmm.eq.3.or.pot.eq.'mm') write(*,nml=qmmm)

      if (pot.eq."nab".or.iqmmm.eq.2)then
       if (alpha_pme.lt.0) alpha_pme = pi/cutoff
       if (kappa_pme.lt.0) kappa_pme = alpha_pme
       call nab_init(alpha_pme,cutoff,nsnb,ipbc,ips,iqmmm) !C function...see nabinit.c

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
!       write(*,*)x(iat,1)/ang,y(iat,1)/ang,z(iat,1)/ang,charges(iat)
        call ewald(VECPTR,VECPTR,charges,REALPTR,boxx,boxy,boxz,cutoff,alpha_pme,kappa_pme,epsinf,natom,ipom2)
        boxx=boxx*ang
        boxy=boxy*ang
        boxz=boxz*ang
        boxx2=0.5d0*boxx
        boxy2=0.5d0*boxy
        boxz2=0.5d0*boxz
       endif
      endif


!--------END OF INITIALIZATION-------------------

   call flush()

   contains

   subroutine check_inputsanity()
   character(len=*),parameter :: chknow='If you know what you are doing, &
    &  set  iknow=1 (namelist general) to proceed.'

      !we should exclude all non-abinitio options, but whatever....
!$    if(nthreads.gt.1.and.(ipimd.ne.1.or.pot.eq.'nab'))then
!$     write(*,*)'ERROR: Parallel execution is currently only supported with ab initio PIMD (ipimd=1)'
!$     call abinerror('init')
!$    endif

      if(nproc.gt.1)then
!$       if(.false.)then
         write(*,*)'ERROR: This executable was not compiled with parallel support.'
         error=1
!$       end if
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
      if(pot.eq.'none')then
       write(*,*)'FATAL: Variable "pot" not specified.Exiting now...'
       stop 1
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

      if(ipimd.eq.1.and.inose.ne.1.and.inose.ne.2)then
       write(*,*)'You have to use Nosé-Hoover or quantum thermostat with PIMD!(inose=1 or 2)'
       error=1
      endif
      if(ipimd.lt.0.or.ipimd.gt.3)then
       write(*,*)'ipimd has to be 0,1,2 or 3.'
       error=1
      endif
      if(istage.ne.1.and.istage.ne.0.and.istage.ne.2)then
       write(*,*)'istage has to be 0,1 or 2'
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
      if(inose.eq.3)then
       write(*,*)'Langevin thermostat was not extensively tested.'
       write(*,*)chknow
       if (iknow.ne.1) error=1
      endif
      if(istage.eq.1.and.ipimd.ne.1)then
      write(*,*)'The staging transformation is only meaningful for PIMD'
       error=1
      endif
      if(istage.eq.2.and.ipimd.ne.1)then
      write(*,*)'The normal mode transformation is only meaningful for PIMD. Exiting...'
       error=1
      endif
      if(istage.eq.0.and.ipimd.eq.1.and.inose.ne.2)then
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
       write(*,*)chknow
       if (iknow.ne.1)error=1
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

      INQUIRE(FILE="movie.xyz", EXIST=file_exists)
      if(file_exists)then
       if(irest.eq.0)then
        write(*,*)'File "movie.xyz" exists.Please (re)move it or set irest=1.'
        error=1
       else
        write(*,*)'File "movie.xyz" exists and irest=1.Trajectory will be appended.'
       endif
      endif

      INQUIRE(FILE="restart.xyz", EXIST=file_exists)
      if(file_exists)then
       if(irest.eq.0)then
        write(*,*)'File "restart.xyz" exists. Please (re)move it or set irest=1.'
        error=1
       endif
      else
       if(irest.eq.1)then
         write(*,*)'File restart.xyz not found.' 
         error=1
       endif 
      endif

      if(error.eq.1)then
         write(*,*)'Input errors were found! Exiting now...'
         call abinerror('init')
      endif


   end subroutine check_inputsanity

end subroutine init


