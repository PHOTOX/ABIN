!--Input and Initialization                         by Daniel Hollas,9.2.2012
!-- some comment
      subroutine init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt)
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
      use mod_sbc
      use mod_fftw3
      implicit none
      real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(out)   :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
      real*8,intent(out)   :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
      real*8,intent(inout) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8,intent(inout) :: dt
      real*8  :: ran1,rans(10)
      real*8  :: vx1(npartmax),vy1(npartmax),vz1(npartmax)
      integer :: irnd(nwalkmax),scaleveloc=1,readNHC
      integer :: iw,iat,inh,natom1,ifirst,itrj,ist1,imol,shiftdihed=1
      integer :: error=0,getpid,nproc=1,iknow=0,ipom,ipom2=0,is
      character(len=2)  :: shit
      character(len=10) :: chaccess
      LOGICAL :: file_exists
      real*8  :: wnw=5.0e-5,pom=0.0d0,ekin_mom=0.0d0,temp_mom=0.0d0,scal
!$    integer :: nthreads,omp_get_max_threads
! wnw "optimal" frequency for langevin (inose=3) 
      REAL, POINTER, DIMENSION(:) :: VECPTR => NULL ()  !null pointer
!      REAL, POINTER, DIMENSION(:,:) :: VECPTR2 => NULL ()
      REAL, POINTER  :: REALPTR => NULL ()

      namelist /general/pot,ipimd,istage,nwalk,nstep,icv,ihess,imini,nproc,iqmmm, &
               nwrite,nwritex,nwritev,dt,irandom,nabin,irest,nrest,anal_ext,isbc,ibag,rb_sbc,kb_sbc,gamm,gammthr,conatom, &
               parrespa,dime,ncalc,idebug,enmini,rho,iknow

      namelist /nhcopt/ inose,temp,nchain,ams,tau0,imasst,wnw,nrespnose,nyosh,scaleveloc,readNHC,initNHC,nmolt,natmolt
      namelist /system/ natom,am,imass_init,nbin,nbin_ang,ndist,dist1,dist2,xmin,xmax, &
                        nang,ang1,ang2,ang3,ndih,dih1,dih2,dih3,dih4,shiftdihed, &
                        k,r0,k1,k2,k3,De,a, &
                        Nshake,ishake1,ishake2,shake_tol,nmol,nshakemol,names,natmol
      namelist /sh/     ntraj,istate_init,nstate,substep,deltae,integ,inac,nohop,alpha,popthr,nac_accu1,nac_accu2 !TODO: some checking for accu1/2
      namelist /qmmm/   natqm,natmm,LJcomb,q,rmin,eps,attypes,inames,qmmmtype
      namelist /nab/    ipbc,alpha_pme,kappa_pme,cutoff,nsnb,ips,epsinf


      do iat=1,npartmax
       am(iat)=0.0d0
       names(iat)=''
       attypes(iat)=''
      enddo
      dt=-1  

      xmin=0.5
      xmax=3.5

      open(150,file="input.in", status='OLD', delim='APOSTROPHE', action = "READ") !here ifort has some troubles
      read(150,general)
      if(irest.eq.1)then
       readnhc=1   !readnhc has precedence before initNHC
       initNHC=0   !i.e. if(readnhc.eq.1.and.initNHC.eq.1)
       scaleveloc=0  !do not scale velocities when restarting a job
      else         !then nhc momenta from restart.xyz will be used
       readnhc=0
       initNHC=1
       scaleveloc=1
      endif
      read(150,system)
      read(150,nhcopt)
      if(ipimd.eq.2) read(150,sh)
      if(iqmmm.eq.1.or.pot.eq.'mm') read(150,qmmm)
      if(qmmmtype.eq."nab".or.pot.eq.'nab') read(150,nab)
      close(150)

!$    call OMP_set_num_threads(nproc)
!$    nthreads=omp_get_max_threads()
!$    write(*,*)nthreads

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
      endif

      if(ipimd.eq.1)then
              write(*,*)'ipimd=1,using RESPA integrator'
              md=1
      endif

!-------------------------INITIALIZATION OF FFTW ROUTINES-------
      if(istage.eq.2) call fftw_init(nwalk)

!-----HERE WE CHECK FOR ERRORS IN INPUT-----------------------------------------------
!$    if(nthreads.gt.1.and.(pot.ne.'nab'.and.pot.ne.'harm'.and.pot.ne.'guillot'.and.pot.ne.'g09'))then
!$     write(*,*)'Parallel execution is currently only supported with NAB,G09 or model potentials.Exiting...'
!$     stop
!$    endif

      !-----Check,whether input variables don't exceeds array limits
      if(ntraj.gt.ntrajmax)then
       write(*,*)'Maximum number of trajectories is:'
       write(*,*)ntrajmax
       write(*,*)'Adjust variable ntrajmax in modules.f90'
       stop
      endif
      if(nstate.gt.nstmax)then
       write(*,*)'Maximum number of states is:'
       write(*,*)nstmax
       write(*,*)'Adjust variable nstmax in modules.f90'
       stop
      endif
      if(nwalk.ge.nwalkmax)then
       write(*,*)'Maximum number of random walkers is:'
       write(*,*)nwalkmax
       write(*,*)'Adjust variable nwalkmax in modules.f90'
       stop
      endif
      if(natom.gt.npartmax)then
       write(*,*)'Maximum number of atoms is:'
       write(*,*)npartmax
       write(*,*)'Adjust variable npartmax in modules.f90'
       stop
      endif
      if(nchain.gt.maxchain)then
       write(*,*)'Maximum number of Nose-Hoover chains is:'
       write(*,*)maxchain
       write(*,*)'Adjust variable maxchain in modules.f90'
       stop
      endif
      if(nshake.gt.nshakemax)then
       write(*,*)'Maximum number of bonds with SHAKE is:'
       write(*,*)nshakemax
       write(*,*)'Adjust variable nshakemax in modules.f90'
       stop
      endif
      if(ndist.ge.ndistmax)then
       write(*,*)'Maximum number of bonds for binning is:'
       write(*,*)ndistmax
       write(*,*)'Adjust variable ndistmax in modules.f90'
       stop
      endif
      if(nbin.ge.nbinmax)then
       write(*,*)'Maximum number of bins for densities is:'
       write(*,*)nbinmax
       write(*,*)'Adjust variable nbinmax in modules.f90'
              error=1
       stop
      endif
!----------HERE we check for errors in input.      
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
      if(iqmmm.eq.1.and.qmmmtype.ne.'nab'.and.qmmmtype.ne.'abin')then
       write(*,*)'Set qmmmtype to "abin" or "nab"(using Amber ff)'
       error=1
      endif
      if(iqmmm.eq.1.and.ipbc.eq.1)then
       write(*,*)'QM/MM with PBC not supported !'
       error=1
      endif
      if(integ.ne.'euler'.and.integ.ne.'rk4'.and.integ.ne.'butcher')then
       write(*,*)'integ must be "euler" or "rk4".'
       error=1
      endif
      if(deltae.lt.0)then
       write(*,*)'Parameter deltae must be non-negative number.'
       error=1
      endif
!      if(irandom.le.9999)then
!              write(*,*)'Error in getting irandom from input.Possibly too short&
!              or too long'
!              error=1
!      endif
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
              write(*,*)'Input error: nwritex must be positive.'
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
       write(*,*)'Parameter INAC must be 0,1 or 2.'   !be very carefull if you change this!
       if(iknow.ne.1) error=1
      endif
      if(irest.eq.1.and.scaleveloc.eq.1)then
       write(*,*)'irest=1 AND scaleveloc=1.'
       write(*,*)'You are trying to scale to velocities read from restart.xyz.'
       write(*,*)'I assume this is an error in input. Exiting...'
       write(*,*)'If you know, what you are doing, set  iknow=1 (section general) to proceed.'
       if(iknow.ne.1) error=1
      endif
      if(inose.eq.1.and.ipimd.eq.2)then
       write(*,*)'Thermostating not meaningful for surface hopping simulation.Exiting.'
       write(*,*)'If you know, what you are doing, set  iknow=1 (section general) to proceed.'
       if(iknow.ne.1) error=1
      endif
      if(istate_init.gt.nstate)then
       write(*,*)'Error:Initial state > number of computed states. Exiting...'
       error=1
      endif
      if(nac_accu1.le.0.or.nac_accu2.lt.0)then
       write(*,*)'Input error:NAC precision must be a positive integer.'
       write(*,*)'The treshold is then 10^-(nac_accu)'
       error=1
      endif
      if(nac_accu1.le.nac_accu2)then
       write(*,*)'nac_accu1 < nac_accu2'
       write(*,*)'I will compute NAC only with default accuracy:',nac_accu1
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
      if(imasst.eq.0)then
       do imol=1,nmolt
        if(natmolt(imol).le.0)then
         write(*,*)'Number of atoms in molecules not specified!Set array natmolt properly.'
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
      if(istage.eq.1.and.ipimd.ne.1)then
      write(*,*)'The Normal Mode transformation is only meaningful for PIMD. Exiting...'
       error=1
      endif
      if(istage.eq.0.and.ipimd.eq.1)then
       write(*,*)'PIMD should be done with staging or normal mode transformation! Exiting...'
       write(*,*)'If you know, what you are doing, set iknow=1 (section general) to proceed.'
       if (iknow.ne.1) error=1
      endif
      if(istage.eq.1.and.inose.eq.2)then
       write(*,*)'The staging transformation is not compatible with GLE thermostat'
       error=1
      endif
      if(nyosh.ne.1.and.nyosh.ne.3.and.nyosh.ne.7)then
       write(*,*)'Variable nyosh(order of Suzuki-Yoshiga scheme) must be 1,3 or 7'
       error=1
      endif
      if(nyosh.le.1.and.inose.eq.1)then
       write(*,*)'It is strongly reccomended to use Suzuki-Yoshida scheme when using Nose-Hoover thermostat (nyosh 3 or 7).'
       write(*,*)'If you know, what you are doing, set iknow=1 (section general) to proceed.'
       if (iknow.ne.1)error=1
      endif
      if(nrespnose.lt.3.and.inose.eq.1)then
       write(*,*)'Variable nrespnose < 3! Assuming this is an error in input and exiting.'
       write(*,*)'If you know, what you are doing, set iknow=1 (section general) to proceed.'
       if (iknow.ne.1)error=1
      endif
      if(nrespnose.le.0)then
       write(*,*)'Variable nrespnose must be positive integer'
       error=1
      endif
      if(irest.eq.1.and.irest.eq.0)then
      write(*,*)'irest has to be 1 or zero'
       error=1
      endif
      if(nshake.gt.0.and.ipimd.eq.1)then
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
       write(*,*)'PIMD with SHAKE cannot use massive thermostating!Exiting... !'
       error=1
      endif

      if(pot.eq.'2dho'.and.natom.gt.1)then
       write(*,*)'Only 1 particle is allowed for 2D harmonic oscillator!'
       error=1
      endif
      if(pot.eq.'mm'.and.iqmmm.gt.1)then
       write(*,*)'Pot="mm"is not compatible with iqmmm=1!'
       error=1
      endif
      if((natmm+natqm.ne.natom).and.iqmmm.eq.1)then
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
       write(*,*)'If you know, what you are doing, set  iknow=1 (section general) to proceed.'
        if(iknow.ne.1) error=1
       endif
      endif

      if(temp.lt.1.and.inose.ge.1)then
       write(*,*)'WARNING!:Temperature below 1K. Are you sure?'
       write(*,*)'This is probably mistake in input.Exiting...'
       write(*,*)'If you know, what you are doing, set  iknow=1 (section general) to proceed.'
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
              stop 1
      endif
!--END OF ERROR CHECKING

!------------READING GEOMETRY
      open(111,file='mini.dat',status = "old", action = "read") 
      read(111,*)natom1
      if(natom1.ne.natom)then
        write(*,*)'No. of atoms in input.in and in mini.dat do not match.'
        stop 1
      endif
      read(111,*)
      do iat=1,natom

        read(111,*)shit,x(iat,1),y(iat,1),z(iat,1)
        
        if(imass_init.eq.0)then
         if(am(iat).le.0)then
          write(*,*)'Error on input:Mass array do not match names array.'
         endif
         if(shit.ne.names(iat))then
          write(*,*)'Names of atoms in mini.dat and input.in do not match!'
          stop 1
         endif
        else
         names(iat)=shit
        endif

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

      it=1


!conversion of temperature from K to au
      write(*,*)'Target temperature in Kelvins =',temp
      temp=temp/autok
!-----ZEROING some arrays
      do iw=1,nwalk
       cvhess_cumul(iw)=0.0d0
       do iat=1,natom
        fxq(iat,iw)=0.0d0 
        fyq(iat,iw)=0.0d0 
        fzq(iat,iw)=0.0d0 
        fxc(iat,iw)=0.0d0 
        fyc(iat,iw)=0.0d0 
        fzc(iat,iw)=0.0d0 
        vx(iat,iw)=0.0d0 
        vy(iat,iw)=0.0d0 
        vz(iat,iw)=0.0d0 
        enddo
       enddo
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

!-----THERMOSTAT INITIALIZATION------------------
     if (inose.gt.0)then
      call gautrg(rans,0,IRandom,6)  !inicializace prng
     endif
     if (inose.eq.1) call nhc_init()
     if (inose.eq.2.and.md.eq.2) call gle_init(dt*0.5,irandom)
     if (inose.eq.2.and.md.eq.1) call gle_init(dt*0.5/nabin,irandom)  
     if (inose.eq.3) call wn_init(dt*0.5,wnw,irandom)

!in case of big systems, we dont want to manually set am and names arrays
!----currently supported only for some elements
     if (imass_init.eq.1) call mass_init()

!-----performing RESTART from restart.xyz
     if(irest.eq.1)then
      write(*,*)'irest=1,Reading geometry,velocities and NHC momenta from restart.xyz'

      open(111,file='restart.xyz',status = "OLD", action = "READ")
      read(111,*)it
      it=it+1
      read(111,*)
      do iw=1,nwalk
       do iat=1,natom
        read(111,*)x(iat,iw),y(iat,iw),z(iat,iw)
       enddo
      enddo

      read(111,*)
      do iw=1,nwalk
       do iat=1,natom
        read(111,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
       enddo
      enddo

      if(ipimd.eq.2)then
       read(111,*)
       do itrj=1,ntraj
        read(111,*)istate(itrj)
        do ist1=1,nstate
         read(111,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
        enddo
       enddo
      endif

      if(inose.eq.1.and.readNHC.eq.1)then
       read(111,*)
      if(imasst.eq.1)then
       do inh=1,nchain
        do iw=1,nwalk
         do iat=1,natom
          read(111,*)pnhx(iat,iw,inh),pnhy(iat,iw,inh),pnhz(iat,iw,inh)
         enddo
        enddo
       enddo

     else

      do inh=1,nchain
       do iw=1,nwalk
        do iat=1,nmolt
          read(111,*)pnhx(iat,iw,inh)
        enddo
       enddo
      enddo
     endif

     endif

     if(inose.eq.2)then
      read(111,*)
      if(nwalk.eq.1)then
       do iat=1,natom
        read(111,*)(gp(iat,is),is=1,ns)
       enddo
      else
       do iw=1,nwalk
        do iat=1,natom
         read(111,*)(ps(iat,is,iw),is=1,ns)
        enddo
       enddo
     endif
      read(111,*)langham
     endif

!reading cumulative averages of various estimators
       read(111,*)
       read(111,*)est_temp_cumul
       read(111,*)est_prim_cumul,est_vir_cumul
       read(111,*)entot_cumul
      
       if(icv.eq.1)then
        read(111,*)est_prim2_cumul,est_prim_vir,est_vir2_cumul
        read(111,*)cv_prim_cumul,cv_vir_cumul
        if(ihess.eq.1)then 
         read(111,*)cv_dcv_cumul
         do iw=1,nwalk
          read(111,*)cvhess_cumul(iw)
         enddo
        endif
       endif

      
      close(111)
     endif

!!!!! END OF RESTART 


      write(*,*)'-------SIMULATION PARAMETERS---------------'
      write(*,nml=general)
      if ( inose.ge.1 ) write(*,nml=nhcopt)
      if ( ipimd.eq.2 ) write(*,nml=sh)
      write(*,nml=system)

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

open(1,file='energies.dat',access=chaccess,action='write')
write(1,*)'#        Time[fs] E-potential           E-kinetic     E-Total    E-Total-Avg'
if(ipimd.eq.1.or.icv.eq.1)then
  open(7,file='est_energy.dat',access=chaccess,action='write')
  write(7,*)'#     Time[fs] E-potential  E-primitive   E-virial  CumulAvg_prim  CumulAvg_vir'
endif

if(inose.gt.0)then
  open(2,file='temper.dat',access=chaccess,action='write')
  write(2,*)'#        Time[fs] Temperature           Temp-Average'
endif

if(ipimd.eq.2)then
  open(3,file='pop.dat',access=chaccess,action='write')
  write(3,*)'#    Time[fs] CurrentState   Populations Sum-of-Populations'
  open(4,file='prob.dat')
  write(4,*)'#    Time[fs] CurrentState   Probabilities'
  open(8,file='PES.dat')
  write(8,*)'#    Time[fs] Potential energies'
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

      do iat=1,natom
       am(iat)=am(iat)*amu
      enddo

!-----PID of the proccess
      pid=GetPID()
      write(*,*)'Pid of the current proccess is:',pid

!-----inames initialization, currently only for guillot. 
!-----we do this because string comparison is very costly
      if(pot.eq.'guillot') call inames_guillot()


      if(pot.eq.'2dho')then
       f=0 !temporary hack
      endif
      if(nchain.gt.1)then
       f=0 !what about nchain=1?
       !what about massive therm?
      endif
!-----SETTING initial velocities according to the Maxwell-Boltzmann distribution
      if(irest.eq.0)then
!-----Random number generation. Used for Maxwell-Boltzmann sampling      
      ifirst=-irandom
      do iw=1,nwalk
       irnd(iw)=ran1(ifirst)*irandom
       ifirst=-irnd(iw)
      enddo

!  We use vx1 variables because vinit needs one-dimensional array.
       do iw=1,nwalk
        ifirst=irnd(iw)
        call vinit(TEMP,am,vx1,vy1,vz1,110,ifirst)

        do iat=1,natom
         vx(iat,iw)=vx1(iat)
         vy(iat,iw)=vy1(iat)
         vz(iat,iw)=vz1(iat)
        enddo
       enddo

        write(*,*)'Removing center of mass velocity'


! TODO: odstraneni rotace
!       irest endif    
       endif

!Nulovani rychlosti pro constrainovane atomy
       if(conatom.gt.0)then
       write(*,*)'Removing initial velocity of constrained atoms.'
        do iw=1,nwalk
         do iat=1,conatom
          vx(iat,iw)=0.0
          vy(iat,iw)=0.0
          vz(iat,iw)=0.0
         enddo
        enddo
       endif

! Scaling of velocities to correct temperature after COM velocity removal           
! or if restarting to different temperature
       ekin_mom=0.0d0
       do iw=1,nwalk
       do iat=1,natom
        pom=vx(iat,iw)**2+vy(iat,iw)**2+vz(iat,iw)**2
        pom=0.5*pom*am(iat)
        ekin_mom=ekin_mom+pom
       enddo
       enddo
       if(ekin_mom.gt.0.1d-14)then
        temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
       else
        temp_mom=0.0d0
       endif
       write(*,*)'Initial temperature [K]:',temp_mom*autok

! TODO: pro normal modes nemusi nutne fungovat!
       if(scaleveloc.eq.1.and.temp_mom.gt.0.1e-10)then
        ekin_mom=0.0d0
        write(*,*)'Scaling velocities to correct temperature.'
        scal=dsqrt(temp/temp_mom)
        do iw=1,nwalk
         do iat=1,natom
         vx(iat,iw) = vx(iat,iw)*scal 
         vy(iat,iw) = vy(iat,iw)*scal
         vz(iat,iw) = vz(iat,iw)*scal
         pom=vx(iat,iw)**2+vy(iat,iw)**2+vz(iat,iw)**2
         pom=0.5*pom*am(iat)
         ekin_mom=ekin_mom+pom
        enddo
       enddo
       temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
       write(*,*)'Temperature after scaling [K]:',temp_mom*autok
       endif


       
!----some stuff for spherical boundary onditions
      if(isbc.eq.1)then
       call sbc_init(x,y,z)
      endif

      if(iqmmm.eq.1.or.pot.eq.'mm')then
       if(qmmmtype.ne.'nab')then
        call inames_init()
        call ABr_init()
       endif
       write(*,nml=qmmm)
      endif

      if (ihess.eq.1)then 
       allocate ( hess(natom*3,natom*3,nwalk) )
       if (pot.eq.'nab')then
!!$OMP PARALLEL
        allocate ( h(natom*natom*9) )
!!$OMP END PARALLEL
       endif
      endif

      if (pot.eq."nab".or.qmmmtype.eq."nab")then
       if (alpha_pme.lt.0) alpha_pme = pi/cutoff
       if (kappa_pme.lt.0) kappa_pme = alpha_pme
       call nab_init(alpha_pme,cutoff,nsnb,ipbc,ips,iqmmm) !C function...see nabinit.c

       if(ipbc.eq.1)then
        if (nchain.eq.1) f=3
        call nab_getbox(boxx,boxy,boxz) !see nabinit.c
        call nab_getcharges(charges) !see nabinit.c

        do iat=1,natom
         charges(iat)=charges(iat)*ambtoau*sqrt(167100.75)  !for macsimus units
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


      end



      subroutine sh_init(nacx_old,nacy_old,nacz_old,vx_old,vy_old,vz_old,en_array_old,dt)
      use mod_array_size
      use mod_general,ONLY:natom,irandom,irest
      use mod_sh
      implicit none
      real*8 nacx_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
      real*8 en_array_old(nstmax,ntrajmax),dt,pop,pop2
      integer :: itrj,ist1,ist2,iat
      real*8  :: ran_test(10),pom=0.0d0,maxosc=0.0d0

     deltaE=deltaE/27.2114
     dtp=dt/substep
!-----random number initialization for surface hopping
!      call gautrg(ran_test,0,IRandom,6)
!      call vranf(ran_test,ntraj,0,6)
       call vranf(ran_test,0,irandom,6)

      do itrj=1,ntraj
       do ist1=1,nstate
        do ist2=1,nstate
         tocalc(ist1,ist2)=1
         dotproduct_old(ist1,ist2,itrj)=0.0d0
         do iat=1,natom
          nacx(iat,itrj,ist1,ist2)=0.0d0
          nacy(iat,itrj,ist1,ist2)=0.0d0
          nacz(iat,itrj,ist1,ist2)=0.0d0
         enddo
        enddo
       enddo   
      enddo

      if(inac.eq.2)then  ! for Adiab. dynamics
       do ist1=1,nstate
        do ist2=1,nstate
         tocalc(ist1,ist2)=0
        enddo
       enddo
      endif

     if(irest.ne.1)then

!automatic determination of initial state based on osc. streght
! NOT TESTED yet
        if(istate_init.eq.-1)then
         open(151,file='oscil.dat')
        endif
       do itrj=1,ntraj

        if(istate_init.eq.-1)then
         do ist1=1,nstate
          read(151,*)pom
          if(pom.gt.maxosc)then
           istate(itrj)=ist1
           maxosc=pom
          endif
         enddo
        else
         istate(itrj)=istate_init
        endif

        do ist1=1,nstate
         cel_re(ist1,itrj)=0.0d0
         cel_im(ist1,itrj)=0.0d0
        enddo
        cel_re(istate(itrj),itrj)=1.0d0
       enddo

        if(istate_init.eq.-1)then
         close(151)
        endif
!irest endif
      endif

      if(popthr.gt.0)then  
!--COMPUTE NACs only if population of the states is gt.popthr
      do ist1=1,nstate-1
       pop=cel_re(ist1,itrj)**2+cel_im(ist1,itrj)**2
       do ist2=ist1+1,nstate
        pop2=cel_re(ist2,itrj)**2+cel_im(ist2,itrj)**2
        if(pop.lt.popthr.and.pop2.lt.popthr.and.ist1.ne.istate(itrj).and.ist2.ne.istate(itrj)) tocalc(ist1,ist2)=0
       enddo
      enddo
      endif

      do itrj=1,ntraj
       do ist1=1,nstate
        en_array_old(ist1,itrj)=0.0d0

        do ist2=1,nstate
         do iat=1,natom
          nacx_old(iat,itrj,ist1,ist2)=0.0d0
          nacy_old(iat,itrj,ist1,ist2)=0.0d0
          nacz_old(iat,itrj,ist1,ist2)=0.0d0
          vx_old(iat,itrj)=0.0d0
          vy_old(iat,itrj)=0.0d0
          vz_old(iat,itrj)=0.0d0
         enddo
        enddo
       enddo   
      enddo

   end

   subroutine finish()
   use mod_general
   use mod_nhc
   use mod_estimators, only: hess,h
   use mod_fftw3
   implicit none
!   integer :: iter=-3
!   CHARACTER(len=150) :: cleanup
!   CHARACTER(len=20)  :: user,chpid
   close(1)
   close(2)
   if(ipimd.eq.1.or.icv.eq.1)then
    close(7)
   endif
   if(ipimd.eq.2)then
    close(3)
    close(4)
    close(5)
   endif
   if(nwritev.gt.0) close(13)
!--------------CLEANING-------------------------
   if (ihess.eq.1) deallocate ( hess )
   if (ihess.eq.1.and.pot.eq.'nab') deallocate ( h )
   if (istage.eq.2) call fftw_end()

   if(inose.eq.1)then
    deallocate( w )
    deallocate( Qm )
    deallocate( ms )
    if (imasst.eq.1)then
      deallocate( pnhx )
      deallocate( pnhy )
      deallocate( pnhz )
      deallocate( xi_x )
      deallocate( xi_y )
      deallocate( xi_z )
     else
      deallocate( pnhx )
      deallocate( xi_x )
    endif
   endif
!TODO dealokovat pole v NABU ...tj zavolet mme rutinu s iter=-3 nebo tak neco
!   if(pot.eq.'nab') call mme(NULL,NULL,iter)

!--following is needed for cleaning up after ab initio jobs      
!      if(pot.ne.'nab'.and.pot.ne.'harm'.and.pot.ne.'2dho'.and.pot.ne.'morse'.and.pot.ne.'guillot'.and.pot.ne.'nab')then
!       CALL get_environment_variable("USER", user)
!       write (chpid,*)pid
!       cleanup='rm -r /scratch/'//trim(user)//'/ABIN_SCRATCH_'//adjustl(chpid)
!       call system(trim(cleanup))
!      endif

   end

