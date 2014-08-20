! modules and routines for Surface hopping dynamics
! D. Hollas, O. Svoboda, P. Slavíček, M. Ončák

module mod_sh
   use mod_const, only: DP
   use mod_array_size, only: nstmax, ntrajmax
   use mod_utils, only: abinerror, printf
   implicit none
   private
   public :: istate_init,substep,deltae,integ,inac,nohop,alpha,popthr,nac_accu1,nac_accu2 
   public :: surfacehop, sh_init, istate, ntraj, nstate, cel_re, cel_im, tocalc, en_array
   public :: move_vars, get_nacm, write_nacmrest, read_nacmrest
   public :: energydifthr, energydriftthr, popsumthr, phase, adjmom, revmom

   integer,parameter :: ntraj = ntrajmax
   integer   :: istate_init=1, nstate=1, substep=100
   integer   :: inac=0, nohop=0, phase=0, adjmom=0, revmom=0
   integer   :: nac_accu1=7, nac_accu2=5 !7 is MOLPRO default
   real(DP)  :: dtp, alpha=0.1d0
   real(DP)  :: deltae=5.0d0, popthr=0.001d0
   real(DP)  :: popsumthr=0.001d0, energydifthr=1.0d0, energydriftthr=1.0d0 !eV
   character(len=10)    :: integ='butcher'
   real(DP),allocatable :: nacx(:,:,:,:), nacy(:,:,:,:), nacz(:,:,:,:)
   real(DP),allocatable :: nacx_old(:,:,:,:), nacy_old(:,:,:,:), nacz_old(:,:,:,:)
   real(DP),allocatable :: dotproduct(:,:,:), dotproduct_old(:,:,:) !for inac=1
   real(DP),allocatable :: en_array(:,:), en_array_old(:,:)
   real(DP)  :: cel_re(nstmax,ntrajmax), cel_im(nstmax,ntrajmax)
   real(DP)  :: eshift, entot0, gama(nstmax,nstmax,ntrajmax)
   integer, allocatable :: tocalc(:,:)
   integer   :: istate(ntrajmax)
   save
   contains

   subroutine sh_init(x, y, z, vx, vy, vz, dt)
   use mod_const,    only: AUtoEV
   use mod_general,  only: irest, natom, it
   use mod_forces,   only: force_clas
   use mod_kinetic,  only: ekin_v
   real(DP),intent(inout)  :: x(:,:), y(:,:), z(:,:)
   real(DP),intent(in)     :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: dt
   real(DP)  :: dum_fx(size(x,1),size(x,2))
   real(DP)  :: dum_fy(size(y,1),size(y,2))
   real(DP)  :: dum_fz(size(z,1),size(z,2))
   real(DP)  :: dum_eclas
   integer   :: itrj,ist1
   real(DP)  :: pom=0.0d0,maxosc=0.0d0

   deltaE=deltaE/AUtoEV
   dtp=dt/substep
   gama=0.0d0

   allocate( nacx(natom, ntrajmax, nstate, nstate) )
   nacx=0.0d0
   !automatic allocation and data copying
   nacy=nacx
   nacz=nacx
   nacx_old=nacx
   nacy_old=nacx
   nacz_old=nacx
   if (inac.eq.1)then
      allocate( dotproduct_old(nstate, nstate, ntrajmax) )
      dotproduct_old=0.0d0
      dotproduct=dotproduct_old
      adjmom=1
   end if

   allocate( en_array(nstate, ntraj) )
   en_array=0.0d0
   en_array_old=en_array

!--Determining the initial WF coefficients
   if(irest.ne.1)then

!--automatic determination of initial state based on osc. streght
!  NOT TESTED yet
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

!--computing only energies, used for subsequent determination of tocalc according to deltae
   allocate( tocalc(nstate, nstate) )
   tocalc=0
   dum_eclas=0.0d0
   dum_fx=0.0d0   ; dum_fy=0.0d0    ;dum_fz=0.0d0

   !restarting the SH, reading NACM

   call force_clas(dum_fx,dum_fy,dum_fz,x,y,z,dum_eclas)
   !WARNING: itrj hack
   itrj=1
   eshift = -en_array(1, itrj)
   entot0 = en_array(istate(itrj), itrj)+ekin_v(vx, vy, vz)
   entot0 = entot0 * AUtoEV

   do itrj=1,ntraj
      call set_tocalc(itrj)
   end do
   if(it.ne.0.and.irest.eq.1) call read_nacmrest()

   end subroutine sh_init

   subroutine get_nacm(itrj)
   integer, intent(in) :: itrj
   integer :: iost
   if(inac.eq.0)then
      iost=readnacm(itrj)
      if(iost.ne.0.and.nac_accu1.gt.nac_accu2)then
!----------if NACME NOT COMPUTED: TRY TO DECREASE ACCURACY--------------
         call calcnacm(itrj)
         iost=readnacm(itrj)
      endif
      if(iost.ne.0)then
         write(*,*)'Some NACMEs not read. Exiting...'
         call abinerror('main program')
      endif
      !we always have to set tocalc, as we change it in readnacm
      call set_tocalc(itrj)
   endif
   end subroutine Get_Nacm

   subroutine set_tocalc(itrj)
   integer,intent(in) :: itrj 
   integer   :: ist1,ist2
   real(DP)  :: pop,pop2
   
   do ist1=1,nstate-1
      do ist2=ist1+1,nstate
         if(abs(en_array(ist1,itrj)-en_array(ist2,itrj)).lt.deltae) then
            tocalc(ist1,ist2)=1 
         else
            tocalc(ist1,ist2)=0 
         endif
      enddo
   enddo
   
   if(inac.eq.2)then  ! for ADIABATIC dynamics
      do ist1=1,nstate-1
         do ist2=ist1+1,nstate
         tocalc(ist1,ist2)=0
         enddo
      enddo
   endif
   
   if(popthr.gt.0)then  
      !COMPUTE NACM only if population of the states is gt.popthr
      do ist1=1,nstate-1
         pop=cel_re(ist1,itrj)**2+cel_im(ist1,itrj)**2
         do ist2=ist1+1,nstate
            pop2=cel_re(ist2,itrj)**2+cel_im(ist2,itrj)**2
            if(pop.lt.popthr.and.pop2.lt.popthr.and.ist1.ne.istate(itrj).and.ist2.ne.istate(itrj)) tocalc(ist1,ist2)=0
         enddo
      enddo
   endif
   end subroutine set_tocalc

   subroutine Write_nacmrest()
   use mod_qmmm, only:natqm
   integer :: ist1,ist2,iat,itrj
   integer :: iunit1, iunit2
   open(newunit=iunit1,file='nacmrest.dat',action="write")
   if (phase.eq.1) open(newunit=iunit2,file='phaserest.dat',action="write")
   do itrj=1, ntraj

   do ist1=1,nstate-1

      do ist2=ist1+1,nstate
   
         if(tocalc(ist1,ist2).eq.1)then
   
            write(iunit1,*)'NACM between states',ist1,ist2
            do iat=1,natqm              ! reading only for QM atoms
               write(iunit1,*)nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
            enddo
   
   !--------if tocalc 
         endif

   
      enddo
   enddo

   if(phase.eq.1)then
      do ist1=1,nstate
         write(iunit2,*)(gama(ist1,ist2,itrj),ist2=1,nstate)
      end do
   end if

   end do

   close(iunit1)
   if(phase.eq.1) close(iunit2)

   end subroutine write_nacmrest

   subroutine read_nacmrest()
   use mod_general, only: it
   use mod_qmmm, only:natqm
   integer :: iost,ist1,ist2,iat,itrj
   integer :: iunit1, iunit2
   character(len=200) :: chmsg
   character(len=20) :: chit
   character(len=60) :: chrestart

   write(*,*)'Reading NACM from nacmrest.dat'
   open(newunit=iunit1,file='nacmrest.dat',action="read")
   if (phase.eq.1)then
      write(*,*)'Reading phase from phaserest.dat'
      open(newunit=iunit2,file='phaserest.dat',action="read")
   end if

   do itrj=1, ntraj

   do ist1=1,nstate-1

      do ist2=ist1+1,nstate
   
         if(tocalc(ist1,ist2).eq.1)then
   
            read(iunit1,*, iostat=iost)
            do iat=1,natqm              ! reading only for QM atoms
               read(iunit1,*,iomsg=chmsg)nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
               if (iost.ne.0)then
                  write(*,*)'Error reading NACM from file nacmrest.'
                  write(*,*)chmsg
                  call abinerror('write_nacm_rest')
               end if
               nacx(iat,itrj,ist2,ist1)=-nacx(iat,itrj,ist1,ist2)
               nacy(iat,itrj,ist2,ist1)=-nacy(iat,itrj,ist1,ist2)
               nacz(iat,itrj,ist2,ist1)=-nacz(iat,itrj,ist1,ist2)
            enddo
   
   !--------if tocalc 
         endif
   
      enddo
   enddo

   if(phase.eq.1)then
      do ist1=1,nstate
         read(iunit2,*)(gama(ist1,ist2,itrj),ist2=1,nstate)
      end do
   end if

   end do

   close(iunit1)

   write (chit,*)it
   chrestart='cp nacmrest.dat nacmrest.dat.'//adjustl(chit)
   write(*,*)'Making backup of the nacmrest.dat'
   write(*,*)chrestart
   call system(chrestart)  

   if (phase.eq.1)then
      close(iunit2)
      chrestart='cp phaserest.dat phaserest.dat.'//adjustl(chit)
      write(*,*)'Making backup of the phaserest.dat'
      write(*,*)chrestart
      call system(chrestart)  
   end if

   end subroutine read_nacmrest


   integer function readnacm(itrj)
   use mod_qmmm,only:natqm
   integer :: iost,ist1,ist2,iat,itrj
   iost=0  ! needed if each tocalc=0
   open(127,file='nacm.dat')
   do ist1=1,nstate-1
      do ist2=ist1+1,nstate
   
         if(tocalc(ist1,ist2).eq.1)then
   
            do iat=1,natqm              ! reading only for QM atoms
               read(127,*,IOSTAT=iost)nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
               if(iost.eq.0)then
                  tocalc(ist1,ist2)=0   !marking as read, useful if we do decreased accuracy
                  nacx(iat,itrj,ist2,ist1)=-nacx(iat,itrj,ist1,ist2)
                  nacy(iat,itrj,ist2,ist1)=-nacy(iat,itrj,ist1,ist2)
                  nacz(iat,itrj,ist2,ist1)=-nacz(iat,itrj,ist1,ist2)
               else
                  close(127,status='delete')
                  write(*,*)'WARNING:NACM between states',ist1,ist2,'not read.'
                  readnacm=iost
                  return
               endif
            enddo
   
   !--------if tocalc 
         endif
   
      enddo
   enddo

   close(127,status='delete')
   readnacm=iost
   return
   end function readnacm

   subroutine calcnacm(itrj)
   use mod_general, only: it, pot
   integer, intent(in) :: itrj
   integer :: ist1,ist2
   character(len=100) :: chsystem
   open(unit=510,file='state.dat')
   write(510,'(I2)')istate(itrj)
   write(510,'(I2)')nstate
   ! tocalc is upper triangular part of a matrix without diagonal
   ! tocalc(,)=1 -> calculate NACME
   ! tocalc(,)=0 -> don't calculate NACME
   do ist1=1,nstate-1
      do ist2=ist1+1,nstate
         write(510,'(I1,A1)',advance='no')tocalc(ist1,ist2),' ' 
      enddo
   enddo
   close(510) 

   if(pot.eq.'molpro')then
      write(*,*)'WARNING: Some NACME not computed.Trying with decreased accuracy.'
      write(*,*)'Calling script r.molpro with accuracy:',nac_accu2
      write(chsystem,'(A20,I13,I4.3,I3,A12)')'./MOLPRO/r.molpro ',it,itrj,nac_accu2,' < state.dat'
   else
      write(*,*)'Different accuracy for NACME is currently supported only by MOLPRO.'
      call abinerror('calcnacm')
   endif

   call system(chsystem)
   end subroutine calcnacm

   subroutine move_vars(vx,vy,vz,vx_old,vy_old,vz_old,itrj)
   use mod_general,only:natom
   real(DP),intent(in)   :: vx(:,:),vy(:,:),vz(:,:)
   real(DP),intent(out)  :: vx_old(:,:),vy_old(:,:),vz_old(:,:)
   integer, intent(in) :: itrj
   integer :: ist1,ist2,iat
   !---moving new to old variables
   do ist1=1,nstate
      en_array_old(ist1,itrj)=en_array(ist1,itrj)
   
      if(inac.eq.0)then   ! nedelame, pokud nacitame rovnou dotprodukt
         do ist2=1,nstate
            do iat=1,natom
               nacx_old(iat,itrj,ist1,ist2)=nacx(iat,itrj,ist1,ist2)
               nacy_old(iat,itrj,ist1,ist2)=nacy(iat,itrj,ist1,ist2)
               nacz_old(iat,itrj,ist1,ist2)=nacz(iat,itrj,ist1,ist2)
            enddo
         enddo
      endif
   
   enddo
   
   do iat=1,natom
      vx_old(iat,itrj)=vx(iat,itrj)
      vy_old(iat,itrj)=vy(iat,itrj)
      vz_old(iat,itrj)=vz(iat,itrj)
   enddo

   if(inac.eq.1)then
      do ist1=1,nstate
         do ist2=1,nstate
            dotproduct_old(ist1,ist2,itrj)=dotproduct(ist1,ist2,itrj)
         enddo
      enddo
   endif

   end subroutine move_vars
   
   subroutine check_popsum(itrj,popsumin)
      real(DP),intent(in),optional :: popsumin
      integer,intent(in)           :: itrj
      real(DP)                     :: popsum
      integer                      :: ist1

      popsum=0.0d0
      if(present(popsumin))then
         popsum=popsumin
      else
         do ist1=1,nstate 
            popsum=popsum+cel_re(ist1,itrj)**2+cel_im(ist1,itrj)**2
         end do
      end if

      if (abs(popsum-1.0d0).gt.popsumthr)then
         write(*,*)'ERROR:Sum of populations=',popsum, &
          'which differs from 1.0 by more than popsumthr=',popsumthr
         write(*,*)'Increase number of substeps or use more accurate integrator.'
         call abinerror('surfacehop')
      end if
   end subroutine check_popsum


   subroutine surfacehop(x,y,z,vx,vy,vz,vx_old,vy_old,vz_old,dt,eclas)
      use mod_const,   only: ANG, AUTOFS
      use mod_general, only: natom, nwrite, idebug, it
      use mod_system,  ONLY: names
      use mod_qmmm,    ONLY: natqm
      use mod_random,  ONLY: vranf
      use mod_kinetic, ONLY: ekin_v
      real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(inout) :: vx_old(:,:),vy_old(:,:),vz_old(:,:)
      real(DP),intent(inout) :: eclas
      real(DP),intent(in)    :: dt
      real(DP)  :: vx_int(size(vx,1),size(vx,2))
      real(DP)  :: vy_int(size(vy,1),size(vy,2))
      real(DP)  :: vz_int(size(vz,1),size(vz,2))
      real(DP)  :: vx_newint(size(vx,1),size(vx,2))
      real(DP)  :: vy_newint(size(vy,1),size(vy,2))
      real(DP)  :: vz_newint(size(vz,1),size(vz,2))
      real(DP)  :: en_array_int(nstmax,ntrajmax),en_array_newint(nstmax,ntrajmax)
      real(DP)  :: nacx_int(size(vx,1),ntrajmax,nstmax,nstmax)  
      real(DP)  :: nacy_int(size(vx,1),ntrajmax,nstmax,nstmax)
      real(DP)  :: nacz_int(size(vx,1),ntrajmax,nstmax,nstmax)
      real(DP)  :: nacx_newint(size(vx,1),ntrajmax,nstmax,nstmax)
      real(DP)  :: nacy_newint(size(vx,1),ntrajmax,nstmax,nstmax)
      real(DP)  :: nacz_newint(size(vx,1),ntrajmax,nstmax,nstmax)
      real(DP)  :: dotproduct_int(nstmax,nstmax,ntrajmax),dotproduct_newint(nstmax,nstmax,ntrajmax) !rename dotproduct_int
      real(DP)  :: t(nstmax,nstmax)           !switching probabilities
      real(DP)  :: t_tot(nstmax,nstmax)       !cumulative switching probabilities
      real(DP)  :: ran(10)
      real(DP)  :: pop(size(tocalc,1),ntrajmax),popsum !populations
      integer :: iat,ist1,ist2,itrj,itp     !iteration counters
      integer :: ist                    ! =istate(itrj)
      real(DP)  :: vect_olap,fr,frd
      real(DP)  :: ekin_mom,apom,edif,tau,fact,sum_norm
      integer :: ihop,ijunk
      real(DP)  :: pop0, a_re, a_im, prob(nstmax), cn, stepfs
      character(len=500) :: formt
      character(len=20) :: chist,chihop,chit
      integer :: iost, iunit=1


      do itrj=1,ntraj

         call check_energy(vx_old, vy_old, vz_old, vx, vy, vz, itrj)
         call check_energydrift(vx, vy, vz, itrj)

         t_tot=1.0d0

         popsum=0.0d0

!-------READING NACM-----------------------     
     if(inac.eq.0)then

      do ist1=1,nstate-1
         do ist2=ist1+1,nstate
            if(tocalc(ist1,ist2).eq.0)then
               write(*,*)'Not computing NACME between states',ist1,ist2
               do iat=1,natqm        ! MUSIME NULOVAT UZ TADY,bo pak menime tocalc behem cteni
                  nacx(iat,itrj,ist1,ist2)=0.0d0
                  nacy(iat,itrj,ist1,ist2)=0.0d0
                  nacz(iat,itrj,ist1,ist2)=0.0d0
                  nacx(iat,itrj,ist2,ist1)=0.0d0
                  nacy(iat,itrj,ist2,ist1)=0.0d0
                  nacz(iat,itrj,ist2,ist1)=0.0d0
               enddo
            endif
         enddo
      enddo

       iost=readnacm(itrj)
!------------if NACM NOT COMPUTED: TRY TO DECREASE ACCURACY--------------
       if(iost.ne.0.and.nac_accu1.gt.nac_accu2)then
       call calcnacm(itrj)

       iost=readnacm(itrj)
       endif
!------------if NACM STILL NOT COMPUTED: USE OLD NACM--------------
       if(iost.ne.0)then
        write(*,*)'ERROR:Some NACMEs not read.'
        call abinerror('surfacehop')
!        do ist1=1,nstate-1
!         do ist2=ist1+1,nstate
!
!         if(tocalc(ist1,ist2).eq.1)then !!po uspesnem precteni nulujeme tocalc
!          write(*,*)'Warning! NACM between states',ist1,'and',ist2,'not computed.'
!          write(*,*)'Using NACM from previous step.'
!          do iat=1,natom
!           nacx(iat,itrj,ist1,ist2)=nacx_old(iat,itrj,ist1,ist2)
!           nacy(iat,itrj,ist1,ist2)=nacy_old(iat,itrj,ist1,ist2)
!           nacz(iat,itrj,ist1,ist2)=nacz_old(iat,itrj,ist1,ist2)
!           nacx(iat,itrj,ist2,ist1)=-nacx(iat,itrj,ist1,ist2)
!           nacy(iat,itrj,ist2,ist1)=-nacy(iat,itrj,ist1,ist2)
!           nacz(iat,itrj,ist2,ist1)=-nacz(iat,itrj,ist1,ist2)
!          enddo
!         endif

!         enddo
!        enddo
       endif

!--------------calculating overlap between nacmes-------------------
      do ist1=1,nstate
       do ist2=1,nstate
        vect_olap=0.0d0
        do iat=1,natom
         vect_olap=vect_olap+nacx_old(iat,itrj,ist1,ist2)*nacx(iat,itrj,ist1,ist2)
         vect_olap=vect_olap+nacy_old(iat,itrj,ist1,ist2)*nacy(iat,itrj,ist1,ist2)
         vect_olap=vect_olap+nacz_old(iat,itrj,ist1,ist2)*nacz(iat,itrj,ist1,ist2)
        enddo
        if(vect_olap.lt.0)then
         do iat=1,natom
          nacx(iat,itrj,ist1,ist2)=-nacx(iat,itrj,ist1,ist2)
          nacy(iat,itrj,ist1,ist2)=-nacy(iat,itrj,ist1,ist2)
          nacz(iat,itrj,ist1,ist2)=-nacz(iat,itrj,ist1,ist2)
         enddo 
        endif

        enddo
       enddo

!----------- INAC=1  endif
     endif

!------READING time-derivative couplings----------------------------------
   if ( inac.eq.1 ) then

      open(100,file='tdcoups.dat')
      read(100,*)
      read(100,*)
      do ist1=1,nstate
         read(100,*)ijunk,(dotproduct(ist1,ist2,itrj),ist2=1,nstate)
      enddo
      close(100)

      do ist1=1,nstate
         do ist2=1,nstate
            dotproduct(ist1,ist2,itrj)=-dotproduct(ist1,ist2,itrj)/dt
            if(ist1.eq.ist2)  dotproduct(ist1,ist2,itrj)=0.0d0
         enddo
      enddo

      do ist1=1,nstate-1
         do ist2=ist1+1,nstate
            if(tocalc(ist1,ist2).eq.0)then
               write(*,*)'Not computing NACM(tocalc=0) for states',ist1,ist2
               dotproduct(ist1,ist2,itrj)=0.0d0
               dotproduct(ist2,ist1,itrj)=0.0d0
            endif
         enddo
      enddo


!v prvnim korku nemame mezi cim interpolovat..bereme nulty krok jako prvni
      if( it.eq.1) dotproduct_old=dotproduct
       !TDC end if
   endif
!------------------END-OF-TDC-------------------------------

! Smaller time step for electron population transfer
   do itp=1,substep      

      ist=istate(itrj)
      !pop0 is later used for tully's fewest switches
      pop0=cel_re(ist,itrj)**2+cel_im(ist,itrj)**2

!-----------------INTERPOLACE--------------------------------------------------


       fr=real(itp,DP)/real(substep,DP)
       frd=1.0d0-fr

       call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_newint,vy_newint,vz_newint, &
                      nacx_newint,nacy_newint,nacz_newint,en_array_newint, &
                      dotproduct_newint,fr,frd,itrj)

       if(inac.eq.1)then
        call interpolate_dot(dotproduct_newint,fr,frd,itrj)
       endif


       fr=real(itp-1,DP)/real(substep,DP)
       frd=1.0d0-fr

       call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_int,vy_int,vz_int, &
                      nacx_int,nacy_int,nacz_int,en_array_int, &
                      dotproduct_int,fr,frd,itrj)

       if(inac.eq.1)then
        call interpolate_dot(dotproduct_int,fr,frd,itrj)
       endif


      !------END-OF-INTERPOLATIONS-------------------------------------
       
      if(integ.eq.'rk4') call rk4step(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)

      if(integ.eq.'butcher') call butcherstep(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)

      if(integ.eq.'euler') call eulerstep(en_array_int, en_array_newint, dotproduct_int, itrj)

      call check_popsum(itrj)

!--calculation of switching probabilities
!- do not calculated the whole transition matrix
!- t array could be one dimensional but whatever
      do ist2=1,nstate 
        pop(ist2,itrj)=cel_re(ist2,itrj)**2+cel_im(ist2,itrj)**2
        a_re=( cel_re(ist,itrj)*cel_re(ist2,itrj)+cel_im(ist,itrj)*cel_im(ist2,itrj) )
        if (phase.eq.1)then
         a_im=(- cel_im(ist,itrj)*cel_re(ist2,itrj) + cel_re(ist, itrj)*cel_im(ist2, itrj) )
         t(ist,ist2)=a_re*cos(gama(ist, ist2, itrj))-sin(gama(ist, ist2, itrj))*a_im
         t(ist,ist2)=t(ist,ist2)*(dotproduct_int(ist,ist2,itrj) + dotproduct_newint(ist,ist2,itrj) )
        else
!         t(ist,ist2)=2*a_re*dotproduct_int(ist,ist2,itrj)
          t(ist,ist2)=a_re*(dotproduct_int(ist,ist2,itrj)+dotproduct_newint(ist,ist2,itrj) )
        endif
        t(ist,ist2)=t(ist,ist2)*dtp/(pop0+1d-20)
      enddo        

      if(idebug.gt.1)then
         stepfs=(it*substep+itp-substep)*dt*AUtoFS/substep
         write(formt,'(A7,I3,A7)')'(F15.2,',nstate,'E20.10)'
         write(16,fmt=formt)stepfs,(t(ist,ist1),ist1=1,nstate)

         write(17,fmt=formt,advance="no")stepfs,(cel_re(ist1,itrj),ist1=1,nstate)
         write(formt,'(A1,I3,A7)')'(',nstate,'E20.10)'
         write(17,fmt=formt,advance="no")(cel_im(ist1,itrj),ist1=1,nstate)
         write(17,*)''

         if(phase.eq.1)then
            write(18,'(F15.2,E20.10)',advance="no")stepfs,gama(2,1,itrj)
            do ist1=3,nstate
               write(formt,'(A1,I3,A7)')'(',ist1-1,'E20.10)'
               write(18,fmt=formt,advance="no")(gama(ist1,ist2,itrj),ist2=1,ist1-1)
            end do
            write(18,*)''
         end if
      end if

      do ist2=1,nstate 
         if(t(ist,ist2).gt.1.0_DP)then
            write(*,*)'ERROR: Hopping probability greater than 1.'
            call abinerror('surfacehop')
         end if
         if(t(ist,ist2).lt.0.0d0)  t(ist,ist2)=0.0d0
         !cumulative probability over whole big step
         t_tot(ist,ist2)=t_tot(ist,ist2)*(1-t(ist,ist2))
      enddo


      prob=0.0d0

! Auxiliary calculations of probabilities
      do ist1=1,nstate
       if(ist1.eq.ist.and.ist1.eq.1)then
        prob(ist1)=0.0d0
       endif
       if(ist1.ne.ist.and.ist1.eq.1)then
        prob(ist1)=t(ist,ist1)
       endif
       if(ist1.ne.ist.and.ist1.ne.1)then
        prob(ist1)=prob(ist1-1)+t(ist,ist1)
       endif
       if(ist1.eq.ist.and.ist1.ne.1)then
        prob(ist1)=prob(ist1-1)
       endif
      enddo

      !should we hop before decoherence?????
      ! every article says something different
      ! newton-X hops before decoherence
! HOPPING      
   if (nohop.ne.1)then

      ihop=0
      call vranf(ran,1,0,6)
      cn=ran(1)
      if(ist.ne.1)then
       if(cn.ge.0.0d0.and.cn.lt.prob(1))then
        ihop=1
       endif
      endif

      do ist1=2,nstate
      if(ist.ne.ist1)then
       if(cn.ge.prob(ist1-1).and.cn.lt.prob(ist1))then 
        ihop=ist1
       endif
      endif
      enddo

!-------does HOP occured???-----------------
      if(ihop.ne.0)then
         if(adjmom.eq.0) call hop(vx,vy,vz,ist,ihop,itrj,eclas)
         if(adjmom.eq.1) call hop_dot(vx,vy,vz,ist,ihop,itrj,eclas)
         write(formt,'(A8,I3,A7)')'(A1,I10,',nstate+1,'E20.10)'
         write(3,*)'#Substep   RandomNum   Probabilities'
         write(3,fmt=formt)'#',itp,cn,(t(ist,ist1),ist1=1,nstate)
         write(formt,'(A5,I10,A1,I2,A1,I2)')'geom.',it,'.',ist,'.',ihop
         write(chist,*)ist
         write(chihop,*)ihop
         write(chit,*)it
         formt='geom.'//trim(adjustl(chist))//'.'//trim(adjustl(chihop))//'.'//adjustl(chit)
         open(100,file=trim(formt))
         write(100,*)natom
         write(100,*)''
         do iat=1,natom
            write(100,*)names(iat),x(iat,itrj)/ang,y(iat,itrj)/ang,z(iat,itrj)/ang
         enddo
         close(100)
      endif   

      !nohop endif
   endif

      if(alpha.gt.0)then
! Quantum decoherence part----------------------------------
      
      ekin_mom=ekin_v(vx_int,vy_int,vz_int)

      if(ekin_mom.gt.1.0d-4)then !Why this number?In almost all cases it won't matter

      do ist1=1,nstate
       if(ist1.ne.istate(itrj)) then

!Calculation of exponential factor               
       edif=abs( en_array_int(ist1,itrj)-en_array_int(istate(itrj),itrj) )
       tau=1/edif
       apom=alpha/ekin_mom
       tau=tau*(1.0d0+apom)
       fact=dexp(-dtp/tau)

       cel_re(ist1,itrj)=cel_re(ist1,itrj)*fact
       cel_im(ist1,itrj)=cel_im(ist1,itrj)*fact
       endif
      enddo

! RENORMALIZATION OF ISTATE     
      sum_norm=1.0d0
      do ist1=1,nstate
         if(ist1.ne.istate(itrj)) then
            sum_norm=sum_norm-cel_re(ist1,itrj)**2-cel_im(ist1,itrj)**2
         endif
      enddo
      fact=sum_norm/(cel_re(istate(itrj),itrj)**2+cel_im(istate(itrj),itrj)**2+1.0d-7) !TODO:smaller number
      !Following should never happen as we check for popsumthr later in this subroutine
      if(fact.lt.0.0d0)then
         write(*,*)'Fatal error in surfacehop during decoherence renormalization.'
         write(*,*)'fact=',fact,'but should be > 0'
         write(*,*)'This usually means inaccurate integration of electronic SE.'
         write(*,*)'Increase number of substeps or use more accurate integrator.'
         call abinerror('surfacehop')
      end if
      fact=sqrt(fact)
!     write(*,*)'renomr',fact,istate(itrj)

      cel_re(istate(itrj),itrj)=cel_re(istate(itrj),itrj)*fact
      cel_im(istate(itrj),itrj)=cel_im(istate(itrj),itrj)*fact

      !ekin endif
      endif

!--------END-OF-DECOHERENCE------------------------      
      endif




!itp loop
      enddo

!    set tocalc array for the next step
     call set_tocalc(itrj)

     if(idebug.eq.1)then
       if(it.eq.1)then
         open(newunit=iunit,file='dotprodmatrix.dat')
       else
         open(newunit=iunit,file='dotprodmatrix.dat',access='append')
       endif
      write(iunit,*)'Step: ',it
      do ist1=1,nstate
       write(iunit,*)(dotproduct_int(ist1,ist2,itrj),ist2=1,nstate)
      enddo
      close(iunit)
     endif


      do ist1=1,nstate
       popsum=popsum+pop(ist1,itrj)
      enddo

     call check_popsum(itrj,popsum)
     call move_vars(vx,vy,vz,vx_old,vy_old,vz_old,itrj)

   if(modulo(it,nwrite).eq.0)then
      stepfs=it*dt*AUtoFS
      write(formt,'(A10,I3,A13)')'(F15.2,I3,',nstate,'F10.5,1F10.7)'
      write(3,fmt=formt)stepfs,istate(itrj),(pop(ist1,itrj), ist1=1,nstate),popsum

      t_tot=1-t_tot !up to know, t_tot is the probability of not hopping
      write(formt,'(A10,I3,A6)')'(F15.2,I3,',nstate,'F10.5)'
      write(4,fmt=formt)stepfs,istate(itrj),(t_tot(ist,ist1),ist1=1,nstate)
      write(formt,'(A7,I3,A7)')'(F15.2,',nstate,'E20.10)'
      write(8,fmt=formt)stepfs,(en_array(ist1,itrj),ist1=1,nstate)
      if (inac.eq.0)  write(14,*)'Time step:',it
      do ist1=1,nstate-1
         do ist2=ist1+1,nstate

            if (ist1.eq.1.and.ist2.eq.2)then
               write(15,'(F15.2,E20.10)',advance='no')stepfs,dotproduct_int(ist1,ist2,itrj)
            else
               write(15,'(E20.10)',advance='no')dotproduct_int(ist1,ist2,itrj)
            end if

            if (inac.eq.0)then
               write(14,*)'NACME between states:',ist1,ist2
               do iat=1,natom
                  write(14,'(3E20.10)')nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
               enddo
            endif

         end do
      end do
      write(15,*)''

   endif

   ! ntraj enddo       
   enddo

   end subroutine surfacehop

   subroutine hop(vx,vy,vz,instate,outstate,itrj,eclas)
      use mod_general, ONLY:natom
      use mod_system, ONLY: am
      use mod_kinetic, only: ekin_v
      real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(inout) :: eclas
      integer, intent(in)    :: itrj,instate,outstate
      real(DP)  :: a_temp,b_temp,c_temp,g_temp
      real(DP)  :: ekin, ekin_new
      integer   :: iat

!  Checking for frustrated hop

      a_temp=0.d0
      b_temp=0.d0

      ekin=ekin_v(vx,vy,vz)

      do iat=1,natom
         a_temp=a_temp+nacx(iat,itrj,instate,outstate)**2/am(iat)
         a_temp=a_temp+nacy(iat,itrj,instate,outstate)**2/am(iat)
         a_temp=a_temp+nacz(iat,itrj,instate,outstate)**2/am(iat)
         b_temp=b_temp+nacx(iat,itrj,instate,outstate)*vx(iat,itrj)
         b_temp=b_temp+nacy(iat,itrj,instate,outstate)*vy(iat,itrj)
         b_temp=b_temp+nacz(iat,itrj,instate,outstate)*vz(iat,itrj)
      enddo
      a_temp=0.5d0*a_temp
      c_temp=b_temp**2+4*a_temp*(en_array(instate,itrj)-en_array(outstate,itrj))

      if(c_temp.lt.0)then
         write(3,*)'#Not enough momentum in the direction of NAC vector.'
         ! Try, whether there is enough total kinetic energy and scale velocities.
         call hop_dot(vx,vy,vz,instate,outstate,itrj,eclas)
         return
      endif

      istate(itrj)=outstate
      eclas=en_array(outstate,itrj)
      write(3,'(A24,I3,A10,I3)')'#Hop occured from state ',instate,' to state ',outstate


!------------- Rescaling the velocities------------------------

      if(b_temp.lt.0) then
         g_temp=( b_temp + sqrt(c_temp) )/2.0d0/a_temp
      endif

      if(b_temp.ge.0) then
         g_temp=(b_temp-sqrt(c_temp))/2.0d0/a_temp
      endif


      do iat=1,natom
         vx(iat,itrj)=vx(iat,itrj)-g_temp*nacx(iat,itrj,instate,outstate)/am(iat)
         vy(iat,itrj)=vy(iat,itrj)-g_temp*nacy(iat,itrj,instate,outstate)/am(iat)
         vz(iat,itrj)=vz(iat,itrj)-g_temp*nacz(iat,itrj,instate,outstate)/am(iat)
      enddo
      ekin_new=ekin_v(vx,vy,vz)

      write(3,'(A31,2E20.10)')'#deltaE_pot     E_kin-total',en_array(outstate,itrj)-en_array(instate,itrj),ekin
      write(3,'(A,2E20.10)')'#Total_energy_old   Total_energy_new :',ekin+en_array(instate,itrj),ekin_new+en_array(outstate,itrj)

   end subroutine hop

   subroutine interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_int,vy_int,vz_int, &
                      nacx_int,nacy_int,nacz_int,en_array_int, &
                      dotproduct_int,fr,frd,itrj)
      use mod_general, only: natom
      real(DP),intent(out) :: dotproduct_int(:,:,:)
      real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(in) :: vx_old(:,:),vy_old(:,:),vz_old(:,:)
      real(DP),intent(out) :: vx_int(:,:),vy_int(:,:),vz_int(:,:)
      real(DP),intent(out) :: nacx_int(:,:,:,:)
      real(DP),intent(out) :: nacy_int(:,:,:,:)
      real(DP),intent(out) :: nacz_int(:,:,:,:)
      real(DP),intent(out) :: en_array_int(:,:)
      real(DP) :: fr,frd
      integer  :: iat,ist1,ist2,itrj     !iteration counters
      
      do ist1=1,nstate

        en_array_int(ist1,itrj)=en_array(ist1,itrj)*fr+en_array_old(ist1,itrj)*frd
        do ist2=1,nstate
         dotproduct_int(ist1,ist2,itrj)=0.0d0
         do iat=1,natom
         nacx_int(iat,itrj,ist1,ist2)=nacx(iat,itrj,ist1,ist2)*fr+nacx_old(iat,itrj,ist1,ist2)*frd
         nacy_int(iat,itrj,ist1,ist2)=nacy(iat,itrj,ist1,ist2)*fr+nacy_old(iat,itrj,ist1,ist2)*frd
         nacz_int(iat,itrj,ist1,ist2)=nacz(iat,itrj,ist1,ist2)*fr+nacz_old(iat,itrj,ist1,ist2)*frd
         vx_int(iat,itrj)=vx(iat,itrj)*fr+vx_old(iat,itrj)*frd
         vy_int(iat,itrj)=vy(iat,itrj)*fr+vy_old(iat,itrj)*frd
         vz_int(iat,itrj)=vz(iat,itrj)*fr+vz_old(iat,itrj)*frd
         dotproduct_int(ist1,ist2,itrj)=dotproduct_int(ist1,ist2,itrj)+vx_int(iat,itrj)*nacx_int(iat,itrj,ist1,ist2)
         dotproduct_int(ist1,ist2,itrj)=dotproduct_int(ist1,ist2,itrj)+vy_int(iat,itrj)*nacy_int(iat,itrj,ist1,ist2)
         dotproduct_int(ist1,ist2,itrj)=dotproduct_int(ist1,ist2,itrj)+vz_int(iat,itrj)*nacz_int(iat,itrj,ist1,ist2)
         enddo
        enddo
       enddo

    end subroutine interpolate


   subroutine hop_dot(vx,vy,vz,instate,outstate,itrj,eclas)
      use mod_general, ONLY:natom,idebug
      use mod_kinetic ,ONLY: ekin_v
      real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(inout) :: eclas
      integer,intent(in)     :: itrj,instate,outstate
      integer   :: iat
      real(DP)  :: de,ekin,alfa,ekin_new

      ekin=0.0d0
      ekin_new=0.0d0

      dE=en_array(outstate,itrj)-en_array(instate,itrj)
      ekin=ekin_v(vx,vy,vz)

      if(ekin.ge.de)then

        alfa=sqrt(1-de/ekin)

        do iat=1,natom
         vx(iat,itrj)=alfa*vx(iat,itrj) 
         vy(iat,itrj)=alfa*vy(iat,itrj) 
         vz(iat,itrj)=alfa*vz(iat,itrj) 
        enddo
        istate(itrj)=outstate
        eclas=en_array(outstate,itrj)
        ekin_new=ekin_v(vx,vy,vz)

        write(3,'(A24,I3,A10,I3)')'#Hop occured from state ', instate, ' to state ', outstate
        write(3,*)'#Adjusting velocities by simple scaling.'
        write(3,'(A,2E20.10)')'#Total_energy_old   Total_energy_new :',ekin+en_array(instate,itrj),ekin_new+en_array(outstate,itrj)

      else

       write(3,'(A35,I3,A10,I3)')'#Frustrated Hop occured from state ',instate,' to state ',outstate
       if(revmom.eq.1)then
          write(3,*)'#Reversing momentum direction.'
          vx=-vx
          vy=-vy
          vz=-vz
       end if

      endif

      write(3,'(A31,2E20.10)')'#deltaE_pot  E_kin-total',dE,ekin

   end subroutine hop_dot

   subroutine integstep(k_re,k_im,en,y_re,y_im,dotprod,gam)
   real(DP),intent(out) :: k_re(nstmax),k_im(nstmax)
   real(DP),intent(in)  :: dotprod(nstmax,nstmax), gam(nstmax,nstmax)
   real(DP),intent(in)  :: en(nstmax),y_im(nstmax),y_re(nstmax)
   real(DP) :: g
   integer  :: ist1,ist2


   if (phase.eq.0)then
      do ist1=1,nstate
         k_re(ist1)=en(ist1)*y_im(ist1)
         k_im(ist1)=-en(ist1)*y_re(ist1)
         do ist2=1,nstate
            if(ist1.ne.ist2)then
               k_re(ist1)=k_re(ist1)-y_re(ist2)*dotprod(ist1,ist2)
               k_im(ist1)=k_im(ist1)-y_im(ist2)*dotprod(ist1,ist2)
            end if
         enddo
         k_re(ist1)=dtp*k_re(ist1)
         k_im(ist1)=dtp*k_im(ist1)
      enddo
   end if

   if (phase.eq.1)then
      do ist1=1,nstate
         k_re(ist1)=0.0_DP
         k_im(ist1)=0.0_DP
         do ist2=1,nstate
            if(ist1.ne.ist2)then
               g=gam(ist1,ist2)
               k_re(ist1)=k_re(ist1)-(y_re(ist2)*cos(g)-y_im(ist2)*sin(g))*dotprod(ist1,ist2)
               k_im(ist1)=k_im(ist1)-(y_im(ist2)*cos(g)+y_re(ist2)*sin(g))*dotprod(ist1,ist2)
            end if
         end do
         k_re(ist1)=dtp*k_re(ist1)
         k_im(ist1)=dtp*k_im(ist1)
      end do
   end if
   end subroutine integstep

   subroutine integ_gama(en_array_int, en_array_newint, itrj)
   real(DP),intent(in) :: en_array_int(:,:)
   real(DP),intent(in) :: en_array_newint(:,:)
   integer, intent(in) :: itrj
   real(DP) :: g
   integer  :: ist1, ist2

   do ist1=1,nstate
      do ist2=1,nstate
      if (ist1.ne.ist2)then
         g=(en_array_int(ist1,itrj)+en_array_newint(ist1,itrj))/2
         g=g-(en_array_int(ist2,itrj)+en_array_newint(ist2,itrj))/2
         gama(ist1,ist2,itrj)=gama(ist1,ist2,itrj)+g*dtp
      end if
      end do
   end do
   end subroutine integ_gama

   subroutine eulerstep(en_array_int, en_array_newint, dotproduct_int,itrj)
      real(DP),intent(in) :: en_array_int(:,:)
      real(DP),intent(in) :: en_array_newint(:,:)
      real(DP),intent(in) :: dotproduct_int(:,:,:)
      integer,intent(in)  :: itrj
      real(DP) :: dotprod0(nstmax,nstmax), gam0(nstmax, nstmax)
      real(DP) :: k1_re(nstmax),k1_im(nstmax)
      real(DP) :: y_im(nstmax),y_re(nstmax)
      real(DP) :: en0(nstmax)
      integer   :: ist1, ist2

      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
         dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
         gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      end do

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0, gam0)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)
      enddo

      call integ_gama(en_array_int, en_array_newint, itrj)
   end subroutine eulerstep

   subroutine rk4step(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)
      real(DP),intent(in) :: en_array_int(nstmax,ntrajmax)
      real(DP),intent(in) :: en_array_newint(nstmax,ntrajmax)
      real(DP),intent(in) :: dotproduct_int(nstmax,nstmax,ntrajmax)
      real(DP),intent(in) :: dotproduct_newint(nstmax,nstmax,ntrajmax)
      integer, intent(in) :: itrj
      real(DP) :: dotprod2(nstmax,nstmax)
      real(DP) :: dotprod0(nstmax,nstmax),dotprod1(nstmax,nstmax)
      real(DP) :: k1_re(nstmax),k1_im(nstmax)
      real(DP) :: k2_re(nstmax),k2_im(nstmax)
      real(DP) :: k3_re(nstmax),k3_im(nstmax)
      real(DP) :: k4_re(nstmax),k4_im(nstmax)
      real(DP) :: y_im(nstmax),y_re(nstmax)
      real(DP) :: en0(nstmax),en1(nstmax),en2(nstmax)
      real(DP) :: gam0(nstmax,nstmax),gam1(nstmax,nstmax),gam2(nstmax,nstmax)
      integer :: ist1,ist2     

!     initial interpolations
      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       en1(ist1)=en_array_newint(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
       dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
       dotprod1(ist1,ist2)=dotproduct_newint(ist1,ist2,itrj)
       dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
       if(phase.eq.1) gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      enddo

      if (phase.eq.1)then
         call integ_gama(en_array_int, en_array_newint, itrj)
         do ist1=1,nstate
            do ist2=1,nstate
               gam1(ist1,ist2)=gama(ist1,ist2,itrj)
               gam2(ist1,ist2)=(gam0(ist1,ist2)+gam1(ist1,ist2))/2
            end do
         end do 
      end if

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0,gam0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/2
      enddo
      call integstep(k2_re,k2_im,en2,y_re,y_im,dotprod2,gam2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k2_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k2_im(ist1)/2
      enddo
      call integstep(k3_re,k3_im,en2,y_re,y_im,dotprod2,gam2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en1,y_re,y_im,dotprod1,gam1)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)/6+k2_re(ist1)/3+k3_re(ist1)/3+k4_re(ist1)/6
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)/6+k2_im(ist1)/3+k3_im(ist1)/3+k4_im(ist1)/6
      enddo

   end subroutine rk4step


   subroutine butcherstep(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)
      real(DP),intent(in) :: en_array_int(nstmax,ntrajmax)
      real(DP),intent(in) :: en_array_newint(nstmax,ntrajmax)
      real(DP),intent(in) :: dotproduct_int(nstmax,nstmax,ntrajmax)
      real(DP),intent(in) :: dotproduct_newint(nstmax,nstmax,ntrajmax)
      real(DP) :: dotprod2(nstmax,nstmax),dotprod4(nstmax,nstmax),dotprod34(nstmax,nstmax)
      real(DP) :: dotprod0(nstmax,nstmax),dotprod1(nstmax,nstmax)
      real(DP) :: k1_re(nstmax),k1_im(nstmax)
      real(DP) :: k2_re(nstmax),k2_im(nstmax)
      real(DP) :: k3_re(nstmax),k3_im(nstmax)
      real(DP) :: k4_re(nstmax),k4_im(nstmax)
      real(DP) :: k5_re(nstmax),k5_im(nstmax)
      real(DP) :: k6_re(nstmax),k6_im(nstmax)
      real(DP) :: y_im(nstmax),y_re(nstmax)
      real(DP) :: en0(nstmax),en1(nstmax),en2(nstmax),en4(nstmax),en34(nstmax)
      real(DP) :: gam0(nstmax,nstmax),gam1(nstmax,nstmax),gam2(nstmax,nstmax)
      real(DP) :: gam4(nstmax,nstmax),gam34(nstmax,nstmax)
      integer :: ist1,ist2,itrj     !iteration counters

!pripravne interpolace....
      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       en1(ist1)=en_array_newint(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       en4(ist1)=(en0(ist1)+en2(ist1))/2
       en34(ist1)=(en2(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
        dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
        dotprod1(ist1,ist2)=dotproduct_newint(ist1,ist2,itrj)
        dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
        dotprod4(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod2(ist1,ist2))/2
        dotprod34(ist1,ist2)=(dotprod2(ist1,ist2)+dotprod1(ist1,ist2))/2
        if (phase.eq.1) gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      enddo

      if (phase.eq.1)then
         call integ_gama(en_array_int, en_array_newint, itrj)
         do ist1=1,nstate
            do ist2=1,nstate
               gam1(ist1,ist2)=gama(ist1,ist2,itrj)
               gam2(ist1,ist2)=(gam0(ist1,ist2)+gam1(ist1,ist2))/2
               gam4(ist1,ist2)=(gam0(ist1,ist2)+gam2(ist1,ist2))/2
               gam34(ist1,ist2)=(gam2(ist1,ist2)+gam1(ist1,ist2))/2
            end do
         end do 
      end if

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0,gam0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/4
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/4
      enddo
      call integstep(k2_re,k2_im,en4,y_re,y_im,dotprod4,gam4)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/8+k2_re(ist1)/8
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/8+k2_im(ist1)/8
      enddo
      call integstep(k3_re,k3_im,en4,y_re,y_im,dotprod4,gam4)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-k2_re(ist1)/2+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)-k2_im(ist1)/2+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en2,y_re,y_im,dotprod2,gam2)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+3*k1_re(ist1)/16+9*k4_re(ist1)/16
       y_im(ist1)=cel_im(ist1,itrj)+3*k1_im(ist1)/16+9*k4_im(ist1)/16
      enddo
      call integstep(k5_re,k5_im,en34,y_re,y_im,dotprod34,gam34)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-3*k1_re(ist1)/7+2*k2_re(ist1)/7+12*k3_re(ist1)/7 &
                  -12*k4_re(ist1)/7+8*k5_re(ist1)/7
       y_im(ist1)=cel_im(ist1,itrj)-3*k1_im(ist1)/7+2*k2_im(ist1)/7+12*k3_im(ist1)/7 &
                  -12*k4_im(ist1)/7+8*k5_im(ist1)/7
      enddo
      call integstep(k6_re,k6_im,en1,y_re,y_im,dotprod1,gam1)


      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+7*k1_re(ist1)/90+32*k3_re(ist1)/90+12*k4_re(ist1)/90 &
                        +32*k5_re(ist1)/90+7*k6_re(ist1)/90
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+7*k1_im(ist1)/90+32*k3_im(ist1)/90+12*k4_im(ist1)/90 &
                        +32*k5_im(ist1)/90+7*k6_im(ist1)/90
      enddo

   end subroutine butcherstep


   subroutine interpolate_dot(dotproduct_int,fr,frd,itrj)
      real(DP) dotproduct_int(nstmax,nstmax,ntrajmax)
      integer :: ist1,ist2,itrj     !iteration counters
      real(DP) :: fr,frd

      do ist1=1,nstate
       do ist2=1,nstate
        dotproduct_int(ist1,ist2,itrj)=dotproduct(ist1,ist2,itrj)*fr + &
                dotproduct_old(ist1,ist2,itrj)*frd
       enddo
      enddo

   end subroutine interpolate_dot

   subroutine check_energy(vx_old, vy_old, vz_old, vx, vy, vz, itrj)
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(in) :: vx_old(:,:),vy_old(:,:),vz_old(:,:)
      integer, intent(in) :: itrj
      real(DP)            :: ekin, ekin_old, entot, entot_old

      ekin=ekin_v(vx, vy, vz)
      ekin_old=ekin_v(vx_old, vy_old, vz_old)

      entot=(ekin+en_array(istate(itrj), itrj) )*AUtoEV
      entot_old=(ekin_old+en_array_old(istate(itrj), itrj) )*AUtoEV

      if (abs(entot-entot_old).gt.energydifthr)then
         write(*,*)'ERROR:Poor energy conservation. Exiting...'
         write(*,*)'Total energy difference [eV] is:', entot-entot_old
         write(*,*)'The threshold was:',energydifthr
         call abinerror('check_energy')
      end if

   end subroutine check_energy

   subroutine check_energydrift(vx, vy, vz, itrj)
      use mod_const, only: AUtoEV
      use mod_kinetic, only: ekin_v
      real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
      integer, intent(in) :: itrj
      real(DP)            :: ekin, entot

      ekin=ekin_v(vx, vy, vz)

      entot=(ekin+en_array(istate(itrj), itrj) )*AUtoEV

      if (abs(entot-entot0).gt.energydriftthr)then
         write(*,*)'ERROR: Energy drift exceeded threshold value. Exiting...'
         write(*,*)'Total energy difference [eV] is:', entot-entot0
         write(*,*)'The threshold was:',energydriftthr
         call abinerror('check_energy')
      end if

   end subroutine check_energydrift

end module mod_sh
