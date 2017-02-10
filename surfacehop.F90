! Driver routines for Surface hopping dynamics
! D. Hollas, O. Svoboda, P. Slavíček, M. Ončák

module mod_sh
   use mod_const, only: DP
   use mod_array_size, only: NSTMAX, NTRAJMAX
   use mod_utils, only: abinerror, printf
   use mod_sh_integ
   implicit none
   private
   public :: istate_init, substep, deltaE, inac,nohop, alpha,popthr, nac_accu1, nac_accu2
   public :: surfacehop, sh_init, istate, ntraj, tocalc, en_array
   public :: nacx, nacy, nacz
   public :: move_vars, get_nacm, write_nacmrest, read_nacmrest
   public :: energydifthr, energydriftthr, popsumthr, adjmom, revmom
   public :: check_CIVector
   public :: ignore_state

   integer,parameter :: ntraj = NTRAJMAX
   integer   :: istate_init=1, substep=100
   integer   :: inac=0, nohop=0, adjmom=0, revmom=0
   integer   :: nac_accu1=7, nac_accu2=5 !7 is MOLPRO default
   real(DP)  :: alpha=0.1d0
   real(DP)  :: deltae=5.0d0, popthr=0.001d0
   real(DP)  :: popsumthr=0.001d0, energydifthr=1.0d0, energydriftthr=1.0d0 !eV
   real(DP),allocatable :: nacx(:,:,:,:), nacy(:,:,:,:), nacz(:,:,:,:)
   real(DP),allocatable :: nacx_old(:,:,:,:), nacy_old(:,:,:,:), nacz_old(:,:,:,:)
   real(DP),allocatable :: dotproduct(:,:,:), dotproduct_old(:,:,:) !for inac=1
   real(DP),allocatable :: en_array(:,:), en_array_old(:,:)
   real(DP) :: entot0
   integer, allocatable :: tocalc(:,:)
   integer  :: istate(NTRAJMAX)
   ! for ethylene 2-state SA3 dynamics
   integer  :: ignore_state = 0
   save

   CONTAINS

   subroutine sh_init(x, y, z, vx, vy, vz, dt)
   use mod_const,      only: AUtoEV
   use mod_general,    only: irest, natom, it, pot
   use mod_interfaces, only: force_clas
   use mod_kinetic,    only: ekin_v
   real(DP),intent(inout)  :: x(:,:), y(:,:), z(:,:)
   real(DP),intent(in)     :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: dt
   real(DP)  :: dum_fx(size(x,1),size(x,2))
   real(DP)  :: dum_fy(size(y,1),size(y,2))
   real(DP)  :: dum_fz(size(z,1),size(z,2))
   real(DP)  :: dum_eclas
   integer   :: itrj,ist1

   deltaE = deltaE / AUtoEV
   dtp = dt / substep

   allocate( nacx(natom, NTRAJMAX, nstate, nstate) )
   allocate( nacy(natom, NTRAJMAX, nstate, nstate) )
   allocate( nacz(natom, NTRAJMAX, nstate, nstate) )
   allocate( nacx_old(natom, NTRAJMAX, nstate, nstate) )
   allocate( nacy_old(natom, NTRAJMAX, nstate, nstate) )
   allocate( nacz_old(natom, NTRAJMAX, nstate, nstate) )
   nacx=0.0d0
   nacy=0.0d0
   nacz=0.0d0
   nacx_old=nacx
   nacy_old=nacx
   nacz_old=nacx
   if (inac.eq.1)then
      allocate( dotproduct_old(nstate, nstate, NTRAJMAX) )
      allocate( dotproduct(nstate, nstate, NTRAJMAX) )
      dotproduct_old=0.0d0
      dotproduct=dotproduct_old
      adjmom=1
   end if

   allocate( en_array(nstate, ntraj) )
   allocate( en_array_old(nstate, ntraj) )
   en_array=0.0d0
   en_array_old=en_array

!  Determining the initial state
   if(irest.ne.1)then

      do itrj=1,ntraj
 
!        Automatic determination of initial state based on osc. strenght
         if(istate_init.eq.-1)then
            call choose_initial_state(itrj)
         else
            istate(itrj) = istate_init
         endif
 
      enddo
 
   endif

!--computing only energies, used for subsequent determination of TOCALC according to deltaE threshold
   allocate( tocalc(nstate, nstate) )
   tocalc = 0
   dum_eclas = 0.0d0
   dum_fx=0.0d0   ; dum_fy=0.0d0    ;dum_fz=0.0d0

   call force_clas(dum_fx, dum_fy, dum_fz, x, y , z, dum_eclas, pot)

   do itrj = 1, ntraj

      call sh_set_initialwf(istate(itrj), en_array(1, itrj), itrj)

      ! WARNING: entot0 does not respect itrj, same for tocalc
      entot0 = en_array(istate(itrj), itrj)+ekin_v(vx, vy, vz)
      entot0 = entot0 * AUtoEV

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
   ! WARNING: tocalc array is currently the same for all trajs
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
      !COMPUTE NACME only if population of the states is gt.popthr
      do ist1=1,nstate-1
         pop = cel_re(ist1,itrj)**2 + cel_im(ist1,itrj)**2
         do ist2=ist1+1,nstate
            pop2 = cel_re(ist2,itrj)**2 + cel_im(ist2,itrj)**2
            if(pop.lt.popthr.and.pop2.lt.popthr.and.ist1.ne.istate(itrj).and.ist2.ne.istate(itrj)) tocalc(ist1,ist2)=0
         enddo
      enddo
   endif
   end subroutine set_tocalc

   subroutine Write_nacmrest()
   use mod_general, only: narchive, it
   use mod_qmmm, only:natqm
   use mod_utils, only: archive_file
   integer :: ist1,ist2,iat,itrj
   integer :: iunit1
   iunit1 = 600
   open(iunit1,file='nacmrest.dat',action="write")

   do itrj=1, ntraj

   do ist1=1,nstate-1

      do ist2=ist1+1,nstate
   
         if(tocalc(ist1,ist2).eq.1.and.inac.eq.0)then
   
            write(iunit1,*)'NACME between states',ist1,ist2
            do iat=1,natqm              ! reading only for QM atoms
               write(iunit1,*)nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
            enddo
   
         endif
   
      enddo
   enddo

   ! Printing total energy at t=0, so that we can safely restart
   ! and we do not break checking for energy drift
   ! For now, energy jump check is not handled well.
   write(iunit1, *)'Total Energy(t=0) [eV]'
   write(iunit1, *)entot0

   if(phase.eq.1)then
      call write_phaserest(iunit1, itrj)
   end if
   
   ! ntraj enddo
   end do

   close(iunit1)


   if(modulo(it,narchive).eq.0)then
      call archive_file('nacmrest.dat',it)
   end if

   end subroutine write_nacmrest

   ! TODO: rename to sh_restart
   ! save everything SH related to file sh_restart.dat
   subroutine read_nacmrest()
   use mod_general, only: it
   use mod_qmmm,  only: natqm
   use mod_utils, only: archive_file
   integer  :: iost,ist1,ist2,iat,itrj
   integer  :: iunit1, iunit2
   character(len=200)   :: chmsg
   character(len=20)    :: chit
   character(len=60)    :: chrestart
   iunit1=600; iunit2=601

   if (inac.eq.0) write(*,*)'Reading NACME from nacmrest.dat'
   open(iunit1,file='nacmrest.dat',action="read")
   if (phase.eq.1)then
      write(*,*)'Reading phase from phaserest.dat'
      open(iunit2,file='phaserest.dat',action="read")
   end if

   do itrj=1, ntraj

   do ist1=1,nstate-1

      do ist2=ist1+1,nstate
   
         if(tocalc(ist1,ist2).eq.1.and.inac.eq.0)then
   
            read(iunit1,*, iostat=iost)
            do iat=1,natqm              ! reading only for QM atoms
               read(iunit1,*,iomsg=chmsg)nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
               if (iost.ne.0)then
                  write(*,*)'Error reading NACME from file nacmrest.'
                  write(*,*)chmsg
                  call abinerror('read_nacmrest')
               end if
               nacx(iat,itrj,ist2,ist1)=-nacx(iat,itrj,ist1,ist2)
               nacy(iat,itrj,ist2,ist1)=-nacy(iat,itrj,ist1,ist2)
               nacz(iat,itrj,ist2,ist1)=-nacz(iat,itrj,ist1,ist2)
            enddo
   
   !--------if tocalc 
         endif
   
      enddo
   enddo

   ! Reading total energy at t=0, so that we can safely restart
   ! and we do not break checking for energy drift
   ! For now, energy jump check is not handled well.
   read(iunit1, *)
   read(iunit1, *)entot0

   if(phase.eq.1)then
      do ist1=1,nstate
         read(iunit2,*)(gama(ist1,ist2,itrj),ist2=1,nstate)
      end do
   end if

   end do

   close(iunit1)
   call archive_file('nacmrest.dat',it)

   if (phase.eq.1)then
      close(iunit2)
      call archive_file('phaserest.dat',it)
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
                  write(*,*)'WARNING: NACME between states',ist1,ist2,'not read.'
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
   use mod_utils, only: LowerToUpper
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

   chsystem='./'//trim(LowerToUpper(pot))//'/r.'//pot
   ! We do not check whether the script is sensitive to the 3rd parameter
   write(*,*)'WARNING: Some NACMs not computed. Trying with decreased accuracy...'
   write(*,*)'Calling script r.'//pot//'with accuracy:',nac_accu2
   write(chsystem,'(A20,I13,I4.3,I3,A12)')chsystem,it,itrj,nac_accu2,' < state.dat'

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
         write(*,*)'ERROR:Sum of electronic populations = ',popsum
         write(*,*)'which differs from 1.0 by more than popsumthr = ',popsumthr
         write(*,*)'Increase the number of substeps or use more accurate integrator.'
         call abinerror('surfacehop')
      end if
   end subroutine check_popsum


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This is the main routine !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine surfacehop(x,y,z,vx,vy,vz,vx_old,vy_old,vz_old,dt,eclas)
   use mod_const,    only: ANG, AUTOFS
   use mod_general,  only: natom, pot, nwrite, idebug, it, sim_time
   use mod_system,   ONLY: names
   use mod_files,    ONLY: UPOP,UPROB,UPES,UWFCOEF,UWFCOEF,UNACME, UBKL,UPHASE,UDOTPROD
   use mod_qmmm,     ONLY: natqm
   use mod_random,   ONLY: vranf
   use mod_kinetic,  ONLY: ekin_v
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
   real(DP)  :: en_array_int(NSTMAX,NTRAJMAX),en_array_newint(NSTMAX,NTRAJMAX)
   real(DP)  :: nacx_int(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)  
   real(DP)  :: nacy_int(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)
   real(DP)  :: nacz_int(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)
   real(DP)  :: nacx_newint(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)
   real(DP)  :: nacy_newint(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)
   real(DP)  :: nacz_newint(size(vx,1),NTRAJMAX,NSTMAX,NSTMAX)
   real(DP)  :: dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX),dotproduct_newint(NSTMAX,NSTMAX,NTRAJMAX) !rename dotproduct_int
   real(DP)  :: t(NSTMAX,NSTMAX)           !switching probabilities
   real(DP)  :: t_tot(NSTMAX,NSTMAX)       !cumulative switching probabilities
   real(DP)  :: ran(10)
   real(DP)  :: pop(size(tocalc,1),NTRAJMAX),popsum !populations
   integer :: iat,ist1,ist2,itrj,itp     !iteration counters
   integer :: ist                    ! =istate(itrj)
   real(DP)  :: vect_olap,fr,frd
   real(DP)  :: ekin_mom,apom,edif,tau,fact,sum_norm
   integer :: ihop,ijunk
   real(DP)  :: pop0, a_re, a_im, prob(NSTMAX), hop_rdnum, stepfs
   character(len=500) :: formt
   character(len=20) :: chist,chihop,chit
   integer :: iost, iunit=600


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

      if(pot.ne.'_tera_')then
         iost=readnacm(itrj)
!------------if NACME NOT COMPUTED: TRY TO DECREASE ACCURACY--------------
         if(iost.ne.0.and.nac_accu1.gt.nac_accu2)then
            call calcnacm(itrj)

            iost=readnacm(itrj)
         endif
!------------if NACME STILL NOT COMPUTED: USE OLD NACM--------------
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

      end if

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
               write(*,*)'Not computing NACM for states',ist1,ist2
               dotproduct(ist1,ist2,itrj)=0.0d0
               dotproduct(ist2,ist1,itrj)=0.0d0
            endif
         enddo
      enddo


      !  we don't have interpolation in the zeroth step
      if( it.eq.1) dotproduct_old=dotproduct
   ! TDC end if
   endif
!------------------END-OF-TDC-------------------------------

!#### MAIN LOOP

!  Smaller time step for electronic population transfer
   do itp=1,substep      

      ist=istate(itrj)
      ! pop0 is later used for Tully's fewest switches
      pop0=cel_re(ist,itrj)**2+cel_im(ist,itrj)**2


!-----INTERPOLATION

      fr = real(itp,DP) / real(substep,DP)
      frd = 1.0d0 - fr

      call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_newint,vy_newint,vz_newint, &
                     nacx_newint,nacy_newint,nacz_newint,en_array_newint, &
                     dotproduct_newint,fr,frd,itrj)

      ! In HST model, we do not have NACME, only dotproduct
      if(inac.eq.1)then
         call interpolate_dot(dotproduct_newint, fr, frd, itrj)
      endif

      fr = real(itp-1,DP) / real(substep,DP)
      frd = 1.0d0 - fr

      call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_int,vy_int,vz_int, &
                      nacx_int,nacy_int,nacz_int,en_array_int, &
                      dotproduct_int,fr,frd,itrj)

      if(inac.eq.1)then
        call interpolate_dot(dotproduct_int,fr,frd,itrj)
      endif

!-----END-OF-INTERPOLATIONS
       

      ! Integrate electronic wavefunction for one dtp time step
      call sh_integrate_wf(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)

      ! Check whether total population is 1.0
      call check_popsum(itrj)

      !--calculation of switching probabilities
      !- do not calculated the whole transition matrix
      !- t array could be one dimensional but whatever

      ! TODO: make this a separate function
      ! Prepare it for other schemes, i.e. for autoionization

      ! call sh_calc_hopprob(pop0, trans_matrix)
      do ist2=1,nstate 
         pop(ist2,itrj)=cel_re(ist2,itrj)**2+cel_im(ist2,itrj)**2
         a_re=( cel_re(ist,itrj)*cel_re(ist2,itrj)+cel_im(ist,itrj)*cel_im(ist2,itrj) )
         if (phase.eq.1)then
             a_im=(- cel_im(ist,itrj)*cel_re(ist2,itrj) + cel_re(ist, itrj)*cel_im(ist2, itrj) )
             t(ist,ist2)=a_re*cos(gama(ist, ist2, itrj))-sin(gama(ist, ist2, itrj))*a_im
             t(ist,ist2)=t(ist,ist2)*(dotproduct_int(ist,ist2,itrj) + dotproduct_newint(ist,ist2,itrj) )
         else
!            t(ist,ist2)=2*a_re*dotproduct_int(ist,ist2,itrj)
             t(ist,ist2)=a_re*(dotproduct_int(ist,ist2,itrj)+dotproduct_newint(ist,ist2,itrj) )
         endif
         t(ist,ist2)=t(ist,ist2)*dtp/(pop0+1d-20)
      enddo

      if(idebug.gt.1)then
         ! WaRNING: this will not work for adaptive time step
         stepfs=(it*substep+itp-substep)*dt*AUtoFS/substep
         write(formt,'(A7,I3,A7)')'(F15.2,',nstate,'E20.10)'
         write(UBKL,fmt=formt)stepfs,(t(ist,ist1),ist1=1,nstate)

         write(UWFCOEF,fmt=formt,advance="no")stepfs,(cel_re(ist1,itrj),ist1=1,nstate)
         write(formt,'(A1,I3,A7)')'(',nstate,'E20.10)'
         write(UWFCOEF,fmt=formt,advance="no")(cel_im(ist1,itrj),ist1=1,nstate)
         write(UWFCOEF,*)''

         if(phase.eq.1)then
            write(UPHASE,'(F15.2,E20.10)',advance="no")stepfs,gama(2,1,itrj)
            do ist1=3,nstate
               write(formt,'(A1,I3,A7)')'(',ist1-1,'E20.10)'
               write(UPHASE,fmt=formt,advance="no")(gama(ist1,ist2,itrj),ist2=1,ist1-1)
            end do
            write(UPHASE,*)''
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

!     Auxiliary calculations of probabilities
      do ist1=1,nstate
         if(ist1.eq.ist.and.ist1.eq.1)  prob(ist1) = 0.0d0
         if(ist1.ne.ist.and.ist1.eq.1)  prob(ist1) = t(ist,ist1)
         if(ist1.ne.ist.and.ist1.ne.1)  prob(ist1) = prob(ist1-1)+t(ist,ist1)
         if(ist1.eq.ist.and.ist1.ne.1)  prob(ist1) = prob(ist1-1)
      enddo

      ! Should we hop before decoherence?
      ! Every article says something different.
      ! Newton-X hops before decoherence.

!     HOPPING SECTION
      if (nohop.ne.1)then

         ihop=0
         call vranf(ran,1,0,6)
         hop_rdnum = ran(1)
         ! TODO streamline this
         !do ist1=1, nstate
         !   if (ist1.eq.ist) cycle
         !   if(hop_rdum.lt.prob(ist1))then
         !      ihop = ist1
         !      exit
         !   end if
         !end do

         if(ist.ne.1)then
            if(hop_rdnum.ge.0.0d0.and.hop_rdnum.lt.prob(1))then
               ihop=1
            endif
         endif
       
         do ist1=2,nstate
            if(ist.ne.ist1)then
               if(hop_rdnum.ge.prob(ist1-1).and.hop_rdnum.lt.prob(ist1))then
                  ihop=ist1
               endif
            endif
         enddo
       
!        Did HOP occur?
         if(ihop.ne.0)then
            if(adjmom.eq.0) call hop(vx,vy,vz,ist,ihop,itrj,eclas)
            if(adjmom.eq.1) call hop_dot(vx,vy,vz,ist,ihop,itrj,eclas)
            write(formt,'(A8,I3,A7)')'(A1,I10,',nstate+1,'E20.10)'
            if(idebug.gt.0)then
              write(UPOP,*)'# Substep RandomNum   Probabilities'
              write(UPOP,fmt=formt)'#',itp, hop_rdnum, (t(ist,ist1),ist1=1,nstate)
            end if
       
            ! write current geometry
            write(formt,'(A5,I10,A1,I2,A1,I2)')'geom.',it,'.',ist,'.',ihop
            write(chist,*)ist
            write(chihop,*)ihop
            write(chit,*)it
            formt='geom.'//trim(adjustl(chist))//'.'//trim(adjustl(chihop))//'.'//adjustl(chit)
            open(500,file=trim(formt))
            write(500,*)natom
            write(500,*)''
            do iat=1,natom
               write(500,*)names(iat),x(iat,itrj)/ANG,y(iat,itrj)/ANG,z(iat,itrj)/ANG
            enddo
            close(500)
         endif   
       
         !nohop endif
      endif

!     Quantum decoherence part----------------------------------
      if(alpha.gt.0)then
      
         ekin_mom=ekin_v(vx_int,vy_int,vz_int)

         if(ekin_mom.gt.1.0d-4)then !Why this number? In almost all cases it won't matter

            do ist1=1,nstate
               if(ist1.ne.istate(itrj)) then

!                 Calculation of exponential factor               
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
            fact=sum_norm/(cel_re(istate(itrj),itrj)**2+cel_im(istate(itrj),itrj)**2+1.0d-7)

            !Following should never happen as we check for popsumthr later in this subroutine
            if(fact.lt.0.0d0)then
               write(*,*)'Fatal error in surfacehop during decoherence renormalization.'
               write(*,*)'fact=',fact,'but should be > 0'
               write(*,*)'This usually means inaccurate integration of electronic SE.'
               write(*,*)'Increase number of substeps or use more accurate integrator.'
               call abinerror('surfacehop')
            end if

            fact=sqrt(fact)

            cel_re(istate(itrj),itrj)=cel_re(istate(itrj),itrj)*fact
            cel_im(istate(itrj),itrj)=cel_im(istate(itrj),itrj)*fact

         ! ekin endif
         endif

!-----END-OF-DECOHERENCE------------------------      
      endif



!itp loop
   enddo

!  set tocalc array for the next step
   call set_tocalc(itrj)

   if(idebug.eq.1)then
      ! this is rarely needed, since we always print dotproduct into
      ! dotprod.dat
      ! this is only to check that the matrix is really symmetrical
      if(it.eq.1)then
         open(iunit,file='dotprodmatrix.dat')
      else
         open(iunit,file='dotprodmatrix.dat',access='append')
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
      stepfs = sim_time * AUtoFS
      write(formt,'(A10,I3,A13)')'(F15.2,I3,',nstate,'F10.5,1F10.7)'
      write(UPOP,fmt=formt)stepfs,istate(itrj),(pop(ist1,itrj), ist1=1,nstate),popsum

      t_tot=1-t_tot !up to know, t_tot is the probability of not hopping
      write(formt,'(A10,I3,A6)')'(F15.2,I3,',nstate,'F10.5)'
      write(UPROB,fmt=formt)stepfs,istate(itrj),(t_tot(ist,ist1),ist1=1,nstate)
      write(formt,'(A7,I3,A7)')'(F15.2,',nstate,'E20.10)'
      write(UPES,fmt=formt)stepfs,(en_array(ist1,itrj),ist1=1,nstate)
      if (inac.eq.0)  write(UNACME,*)'Time step:',it
      do ist1=1,nstate-1
         do ist2=ist1+1,nstate

            if (ist1.eq.1.and.ist2.eq.2)then
               write(UDOTPROD,'(F15.2,E20.10)',advance='no')stepfs,dotproduct_int(ist1,ist2,itrj)
            else
               write(UDOTPROD,'(E20.10)',advance='no')dotproduct_int(ist1,ist2,itrj)
            end if

            if (inac.eq.0)then
               write(UNACME,*)'NACME between states:', ist1, ist2
               do iat=1,natom
                  write(UNACME,'(3E20.10)')nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
               enddo
            endif

         end do
      end do
      write(UDOTPROD,*)''

   endif

! ntraj enddo       
enddo

   end subroutine surfacehop

   subroutine hop(vx,vy,vz,instate,outstate,itrj,eclas)
      use mod_general,  ONLY: natom, pot
      use mod_system,   ONLY: am
      use mod_files,    ONLY: UPOP
      use mod_kinetic,  only: ekin_v
      use mod_arrays,   ONLY: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(inout) :: eclas
      integer, intent(in)    :: itrj,instate,outstate
      real(DP)  :: a_temp,b_temp,c_temp,g_temp
      real(DP)  :: ekin, ekin_new
      integer   :: iat
      integer, allocatable :: tocalc_temp(:,:)

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
         write(*,*)'# Not enough momentum in the direction of NAC vector.'
         ! Try, whether there is enough total kinetic energy and scale velocities.
         call hop_dot(vx,vy,vz,instate,outstate,itrj,eclas)
         return
      endif

      istate(itrj)=outstate
      eclas=en_array(outstate,itrj)
      write(*,'(A24,I3,A10,I3)')'# Hop occured from state ',instate,' to state ',outstate


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

      write(*,'(A31,2E20.10)')'# deltaE_pot     E_kin-total',en_array(outstate,itrj)-en_array(instate,itrj),ekin
      write(*,'(A,2E20.10)')'# Total_energy_old   Total_energy_new :',ekin+en_array(instate,itrj),ekin_new+en_array(outstate,itrj)

      allocate(tocalc_temp(nstate, nstate))
      tocalc_temp = tocalc
      tocalc = 0
      write(*,*)'Calculating correct forces for the new state.'
      call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
      tocalc = tocalc_temp
      deallocate(tocalc_temp)

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
      use mod_general, ONLY:natom, idebug, pot
      use mod_kinetic ,ONLY: ekin_v
      use mod_arrays,   ONLY: fxc, fyc, fzc, x, y, z
      use mod_interfaces, only: force_clas
      real(DP),intent(inout)  :: vx(:,:),vy(:,:),vz(:,:)
      real(DP),intent(inout)  :: eclas
      integer,intent(in)      :: itrj,instate,outstate
      integer   :: iat
      real(DP)  :: de,ekin,alfa,ekin_new
      integer, allocatable    :: tocalc_temp(:,:)

      ekin=0.0d0
      ekin_new=0.0d0

      dE=en_array(outstate,itrj)-en_array(instate,itrj)
      ekin=ekin_v(vx,vy,vz)

      if(ekin.ge.de)then

         alfa=sqrt(1-de/ekin)

         do iat=1,natom
            vx(iat,itrj) = alfa * vx(iat,itrj) 
            vy(iat,itrj) = alfa * vy(iat,itrj) 
            vz(iat,itrj) = alfa * vz(iat,itrj) 
         enddo
         istate(itrj) = outstate
         eclas = en_array(outstate,itrj)
         ekin_new = ekin_v(vx,vy,vz)

         write(*,'(A24,I3,A10,I3)')'# Hop occured from state ', instate, ' to state ', outstate
         write(*,*)'# Adjusting velocities by simple scaling.'
         write(*,'(A,2E20.10)')'#Total_energy_old   Total_energy_new :',ekin+en_array(instate,itrj),ekin_new+en_array(outstate,itrj)

         allocate(tocalc_temp(nstate, nstate))
         tocalc_temp = tocalc
         tocalc = 0
         write(*,*)'Calculating correct forces for the new state.'
         call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
         tocalc = tocalc_temp
         deallocate(tocalc_temp)

      else

         write(*,'(A35,I3,A10,I3)')'# Frustrated Hop occured from state ',instate,' to state ',outstate
         if(revmom.eq.1)then
            write(*,*)'# Reversing momentum direction.'
            vx = -vx
            vy = -vy
            vz = -vz
         end if

      endif

      write(*,'(A31,2E20.10)')'# deltaE_pot  E_kin-total',dE,ekin


   end subroutine hop_dot



   subroutine interpolate_dot(dotproduct_int,fr,frd,itrj)
      real(DP) dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX)
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

   integer function check_CIVector(CIvecs, CIvecs_old, ci_len, nstates)
   use mod_const, only:AUtoFS
   use mod_files, only: UDOTPRODCI
   use mod_general, only: sim_time
   real(DP), allocatable, intent(in) :: CIvecs(:,:), CIvecs_old(:,:)
   integer, intent(in) :: ci_len, nstates
   real(DP)  :: cidotprod(NSTMAX)
   integer   :: ist1, ist2, i
   character(len=20) :: formt

   do ist1=1, nstates
      cidotprod(ist1)=0.0d0
      do i=1,ci_len
         cidotprod(ist1) = cidotprod(ist1) + CIvecs(i,ist1)*CIvecs_old(i,ist1)
      end do
      if(cidotprod(ist1).lt.1/dsqrt(2.0d0))then
         write(*,*)"Warning: cidotprod too low."
         ! call abinerror("surfacehop.f90")
      end if
   end do

   write(formt,'(A7,I3,A7)')'(F15.2,',nstates,'F15.10)'
   write(UDOTPRODCI,fmt=formt)sim_time * AUtoFS,(cidotprod(ist1),ist1=1,nstates)

   check_CIVector=0
   return
   end function check_CIVector


   ! Choose initial state according to oscillator stretgth
   subroutine choose_initial_state(itrj)
      integer, intent(in)  :: itrj
      real(DP) :: pom, maxosc
      integer  :: ist1
      open(500,file='oscil.dat')

      pom = 0.0d0
      maxosc = 0.0d0

      do ist1=1,nstate
         read(500,*)pom
         if(pom.gt.maxosc)then
            istate(itrj) = ist1
            maxosc = pom
         endif
      enddo
 
      close(500)
   end subroutine choose_initial_state

end module mod_sh
