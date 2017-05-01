!----Initial version                    by Daniel Hollas,9.2.2012
!- This module contains some routines that do analyses and I/O operations.
!- It also contains routines performing restart.

module mod_analysis
   use mod_const, only: DP
   use mod_files
   implicit none
   private 
   public :: trajout, restout, analysis, restin
   !These are the character string that we check for during restart.
   character(len=*),parameter :: chnose='NHC momenta',chQT='Quantum Thermostat',chLT='Langevin Thermostat', &
            chSH='Coefficients for SH', chAVG='Cumulative averages of various estimators', &
            chcoords='Cartesian Coordinates [au]', chvel='Cartesian Velocities [au]'

   CONTAINS

!  Contains all analysis stuff
   SUBROUTINE analysis(x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
   use mod_analyze_ext, only: analyze_ext
   use mod_estimators , only: estimators
   use mod_general,     only: it, ipimd, icv,nwrite, nwritef, nwritev, &
                              nrest, nwritex, imini, nstep, anal_ext, idebug
   use mod_system
   use mod_density
   use mod_io 
   use mod_vinit,    only: remove_rotations
   use mod_system,   only: am
   implicit none
   !intent inout because of estimators, writing to nwalk+1
   real(DP),intent(inout) :: x(:,:),   y(:,:),   z(:,:)
   real(DP),intent(in)    :: fxc(:,:), fyc(:,:), fzc(:,:)
   real(DP),intent(inout)    :: vx(:,:),  vy(:,:),  vz(:,:)
   real(DP),intent(in)    :: amt(:,:)
   real(DP),intent(in)    :: eclas, equant
   real(DP) :: dt  !,energy

!  eclas comes from force_clas,equant from force_quantum
!  energy=eclas+equant

   if(modulo(it,nwrite).eq.0.and.idebug.gt.0) call remove_rotations(x, y, z, vx, vy, vz, am, .false.)

   if (ipimd.eq.1.or.icv.eq.1)then
      call estimators(x, y, z, fxc, fyc, fzc, eclas, dt)
   endif

   if(modulo(it,nwritex).eq.0)then
      call trajout(x, y, z, it)
   endif

   if(nwritev.gt.0)then
      if(modulo(it,nwritev).eq.0)then
         call velout(vx, vy, vz)
      endif
   endif

   if(nwritef.gt.0)then
      if(modulo(it,nwritef).eq.0)then
         call forceout(x, y, z, fxc, fyc, fzc, UFORCE )
      endif
   endif

   if(ndist.ge.1.and.it.gt.imini)then
      call density(x, y, z)
   endif

   if(nang.ge.1.and.it.gt.imini)then
      call density_ang(x, y, z)
   endif

   if(ndih.ge.1.and.it.gt.imini)then
      call density_dih(x, y, z)
   endif

   if((modulo(it,nrest).eq.0).or.it.eq.nstep)then
      call restout(x,y,z,vx,vy,vz,it)
   endif

   if (anal_ext.eq.1)then
      call analyze_ext(x, y, z, vx, vy, vz, amt)
   endif

   end subroutine analysis

   subroutine trajout(x,y,z,time_step)
   use mod_const, only: ANG
   use mod_files, only: UMOVIE
   use mod_general,  only: imini, nwalk, natom, iremd, my_rank, sim_time
   use mod_system,   only: names
   implicit none
   real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
   integer,intent(in) :: time_step
   integer            :: iat,iw
   character(len=20)  :: fgeom, chout
   logical            :: lopened

   ! close has no effect if UMOVIE is not open according to the standard
   if (imini.eq.time_step) close(UMOVIE) 

   INQUIRE(UMOVIE,opened=lopened)

   if(.not.lopened)then
      if (imini.gt.time_step)then
         chout='movie_mini.xyz'
      else
         chout='movie.xyz'
      end if
      if(iremd.eq.1) write(chout, '(A,I2.2)')trim(chout)//'.', my_rank
      open(UMOVIE, file=chout, access='append', action="write")
   end if

!  printing with slightly lower precision for saving space
!  could be probably much lower
10 format(A2,3E18.8E2)

   do iw=1,nwalk
       write(UMOVIE,*)natom
       ! In the future, we should get rid of time step
       write(UMOVIE,'(A10,I20,A15,F15.2)')'Time step:',time_step,' Sim. Time [au]',sim_time
       do iat=1,natom
          write(UMOVIE,10)names(iat), x(iat,iw)/ANG, y(iat,iw)/ANG, z(iat,iw)/ANG
       enddo
   enddo

   end subroutine trajout

   subroutine forceout(x, y, z, fx,fy,fz,fUNIT)
   use mod_general, only: nwalk, natom, it
   use mod_system, ONLY: names
   real(DP),intent(in) :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(in) :: fx(:,:),fy(:,:),fz(:,:)
   integer, intent(in) :: funit
   real(DP)  :: fx_tot, fy_tot, fz_tot
   real(DP)  :: fx_rot, fy_rot, fz_rot
   integer   :: iat, iw
   character(len=40)  :: fgeom, fkom

   ! Calculate net translational and rotational gradient (should be close to zero)
   fx_tot = 0.0d0; fy_tot=0.0d0; fz_tot=0.0d0
   fx_rot = 0.0d0; fy_rot=0.0d0; fz_rot=0.0d0

   do iw=1,nwalk
      do iat=1,natom
         fx_tot = fx_tot + fx(iat, iw)
         fy_tot = fy_tot + fy(iat, iw)
         fz_tot = fz_tot + fz(iat, iw)
         fx_rot = fx_rot + fy(iat, iw)*z(iat, iw) - fz(iat,iw)*y(iat,iw)
         fy_rot = fy_rot + fz(iat, iw)*x(iat, iw) - fx(iat,iw)*z(iat,iw)
         fz_rot = fz_rot + fx(iat, iw)*y(iat, iw) - fy(iat,iw)*x(iat,iw)
      enddo
   enddo
   
   fgeom='(A2,3E18.10E2)'
   fkom='(A10,3E13.5E2,A14,3E13.5E2)'

   write(funit,*)natom
   write(funit,fkom)'net force:',fx_tot,fy_tot,fz_tot,' torque force:',fx_rot,fy_rot,fz_rot
   do iw=1,nwalk
      do iat=1,natom
         ! Printing in atomic units
         write(funit,fgeom)names(iat),fx(iat,iw),fy(iat,iw),fz(iat,iw)
      enddo
   enddo

   end subroutine forceout


   subroutine velout(vx,vy,vz)
   use mod_general, only: nwalk, natom, it
   use mod_system,  only: names
   use mod_files,   only: UVELOC
   real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
   integer :: iat,iw
   character(len=20)  :: fgeom
   
   
   fgeom='(A2,3E18.10E2)'
   write(UVELOC,*)natom
   write(UVELOC,*)'Time step:',it

   do iw=1,nwalk
      do iat=1,natom
         ! Printing in atomic units
         write(UVELOC,fgeom)names(iat),vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
   enddo

   end subroutine velout

   subroutine restout(x,y,z,vx,vy,vz,time_step)
   use mod_general,only: icv, ihess, nwalk, ipimd, natom, &
                         iremd, my_rank, pot, narchive, sim_time
   use mod_utils,  only: archive_file
   use mod_nhc,    only: inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
   use mod_estimators
   use mod_kinetic,only: entot_cumul, est_temp_cumul
   use mod_sh,     only: write_nacmrest,cel_re,cel_im,ntraj,nstate,istate
   use mod_gle
   use mod_random
   use mod_terampi_sh, only: write_wfn
   real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(in)  :: vx(:,:),vy(:,:),vz(:,:)
   integer,intent(in) :: time_step
   integer :: iat,iw,inh,itrj,ist1,is
   LOGICAL :: file_exists
   character(len=200)    :: chout, chsystem

   if(pot.eq.'_tera_'.and.ipimd.eq.2) call write_wfn()

   if(iremd.eq.1)then
      write(chout, '(A,I2.2)')'restart.xyz.', my_rank
   else
      chout='restart.xyz'
   end if

   INQUIRE(FILE=chout, EXIST=file_exists)
   chsystem='cp '//trim(chout)//'  '//trim(chout)//'.old'
   if(file_exists) call system(chsystem)

!  open(102, file=chout, action='WRITE',recl=250)
!  intel compilers don't write too many columns on single line
   open(102, file=chout, action='WRITE')

   write(102,*)time_step, sim_time

   write(102,*)chcoords
   do iw=1,nwalk
      do iat=1,natom
         write(102,*)x(iat,iw),y(iat,iw),z(iat,iw)
      enddo
   enddo

   write(102,*)chvel

   do iw=1,nwalk
      do iat=1,natom
         write(102,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
   enddo

   if(ipimd.eq.2)then
      call write_nacmrest()
      write(102,*)chSH
      do itrj=1,ntraj
         write(102,*)istate(itrj)
         do ist1=1,nstate
            write(102,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
         enddo
      enddo
   endif

   if(inose.eq.1)then
      write(102,*)chnose

      if(imasst.eq.1)then
         do inh=1,nchain
            do iw=1,nwalk
               do iat=1,natom
                  write(102,*)pnhx(iat,iw,inh),pnhy(iat,iw,inh),pnhz(iat,iw,inh)
               enddo
            enddo
         enddo

      else

         do inh=1,nchain
            do iw=1,nwalk
               do iat=1,nmolt
                  write(102,*)pnhx(iat,iw,inh)
               enddo
            enddo
         enddo
      endif

   endif

   if(inose.eq.2)then
      write(102,*)chQT
      do iw=1,nwalk
         do iat=1,natom*3
            write(102,*)(ps(iat,is,iw),is=1,ns)
         enddo
      enddo
      write(102,*)langham
   endif
   if(inose.eq.3)then
      write(102,*)chLT
      write(102,*)langham
   endif

   write(102,*)chAVG
   write(102,*)est_temp_cumul
   write(102,*)est_prim_cumul,est_vir_cumul
   write(102,*)entot_cumul

   if(icv.eq.1)then
      write(102,*)est_prim2_cumul,est_prim_vir,est_vir2_cumul
      write(102,*)cv_prim_cumul,cv_vir_cumul
      if(ihess.eq.1)then 
         write(102,*)cv_dcv_cumul
         do iw=1,nwalk
            write(102,*)cvhess_cumul(iw)
         enddo
      endif
   endif

   ! write current state of PRNG
   call rsavef(102)

   close(102)

   if(modulo(time_step,narchive).eq.0)then
      call archive_file('restart.xyz', time_step)
   end if

   end subroutine restout


   ! Subroutine that reads from restart.xyz during restart
   ! It is called from subroutine init.
   subroutine restin(x, y, z, vx, vy, vz, it)
   use mod_general,  only: icv, ihess, nwalk, ipimd, natom, &
                         iremd, my_rank, pot,sim_time
   use mod_nhc,      only: readNHC,inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
   use mod_estimators
   use mod_kinetic,  only: entot_cumul, est_temp_cumul
   use mod_sh,       only: write_nacmrest,cel_re,cel_im,ntraj,nstate,istate
   use mod_gle
   use mod_random
   use mod_terampi_sh, only: read_wfn
   real(DP),intent(out)  :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out)  :: vx(:,:),vy(:,:),vz(:,:)
   integer,intent(out)   :: it
   integer :: iat,iw,inh,itrj,ist1,is
   character(len=100) :: chtemp
   logical :: prngread
   character(len=20)  :: chin


   prngread=.false. 

   if(iremd.eq.1)then
      write(chin, '(A,I2.2)')'restart.xyz.', my_rank
   else
      chin='restart.xyz'
   end if

   if (my_rank.eq.0)then
      write(*,*)'irest=1, Reading geometry, velocities and other information from restart.xyz.'
   end if

   open(111,file=chin,status = "OLD", action = "READ")
   read(111,*)it, sim_time
   read(111,'(A)')chtemp
   call checkchar(chtemp, chcoords)
   do iw=1,nwalk
      do iat=1,natom
         read(111,*)x(iat,iw),y(iat,iw),z(iat,iw)
      enddo
   enddo

   read(111,'(A)')chtemp
   call checkchar(chtemp, chvel)
   do iw=1,nwalk
      do iat=1,natom
         read(111,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
   enddo

   if(ipimd.eq.2)then
      read(111,'(A)')chtemp
      call checkchar(chtemp, chsh)
      do itrj=1,ntraj
         read(111,*)istate(itrj)
         do ist1=1,nstate
            read(111,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
         enddo
      enddo
   endif

   if(inose.eq.1.and.readNHC.eq.1)then
      read(111,'(A)')chtemp
      call checkchar(chtemp, chnose)
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

   if(inose.eq.2.and.readQT.eq.1)then
      read(111,'(A)')chtemp
      call checkchar(chtemp, chqt)
      do iw=1,nwalk
         do iat=1,natom*3
            read(111,*)(ps(iat, is, iw),is=1, ns)
         enddo
      enddo
      read(111,*)langham
   endif

   if(inose.eq.3.and.readQT.eq.1)then
      read(111,'(A)')chtemp
      call checkchar(chtemp, chLT)
      read(111,*)langham
   endif

!- reading cumulative averages of various estimators
   read(111,'(A)')chtemp
   call checkchar(chtemp, chavg)
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

!- Trying to restart PRNG
!- prngread is optional argument determining, whether we write or read
!- currently,prngread is not used, since vranf is initialize BEFORE restart
!  and is possibly rewritten here
   call rsavef(111,prngread) 
    
   close(111)

   if(pot.eq.'_tera_'.and.ipimd.eq.2) call read_wfn()

   contains

   subroutine checkchar(chin, chref)
      use mod_utils,  only: abinerror
      character(len=*) :: chin, chref

      if(trim(adjustl(chin)).ne.trim(chref))then
         write(*,*)'ERROR while reading from restart.xyz.'
         write(*,*)'I read ',chin
         write(*,*)'but expected ',chref
         call abinerror('restin')
      end if
   end subroutine checkchar

   end subroutine restin

end module mod_analysis


