!----Initial version                    by Daniel Hollas,9.2.2012
!- This module contains some routines that do analyses and do I/O operations.
!- Crucially, it also contains routines performing restart.

module mod_analysis
   use mod_const, only: DP
   implicit none
   private 
   public :: trajout, restout, analysis, restin
   !These are the character string that we check for during restart.
   character(len=*),parameter :: chnose='NHC momenta',chQT='Quantum Thermostat', &
            chSH='Coefficients for SH', chAVG='Cumulative averages of various estimators', &
            chcoords='Cartesian Coordinates [au]', chvel='Cartesian Velocities [au]'

   contains

!----Contains all analysis stuff
   subroutine analysis(x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
     use mod_analyze_ext, only:analyze_ext
     use mod_estimators ,only:estimators
     use mod_general
     use mod_system
     use mod_density
     implicit none
     !intent inout because of estimators, writing to nwalk+1
     real(DP),intent(inout) :: x(:,:),y(:,:),z(:,:)
     real(DP),intent(in) :: fxc(:,:),fyc(:,:),fzc(:,:)
     real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
     real(DP),intent(in) :: amt(:,:)
     real(DP),intent(in) :: eclas,equant
     real(DP) :: dt  !,energy

!     eclas comes from force_clas,equant from force_quantum
!     energy=eclas+equant

      if (ipimd.eq.1.or.icv.eq.1)then
       call estimators(x,y,z,fxc,fyc,fzc,eclas,dt)
      endif

      if(modulo(it,nwritex).eq.0)then
       call trajout(x,y,z,it)
      endif

      if(nwritev.gt.0)then
       if(modulo(it,nwritev).eq.0)then
        call velout(vx,vy,vz)
       endif
      endif

      if(ndist.ge.1.and.it.gt.imini)then
       call density(x,y,z)
      endif

      if(nang.ge.1.and.it.gt.imini)then
       call density_ang(x,y,z)
      endif

      if(ndih.ge.1.and.it.gt.imini)then
       call density_dih(x,y,z)
      endif

      if((modulo(it,nrest).eq.0).or.it.eq.nstep)then
       call restout(x,y,z,vx,vy,vz,it)
      endif

      if (anal_ext.eq.1)then
       call analyze_ext(x,y,z,vx,vy,vz,amt)
      endif
      

   end subroutine analysis

   subroutine trajout(x,y,z,it)
     use mod_const, only: ANG
     use mod_general, only: imini,nwalk,natom
     use mod_system, ONLY: names
     implicit none
     real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
     integer,intent(in) :: it
     integer            :: iat,iw
     character(len=20)  :: fgeom
     
     if(it.le.imini)then
           open(101,file='movie_mini.xyz',access='append')
     else
           open(101,file='movie.xyz',access='append')
     endif
!printing with slightly lower precision for saving space
     fgeom='(A2,3E18.10E2)'

     do iw=1,nwalk
      write(101,*)natom
      write(101,*)'Time step:',it
      do iat=1,natom
       write(101,fgeom)names(iat),x(iat,iw)/ang,y(iat,iw)/ang,z(iat,iw)/ang
      enddo
     enddo

     close(101)


   end subroutine trajout

   subroutine velout(vx,vy,vz)
     use mod_general, only: nwalk, natom, it
     real(DP),intent(in) :: vx(:,:),vy(:,:),vz(:,:)
     integer :: iat,iw
     

     write(13,*)'Time step:',it
     do iw=1,nwalk
      do iat=1,natom
       write(13,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
     enddo


   end subroutine velout

   subroutine restout(x,y,z,vx,vy,vz,it)
     use mod_general,only:icv, ihess, nwalk, ipimd, natom
     use mod_nhc,    only: inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
     use mod_estimators
     use mod_kinetic,only: entot_cumul, est_temp_cumul
     use mod_sh,     only: write_nacmrest,cel_re,cel_im,ntraj,nstate,istate
     use mod_gle
     use mod_random
     real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
     real(DP),intent(in)  :: vx(:,:),vy(:,:),vz(:,:)
     integer,intent(in) :: it
     integer :: iat,iw,inh,itrj,ist1,is
     LOGICAL :: file_exists


     INQUIRE(FILE='restart.xyz', EXIST=file_exists)
     if(file_exists) call system('cp restart.xyz restart.xyz.old')

     open(102,file='restart.xyz',action='WRITE')

     write(102,*)it

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
      if(nwalk.eq.1)then
       do iat=1,natom*3
        write(102,*)(gp(iat,is),is=2,ns+1)
       enddo
      else
       do iw=1,nwalk
        do iat=1,natom*3
         write(102,*)(ps(iat,is,iw),is=1,ns)
        enddo
       enddo
     endif
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

     call rsavef(102)

     close(102)


   end subroutine restout


   ! Subroutine that reads from restart.xyz during restart
   ! It is called from subroutine init.
   subroutine restin(x,y,z,vx,vy,vz,it)
     use mod_general,only: icv, ihess, nwalk, ipimd, natom
     use mod_nhc,    only: readNHC,inose, pnhx, pnhy, pnhz, imasst, nmolt, nchain
     use mod_estimators
     use mod_kinetic,only: entot_cumul, est_temp_cumul
     use mod_sh,     only: write_nacmrest,cel_re,cel_im,ntraj,nstate,istate
     use mod_gle
     use mod_random
     real(DP),intent(out)  :: x(:,:),y(:,:),z(:,:)
     real(DP),intent(out)  :: vx(:,:),vy(:,:),vz(:,:)
     integer,intent(out)   :: it
     integer :: iat,iw,inh,itrj,ist1,is
     character(len=100) :: chtemp
     logical :: prngread

     prngread=.false. 

     write(*,*)'irest=1, Reading geometry, velocities and other information from restart.xyz.'

     open(111,file='restart.xyz',status = "OLD", action = "READ")
     read(111,*)it
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
      if(nwalk.eq.1)then
       do iat=1,natom*3
        read(111,*)(gp(iat,is),is=2,ns+1)
       enddo
      else
       do iw=1,nwalk
        do iat=1,natom*3
         read(111,*)(ps(iat,is,iw),is=1,ns)
        enddo
       enddo
     endif
      read(111,*)langham
     endif

!-   reading cumulative averages of various estimators
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

!-   Trying to restart PRNG
!-   prngread is optional argument determining, whether we write or read
!-   currently,prngread is not used, since vranf is initialize BEFORE restart
!    and is possibly rewritten here
     call rsavef(111,prngread) 
      
     close(111)

     contains

     subroutine checkchar(chin,chref)
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


