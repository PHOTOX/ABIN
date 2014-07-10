!----Initial version                    by Daniel Hollas,9.2.2012

!----Contains all analysis stuff
    subroutine analysis(x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
     use mod_array_size
     use mod_analyze_ext, only:analyze_ext
     use mod_estimators ,only:estimators
     use mod_general
     use mod_system
     use mod_density
     implicit none
     !intent inout because of estimators, writing to nwalk+1
     real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
     real*8,intent(in) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
     real*8,intent(in) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
     real*8,intent(in) :: amt(npartmax,nwalkmax)
     real*8,intent(in) :: eclas,equant
     real*8 :: dt  !,energy

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
      

       end



     subroutine trajout(x,y,z,it)
     use mod_array_size
     use mod_general, only: imini,nwalk,natom
     use mod_system, ONLY: names
     implicit none
     real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
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


     end

     subroutine velout(vx,vy,vz)
     use mod_array_size
     use mod_general
     implicit none
     real*8,intent(in) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
     integer :: iat,iw
     

     write(13,*)'Time step:',it
     do iw=1,nwalk
      do iat=1,natom
       write(13,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
     enddo


     end

     subroutine restout(x,y,z,vx,vy,vz,it)
     use mod_array_size
     use mod_general,only:icv,ihess,nwalk,ipimd,natom
     use mod_nhc
     use mod_estimators
     use mod_kinetic, only: entot_cumul, est_temp_cumul
     use mod_sh,only:cel_re,cel_im,ntraj,nstate,istate
     use mod_gle
     use mod_random
     implicit none
     real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
     real*8,intent(in)  :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
     integer,intent(in) :: it
     integer :: iat,iw,inh,itrj,ist1,is
     LOGICAL :: file_exists

!---NOTE: trajectories after restart are not precisely the same as without
!restart, since we don't pass fxc and fxq. Therefore,first step after
!restart is incorrect (see respa.f90).
! I think the above argument is invalid, as we now always do the zero step forces

     INQUIRE(FILE='restart.xyz', EXIST=file_exists)
     if(file_exists) call system('cp restart.xyz restart.xyz.old')

     open(102,file='restart.xyz',action='WRITE')

     write(102,*)it

     write(102,*)'Cartesian Coordinates [au]'
     do iw=1,nwalk
      do iat=1,natom
       write(102,*)x(iat,iw),y(iat,iw),z(iat,iw)
      enddo
     enddo

     write(102,*)'Cartesian Velocities [au]'

     do iw=1,nwalk
      do iat=1,natom
       write(102,*)vx(iat,iw),vy(iat,iw),vz(iat,iw)
      enddo
     enddo

     if(ipimd.eq.2)then
     write(102,*)'Coefficients for SH'
      do itrj=1,ntraj
       write(102,*)istate(itrj)
       do ist1=1,nstate
        write(102,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
       enddo
      enddo
     endif

     if(inose.eq.1)then
      write(102,*)'NHC momenta'
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
      write(102,*)'Quantum Thermostat'
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

     write(102,*)'Cumulative averages of various estimators'
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

     end


