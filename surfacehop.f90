
module mod_sh
use mod_array_size
implicit none
integer :: istate_init=1,nstate=1,ntraj=1,substep=10000
integer :: inac=0,nohop=0,nac_accu1=7,nac_accu2=5 !7 is MOLPRO default
real*8  :: dtp,alpha=0.1d0,eshift,deltae=100.,popthr=-1
integer :: istate(ntrajmax)
real*8  :: nacx(npartmax,ntrajmax,nstmax,nstmax)
real*8  :: nacy(npartmax,ntrajmax,nstmax,nstmax)
real*8  :: nacz(npartmax,ntrajmax,nstmax,nstmax)
real*8  :: dotproduct_old(nstmax,nstmax,ntrajmax)=0.0d0 !for inac=1
real*8  :: dotproduct_new(nstmax,nstmax,ntrajmax)=0.0d0 !for inac=1
real*8  :: en_array(nstmax,ntrajmax)
real*8  :: cel_re(nstmax,ntrajmax),cel_im(nstmax,ntrajmax)
integer :: tocalc(nstmax,nstmax)
character(len=10) :: integ='rk4'
save
contains
subroutine set_tocalc()
integer :: ist1,ist2,itrj=1
real*8  :: pop,pop2

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
end subroutine

integer function readnacm(itrj)
use mod_qmmm,only:natqm
implicit none
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
end function

subroutine calcnacm(itrj)
use mod_general, only: it,pot
implicit none
integer :: ist1,ist2,itrj
character(len=100) :: chsystem
open(unit=510,file='state.dat')
write(510,'(I2)')istate(itrj)
write(510,'(I2)')nstate
! horni troj. matice bez diagonaly
! tocalc(,)=1 -> pocitame couplingy
! tocalc(,)=0 -> nepocitame couplingy
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
   write(*,*)'Different accuracy for NACME is currently supported only by molpro.'
   write(*,*)'Exiting...'
   stop 1
endif

call system(chsystem)
end subroutine

subroutine move_vars(en_array_old,nacx_old,nacy_old,nacz_old,vx,vy,vz,vx_old,vy_old,vz_old,itrj)
use mod_general,only:natom
implicit none
real*8,intent(out) :: en_array_old(nstmax,ntrajmax)
real*8,intent(out) :: vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
real*8,intent(out) :: nacx_old(npartmax,ntrajmax,nstmax,nstmax)
real*8,intent(out) :: nacy_old(npartmax,ntrajmax,nstmax,nstmax)
real*8,intent(out) :: nacz_old(npartmax,ntrajmax,nstmax,nstmax)
real*8,intent(in)  :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
integer :: itrj,ist1,ist2,iat
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
         dotproduct_old(ist1,ist2,itrj)=dotproduct_new(ist1,ist2,itrj)
      enddo
   enddo
endif

end subroutine

end module

      subroutine surfacehop(x,y,z,vx,vy,vz,nacx_old,nacy_old,nacz_old,vx_old,vy_old,vz_old,en_array_old,dt)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: am,names
      use mod_sh
      use mod_qmmm, ONLY:natqm
      implicit none
      real*8,intent(in)    :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(inout) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8,intent(inout) :: vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
      real*8,intent(inout) :: en_array_old(nstmax,ntrajmax)
      real*8,intent(inout) :: nacx_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8,intent(inout) :: nacy_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8,intent(inout) :: nacz_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: vx_int(npartmax,nwalkmax),vy_int(npartmax,nwalkmax),vz_int(npartmax,nwalkmax)
      real*8  :: vx_newint(npartmax,nwalkmax),vy_newint(npartmax,nwalkmax),vz_newint(npartmax,nwalkmax)
      real*8  :: en_array_int(nstmax,ntrajmax),en_array_newint(nstmax,ntrajmax)
      real*8  :: dcel_re(nstmax,ntrajmax),dcel_im(nstmax,ntrajmax)
      real*8  :: ancx(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: ancy(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: ancz(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: ancx_newint(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: ancy_newint(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: ancz_newint(npartmax,ntrajmax,nstmax,nstmax)
      real*8  :: dotproduct(nstmax,nstmax,ntrajmax),dotproduct_newint(nstmax,nstmax,ntrajmax)
      real*8  :: t(nstmax,nstmax)           !switching probabilities
      real*8  :: t_tot(nstmax,nstmax)       !switching probabilities
      real*8  :: ran(10)
      real*8  ::  pop(nstmax,ntrajmax),popsum
      integer :: itp
      integer :: iat,ist1,ist2,itrj     !iteration counters
      integer :: ist                    ! =istate(itrj)
      real*8  :: vect_olap,fr,frd,dt,temp1
      real*8  :: ekin_mom,apom,edif,tau,fact,sum_norm
      integer :: ihop,ijunk
      real*8  :: a_re,prob(nstmax),cn
      character(len=500) :: formt
      character(len=20) :: chist,chihop,chit
      integer :: iost


     do itrj=1,ntraj


      do ist1=1,nstate
       do ist2=1,nstate
        t_tot(ist1,ist2)=0.0
       enddo
      enddo

      popsum=0.0d0
!DEBUG
!      if(idebug.eq.1)then
!       call printf(vx_old,vy_old,vz_old)
!       write(*,*)(en_array_old(ist1,itrj),ist1=1,nstate)
!      endif
!-------READING NACM-----------------------     
     if(inac.eq.0)then

      do ist1=1,nstate-1
       do ist2=ist1+1,nstate
       if(tocalc(ist1,ist2).eq.0)then
          write(*,*)'Not computing NACM(tocalc=0) for states',ist1,ist2
          do iat=1,natqm        ! MUSIME NULOVAT UZ TADY,bo pak menime tocalc behem cteni
           nacx(iat,itrj,ist1,ist2)=0.0
           nacy(iat,itrj,ist1,ist2)=0.0
           nacz(iat,itrj,ist1,ist2)=0.0
           nacx(iat,itrj,ist2,ist1)=0.0
           nacy(iat,itrj,ist2,ist1)=0.0
           nacz(iat,itrj,ist2,ist1)=0.0
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
        write(*,*)'ERROR:Some NACs not read. Exiting...'
        stop 1
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

      if(idebug.eq.1)then
      open(17,file='debug.nacm',access='append')
      write(*,*)'Time step:',it
        do ist1=1,nstate
         do ist2=1,nstate
          do iat=1,natom
          write(17,'(3I3,3E20.8)')iat,ist1,ist2,nacx(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2)
         enddo
         enddo
        enddo
      close(17)
      endif

!----------- INAC=1  endif
     endif

!------READING time-derivative couplings----------------------------------
      if ( inac.eq.1 ) then

      open(100,file='tdcoups.dat')
      read(100,*)
      read(100,*)
       do ist1=1,nstate
        read(100,*)ijunk,(dotproduct_new(ist1,ist2,itrj),ist2=1,nstate)
       enddo
      close(100)

       do ist1=1,nstate
        do ist2=1,nstate
         dotproduct_new(ist1,ist2,itrj)=-dotproduct_new(ist1,ist2,itrj)/dt
         if(ist1.eq.ist2)  dotproduct_new(ist1,ist2,itrj)=0.0d0
        enddo
       enddo

      do ist1=1,nstate-1
       do ist2=ist1+1,nstate
        if(tocalc(ist1,ist2).eq.0)then
         write(*,*)'Not computing NACM(tocalc=0) for states',ist1,ist2
         dotproduct_new(ist1,ist2,itrj)=0.0d0
         dotproduct_new(ist2,ist1,itrj)=0.0d0
        endif
       enddo
      enddo


!v prvnim korku nemame mezi cim interpolovat..bereme nulty krok jako prvni
       if( it.eq.1) dotproduct_old=dotproduct_new

       endif
!------------------END-OF-TDC-------------------------------

! Smaller time step for electron population transfer
      do itp=1,substep      

      ist=istate(itrj)

!-----------------INTERPOLACE--------------------------------------------------


      if(integ.ne.'euler')then 
        fr=float(itp)/float(substep)
        frd=1.0d0-fr

        call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_newint,vy_newint,vz_newint, &
                      ancx_newint,ancy_newint,ancz_newint,nacx_old,nacy_old,nacz_old,en_array_newint,en_array_old, &
                      dotproduct_newint,fr,frd,itrj)

       if(inac.eq.1)then
        call interpolate_dot(dotproduct_newint,fr,frd,itrj)
       endif

      endif

       fr=float(itp-1)/float(substep)
       frd=1.0d0-fr

       call interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_int,vy_int,vz_int, &
                      ancx,ancy,ancz,nacx_old,nacy_old,nacz_old,en_array_int,en_array_old, &
                      dotproduct,fr,frd,itrj)

       if(inac.eq.1)then
        call interpolate_dot(dotproduct,fr,frd,itrj)
       endif

!      if(idebug.eq.1)then
!      open(100,file='debug.v',access='append')
!      open(101,file='debug.en',access='append')
!      open(102,file='debug.ancm',access='append')
!      write(100,*)'HALO it inp',it,itp
!      write(101,*)'HALO it inp',it,itp
!      write(102,*)'HALO it inp',it,itp
!      do iat=1,natom
!      write(100,*)vx_old(iat,itrj),vx(iat,itrj),vx_int(iat,itrj)
!      write(100,*)vy_old(iat,itrj),vy(iat,itrj),vy_int(iat,itrj)
!      write(100,*)vz_old(iat,itrj),vz(iat,itrj),vz_int(iat,itrj)
!      enddo
!      do ist1=1,nstate
!      write(101,*)en_array_old(ist1,itrj),en_array(ist1,itrj),en_array_int(ist1,itrj)
!      do ist2=1,nstate
!      write(102,*)'state1 state2',ist1,ist2
!      do iat=1,natom
!      write(102,*)nacx_old(iat,itrj,ist1,ist2),nacx(iat,itrj,ist1,ist2),ancx(iat,itrj,ist1,ist2)
!      write(102,*)nacy_old(iat,itrj,ist1,ist2),nacy(iat,itrj,ist1,ist2),ancy(iat,itrj,ist1,ist2)
!      write(102,*)nacz_old(iat,itrj,ist1,ist2),nacz(iat,itrj,ist1,ist2),ancz(iat,itrj,ist1,ist2)
!      enddo
!      enddo
!      enddo
!      close(100)
!      close(101)
!      close(102)
!      endif

      !------END-OF-INTERPOLATIONS-------------------------------------

       
       if(integ.eq.'rk4')then
        call rk4step(en_array_int,en_array_newint,dotproduct,dotproduct_newint,itrj)
        !interpolujeme dotprodukt poctive i uvnitr integratoru...
!        call rk4step_new(en_array_int,en_array_newint,dotproduct,dotproduct_newint, &
!                      vx_int,vy_int,vz_int,vx_newint,vy_newint,vz_newint,  &
!                      ancx,ancy,ancz,ancx_newint,ancy_newint,ancz_newint,itrj)
       endif

       if(integ.eq.'butcher')then
        call butcherstep(en_array_int,en_array_newint,dotproduct,dotproduct_newint,itrj)
       endif

       if(integ.eq.'euler')then

      do ist1=1,nstate

      dcel_re(ist1,itrj)=(en_array_int(ist1,itrj)+eshift)*cel_im(ist1,itrj)
      dcel_im(ist1,itrj)=-(en_array_int(ist1,itrj)+eshift)*cel_re(ist1,itrj)
      
       do ist2=1,nstate

       ! dotprodukt pro stejne stavy by mel byt nulovy       
!         if(ist1.ne.ist2)then
         dcel_re(ist1,itrj)=dcel_re(ist1,itrj)-cel_re(ist2,itrj)*dotproduct(ist1,ist2,itrj)
         dcel_im(ist1,itrj)=dcel_im(ist1,itrj)-cel_im(ist2,itrj)*dotproduct(ist1,ist2,itrj)
!         endif
  
       enddo
      enddo

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+dcel_re(ist1,itrj)*dtp
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+dcel_im(ist1,itrj)*dtp
      enddo

      endif

!----calculation of switching probabilities(asi nemusime pocitat celou matici)      
      do ist1=1,nstate
       do ist2=1,nstate
        a_re=(cel_re(ist1,itrj)*cel_re(ist2,itrj)+cel_im(ist1,itrj)*cel_im(ist2,itrj))
        t(ist1,ist2)=2*a_re*dotproduct(ist1,ist2,itrj)
       enddo
      enddo

      do ist1=1,nstate
       pop(ist1,itrj)=cel_re(ist1,itrj)**2+cel_im(ist1,itrj)**2
       do ist2=1,nstate
         t(ist1,ist2)=t(ist1,ist2)*dtp/(pop(ist1,itrj))
        if(t(ist1,ist2).lt.0.0d0)then
         t(ist1,ist2)=0.0d0
        endif
        t_tot(ist1,ist2)=t_tot(ist1,ist2)+t(ist1,ist2)
       enddo
      enddo


      do ist1=1,nstate
       prob(ist1)=0.0d0
      enddo

! Auxiliary calculations of prob
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
       if(inac.eq.0) call hop(vx,vy,vz,vx_int,vy_int,vz_int,ancx,ancy,ancz,en_array_int,ist,ihop,itrj)
       if(inac.eq.1) call hop_dot(vx,vy,vz,ist,ihop,itrj)
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
       do iat=1,natom
        temp1=vx_int(iat,itrj)**2+vy_int(iat,itrj)**2+vz_int(iat,itrj)**2
        temp1=0.5*temp1*am(iat)
        ekin_mom=ekin_mom+temp1
       enddo

      if(ekin_mom.gt.1.0d-4)then

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
      fact=sum_norm/(cel_re(istate(itrj),itrj)**2+cel_im(istate(itrj),itrj)**2+1.0d-7)
      fact=dsqrt(fact)
!     write(148,*)'renomr',fact,istate(itrj)

      cel_re(istate(itrj),itrj)=cel_re(istate(itrj),itrj)*fact
      cel_im(istate(itrj),itrj)=cel_im(istate(itrj),itrj)*fact

      !ekin endif
      endif

!--------END-OF-DECOHERENCE------------------------      
      endif




!itp loop
      enddo

!WARNING---------      
      irest=0

! MO for SH: calculate NACM only for given deltaE
     call set_tocalc()
     call move_vars(en_array_old,nacx_old,nacy_old,nacz_old,vx,vy,vz,vx_old,vy_old,vz_old,itrj)

     if(idebug.eq.1)then
      if(it.eq.1)then
       open(150,file='dotprod.dat')
      else
       open(150,file='dotprod.dat',access='append')
      endif
      write(150,*)'Step: ',it
      do ist1=1,nstate
       write(150,*)(dotproduct(ist1,ist2,itrj),ist2=1,nstate)
      enddo
      close(150)
     endif


      do ist1=1,nstate
       popsum=popsum+pop(ist1,itrj)
      enddo

      if(modulo(it,nwrite).eq.0)then
       write(formt,'(A10,I3,A13)')'(F15.2,I3,',nstate,'F10.5,1F10.7)'
       write(3,fmt=formt)it*dt*autofs,istate(itrj),(pop(ist1,itrj), ist1=1,nstate),popsum
       write(formt,'(A10,I3,A6)')'(F15.2,I3,',nstate,'F10.5)'
       write(4,fmt=formt)it*dt*autofs,istate(itrj),(t_tot(ist,ist1),ist1=1,nstate)
       write(formt,'(A7,I3,A7)')'(F15.2,',nstate,'E20.10)'
       write(8,fmt=formt)it*dt*autofs,(en_array(ist1,itrj),ist1=1,nstate)
      endif

       ! ntraj enddo       
       enddo

       end



      subroutine hop(vx,vy,vz,vx_int,vy_int,vz_int,ancx,ancy,ancz,en_array_int,state1,state2,itrj)
      use mod_array_size
      use mod_general, ONLY:natom
      use mod_system, ONLY: am
      use mod_sh
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 vx_int(npartmax,nwalkmax),vy_int(npartmax,nwalkmax),vz_int(npartmax,nwalkmax)
      real*8 en_array_int(nstmax,ntrajmax)
      real*8 ancx(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancy(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancz(npartmax,ntrajmax,nstmax,nstmax)
      integer :: itrj,state1,state2
      integer :: iat
      real*8  :: a_temp,b_temp,c_temp,g_temp

!  Checking for frustrated hop

      a_temp=0.
      b_temp=0.

      do iat=1,natom
        a_temp=a_temp+ancx(iat,itrj,state1,state2)**2/am(iat)
        a_temp=a_temp+ancy(iat,itrj,state1,state2)**2/am(iat)
        a_temp=a_temp+ancz(iat,itrj,state1,state2)**2/am(iat)
        b_temp=b_temp+ancx(iat,itrj,state1,state2)*vx_int(iat,itrj)
        b_temp=b_temp+ancy(iat,itrj,state1,state2)*vy_int(iat,itrj)
        b_temp=b_temp+ancz(iat,itrj,state1,state2)*vz_int(iat,itrj)
      enddo
      a_temp=0.5d0*a_temp
      c_temp=b_temp**2+4*a_temp*(en_array_int(state1,itrj)-en_array_int(state2,itrj))

      if(c_temp.lt.0)then
        write(3,'(A35,I3,A10,I3)')'#Frustrated Hop occured from state ',state1,' to state ',state2
        return
      endif

      istate(itrj)=state2
       write(3,'(A24,I3,A10,I3)')'#Hop occured from state ',state1,' to state ',state2

!------------- Rescaling the velocities------------------------

      if(b_temp.lt.0) then
       g_temp=(b_temp+dsqrt(b_temp**2+4*a_temp*(en_array_int(state1,itrj)-en_array_int(state2,itrj))))/2.0d0/a_temp
      endif

      if(b_temp.ge.0) then
       g_temp=(b_temp-dsqrt(b_temp**2+4*a_temp*(en_array_int(state1,itrj)-en_array_int(state2,itrj))))/2.0d0/a_temp
      endif

      write(*,*)a_temp,b_temp,c_temp,g_temp

      !TODO: tohle asi neni uplne spravne..mame pouzit vx_int nebo vx?
      do iat=1,natom
       vx(iat,itrj)=vx(iat,itrj)-g_temp*ancx(iat,itrj,state1,state2)/am(iat)
       vy(iat,itrj)=vy(iat,itrj)-g_temp*ancy(iat,itrj,state1,state2)/am(iat)
       vz(iat,itrj)=vz(iat,itrj)-g_temp*ancz(iat,itrj,state1,state2)/am(iat)
      enddo

      end

      subroutine hop_dot(vx,vy,vz,state1,state2,itrj)
      use mod_array_size
      use mod_general, ONLY:natom,idebug
      use mod_sh, ONLY:en_array,istate
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      integer :: itrj,state1,state2,iat
      real*8  :: de,ekin,alfa,ekin_new
      real*8  :: ekin_v

      ekin=0.0
      ekin_new=0.0

      de=en_array(state2,itrj)-en_array(state1,itrj)
      ekin=ekin_v(vx,vy,vz)

      if(ekin.ge.de)then

        alfa=dsqrt(1-de/ekin)

        if(idebug.eq.1)then
         write(*,*)'Velocity before rescaling'
         call printf(vx,vy,vz)
         write(*,*)'alpha=',alfa
         write(*,*)'de=',de
         write(*,*)'ekin=',ekin
         write(*,*)'epot1=',en_array(state1,itrj)
         write(*,*)'epot2=',en_array(state2,itrj)
        endif

        do iat=1,natom
         vx(iat,itrj)=alfa*vx(iat,itrj) 
         vy(iat,itrj)=alfa*vy(iat,itrj) 
         vz(iat,itrj)=alfa*vz(iat,itrj) 
        enddo
        istate(itrj)=state2
        ekin_new=ekin_v(vx,vy,vz)

        if(idebug.eq.1)then
         write(*,*)'Velocity after rescaling'
         call printf(vx,vy,vz)
         write(*,*)'ekin_new=',ekin_new
        endif

         write(3,'(A24,I3,A10,I3)')'#Hop occured from state ',state1,' to state ',state2
         write(3,'(A,2E20.10)')'#TOT_Energy_old   TOT_Energy_new :',ekin+en_array(state1,itrj),ekin_new+en_array(state2,itrj)

      else

       write(3,'(A35,I3,A10,I3)')'#Frustrated Hop occured from state ',state1,' to state ',state2
       write(3,'(A31,2E20.10)')'deltaE-potential     Ekin-total',dE,ekin
       return

      endif

      end

      subroutine integstep(k_re,k_im,en,y_re,y_im,dotproduct)
      use mod_array_size
      use mod_sh,only:nstate,dtp
      real*8 k_re(nstmax),k_im(nstmax)
      real*8 dotproduct(nstmax,nstmax)
      real*8 en(nstmax),y_im(nstmax),y_re(nstmax)

      do ist1=1,nstate
       k_re(ist1)=en(ist1)*y_im(ist1)
       k_im(ist1)=-en(ist1)*y_re(ist1)
       do ist2=1,nstate
       k_re(ist1)=k_re(ist1)-y_re(ist2)*dotproduct(ist1,ist2)
       k_im(ist1)=k_im(ist1)-y_im(ist2)*dotproduct(ist1,ist2)
       enddo
       k_re(ist1)=dtp*k_re(ist1)
       k_im(ist1)=dtp*k_im(ist1)
      enddo

      end

      subroutine rk4step_new(en_array,en_array_new,dotproduct,dotproduct_new,vx,vy,vz,vx_new,vy_new,vz_new,&
                      nacx,nacy,nacz,nacx_new,nacy_new,nacz_new,itrj)
      use mod_array_size
      use mod_sh,only:nstate,cel_re,cel_im,eshift
      implicit none
      real*8 vx_new(npartmax,nwalkmax),vy_new(npartmax,nwalkmax),vz_new(npartmax,nwalkmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 nacx_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacx(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz(npartmax,ntrajmax,nstmax,nstmax)
      real*8 en_array(nstmax,ntrajmax)
      real*8 en_array_new(nstmax,ntrajmax)
      real*8 dotproduct(nstmax,nstmax,ntrajmax)
      real*8 dotproduct_new(nstmax,nstmax,ntrajmax)
      real*8 dotprod2(nstmax,nstmax)
      real*8 dotprod0(nstmax,nstmax),dotprod1(nstmax,nstmax)
      real*8 k1_re(nstmax),k1_im(nstmax)
      real*8 k2_re(nstmax),k2_im(nstmax)
      real*8 k3_re(nstmax),k3_im(nstmax)
      real*8 k4_re(nstmax),k4_im(nstmax)
      real*8 :: y_im(nstmax),y_re(nstmax)
      real*8 :: en0(nstmax),en1(nstmax),en2(nstmax)
      integer :: ist1,ist2,itrj     !iteration counters

!pripravne interpolace....
      do ist1=1,nstate
       en0(ist1)=en_array(ist1,itrj)+eshift
       en1(ist1)=en_array_new(ist1,itrj)+eshift
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
       dotprod0(ist1,ist2)=dotproduct(ist1,ist2,itrj)
       dotprod1(ist1,ist2)=dotproduct_new(ist1,ist2,itrj)
       enddo
      enddo

      call interpolate2(vx,vy,vz,vx_new,vy_new,vz_new, &
                      nacx,nacy,nacz,nacx_new,nacy_new,nacz_new,en0,en1,en2, &
                      dotprod2,0.5d0,0.5d0,itrj)

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/2
      enddo
      call integstep(k2_re,k2_im,en2,y_re,y_im,dotprod2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k2_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k2_im(ist1)/2
      enddo
      call integstep(k3_re,k3_im,en2,y_re,y_im,dotprod2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en1,y_re,y_im,dotprod1)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)/6+k2_re(ist1)/3+k3_re(ist1)/3+k4_re(ist1)/6
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)/6+k2_im(ist1)/3+k3_im(ist1)/3+k4_im(ist1)/6
      enddo

      end

      subroutine rk4step(en_array,en_array_new,dotproduct,dotproduct_new,itrj)
      use mod_array_size
      use mod_sh,only:nstate,cel_re,cel_im,eshift
      implicit none
      real*8 en_array(nstmax,ntrajmax)
      real*8 en_array_new(nstmax,ntrajmax)
      real*8 dotproduct(nstmax,nstmax,ntrajmax)
      real*8 dotproduct_new(nstmax,nstmax,ntrajmax)
      real*8 dotprod2(nstmax,nstmax)
      real*8 dotprod0(nstmax,nstmax),dotprod1(nstmax,nstmax)
      real*8 k1_re(nstmax),k1_im(nstmax)
      real*8 k2_re(nstmax),k2_im(nstmax)
      real*8 k3_re(nstmax),k3_im(nstmax)
      real*8 k4_re(nstmax),k4_im(nstmax)
      real*8 :: y_im(nstmax),y_re(nstmax)
      real*8 :: en0(nstmax),en1(nstmax),en2(nstmax)
      integer :: ist1,ist2,itrj     !iteration counters

!pripravne interpolace....
      do ist1=1,nstate
       en0(ist1)=en_array(ist1,itrj)+eshift
       en1(ist1)=en_array_new(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
       dotprod0(ist1,ist2)=dotproduct(ist1,ist2,itrj)
       dotprod1(ist1,ist2)=dotproduct_new(ist1,ist2,itrj)
       dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
       enddo
      enddo

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/2
      enddo
      call integstep(k2_re,k2_im,en2,y_re,y_im,dotprod2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k2_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k2_im(ist1)/2
      enddo
      call integstep(k3_re,k3_im,en2,y_re,y_im,dotprod2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en1,y_re,y_im,dotprod1)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)/6+k2_re(ist1)/3+k3_re(ist1)/3+k4_re(ist1)/6
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)/6+k2_im(ist1)/3+k3_im(ist1)/3+k4_im(ist1)/6
      enddo

      end


      subroutine butcherstep(en_array,en_array_new,dotproduct,dotproduct_new,itrj)
      use mod_array_size
      use mod_sh,only:nstate,cel_re,cel_im,eshift
      implicit none
      real*8 en_array(nstmax,ntrajmax)
      real*8 en_array_new(nstmax,ntrajmax)
      real*8 dotproduct(nstmax,nstmax,ntrajmax)
      real*8 dotproduct_new(nstmax,nstmax,ntrajmax)
      real*8 dotprod2(nstmax,nstmax),dotprod4(nstmax,nstmax),dotprod34(nstmax,nstmax)
      real*8 dotprod0(nstmax,nstmax),dotprod1(nstmax,nstmax)
      real*8 k1_re(nstmax),k1_im(nstmax)
      real*8 k2_re(nstmax),k2_im(nstmax)
      real*8 k3_re(nstmax),k3_im(nstmax)
      real*8 k4_re(nstmax),k4_im(nstmax)
      real*8 k5_re(nstmax),k5_im(nstmax)
      real*8 k6_re(nstmax),k6_im(nstmax)
      real*8 :: y_im(nstmax),y_re(nstmax)
      real*8 :: en0(nstmax),en1(nstmax),en2(nstmax),en4(nstmax),en34(nstmax)
      integer :: ist1,ist2,itrj     !iteration counters

!pripravne interpolace....
      do ist1=1,nstate
       en0(ist1)=en_array(ist1,itrj)+eshift
       en1(ist1)=en_array_new(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       en4(ist1)=(en0(ist1)+en2(ist1))/2
       en34(ist1)=(en2(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
        dotprod0(ist1,ist2)=dotproduct(ist1,ist2,itrj)
        dotprod1(ist1,ist2)=dotproduct_new(ist1,ist2,itrj)
        dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
        dotprod4(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod2(ist1,ist2))/2
        dotprod34(ist1,ist2)=(dotprod2(ist1,ist2)+dotprod1(ist1,ist2))/2
       enddo
      enddo

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/4
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/4
      enddo
      call integstep(k2_re,k2_im,en4,y_re,y_im,dotprod4)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/8+k2_re(ist1)/8
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/8+k2_im(ist1)/8
      enddo
      call integstep(k3_re,k3_im,en4,y_re,y_im,dotprod4)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-k2_re(ist1)/2+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)-k2_im(ist1)/2+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en2,y_re,y_im,dotprod2)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+3*k1_re(ist1)/16+9*k4_re(ist1)/16
       y_im(ist1)=cel_im(ist1,itrj)+3*k1_im(ist1)/16+9*k4_im(ist1)/16
      enddo
      call integstep(k5_re,k5_im,en34,y_re,y_im,dotprod34)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-3*k1_re(ist1)/7+2*k2_re(ist1)/7+12*k3_re(ist1)/7 &
                  -12*k4_re(ist1)/7+8*k5_re(ist1)/7
       y_im(ist1)=cel_im(ist1,itrj)-3*k1_im(ist1)/7+2*k2_im(ist1)/7+12*k3_im(ist1)/7 &
                  -12*k4_im(ist1)/7+8*k5_im(ist1)/7
      enddo
      call integstep(k6_re,k6_im,en1,y_re,y_im,dotprod1)


      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+7*k1_re(ist1)/90+32*k3_re(ist1)/90+12*k4_re(ist1)/90 &
                        +32*k5_re(ist1)/90+7*k6_re(ist1)/90
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+7*k1_im(ist1)/90+32*k3_im(ist1)/90+12*k4_im(ist1)/90 &
                        +32*k5_im(ist1)/90+7*k6_im(ist1)/90
      enddo

      end

      subroutine rk4step_old(en_array_int,dotproduct,dotproduct_newint,itrj)
      use mod_array_size
      use mod_general
      use mod_sh
      implicit none
      real*8 en_array_int(nstmax,ntrajmax)
      real*8 dotproduct(nstmax,nstmax,ntrajmax)
      real*8 dotproduct_newint(nstmax,nstmax,ntrajmax)
      real*8 dotproduct2(nstmax,nstmax,ntrajmax)
      real*8 k1_re(nstmax,ntrajmax),k1_im(nstmax,ntrajmax)
      real*8 k2_re(nstmax,ntrajmax),k2_im(nstmax,ntrajmax)
      real*8 k3_re(nstmax,ntrajmax),k3_im(nstmax,ntrajmax)
      real*8 k4_re(nstmax,ntrajmax),k4_im(nstmax,ntrajmax)
      integer :: ist1,ist2,itrj     !iteration counters


      do ist1=1,nstate
       k1_re(ist1,itrj)=(en_array_int(ist1,itrj)+eshift)*cel_im(ist1,itrj)
       k1_im(ist1,itrj)=-(en_array_int(ist1,itrj)+eshift)*cel_re(ist1,itrj)
       do ist2=1,nstate
         k1_re(ist1,itrj)=k1_re(ist1,itrj)-cel_re(ist2,itrj)*dotproduct(ist1,ist2,itrj)
         k1_im(ist1,itrj)=k1_im(ist1,itrj)-cel_im(ist2,itrj)*dotproduct(ist1,ist2,itrj)
         k1_re(ist1,itrj)=dtp*k1_re(ist1,itrj)
         k1_im(ist1,itrj)=dtp*k1_im(ist1,itrj)
       enddo
      enddo

      do ist1=1,nstate
       k2_re(ist1,itrj)=(en_array_int(ist1,itrj)+eshift)*(cel_im(ist1,itrj)+k1_im(ist1,itrj)/2)
       k2_im(ist1,itrj)=-(en_array_int(ist1,itrj)+eshift)*(cel_re(ist1,itrj)+k1_re(ist1,itrj)/2)
       do ist2=1,nstate
        dotproduct2(ist1,ist2,itrj)=(dotproduct(ist1,ist2,itrj)+dotproduct_newint(ist1,ist2,itrj))/2
         k2_re(ist1,itrj)=k2_re(ist1,itrj)-(cel_re(ist2,itrj)+k1_re(ist2,itrj)/2)*dotproduct2(ist1,ist2,itrj)
         k2_im(ist1,itrj)=k2_im(ist1,itrj)-(cel_im(ist2,itrj)+k1_im(ist2,itrj)/2)*dotproduct2(ist1,ist2,itrj)
         k2_re(ist1,itrj)=dtp*k2_re(ist1,itrj)
         k2_im(ist1,itrj)=dtp*k2_im(ist1,itrj)
       enddo
      enddo

      do ist1=1,nstate
       k3_re(ist1,itrj)=(en_array_int(ist1,itrj)+eshift)*(cel_im(ist1,itrj)+k2_im(ist1,itrj)/2)
       k3_im(ist1,itrj)=-(en_array_int(ist1,itrj)+eshift)*(cel_re(ist1,itrj)+k2_re(ist1,itrj)/2)
       do ist2=1,nstate
         k3_re(ist1,itrj)=k3_re(ist1,itrj)-(cel_re(ist2,itrj)+k2_re(ist2,itrj)/2)*dotproduct2(ist1,ist2,itrj)
         k3_im(ist1,itrj)=k3_im(ist1,itrj)-(cel_im(ist2,itrj)+k2_im(ist2,itrj)/2)*dotproduct2(ist1,ist2,itrj)
         k3_re(ist1,itrj)=dtp*k3_re(ist1,itrj)
         k3_im(ist1,itrj)=dtp*k3_im(ist1,itrj)
       enddo
      enddo

      do ist1=1,nstate
       k4_re(ist1,itrj)=(en_array_int(ist1,itrj)+eshift)*(cel_im(ist1,itrj)+k3_im(ist1,itrj))
       k4_im(ist1,itrj)=-(en_array_int(ist1,itrj)+eshift)*(cel_re(ist1,itrj)+k3_re(ist1,itrj))
       do ist2=1,nstate
         k4_re(ist1,itrj)=k4_re(ist1,itrj)-(cel_re(ist2,itrj)+k3_re(ist2,itrj))*dotproduct_newint(ist1,ist2,itrj)
         k4_im(ist1,itrj)=k4_im(ist1,itrj)-(cel_im(ist2,itrj)+k3_im(ist2,itrj))*dotproduct_newint(ist1,ist2,itrj)
         k4_re(ist1,itrj)=dtp*k4_re(ist1,itrj)
         k4_im(ist1,itrj)=dtp*k4_im(ist1,itrj)
       enddo
      enddo

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1,itrj)/6+k2_re(ist1,itrj)/3+k3_re(ist1,itrj)/3+k4_re(ist1,itrj)/6
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1,itrj)/6+k2_im(ist1,itrj)/3+k3_im(ist1,itrj)/3+k4_im(ist1,itrj)/6
      enddo


      end

      subroutine interpolate(vx,vy,vz,vx_old,vy_old,vz_old,vx_int,vy_int,vz_int, &
                      ancx,ancy,ancz,nacx_old,nacy_old,nacz_old,en_array_int,en_array_old, &
                      dotproduct,fr,frd,itrj)
      use mod_array_size
      use mod_general
      use mod_sh
      implicit none
      real*8 dotproduct(nstmax,nstmax,ntrajmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
      real*8 vx_int(npartmax,nwalkmax),vy_int(npartmax,nwalkmax),vz_int(npartmax,nwalkmax)
      real*8 nacx_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancx(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancy(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancz(npartmax,ntrajmax,nstmax,nstmax)
      real*8 en_array_int(nstmax,ntrajmax)
      real*8 en_array_old(nstmax,ntrajmax)
      integer :: iat,ist1,ist2,itrj     !iteration counters
      real*8 :: fr,frd
      
      do ist1=1,nstate

        en_array_int(ist1,itrj)=en_array(ist1,itrj)*fr+en_array_old(ist1,itrj)*frd
        do ist2=1,nstate
         dotproduct(ist1,ist2,itrj)=0.0
         do iat=1,natom
         ancx(iat,itrj,ist1,ist2)=nacx(iat,itrj,ist1,ist2)*fr+nacx_old(iat,itrj,ist1,ist2)*frd
         ancy(iat,itrj,ist1,ist2)=nacy(iat,itrj,ist1,ist2)*fr+nacy_old(iat,itrj,ist1,ist2)*frd
         ancz(iat,itrj,ist1,ist2)=nacz(iat,itrj,ist1,ist2)*fr+nacz_old(iat,itrj,ist1,ist2)*frd
         vx_int(iat,itrj)=vx(iat,itrj)*fr+vx_old(iat,itrj)*frd
         vy_int(iat,itrj)=vy(iat,itrj)*fr+vy_old(iat,itrj)*frd
         vz_int(iat,itrj)=vz(iat,itrj)*fr+vz_old(iat,itrj)*frd
         dotproduct(ist1,ist2,itrj)=dotproduct(ist1,ist2,itrj)+vx_int(iat,itrj)*ancx(iat,itrj,ist1,ist2)
         dotproduct(ist1,ist2,itrj)=dotproduct(ist1,ist2,itrj)+vy_int(iat,itrj)*ancy(iat,itrj,ist1,ist2)
         dotproduct(ist1,ist2,itrj)=dotproduct(ist1,ist2,itrj)+vz_int(iat,itrj)*ancz(iat,itrj,ist1,ist2)
         enddo
        enddo
       enddo

      end

      subroutine interpolate_dot(dotproduct_int,fr,frd,itrj)
      use mod_array_size
      use mod_sh,ONLY:nstate,dotproduct_old,dotproduct_new
      implicit none
      real*8 dotproduct_int(nstmax,nstmax,ntrajmax)
      integer :: ist1,ist2,itrj     !iteration counters
      real*8 :: fr,frd

      do ist1=1,nstate
       do ist2=1,nstate
        dotproduct_int(ist1,ist2,itrj)=dotproduct_new(ist1,ist2,itrj)*fr + &
                dotproduct_old(ist1,ist2,itrj)*frd
       enddo
      enddo

      end subroutine interpolate_dot

      subroutine interpolate2(vx_old,vy_old,vz_old,vx_new,vy_new,vz_new, &
                      nacx_old,nacy_old,nacz_old,nacx_new,nacy_new,nacz_new,en_array_old,en_array_new,en_array_int, &
                      dotproduct,fr,frd,itrj)
      use mod_array_size
      use mod_general,only:natom
      use mod_sh, only:nstate
      implicit none
      real*8 dotproduct(nstmax,nstmax)
      real*8 vx_new(npartmax,nwalkmax),vy_new(npartmax,nwalkmax),vz_new(npartmax,nwalkmax)
      real*8 vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
      real*8 vx_int(npartmax,nwalkmax),vy_int(npartmax,nwalkmax),vz_int(npartmax,nwalkmax)
      real*8 nacx_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacx_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_new(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancx(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancy(npartmax,ntrajmax,nstmax,nstmax)
      real*8 ancz(npartmax,ntrajmax,nstmax,nstmax)
      real*8 en_array_int(nstmax)
      real*8 en_array_old(nstmax)
      real*8 en_array_new(nstmax)
      integer :: iat,ist1,ist2,itrj     !iteration counters
      real*8 :: fr,frd
      
      do ist1=1,nstate

        en_array_int(ist1)=en_array_new(ist1)*fr+en_array_old(ist1)*frd
        do ist2=1,nstate
         dotproduct(ist1,ist2)=0.0
         do iat=1,natom
         ancx(iat,itrj,ist1,ist2)=nacx_new(iat,itrj,ist1,ist2)*fr+nacx_old(iat,itrj,ist1,ist2)*frd
         ancy(iat,itrj,ist1,ist2)=nacy_new(iat,itrj,ist1,ist2)*fr+nacy_old(iat,itrj,ist1,ist2)*frd
         ancz(iat,itrj,ist1,ist2)=nacz_new(iat,itrj,ist1,ist2)*fr+nacz_old(iat,itrj,ist1,ist2)*frd
         vx_int(iat,itrj)=vx_new(iat,itrj)*fr+vx_old(iat,itrj)*frd
         vy_int(iat,itrj)=vy_new(iat,itrj)*fr+vy_old(iat,itrj)*frd
         vz_int(iat,itrj)=vz_new(iat,itrj)*fr+vz_old(iat,itrj)*frd
         dotproduct(ist1,ist2)=dotproduct(ist1,ist2)+vx_int(iat,itrj)*ancx(iat,itrj,ist1,ist2)
         dotproduct(ist1,ist2)=dotproduct(ist1,ist2)+vy_int(iat,itrj)*ancy(iat,itrj,ist1,ist2)
         dotproduct(ist1,ist2)=dotproduct(ist1,ist2)+vz_int(iat,itrj)*ancz(iat,itrj,ist1,ist2)
         enddo
        enddo
       enddo

      end
