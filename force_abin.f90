      subroutine force_abin(x,y,z,fx,fy,fz,eclas)
      use mod_const,    only: DP, ANG
      use mod_general
      use mod_system,   only: names
      use mod_harmon,   only: hess
      use mod_sh,       only: nac_accu1, tocalc, en_array, istate, nstate
      use mod_qmmm,     only: natqm
      use mod_utils,    only: abinerror
      implicit none
      real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out)   :: fx(:,:),fy(:,:),fz(:,:)
      real(DP),intent(out)   :: eclas
      real(DP)  :: temp1
      character(len=100) :: chsystem
      character(len=20) :: chgeom,chforce,chhess,fgeom
      logical :: file_exists
      integer :: iat,iw,iat1,iat2,itest !,nthreads=1, ithread
      integer :: ist1,ist2,iost
!!$    integer :: omp_get_max_threads,OMP_get_thread_num

!!$    nthreads=omp_get_max_threads()

     eclas=0.0d0
!     ithread=1
!----Format for geom.dat; needed,so that Molpro can read it
     fgeom='(A2,3E25.17E2)'

!$OMP PARALLEL  ! REDUCTION(+:eclas) alternativa k atomic
!$OMP DO PRIVATE(temp1,chsystem,chgeom,chforce,chhess,itest,file_exists,iost)
     do iw=1,nwalk

!!$   ithread=OMP_get_thread_num()
     write(chgeom,'(A,I3.3)')'geom.dat.',iw
     write(chforce,'(A,I3.3)')'engrad.dat.',iw
     write(chhess,'(A,I3.3)')'hessian.dat.',iw

!----WRITING GEOMETRY IN ANGSTROMS
     open(unit=20+iw,file=chgeom)
      do iat=1,natqm
       write(20+iw,fgeom)names(iat),x(iat,iw)/ang,y(iat,iw)/ang,z(iat,iw)/ang
      enddo
     close(unit=20+iw)

!---- SH     
      if(ipimd.eq.2)then
       open(unit=20+iw+2*nwalk,file='state.dat')
       write(20+iw+2*nwalk,'(I2)')istate(iw)
       write(20+iw+2*nwalk,'(I2)')nstate

!horni troj. matice bez diagonaly
! tocalc(,)=1 -> pocitame couplingy
! tocalc(,)=0 -> nepocitame couplingy
       do ist1=1,nstate-1
        do ist2=ist1+1,nstate
         write(20+iw+2*nwalk,'(I1,A1)',advance='no')tocalc(ist1,ist2),' ' 
        enddo
       enddo
       close(20+iw+2*nwalk) 
       endif



!--- HERE we decide which program we use to obtain gradients and energies
     select case (pot)
      case ('g09')
              write(chsystem,*)'./G09/r.g09 '
      case ('tera')
              write(chsystem,*)'./TERA/r.tera '
      case ('orca')
              write(chsystem,*)'./ORCA/r.orca '
      case ('molpro')
              write(chsystem,*)'./MOLPRO/r.molpro '
      case ('turbo')
              write(chsystem,*)'./TURBO/r.turbo '
      case ('gamess')
              write(chsystem,*)'./GAMESS/r.gamess '
      case ('qchem')
              write(chsystem,*)'./QCHEM/r.qchem '
      case ('dyn')
              write(chsystem,*)'./DYN/r.dyn '
      case DEFAULT
              write(*,*)'Error in input parameter "pot"'
              call abinerror('force_abin')
     end SELECT

!-Passing arguments to bash script
!-First argument is time step
!-Second argument is the bead index, used for parallel calculations
   write(chsystem,'(A20,I13,I4.3)')chsystem,it,iw

  !for SH, pass the 4th parameter:precision of NACME as 10^(-nac_accu1)
   if(ipimd.eq.2)then
    write(chsystem,'(A40,I3,A12)')chsystem,nac_accu1,' < state.dat'
   endif
     
   call system(chsystem)

!----make sure that the file exist and flush the disc buffer     
     itest=0
     INQUIRE(FILE=chforce, EXIST=file_exists)
     do while(.not.file_exists.and.itest.lt.10)
     write(*,*)'WARNING:File ',chforce,' does not exist. Waiting..'
     call system('sync')    !mel by zajistit flush diskoveho bufferu
     INQUIRE(FILE=chforce, EXIST=file_exists)
     itest=itest+1
     end do
     
     open(unit=20+iw,file=chforce,status='old',ACTION='READ')

!-------READING ENERGY from engrad.dat
      read(20+iw,*,IOSTAT=iost)temp1
      if(iost.ne.0)then
              write(*,*)'Fatal problem with reading energy from file', chforce
              write(*,*)'This usually means, that the ab initio program failed.'
              write(*,*)'See the appropriate output files.'
              call abinerror('force_abin')
      endif
!$OMP ATOMIC
      eclas=eclas+temp1
! SH             
     if(ipimd.eq.2)then
      en_array(1,iw)=temp1
      do ist1=2,nstate
       read(20+iw,*)en_array(ist1,iw)
      enddo
      eclas=en_array(istate(iw),iw)
     endif


!----READING energy gradients from engrad.dat
     do iat=1,natqm
      read(20+iw,*,IOSTAT=iost)fx(iat,iw),fy(iat,iw),fz(iat,iw)
      if(iost.ne.0)then
              write(*,*)'Fatal problem with reading energy from engrad.dat'
              write(*,*)'This usually means, that the ab initio program failed.'
              write(*,*)'See the appropriate output files.'
              call abinerror('force_abin')
      endif
!---Conversion to forces        
       fx(iat,iw)=-fx(iat,iw)
       fy(iat,iw)=-fy(iat,iw)
       fz(iat,iw)=-fz(iat,iw)
     enddo

!----READING of HESSIAN     
     if (ihess.eq.1)then
      INQUIRE(FILE=chhess, EXIST=file_exists)
      do while(.not.file_exists.and.itest.lt.10)
       write(*,*)'WARNING:File ',chhess,' does not exist. Waiting..'
       call system('sync')    !mel by zajistit flush diskoveho bufferu
       INQUIRE(FILE=chhess, EXIST=file_exists)
       itest=itest+1
      end do

      open(unit=20+iw+nwalk,file=chhess,status='old',ACTION='READ')

      do iat2=1,natqm*3
       do iat1=1,natqm*3,3
       read(20+iw+nwalk,*)hess(iat1,iat2,iw),hess(iat1+1,iat2,iw),hess(iat1+2,iat2,iw)
       hess(iat1,iat2,iw)=hess(iat1,iat2,iw)/nwalk
       hess(iat1+1,iat2,iw)=hess(iat1+1,iat2,iw)/nwalk
       hess(iat1+2,iat2,iw)=hess(iat1+2,iat2,iw)/nwalk
       enddo
      enddo
     endif


     close(unit=20+iw,status='delete')
     if(ihess.eq.1)then
      close(unit=20+iw+nwalk,status='delete')
     endif

    end do
!$OMP END DO 
!$OMP END PARALLEL  

end
                                        

