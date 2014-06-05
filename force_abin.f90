      subroutine force_abin(x,y,z,fx,fy,fz,eclas)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: names
      use mod_estimators, ONLY: hess
      use mod_sh
      use mod_qmmm, ONLY:natqm
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8  :: eclas,temp1
      integer :: iat,iw,iat1,iat2,itest,nthreads=1 !,ithread
      character(len=100) :: chsystem
      character(len=20) :: chgeom,chforce,chhess,fgeom
      logical :: file_exists
      integer :: ist1,ist2,iost
!$    integer :: omp_get_max_threads,OMP_get_thread_num

!$    nthreads=omp_get_max_threads()

      eclas=0.0d0
!      ithread=1

!$OMP PARALLEL  ! REDUCTION(+:eclas) alternativa k atomic
!$OMP DO PRIVATE(temp1,chsystem,chgeom,chforce,chhess)
      do iw=1,nwalk

!!$   ithread=OMP_get_thread_num()
!$    if(nthreads.eq.1)then
       write(chgeom,'(A8)')'geom.dat'
       write(chforce,'(A10)')'engrad.dat'
       write(chhess,'(A11)')'hessian.dat'
!$    else
!$     write(chgeom,'(A,I3.3)')'geom.dat.',iw
!$     write(chforce,'(A,I3.3)')'engrad.dat.',iw
!$     write(chhess,'(A,I3.3)')'hessian.dat.',iw
!$    endif

!----WRITING GEOMETRY IN ANGSTROMS
     fgeom='(A2,3E25.17E2)'
     open(unit=20+iw,file=chgeom)
      do iat=1,natqm
       write(20+iw,fgeom)names(iat),x(iat,iw)/ang,y(iat,iw)/ang,z(iat,iw)/ang
      enddo
     close(unit=20+iw)

!---- SH     
      if(ipimd.eq.2)then
       open(unit=20+iw+2*nwalk,file='state.dat')
       write(20+iw+2*nwalk,'(I1)')istate(iw)
       write(20+iw+2*nwalk,'(I1)')nstate

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



!-----------------------
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
      case ('nab')
              write(chsystem,*)'./NAB/r.nab '
      case ('gamess')
              write(chsystem,*)'./GAMESS/r.gamess '
      case ('qchem')
              write(chsystem,*)'./QCHEM/r.qchem '
      case ('dyn')
              write(chsystem,*)'./DYN/r.dyn '
      case DEFAULT
              write(*,*)'Error in input parameter "pot"'
              write(*,*)'Exiting, NOW!'
              stop 1
     end SELECT

!TODO: always pass time step for the bash script
!   write(chsystem,*)chsystem,it

!TODO: for OpenMPI, add the index of the walker
!do not use different script as implied here
   if(nthreads.gt.1)then
    write(chsystem,'(A30,A1,I3.3)')chsystem,'.',iw
   endif

   if(ipimd.eq.2)then
    write(chsystem,'(A40,A12)')chsystem,' < state.dat'
   endif
     

   call system(chsystem)

   if(nthreads.eq.1)then

!----make sure that the file exist     
     itest=0
     INQUIRE(FILE=chforce, EXIST=file_exists)
     do while(.not.file_exists.and.itest.lt.10)
     write(*,*)'WARNING:File ',chforce,' does not exist. Waiting..'
!     call sleep(1)
     call system('sync')    !mel by zajistit flush diskoveho bufferu
     INQUIRE(FILE=chforce, EXIST=file_exists)
     itest=itest+1
     end do
     
     open(unit=20+iw,file=chforce,status='old',ACTION='READ')

!-------READING ENERGY from engrad.dat
     if(ipimd.ne.2)then
      read(20+iw,*,IOSTAT=iost)temp1
      if(iost.ne.0)then
              write(*,*)'Fatal problem with reading energy from engrad.dat'
              write(*,*)'This usually means, that the ab initio program failed. See the appropriate output files.'
              write(*,*)'Exiting...'
              stop 1
      endif
!$OMP ATOMIC
      eclas=eclas+temp1
     else
! SH             
      do ist1=1,nstate
       read(20+iw,*)en_array(ist1,iw)
      enddo
      eclas=en_array(istate(iw),iw)
     endif


!----READING energy gradients from engrad.dat
     do iat=1,natqm
      read(20+iw,*,IOSTAT=iost)fx(iat,iw),fy(iat,iw),fz(iat,iw)
      if(iost.ne.0)then
              write(*,*)'Fatal problem with reading gradients from engrad.dat'
              write(*,*)'This usually means, that the ab initio program failed. See the appropriate output files.'
              write(*,*)'Exiting...'
              stop 1
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
!      write(*,*)hess(iat1,iat2,iw),hess(iat1+1,iat2,iw),hess(iat1+2,iat2,iw)
       enddo
      enddo
     endif


     close(unit=20+iw,status='delete')
     if(ihess.eq.1)then
      close(unit=20+iw+nwalk,status='delete')
     endif
!   nthreads endif     
    endif     

    enddo
!$OMP END DO 

!$     if(nthreads.gt.1)then
!$OMP DO PRIVATE(temp1,chforce,chhess)
!$     do iw=1,nwalk
!$     write(chforce,'(A,I3.3)')'engrad.dat.',iw
!$     write(chhess,'(A,I3.3)')'hessian.dat.',iw

!----make sure that the file exist     
!$     INQUIRE(FILE=chforce, EXIST=file_exists)
!$     do while(.not.file_exists.and.itest.lt.10)
!$      write(*,*)'WARNING:File ',chforce,' does not exist. Waiting..'
!$      call system('sync')    !mel by zajistit flush diskoveho bufferu
!$      INQUIRE(FILE=chforce, EXIST=file_exists)
!$      itest=itest+1
!$     end do
!$     
!$     open(unit=20+iw,file=chforce,status='old',ACTION='READ')
!-------READING ENERGY from engrad.dat
!$     read(20+iw,*)temp1
!$OMP ATOMIC
!$     eclas=eclas+temp1
!$
!----READING GRADIENTS from engrad.dat
!$     do iat=1,natqm
!$      read(20+iw,*)fx(iat,iw),fy(iat,iw),fz(iat,iw)
!  In most cases, we actually produce energy gradients         
!$       fx(iat,iw)=-fx(iat,iw)
!$       fy(iat,iw)=-fy(iat,iw)
!$       fz(iat,iw)=-fz(iat,iw)
!$     enddo
!$
!----READING of HESSIAN     
!$     if (ihess.eq.1)then
!$      INQUIRE(FILE=chhess, EXIST=file_exists)
!$      do while(.not.file_exists.and.itest.lt.10)
!$       write(*,*)'WARNING:File ',chhess,' does not exist. Waiting..'
!$       call system('sync')    !mel by zajistit flush diskoveho bufferu
!$       INQUIRE(FILE=chhess, EXIST=file_exists)
!$       itest=itest+1
!$      end do
!$
!$      open(unit=20+iw+nwalk,file=chhess,status='old',ACTION='READ')
!$
!$       do iat2=1,natqm*3
!$        do iat1=1,natqm*3,3
!$        read(20+iw+nwalk,*)hess(iat1,iat2,iw),hess(iat1+1,iat2,iw),hess(iat1+2,iat2,iw)
!$        hess(iat1,iat2,iw)=hess(iat1,iat2,iw)/nwalk
!$        hess(iat1+1,iat2,iw)=hess(iat1+1,iat2,iw)/nwalk
!$        hess(iat1+2,iat2,iw)=hess(iat1+2,iat2,iw)/nwalk
!$!        write(*,*)hess(iat1,iat2,iw),hess(iat1+1,iat2,iw),hess(iat1+2,iat2,iw)
!$        enddo
!$       enddo
!$      endif
!$
!$     close(unit=20+iw,status='delete')
!$     if(ihess.eq.1)then
!$      close(unit=20+iw+nwalk,status='delete')
!$     endif
!$
!$     enddo
!$OMP END DO
!     nthreads endif
!$     endif
!$OMP END PARALLEL 


         end
                                        

