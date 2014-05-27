      program abin_dyn
! ---------------------------------------------------------
!  Program for Born-Oppenheimer MD with potential calculated 
!  on-the-fly by an external procedure placed in ./DYN
!-----------------------------------------------------------  
      use mod_array_size
      use mod_general
      use mod_system, ONLY:nshake
      use mod_nhc, ONLY:imasst
      use mod_sh
      use mod_estimators,ONLY:hess
      use mod_fftw3
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 amt(npartmax,nwalkmax),amg(npartmax,nwalkmax)
      real*8 fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
      real*8 fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
      real*8 transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
      real*8 transfxc(npartmax,nwalkmax),transfyc(npartmax,nwalkmax),transfzc(npartmax,nwalkmax)
      real*8 transxv(npartmax,nwalkmax),transyv(npartmax,nwalkmax),transzv(npartmax,nwalkmax)
      real*8 nacx_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacy_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 nacz_old(npartmax,ntrajmax,nstmax,nstmax)
      real*8 vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
      real*8 en_array_old(nstmax,ntrajmax)
      real*8 dt,TIME,eclas,equant
      integer :: itrj,ist1,ist2,iost
      integer :: iat,iw
      integer,dimension(8) :: values2,values1
      LOGICAL :: file_exists
!$    integer :: nthreads,omp_get_max_threads


      call date_and_time(VALUES=values1)
      call system('rm -f forces.dat nacm.dat hessian.dat geom.dat')


!!   INPUT AND INITIALIZATION SECTION      
     call init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt) 
     if(irest.eq.1)  call system('cp restart.xyz restart0.xyz')  !zaloha puvodniho restartu

!-------SH initialization -- 
     if(ipimd.eq.2)then
      call sh_init(nacx_old,nacy_old,nacz_old,vx_old,vy_old,vz_old,en_array_old,dt)
     endif

      write(*,*)''
!$    nthreads=omp_get_max_threads()
!$    write(*,*)'Number of threads used = ',nthreads
!$    write(*,*)''

! Stage transformation
! Masses, velocities and positions are transformed here into a new set of u variables
! See Tuckermann's article in "Quantum Simulations of Complex Many Body Systems'. 
      if(istage.eq.1)then
      call XtoQ(x,y,z,transx,transy,transz)
      do iat=1,natom
       do iw=2,nwalk ! x(iat,1)=transx(iat,1)
       x(iat,iw)=transx(iat,iw)
       y(iat,iw)=transy(iat,iw)
       z(iat,iw)=transz(iat,iw)
       enddo
      enddo
      call XtoQ(vx,vy,vz,transxv,transyv,transzv)
      do iat=1,natom
       do iw=2,nwalk
       vx(iat,iw)=transxv(iat,iw)
       vy(iat,iw)=transyv(iat,iw)
       vz(iat,iw)=transzv(iat,iw)
       enddo
      enddo

      endif
!------NORMAL MODE TRANSFORMATION-------------
      if(istage.eq.2)then
if(idebug.eq.1)then
        write(*,*)'Positions befir transform'
!        call printf(vx,vy,vz)
        call printf(x,y,z)
endif
       call XtoU(x,y,z,x,y,z)
       call XtoU(vx,vy,vz,vx,vy,vz)
if(idebug.eq.1)then
        write(*,*)'Positions after transform'
        call printf(x,y,z)
       call Utox(x,y,z,transx,transy,transz)
        write(*,*)'Positions after back transform'
        call printf(transx,transy,transz)
        stop 1
endif
      endif

      call init_mass(amg,amt)
! End of transformations

! Note that amt equals am if staging is off
      do iw=1,nwalk
       do iat=1,natom
         px(iat,iw)=amt(iat,iw)*vx(iat,iw)   
         py(iat,iw)=amt(iat,iw)*vy(iat,iw)   
         pz(iat,iw)=amt(iat,iw)*vz(iat,iw)  
       enddo
      enddo

      if (ipimd.eq.3)then

       call minimize(x,y,z,fxc,fyc,fzc,eclas)
       
      else


      write(*,*)'#      Step     Time [fs]'
      
!---------------- PROPAGATION-----------------------------------

!----getting initial forces and energies
      call force_clas(fxc,fyc,fzc,x,y,z,eclas)
!----setting initial values for surface hoping
      if(ipimd.eq.2)then
       itrj=1  ! WARNING: nasty hack
       eshift=-en_array(1,itrj)
       if(inac.eq.0)then
        iost=readnacm(itrj)
        if(iost.ne.0)then
!-------------if NACM NOT COMPUTED: TRY TO DECREASE ACCURACY--------------
        write(*,*)'Some NACM not computed.Trying with decreased accuracy.'
        write(*,*)'Calling script r.molpro.nacm'
        call calcnacm(itrj)
        iost=readnacm(itrj)
        endif
        if(iost.ne.0)then
         write(*,*)'Some nac not read. Exiting...'
         stop 1
        endif
       endif
       call set_tocalc()
       call move_vars(en_array_old,nacx_old,nacy_old,nacz_old,vx,vy,vz,vx_old,vy_old,vz_old,itrj)
      endif

!---------LOOP OVER TIME STEPS
!-----it variable is set to 1 or read from restart.xyz in subroutine init
      do it=it,nstep

      INQUIRE(FILE="EXIT", EXIST=file_exists)
      if(file_exists)then
       write(*,*)'Found file EXIT. Exiting...'
       call system('rm EXIT')
       if (istage.gt.0)then

        if(istage.eq.1)then     
         call QtoX(vx,vy,vz,transxv,transyv,transzv)
         call QtoX(x,y,z,transx,transy,transz)
        endif
        if(istage.eq.2)then
         call UtoX(x,y,z,transx,transy,transz)
         call UtoX(vx,vy,vz,transxv,transyv,transzv)
        endif

        call restout(transx,transy,transz,transxv,transyv,transzv,it-1)
       else
        call restout(x,y,z,vx,vy,vz,it-1)
       endif
       exit                                          !break from time loop
      endif

!-----CALL RESPA or VELOCITY VERLET--------------
      if(nshake.eq.0)then
       if (md.eq.1) call respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas,fxc,fyc,fzc,fxq,fyq,fzq)
       if (md.eq.2) call verletstep(x,y,z,px,py,pz,amt,amg,dt,eclas,fxc,fyc,fzc)
      else
       call respashake(x,y,z,px,py,pz,amt,amg,dt,equant,eclas,fxc,fyc,fzc,fxq,fyq,fzq)
      endif

      do iw=1,nwalk
       do iat=1,natom
        vx(iat,iw)=px(iat,iw)/amt(iat,iw)
        vy(iat,iw)=py(iat,iw)/amt(iat,iw)
        vz(iat,iw)=pz(iat,iw)/amt(iat,iw)
       enddo
      enddo

!-------SURFACE HOPPING SECTION----------------------------      
      if(ipimd.eq.2)then

      call surfacehop(x,y,z,vx,vy,vz,nacx_old,nacy_old,nacz_old,vx_old,vy_old,vz_old,en_array_old,dt)
      do itrj=1,ntraj
       do iat=1,natom
        px(iat,itrj)=amt(iat,itrj)*vx(iat,itrj)
        py(iat,itrj)=amt(iat,itrj)*vy(iat,itrj)
        pz(iat,itrj)=amt(iat,itrj)*vz(iat,itrj)
       enddo
      enddo

      endif


!--------------------SECTION of trajectory ANALYSIS
! In order to analyze the output, we have to perform the back transformation
! Transformed (cartesian) coordinates are stored in trans matrices.
! DHmod-21.12.2012 enter this section only every ncalc step

      if(modulo(it,ncalc).ne.0) cycle

      if(istage.eq.1)then     
       call QtoX(vx,vy,vz,transxv,transyv,transzv)
       call QtoX(x,y,z,transx,transy,transz)
       call FQtoFX(fxc,fyc,fzc,transfxc,transfyc,transfzc)
      endif

      if(istage.eq.2)then
       call UtoX(x,y,z,transx,transy,transz)
       call UtoX(vx,vy,vz,transxv,transyv,transzv)
       call UtoX(fxc,fyc,fzc,transfxc,transfyc,transfzc)
       if(idebug.eq.1) then
        write(*,*)'Back transformed forces'
        call printf(transfxc,transfyc,transfzc)
       endif
      endif

      call temperature(px,py,pz,amt,dt,eclas)

      if(istage.eq.1.or.istage.eq.2)then
        call analysis (transx,transy,transz,transxv,transyv,transzv,  &
                      transfxc,transfyc,transfzc,amt,eclas,equant,dt)
       else
        call analysis (x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
      endif
      

      if(modulo(it,nwrite).eq.0)then
        write(*,'(I20,F15.2)')it,it*dt*autofs
        call flush(6)
      endif

      !------------------------------------------------------------------------
!   Time step loop      
      enddo 

!DUMP restart file at the end of a run even if the final step is not compatible with nrest
!because ncalc might be >1, we have to perform transformation to get the most
!recent coordinates and velocities
       if (istage.gt.0)then
        if(istage.eq.1)then     
         call QtoX(vx,vy,vz,transxv,transyv,transzv)
         call QtoX(x,y,z,transx,transy,transz)
        endif
        if(istage.eq.2)then
         call UtoX(x,y,z,transx,transy,transz)
         call UtoX(vx,vy,vz,transxv,transyv,transzv)
        endif
        call restout(transx,transy,transz,transxv,transyv,transzv,it-1)
       else
        call restout(x,y,z,vx,vy,vz,it-1)
       endif

!   minimization endif
      endif

      call finish()

! some cleaning if we do parallel execution
!!$   if(nthreads.gt.1.and.pot.ne.'guillot') call system('rm geom.*')
!---------TIMING-------------------------------
      call cpu_time(TIME)
      write(*,*)'Total cpu time [s]'
      write(*,*)TIME
      write(*,*)'Total cpu time [hours]'
      write(*,*)TIME/3600.

      call date_and_time(VALUES=values2)
      write(*,*)'Job started at:'
      write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")values1(5),':', &
        values1(6),':',values1(7),'  ',values1(3),'.',values1(2),'.',&
        values1(1)
      write(*,*)'Job finished at:'
      write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")values2(5),':',&
        values2(6),':',values2(7),'  ',values2(3),'.',values2(2),'.',&
        values2(1)
     
      end 
