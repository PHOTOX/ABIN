! -----------------------------------------------------------------
!  ABIN: Program for Born-Oppenheimer MD with potential calculated 
!  on-the-fly by an external procedure placed in ./DYN (and lot more)
!------------------------------------------------------------------  
!  Copyright (C) 2014  D.Hollas,M.Oncak, P.Slavicek
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program in the file LICENSE. If not, see <http://www.gnu.org/licenses/>.
program abin_dyn
   use mod_array_size
   use mod_general
   use mod_system, ONLY:nshake
   use mod_sh
   use mod_fftw3
   implicit none
   INTERFACE
   subroutine init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt)
   use mod_array_size
   real*8,intent(out) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   real*8,intent(out) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(out) :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
   real*8,intent(out) :: dt
   end subroutine init
   subroutine analysis(x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
   use mod_array_size
   implicit none
   real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(in) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(in) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   real*8,intent(in) :: amt(npartmax,nwalkmax)
   real*8,intent(in) :: eclas,equant
   real*8 :: dt  
   end subroutine analysis
   subroutine QtoX(x,y,z,transx,transy,transz)
   use mod_array_size
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine QtoX
   subroutine XtoQ(x,y,z,transx,transy,transz)
   use mod_array_size
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine XtoQ
   END INTERFACE
   real*8  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8  :: amt(npartmax,nwalkmax),amg(npartmax,nwalkmax)
   real*8  :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8  :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
   real*8  :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   real*8  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8  :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   real*8  :: transfxc(npartmax,nwalkmax),transfyc(npartmax,nwalkmax),transfzc(npartmax,nwalkmax)
   real*8  :: transxv(npartmax,nwalkmax),transyv(npartmax,nwalkmax),transzv(npartmax,nwalkmax)
   real*8  :: nacx_old(npartmax,ntrajmax,nstmax,nstmax)
   real*8  :: nacy_old(npartmax,ntrajmax,nstmax,nstmax)
   real*8  :: nacz_old(npartmax,ntrajmax,nstmax,nstmax)
   real*8  :: vx_old(npartmax,nwalkmax),vy_old(npartmax,nwalkmax),vz_old(npartmax,nwalkmax)
   real*8  :: en_array_old(nstmax,ntrajmax)
   real*8  :: dt,eclas,equant
   integer :: itrj,iost
   integer :: iat,iw
   integer,dimension(8) :: values2,values1
   LOGICAL :: file_exists
   character(len=20) :: chit
   character(len=40) :: chrestart
!$ integer :: nthreads,omp_get_max_threads

     call PrintLogo(values1)
     call system('rm -f engrad.dat.* nacm.dat hessian.dat.* geom.dat.*')


!-   INPUT AND INITIALIZATION SECTION      
     call init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt) 
     if(irest.eq.1)then
        write (chit,*)it-1
        chrestart='cp restart.xyz restart.xyz.'//adjustl(chit)
        write(*,*)'Making backup of the current restart file.'
        write(*,*)chrestart
        call system(chrestart)  
     end if

!-------SH initialization -- 
     if(ipimd.eq.2)then
      call sh_init(x,y,z,nacx_old,nacy_old,nacz_old,vx_old,vy_old,vz_old,en_array_old,dt)
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
         x=transx
         y=transy
         z=transz
         call XtoQ(vx,vy,vz,transxv,transyv,transzv)
         vx=transxv
         vy=transyv
         vz=transzv
      endif
!------NORMAL MODE TRANSFORMATION-------------
      if(istage.eq.2)then
if(idebug.eq.1)then
        write(*,*)'Positions before transform'
!        call printf(vx,vy,vz)
        call printf(x,y,z)
endif
       call XtoU(x,y,z,transx,transy,transz)
       x=transx
       y=transy
       z=transz
       call XtoU(vx,vy,vz,transxv,transyv,transzv)
       vx=transxv
       vy=transyv
       vz=transzv
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
!-----End of transformations

!-----Note that amt equals am if staging is off
      px=amt*vx   
      py=amt*vy   
      pz=amt*vz  


      if (ipimd.eq.3)then

       call minimize(x,y,z,fxc,fyc,fzc,eclas)
       
      else


      write(*,*)'#      Step     Time [fs]'
!---------------- PROPAGATION-----------------------------------

!----getting initial forces and energies
   call force_clas(fxc,fyc,fzc,x,y,z,eclas)
   if (ipimd.eq.1) call force_quantum(fxq,fyq,fzq,x,y,z,amg,equant)
!----setting initial values for surface hoping
   if(ipimd.eq.2)then
      itrj=1  ! WARNING: nasty hack
      if(inac.eq.0)then
         iost=readnacm(itrj)
         if(iost.ne.0.and.nac_accu1.gt.nac_accu2)then
!-------------if NACME NOT COMPUTED: TRY TO DECREASE ACCURACY--------------
            call calcnacm(itrj)
            iost=readnacm(itrj)
         endif
         if(iost.ne.0)then
            write(*,*)'Some NACMEs not read. Exiting...'
            stop 1
         endif
      endif
      call set_tocalc()
      call move_vars(en_array_old,nacx_old,nacy_old,nacz_old,vx,vy,vz,vx_old,vy_old,vz_old,itrj)
   endif

!---------LOOP OVER TIME STEPS
!-----it variable is set to 1 or read from restart.xyz in subroutine init
      do it=(it),nstep

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
       if (md.eq.2) call verletstep(x,y,z,px,py,pz,amt,dt,eclas,fxc,fyc,fzc)
      else
       call respashake(x,y,z,px,py,pz,amt,amg,dt,equant,eclas,fxc,fyc,fzc,fxq,fyq,fzq)
      endif

      !TODO: call Update_vel()
!      vx=px/amt
!      vy=py/amt
!      vz=pz/amt
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
      !TODO: px=amt*vx
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

      call finish(values1,values2)

end 


subroutine PrintLogo(values1)
include 'date.inc'
integer,dimension(8),intent(out) :: values1
call date_and_time(VALUES=values1)

print '(a)','                      _ _ _       _       _         _ '
print '(a)','         /\          |      \    | |     | |\      | |'
print '(a)','        /  \         |   _   \   | |     | | \     | |'
print '(a)','       / /\ \        |  |_|  |   | |     | |\ \    | |'
print '(a)','      / /  \ \       |       /   | |     | | \ \   | |'
print '(a)','     / /    \ \      |------/    | |     | |  \ \  | |'
print '(a)','    / /------\ \     |------\    | |     | |   \ \ | |'
print '(a)','   / /--------\ \    |   _   \   | |     | |    \ \| |'
print '(a)','  / /          \ \   |  |_|  |   | |     | |     \ | |'
print '(a)',' / /            \ \  |       /   | |     | |      \| |'
print '(a)','/_/              \_\ |______/    |_|     |_|       |_|'
print '(a)',' '
print '(a)','     version 1.0'
print '(a)','D. Hollas, O.Svoboda, M. Oncak, P. Slavicek       2014'
print '(a)',' '

print *,'Compiled at  ',date
print *,commit
print '(a)',' '

write(*,*)'Job started at:'
write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")values1(5),':', &
        values1(6),':',values1(7),'  ',values1(3),'.',values1(2),'.',&
        values1(1)

end subroutine PrintLogo
