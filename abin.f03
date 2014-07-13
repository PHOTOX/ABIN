! -----------------------------------------------------------------
!  ABIN: Program for Born-Oppenheimer MD with potential calculated 
!  on-the-fly by an external procedure placed in ./DYN (and lot more)
!------------------------------------------------------------------  
!  Copyright (C) 2014                     D.Hollas,M.Oncak, P.Slavicek
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
   use mod_const, only: DP, AUtoFS
   use mod_arrays
   use mod_general
   use mod_sh, only: surfacehop, ntraj, sh_init, get_nacm, move_vars
   use mod_interfaces
   use mod_kinetic, ONLY: temperature
   use mod_utils, only: abinerror, printf
   use mod_shake, only: nshake
   use mod_transform
   use mod_mdstep
   use mod_minimize, only: minimize
   use mod_analysis, only: analysis, restout
   use mod_forces,   only: force_clas, force_quantum
   implicit none
   real(DP)    :: dt,eclas,equant
   integer     :: iat, iw, itrj
   LOGICAL     :: file_exists
   character(len=20) :: chit
   character(len=40) :: chrestart
   integer,dimension(8) :: values2,values1
!$ integer :: nthreads,omp_get_max_threads

     call PrintLogo(values1)
     call system('rm -f ERROR engrad.dat.* nacm.dat hessian.dat.* geom.dat.*')


!-   INPUT AND INITIALIZATION SECTION      
     call init(dt) 
     if(irest.eq.1)then
        write (chit,*)it
        chrestart='cp restart.xyz restart.xyz.'//adjustl(chit)
        write(*,*)'Making backup of the current restart file.'
        write(*,*)chrestart
        call system(chrestart)  
     end if

!-------SH initialization -- 
     if(ipimd.eq.2)then
      call sh_init(x,y,z,vx_old,vy_old,vz_old,dt)
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
      do itrj=1, ntraj
         if (it.eq.0) call get_nacm(itrj)
         call move_vars(vx,vy,vz,vx_old,vy_old,vz_old,itrj)
      end do
   end if

!---------LOOP OVER TIME STEPS
!-----it variable is set to 0 or read from restart.xyz in subroutine init
      it=it+1
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

      if(idebug.eq.2)then
         write(*,*)'positions',it
         call printf(x,y,z)
         write(*,*)'momenta',it
         call printf(px,py,pz)
         write(*,*)'forces',it
         call printf(fxc,fyc,fzc)
      end if


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

      call surfacehop(x,y,z,vx,vy,vz,vx_old,vy_old,vz_old,dt)
      px=amt*vx
      py=amt*vy
      pz=amt*vz

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
        write(*,'(I20,F15.2)')it,it*dt*AUtoFS
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
use iso_fortran_env
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
print '(a)',' D. Hollas, O.Svoboda, M. Oncak, P. Slavicek       2014'
print '(a)',' '

print *,'Compiled at  ',date
print *,commit
!$ print *,'Compiled with parallel OpenMP support for PIMD.'
print *,' '
!COMMENT THIS IF YOU DON'T HAVE FORTRAN 2008
print *, 'This file was compiled by ', &
             compiler_version(), ' using the options: '
print *,     compiler_options()

write(*,*)'Job started at:'
write(*,"(I2,A1,I2.2,A1,I2.2,A2,I2,A1,I2,A1,I4)")values1(5),':', &
        values1(6),':',values1(7),'  ',values1(3),'.',values1(2),'.',&
        values1(1)

end subroutine PrintLogo

subroutine finish(values1,values2)
   use mod_arrays, only: deallocate_arrays
   use mod_general
   use mod_nhc
   use mod_estimators, only: h
   use mod_harmon, only: hess
   use mod_fftw3
   implicit none
   integer,dimension(8),intent(in)  :: values1
   integer,dimension(8),intent(out) :: values2
   real(DP) :: TIME
!   integer :: iter=-3

   call deallocate_arrays( )

   close(1)
   close(2)
   if(ipimd.eq.1.or.icv.eq.1)then
    close(7)
   endif
   if(ipimd.eq.2)then
    close(3)
    close(4)
    close(5)
   endif
   if(nwritev.gt.0) close(13)
!--------------CLEANING-------------------------
   if (ihess.eq.1) deallocate ( hess )
   if (ihess.eq.1.and.pot.eq.'nab') deallocate ( h )
   if (istage.eq.2) call fftw_end()

   if(inose.eq.1)then
    deallocate( w )
    deallocate( Qm )
    deallocate( ms )
    if (imasst.eq.1)then
      deallocate( pnhx )
      deallocate( pnhy )
      deallocate( pnhz )
      deallocate( xi_x )
      deallocate( xi_y )
      deallocate( xi_z )
     else
      deallocate( pnhx )
      deallocate( xi_x )
    endif
   endif
!TODO dealokovat pole v NABU ...tj zavolet mme rutinu s iter=-3 nebo tak neco
!   if(pot.eq.'nab') call mme(NULL,NULL,iter)

   write(*,*)''
   write(*,*)' Job finished!'
   write(*,*)''

   write(*,*)''

!---------TIMING-------------------------------
   call cpu_time(TIME)
   write(*,*)' Total cpu time [s] (does not include ab initio calculations)'
   write(*,*)TIME
   write(*,*)' Total cpu time [hours] (does not include ab initio calculations)'
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
end subroutine finish


