! -----------------------------------------------------------------
!  ABIN: Multipurpose ab initio MD program.
!  Potential is calculated on-the-fly by an external program.
!------------------------------------------------------------------  
!  Copyright (C) 2014             D.Hollas, M.Oncak, O.Svoboda and P.Slavicek
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
#ifdef PLUM
   use mod_plumed, only: plumed
#endif
#ifdef MPI
   use mod_remd, only: nswap, remd_swap
   implicit none
   include 'mpif.h'
#else
   implicit none
#endif
   real(DP)    :: dt=20.0d0, eclas, equant
   integer     :: iat, iw, itrj
   LOGICAL     :: file_exists
   character(len=20) :: chit
   character(len=40) :: chrestart
   character(len=200) :: chsystem
   integer,dimension(8) :: values1, values2
   real(DP) :: TIME
   integer  :: ierr
!$ integer  :: nthreads,omp_get_max_threads

!-   INPUT AND INITIALIZATION SECTION      
   call init(dt, values1) 

   ! TODO: move this bit to init
   if(my_rank.eq.0) call system('rm -f ERROR engrad*.dat.* nacm.dat hessian.dat.* geom.dat.*')

   ! TODO: in case of _cp2k_, this needs to be done only by rank0
   if(irest.eq.1.and.(my_rank.eq.0.or.iremd.eq.1))then
      write (chit,*)it
      if(iremd.eq.1)then
         write(chrestart,'(A,I2.2)')'restart.xyz.',my_rank
      else
         chrestart='restart.xyz'
      end if
      chsystem='cp '//trim(chrestart)//'  '//trim(chrestart)//'.'//adjustl(chit)
      write(*,*)'Making backup of the current restart file.'
      write(*,*)chsystem
      call system(chsystem)  
   end if

!-------SH initialization -- 
   if(ipimd.eq.2)then
      call sh_init(x, y, z, vx, vy, vz, dt)
   endif

!$ nthreads=omp_get_max_threads()
   if (my_rank.eq.0)then
!$    write(*,*)'Number of threads used = ',nthreads
      write(*,*)''
   end if

!  Stage transformation
!  Masses, velocities and positions are transformed here into a new set of u variables
!  See Tuckermann's article in "Quantum Simulations of Complex Many Body Systems'. 
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
         call abinerror('Only debug')
      endif
   endif

   call init_mass(amg, amt)
!-----End of transformations

!-----Note that amt equals am if staging is off
   px = amt * vx   
   py = amt * vy   
   pz = amt * vz  


   if (ipimd.eq.3)then

      call minimize(x, y, z, fxc, fyc, fzc, eclas)
       
   else


      if (my_rank.eq.0)then
         write(*,*)
         write(*,*)'#      Step     Time [fs]'
      end if
!---------------- PROPAGATION-----------------------------------

#ifdef MPI
!     Without this Barrier, ranks > 0 do not write geom.dat in force_clas
!     I don't know why the hell not.
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
!----getting initial forces and energies
      call force_clas(fxc, fyc, fzc, x, y, z, eclas)
      if (ipimd.eq.1) call force_quantum(fxq, fyq, fzq, x, y, z, amg, equant)
!----PLUMED section
     
!----setting initial values for surface hoping
      if(ipimd.eq.2)then
         do itrj=1, ntraj
            if (it.eq.0) call get_nacm(itrj)
            call move_vars(vx, vy, vz, vx_old, vy_old, vz_old, itrj)
         end do
      end if

!---------LOOP OVER TIME STEPS
!---- "it" variable is set to 0 or read from restart.xyz in subroutine init
      it=it+1
      do it=(it),nstep

#ifdef MPI
! This is needed, because all ranks need to see EXIT
! Maybe we should get rid of it for performance reasons
! We could still call Abort from rank0,but we would not be sure that we are on
! the same timestep. Does it matter??
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
         INQUIRE(FILE="EXIT", EXIST=file_exists)
         if(file_exists)then
            write(*,*)'Found file EXIT. Exiting...'
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

#ifdef MPI
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
            if (my_rank.eq.0) call system('rm EXIT')

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
         vx = px / amt
         vy = py / amt
         vz = pz / amt

!-------SURFACE HOPPING SECTION----------------------------      
         if(ipimd.eq.2)then
            call surfacehop(x, y, z, vx, vy, vz, vx_old, vy_old, vz_old, dt, eclas)
            px = amt * vx
            py = amt * vy
            pz = amt * vz
         endif

#ifdef MPI
         if (iremd.eq.1.and.modulo(it,nswap).eq.0.and.it.gt.imini)then
            call remd_swap(x, y, z, x, y, z, fxc, fyc, fzc, eclas)
         end if
#endif

!--------------------SECTION of trajectory ANALYSIS
! In order to analyze the output, we have to perform the back transformation
! Transformed (cartesian) coordinates are stored in trans matrices.
! DHmod-21.12.2012 enter this section only every ncalc step

! maybe we should visit this section only when my_rank.eq.0
! this would not be the case for REMD

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
         
      
         if(modulo(it,nwrite).eq.0.and.my_rank.eq.0)then
           write(*,'(I20,F15.2)')it,it*dt*AUtoFS
           call flush(6)
         endif

      !------------------------------------------------------------------------
!   Time step loop      
      enddo
 
!DUMP restart file at the end of a run even if the final step is not compatible with nrest
!Because ncalc might be >1, we have to perform transformation to get the most
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
   call finish(0)

!---------TIMING-------------------------------
   if (my_rank.eq.0)then
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
   end if

end 



subroutine finish(error_code)
   use mod_arrays, only: deallocate_arrays
   use mod_general
   use mod_nhc
   use mod_estimators, only: h
   use mod_harmon, only: hess

#ifdef USEFFTW
   use mod_fftw3,  only: fftw_end
#endif

#ifdef CP2K
   use mod_cp2k,   only: cp2k_finalize
#endif

#ifdef PLUM
   use mod_plumed, only: plumed,plumedfile
#endif

#ifdef MPI
   use mod_terampi, only: finalize_terachem
   implicit none
   include "mpif.h"
   integer :: errmpi
#else
   implicit none
#endif

   real(DP) :: TIME
   integer  :: i, ierr, error_code
   logical  :: lopen
!   integer :: iter=-3

#ifdef MPI
   if (pot.eq.'_tera_') call finalize_terachem()
#endif


   call deallocate_arrays( )

   close(2);close(3)
   do i=7,20
      inquire(unit=i,opened=lopen)
      if (lopen.and.i.ne.5.and.i.ne.6) close(i)
   end do

!--------------CLEANING-------------------------
   if (allocated(hess)) deallocate ( hess )
   if (allocated(h)) deallocate ( h )

#ifdef USEFFTW
   if (istage.eq.2) call fftw_end()
#endif

   if(allocated(w))     deallocate( w )
   if(allocated(Qm))    deallocate( Qm )
   if(allocated(ms))    deallocate( ms )
   if(allocated(pnhx))  deallocate( pnhx )
   if(allocated(pnhy))  deallocate( pnhy )
   if(allocated(pnhz))  deallocate( pnhz )
   if(allocated(xi_x))  deallocate( xi_x )
   if(allocated(xi_y))  deallocate( xi_y )
   if(allocated(xi_z))  deallocate( xi_z )
!TODO dealokovat pole v NABU ...tj zavolet mme rutinu s iter=-3 nebo tak neco
!   if(pot.eq.'nab') call mme(NULL,NULL,iter)

   if (my_rank.eq.0)then
      write(*,*)''
      if (error_code.eq.0)then
         write(*,*)' Job finished!'
      else
         write(*,*)'Error encountered. See file ERROR for more information.'
         write(*,*)''
      end if
   end if

#ifdef MPI
   if (error_code.eq.0.and.pot.ne."_cp2k_")then
      call MPI_FINALIZE ( errmpi )
      if (errmpi.ne.0)then
         write(*,*)'Bad signal from MPI_FINALIZE: ', errmpi
         ! Let's try to continue
      end if
   else if (error_code.gt.0)then
      call MPI_Abort(MPI_COMM_WORLD, ierr)
   end if
#endif
#ifdef CP2K
!  MPI_FINALIZE is called in this routine as well
   if(pot.eq.'_cp2k_') call cp2k_finalize()
#endif

!   PLUMED closing session
#ifdef PLUMED
    if (plumed.eq.1) then
    call plumed_f_gfinalize()
    end if
#endif
end subroutine finish


