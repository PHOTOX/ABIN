!---mod_transform                  P.Slavicek and D.Hollas, 9.2.2012
!---Transformation routines for PIMD and PI+GLE.
!---Various staging transformation subroutines are placed here.
!---NORMAL MODE TRANSFORMATION added during March 2013

module mod_transform
   use mod_const, only: DP
   use mod_general, only:natom, nwalk
   use mod_utils, only: abinerror
   implicit none
   private
   public :: UtoX, XtoU, QtoX, XtoQ, FXtoFQ, FQtoFX
   public :: init_mass
   public :: equant_nm
   contains

!  This routine transforms staging coordinates to cartesian coordinates, which are stored
!  in trans matrices,values in x,y and z matrices are NOT modified!!      
   subroutine QtoX(x,y,z,transx,transy,transz)
   real(DP),intent(inout)  :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out)    :: transx(:,:),transy(:,:),transz(:,:)
   integer :: iat,iw

   do iat=1,natom
      x(iat,nwalk+1) = x(iat,1)
      y(iat,nwalk+1) = y(iat,1)
      z(iat,nwalk+1) = z(iat,1)

      transx(iat,1) = x(iat,1)
      transy(iat,1) = y(iat,1)
      transz(iat,1) = z(iat,1)
      transx(iat,nwalk+1) = x(iat,1)
      transy(iat,nwalk+1) = y(iat,1)
      transz(iat,nwalk+1) = z(iat,1)

      do iw=nwalk,2,-1
         transx(iat,iw) = x(iat,iw) + ((iw-1.0d0)/iw) * &
            transx(iat,iw+1)+x(iat,1)/iw
         transy(iat,iw) = y(iat,iw) + ((iw-1.0d0)/iw) * &
            transy(iat,iw+1)+y(iat,1)/iw
         transz(iat,iw) = z(iat,iw) + ((iw-1.0d0)/iw) * &
            transz(iat,iw+1)+z(iat,1)/iw
      enddo

   enddo

   end subroutine QtoX

!     This routine transforms cartesian to staging coordinates, which are stored
!     in trans matrices,values in x,y and z matrices are NOT modified!!      
   subroutine XtoQ(x,y,z,transx,transy,transz)
   real(DP),intent(inout)  :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out) :: transx(:,:),transy(:,:),transz(:,:)
   integer :: iat,iw

! (nwalk+1)th walker is identical to the first one in the polymer
   do iat=1,natom
      x(iat,nwalk+1) = x(iat,1)
      y(iat,nwalk+1) = y(iat,1)
      z(iat,nwalk+1) = z(iat,1)
   enddo

   do iat=1,natom

      transx(iat,1) = x(iat,1)
      transy(iat,1) = y(iat,1)
      transz(iat,1) = z(iat,1)

      do iw=2,nwalk
         transx(iat,iw) = ((iw-1.0d0) * x(iat,iw+1) + x(iat,1)) / iw
         transx(iat,iw) = x(iat,iw)-transx(iat,iw)
         transy(iat,iw) = ((iw-1.0d0) * y(iat,iw+1) + y(iat,1)) / iw
         transy(iat,iw) = y(iat,iw) - transy(iat,iw)
         transz(iat,iw) = ((iw-1.0d0) * z(iat,iw+1) + z(iat,1)) / iw
         transz(iat,iw) = z(iat,iw) - transz(iat,iw)
      enddo

   enddo

   end subroutine XtoQ


                                        
   subroutine init_mass(amg, amt)
   use mod_const, only: PI
   use mod_general, only: istage, inormalmodes, idebug
   use mod_system, ONLY: am
   real(DP),intent(out) :: amg(:,:), amt(:,:)
   real(DP)  :: lambda(size(amg,2)+1)
   integer :: iat, iw, k

!  NOTE: am was multiplied by amu earlier in the init.F90

!  Setting mass array, which is only used without the stage transformation
   do iw=1,nwalk
      do iat=1,natom
         amg(iat,iw)=am(iat)
         amt(iat,iw)=am(iat)
      enddo
   enddo

! transforming masses
! There are two sets of masses associated with the transformed coordinates
! amg according to eq. 90 in Tuckermann and amt (which has non-zero first mass)
! amt is associated with the kinetic energy, amg with the potential energy.
   if(istage.eq.1)then
      do iat=1,natom
         amg(iat,1)=0.0d0
         amt(iat,1)=am(iat)
         do iw=2,nwalk
            amg(iat,iw)=(iw/(iw-1.0d0))*am(iat)
            amt(iat,iw)=(iw/(iw-1.0d0))*am(iat)
         enddo
      enddo
   endif

!--SETTING MASSES FOR NORMAL MODES-----------  
   if(inormalmodes.eq.2)then ! this one apparently does not work
      do iat=1,natom
         lambda(1)=0
         do iw=2, nwalk
            k = iw -1
            lambda(iw) = 2 * sin(k * PI / nwalk)
            lambda(iw) = lambda(iw)**2 * nwalk
         enddo
!         do iw=2,nwalk
!            lambda(iw) = 2 * nwalk * (1-cos(PI*iw/nwalk))
!            lambda(iw+1)=lambda(iw)
!             lambda = 2 * nwalk * sin(PI*(iw-1)/nwalk)
!         enddo
         do iw=2,nwalk
            amg(iat,iw) = am(iat)*lambda(iw)
            amt(iat,iw) = am(iat)*lambda(iw)
         enddo
         amg(iat,1)=0.0d0
         amt(iat,1)=am(iat)
         if(idebug.gt.0)then
            write(*,*)'lambda',(lambda(iw),iw=1,nwalk)
            write(*,*)'amg',(amg(iat,iw),iw=1,nwalk)
            write(*,*)'amt',(amt(iat,iw),iw=1,nwalk)
!        write(*,*)'amg2',(amg(iat,iw)*2*sin(pi*(iw-1)/nwalk),iw=1,nwalk)
         end if
      enddo

   else if (inormalmodes.eq.1)then
      do iw=1,nwalk
         do iat=1,natom
            amt(iat,iw) = am(iat)
            amg(iat,iw) = am(iat)
         enddo
      enddo
   endif

   end subroutine init_mass

!     This routine transforms staging forces to cartesian forces, which are stored
!     in trans matrices,values in fx,fy and fz matrices are NOT modified!!
!     used in estimators.f90      
   subroutine FQtoFX(fx,fy,fz,transfx,transfy,transfz)
   real(DP),intent(in)    :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout) :: transfx(:,:),transfy(:,:),transfz(:,:)
   real(DP)  :: sumx, sumy, sumz
   integer   :: iat,iw


   do iat=1,natom

      sumx = 0.0d0
      sumy = 0.0d0
      sumz = 0.0d0
       
      do iw=nwalk,2,-1
         transfx(iat,iw) = fx(iat,iw) - (iw-2.0d0)/(iw-1.0d0) * fx(iat,iw-1) 
         transfy(iat,iw) = fy(iat,iw) - (iw-2.0d0)/(iw-1.0d0) * fy(iat,iw-1) 
         transfz(iat,iw) = fz(iat,iw) - (iw-2.0d0)/(iw-1.0d0) * fz(iat,iw-1) 

         sumx = sumx + transfx(iat,iw)
         sumy = sumy + transfy(iat,iw)
         sumz = sumz + transfz(iat,iw)
      enddo

      transfx(iat,1) = fx(iat,1) - sumx
      transfy(iat,1) = fy(iat,1) - sumy
      transfz(iat,1) = fz(iat,1) - sumz

   enddo

   return 

   end subroutine FQtoFX


! This routine transforms cartesian forces to staging forces, which are stored
! in fx matrices,values in fxab,fyab and fzab matrices are NOT modified!!
! used in force_clas     
! The classical force calculated in cartesian coordinates is transformed
! into staging coordinates.
! CAUTION, the formula 94 in Quantum simulations (Tuckermann) is incorrect
! dphi/dui on the output should be divided by P  (nwalk) both for the i=1 and the others
! futhermore, the dphi/dx(i-1) should actually be dphi/du(i-1)
! See Tuckerman's lecture notes. The force is therefore transformed 
! reccursively. sumux-z contains the (i-1) force used for this purpose.
   subroutine FXtoFQ(fxab,fyab,fzab,fx,fy,fz)
   real(DP),intent(inout)  :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(in)     :: fxab(:,:),fyab(:,:),fzab(:,:)
   real(DP)  :: sumux(size(fx,2)),sumuy(size(fx,2)),sumuz(size(fx,2))
   integer :: iat,iw

   do iat=1,natom
    
      do iw=1,nwalk 
         sumux(iw)=0.0d0
         sumuy(iw)=0.0d0
         sumuz(iw)=0.0d0
      enddo 

      do iw=1,nwalk 
         fx(iat,1)=fx(iat,1)+fxab(iat,iw)/nwalk
         fy(iat,1)=fy(iat,1)+fyab(iat,iw)/nwalk
         fz(iat,1)=fz(iat,1)+fzab(iat,iw)/nwalk
         sumux(1)=sumux(1)+fxab(iat,iw)
         sumuy(1)=sumuy(1)+fyab(iat,iw)
         sumuz(1)=sumuz(1)+fzab(iat,iw)
      enddo 


      do iw=2,nwalk
         sumux(iw)=sumux(iw)+fxab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumux(iw-1)
         sumuy(iw)=sumuy(iw)+fyab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumuy(iw-1)
         sumuz(iw)=sumuz(iw)+fzab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumuz(iw-1)

         fx(iat,iw)=fx(iat,iw)+(fxab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumux(iw-1))/nwalk 
         fy(iat,iw)=fy(iat,iw)+(fyab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumuy(iw-1))/nwalk 
         fz(iat,iw)=fz(iat,iw)+(fzab(iat,iw) + (iw-2.0d0)/(iw-1.0d0)*sumuz(iw-1))/nwalk 
      enddo
   enddo

   return
   end subroutine FXtoFQ


!     This routine transforms normal mode coordinates to cartesian coordinates, which are stored
!     in trans matrices,values in x,y and z matrices are NOT modified!!      
!     masses are also modified
   subroutine UtoX(x,y,z,transx,transy,transz)
#ifdef USEFFTW
      use mod_fftw3
#endif
   use mod_general, only: nwalk, inormalmodes
   real(DP), intent(in) :: x(:,:), y(:,:), z(:,:)
   real(DP), intent(out)   :: transx(:,:), transy(:,:), transz(:,:)
   integer  :: iat, iw
   real(DP) :: dnwalk, fac, equant
   integer  :: nmodes, odd

   if(inormalmodes.eq.1)then
      dnwalk = dsqrt(1.0d0 * nwalk)
      fac = dsqrt(2.0d0)
   else
      dnwalk = dsqrt(1.0d0 * nwalk)
     ! dnwalk = 1.0d0
      fac = 1.0d0
   end if

   nmodes = nwalk / 2
   odd = nwalk - 2 * nmodes ! 0 if even, 1 if odd

   call equant_nm(x, y, z, equant)

#ifdef USEFFTW
   do iat = 1, natom

      cx(1) = complex(x(iat,1),0)
      cy(1) = complex(y(iat,1),0)
      cz(1) = complex(z(iat,1),0)
      do iw = 2, nmodes + odd
         cx(iw) = complex(x(iat,iw), x(iat,nwalk+2-iw)) / fac
         cy(iw) = complex(y(iat,iw), y(iat,nwalk+2-iw)) / fac
         cz(iw) = complex(z(iat,iw), z(iat,nwalk+2-iw)) / fac
      end do
      if(odd.ne.1)then
         cx(nmodes+1) = complex(x(iat,nmodes+1),0)
         cy(nmodes+1) = complex(y(iat,nmodes+1),0)
         cz(nmodes+1) = complex(z(iat,nmodes+1),0)
      end if

!      cx((nwalk+2)/2) = complex(x(iat,nwalk), 0)
!      cy((nwalk+2)/2) = complex(y(iat,nwalk), 0)
!      cz((nwalk+2)/2) = complex(z(iat,nwalk), 0)
!      do iw = 2, nwalk / 2
!         cx(iw) = complex(x(iat,2*iw-2), x(iat,2*iw-1)) / fac
!         cy(iw) = complex(y(iat,2*iw-2), y(iat,2*iw-1)) / fac
!         cz(iw) = complex(z(iat,2*iw-2), z(iat,2*iw-1)) / fac
!      enddo
      call fftw_execute_dft_c2r(plan_utox, cx, x_in)
      call fftw_execute_dft_c2r(plan_utox, cy, y_in)
      call fftw_execute_dft_c2r(plan_utox, cz, z_in)

      do iw=1,nwalk
         transx(iat,iw) = x_in(iw) / dnwalk
         transy(iat,iw) = y_in(iw) / dnwalk
         transz(iat,iw) = z_in(iw) / dnwalk
      enddo

   enddo

   call equant_cart(transx, transy, transz, equant)

!  write(*,*)'normal modes back to cartesian'
!  call printf(transx/ang,transy/ang,transz/ang)
#else
   write(*,*)'FATAL ERROR: The program was not compiled with FFTW libraries.'
   write(*,*)'Normal mode transformations cannot be performed.'
   call abinerror('UtoX')
#endif

   end subroutine UtoX

!  This routine transforms cartesian to normal coordinates, which are stored
!  in trans matrices,values in x,y and z matrices are NOT modified!!
!  see https://github.com/i-pi/i-pi/blob/2a09a6d652b1ffe5f485c4c078c1085db6fcf63a/ipi/utils/nmtransform.py
   subroutine XtoU(x,y,z,transx,transy,transz)
#ifdef USEFFTW
   use mod_fftw3
#endif
   use mod_general, only: idebug, inormalmodes
   real(DP) :: x(:,:), y(:,:), z(:,:)
   real(DP) :: transx(:,:), transy(:,:), transz(:,:)
   integer  :: iat, iw
   real(DP) :: dnwalk, equant, fac
   integer  :: nmodes, odd

   if(inormalmodes.eq.1)then
      dnwalk = dsqrt(1.0d0 * nwalk)
      fac = dsqrt(2.0d0)
   else
      !dnwalk = nwalk
      dnwalk = dsqrt(1.0d0 * nwalk)
      fac = 1.0d0
   end if

   nmodes = nwalk / 2
   odd = nwalk - 2 * nmodes ! 0 if even, 1 if odd

!   call equant_cart(x, y, z, equant)

#ifdef USEFFTW
   do iat = 1, natom
      do iw = 1, nwalk
         x_in(iw) = x(iat,iw)
         y_in(iw) = y(iat,iw)
         z_in(iw) = z(iat,iw)
      enddo
      call fftw_execute_dft_r2c(plan_xtou,x_in,cx)
      call fftw_execute_dft_r2c(plan_xtou,y_in,cy)
      call fftw_execute_dft_r2c(plan_xtou,z_in,cz)

if(idebug.gt.1)then
      write(*,*)'complex coefficients'
      write(*,*)(cx(iw),iw=1,nwalk)
      write(*,*)(cy(iw),iw=1,nwalk)
      write(*,*)(cz(iw),iw=1,nwalk)
endif
      transx(iat,1) = realpart(cx(1))
      transy(iat,1) = realpart(cy(1))
      transz(iat,1) = realpart(cz(1))
      do iw=2,nmodes+odd
         transx(iat,iw) = realpart(cx(iw)) * fac
         transx(iat,nwalk+2-iw) = imagpart(cx(iw)) * fac
         transy(iat,iw) = realpart(cy(iw)) * fac
         transy(iat,nwalk+2-iw) = imagpart(cy(iw)) * fac
         transz(iat,iw) = realpart(cz(iw)) * fac
         transz(iat,nwalk+2-iw) = imagpart(cz(iw)) * fac
      enddo

      if(odd.ne.1)then
         transx(iat,nmodes+1) = realpart(cx((nwalk+2)/2))
         transy(iat,nmodes+1) = realpart(cy((nwalk+2)/2))
         transz(iat,nmodes+1) = realpart(cz((nwalk+2)/2))
      end if

   enddo

   transx = transx / dnwalk
   transy = transy / dnwalk
   transz = transz / dnwalk

!   call equant_nm(transx, transy, transz, equant)

!  write(*,*)'original cartesian to normal modes'
!  call printf(x/ang,y/ang,z/ang)
!  write(*,*)'transformed coordinates to normal modes'
!  call printf(transx/ang,transy/ang,transz/ang)

#else
   write(*,*)'FATAL ERROR: The program was not compiled with FFTW libraries.'
   write(*,*)'Normal mode transformations cannot be performed.'
   call abinerror('XtoU')
#endif


   end subroutine XtoU

   subroutine equant_cart(x, y, z, equant)
   use mod_general , only: nwalk, natom, idebug, inormalmodes
   use mod_system, only: am
   use mod_nhc,    only: temp
   use mod_utils, only: printf
   real(DP), intent(inout) :: x(:,:), y(:,:), z(:,:)
   real(DP), intent(out) :: equant
   real(DP) :: equantx, equanty, equantz, omega_n
   integer  :: iat, iw, k

   omega_n = NWALK * TEMP
   if (inormalmodes.eq.2) then
     omega_n = omega_n * sqrt(NWALK*1.d0)
   end if
   equantx = 0.0d0
   equanty = 0.0d0
   equantz = 0.0d0

   do iat = 1, natom
      x(iat, nwalk+1) = x(iat, 1)
      y(iat, nwalk+1) = y(iat, 1)
      z(iat, nwalk+1) = z(iat, 1)
   end do

   do iat = 1, natom
      do iw = 1, nwalk
         equantx = equantx + 0.5d0 * am(iat) * omega_n**2 * (x(iat,iw)-x(iat,iw+1))**2
         equanty = equanty + 0.5d0 * am(iat) * omega_n**2 * (y(iat,iw)-y(iat,iw+1))**2
         equantz = equantz + 0.5d0 * am(iat) * omega_n**2 * (z(iat,iw)-z(iat,iw+1))**2
      end do
!      write(*,*)"Quantum energy per atom in cartesian coordinates"
!      write(*,*)equantx, equanty, equantz, equantx+equanty+equantz
   end do

   equant = equantx + equanty + equantz

   if(idebug.eq.1)then
      write(*,*)"Quantum energy in cartesian coordinates"
      write(*,*)equantx, equanty, equantz, equant
      write(*,*)"Cartesian coordinates"
      call printf(x,y,z)
   end if

   end subroutine equant_cart

   subroutine equant_nm(x, y, z, equant)
   use mod_const, only: PI
   use mod_general , only: nwalk, natom, idebug, inormalmodes
   use mod_system,   only: am
   use mod_nhc,    only: temp
   use mod_utils, only: printf
   use mod_arrays, only: amg
   real(DP), intent(in) :: x(:,:), y(:,:), z(:,:)
   real(DP), intent(out) :: equant
   real(DP) :: equantx, equanty, equantz, omega_n, omega(size(x,2))
   integer  :: iat, iw, k

   omega_n = NWALK * TEMP
   ! Tuckerman Hamiltonian
   if (inormalmodes.eq.2) then
     omega_n = omega_n * dsqrt(NWALK*1.d0)
   end if
   equantx = 0.0d0
   equanty = 0.0d0
   equantz = 0.0d0

   do iw=1, nwalk
      k = iw -1
      omega(iw) = 2 * omega_n * sin(k * PI / nwalk)
   enddo

   do iat = 1, natom
      do iw = 1, nwalk
         if(inormalmodes.eq.1)then
            equantx = equantx + 0.5d0 * am(iat) * omega(iw)**2 * (x(iat,iw))**2
            equanty = equanty + 0.5d0 * am(iat) * omega(iw)**2 * (y(iat,iw))**2
            equantz = equantz + 0.5d0 * am(iat) * omega(iw)**2 * (z(iat,iw))**2
         else
            equantx = equantx + 0.5d0 * amg(iat,iw) * omega_n**2 * (x(iat,iw))**2
            equanty = equanty + 0.5d0 * amg(iat,iw) * omega_n**2 * (y(iat,iw))**2
            equantz = equantz + 0.5d0 * amg(iat,iw) * omega_n**2 * (z(iat,iw))**2
         end if
      end do
!      write(*,*)"Quantum energy in normal modes coordinates"
!      write(*,*)equantx, equanty, equantz, equantx+equanty+equantz
   end do

   equant = equantx + equanty + equantz

   if(idebug.eq.1)then
      write(*,*)'omegas', (omega(iw),iw=1,nwalk)
      write(*,*)"Quantum energy in normal modes coordinates"
      write(*,*)equantx, equanty, equantz, equant
      write(*,*)"Normal mode coordinates"
      call printf(x, y, z)
   end if

   end subroutine equant_nm

end module mod_transform

