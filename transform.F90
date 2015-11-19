!---mod_transform                  P.Slavicek and D.Hollas, 9.2.2012
!--- Transformation routines for PIMD and PI+GLE.
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
   contains

!     This routine transforms staging coordinates to cartesian coordinates, which are stored
!     in trans matrices,values in x,y and z matrices are NOT modified!!      
   subroutine QtoX(x,y,z,transx,transy,transz)
      real(DP),intent(inout)  :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out)    :: transx(:,:),transy(:,:),transz(:,:)
      integer :: iat,iw

      do iat=1,natom
        x(iat,nwalk+1)=x(iat,1)
        y(iat,nwalk+1)=y(iat,1)
        z(iat,nwalk+1)=z(iat,1)

        transx(iat,1)=x(iat,1)
        transy(iat,1)=y(iat,1)
        transz(iat,1)=z(iat,1)
        transx(iat,nwalk+1)=x(iat,1)
        transy(iat,nwalk+1)=y(iat,1)
        transz(iat,nwalk+1)=z(iat,1)

       do iw=nwalk,2,-1
         transx(iat,iw)=x(iat,iw)+((iw-1.0d0)/iw)*transx(iat,iw+1)+x(iat,1)/iw
         transy(iat,iw)=y(iat,iw)+((iw-1.0d0)/iw)*transy(iat,iw+1)+y(iat,1)/iw
         transz(iat,iw)=z(iat,iw)+((iw-1.0d0)/iw)*transz(iat,iw+1)+z(iat,1)/iw

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
       x(iat,nwalk+1)=x(iat,1)
       y(iat,nwalk+1)=y(iat,1)
       z(iat,nwalk+1)=z(iat,1)
      enddo

      do iat=1,natom
      transx(iat,1)=x(iat,1)
      transy(iat,1)=y(iat,1)
      transz(iat,1)=z(iat,1)
       do iw=2,nwalk
       transx(iat,iw)=((iw-1.0d0)*x(iat,iw+1)+x(iat,1))/iw
       transx(iat,iw)=x(iat,iw)-transx(iat,iw)
       transy(iat,iw)=((iw-1.0d0)*y(iat,iw+1)+y(iat,1))/iw
       transy(iat,iw)=y(iat,iw)-transy(iat,iw)
       transz(iat,iw)=((iw-1.0d0)*z(iat,iw+1)+z(iat,1))/iw
       transz(iat,iw)=z(iat,iw)-transz(iat,iw)
       enddo
      enddo

   end subroutine XtoQ


                                        
   subroutine init_mass(amg,amt)
      use mod_const, only: PI
      use mod_general, only: istage
      use mod_system, ONLY: am
      real(DP),intent(out) :: amg(:,:),amt(:,:)
      real(DP)  :: lambda(size(amg,2)+1)
      integer :: iat,iw

! NOTE: am was multiplied by amu earlier in the init.f


! Setting mass array, which is only used without the stage transformation
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

!---SETTING MASSES FOR NORMAL MODES-----------  
      if(istage.eq.2)then   
       do iat=1,natom
        lambda(1)=0
        do iw=2,nwalk,2
         lambda(iw)=2*nwalk*(1-cos(PI*iw/nwalk))
         lambda(iw+1)=lambda(iw)
        enddo
        do iw=2,nwalk
         amg(iat,iw)=am(iat)*lambda(iw)
         amt(iat,iw)=amg(iat,iw)
        enddo
        amg(iat,1)=0.0d0
        amt(iat,1)=am(iat)
        write(*,*)'lambda',(lambda(iw),iw=1,nwalk)
        write(*,*)'amg',(amg(iat,iw),iw=1,nwalk)
        write(*,*)'amt',(amt(iat,iw),iw=1,nwalk)
!        write(*,*)'amg2',(amg(iat,iw)*2*sin(pi*(iw-1)/nwalk),iw=1,nwalk)
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
       sumx=0.0d0
       sumy=0.0d0
       sumz=0.0d0
       
       do iw=nwalk,2,-1
        transfx(iat,iw)=fx(iat,iw)-(iw-2.0d0)/(iw-1.0d0)*fx(iat,iw-1) 
        transfy(iat,iw)=fy(iat,iw)-(iw-2.0d0)/(iw-1.0d0)*fy(iat,iw-1) 
        transfz(iat,iw)=fz(iat,iw)-(iw-2.0d0)/(iw-1.0d0)*fz(iat,iw-1) 

        sumx=sumx+transfx(iat,iw)
        sumy=sumy+transfy(iat,iw)
        sumz=sumz+transfz(iat,iw)
       enddo

       transfx(iat,1)=fx(iat,1)-sumx
       transfy(iat,1)=fy(iat,1)-sumy
       transfz(iat,1)=fz(iat,1)-sumz
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
      real(DP)  :: x(:,:), y(:,:), z(:,:)
      real(DP)  :: transx(:,:), transy(:,:), transz(:,:)
      integer   :: iat, iw

#ifdef USEFFTW
       do iat = 1, natom

       cx(1) = complex(x(iat,1),0)
       cy(1) = complex(y(iat,1),0)
       cz(1) = complex(z(iat,1),0)
       cx((nwalk+2)/2) = complex(x(iat,nwalk), 0)
       cy((nwalk+2)/2) = complex(y(iat,nwalk), 0)
       cz((nwalk+2)/2) = complex(z(iat,nwalk), 0)
       do iw = 2, nwalk/2
        cx(iw) = complex(x(iat,2*iw-2), x(iat,2*iw-1))
        cy(iw) = complex(y(iat,2*iw-2), y(iat,2*iw-1))
        cz(iw) = complex(z(iat,2*iw-2), z(iat,2*iw-1))
       enddo
      call fftw_execute_dft_c2r(plan_utox, cx, x_in)
      call fftw_execute_dft_c2r(plan_utox, cy, y_in)
      call fftw_execute_dft_c2r(plan_utox, cz, z_in)

      do iw=1,nwalk
       transx(iat,iw) = x_in(iw)
       transy(iat,iw) = y_in(iw)
       transz(iat,iw) = z_in(iw)
      enddo

      enddo

!      write(*,*)'normal modes back to cartesian'
!      call printf(transx/ang,transy/ang,transz/ang)
#else
      write(*,*)'FATAL ERROR: The program was not compiled with FFTW libraries.'
      write(*,*)'Normal mode transformations cannot be performed.'
      call abinerror('UtoX')
#endif

   end subroutine UtoX

!  This routine transforms cartesian to normal coordinates, which are stored
!  in trans matrices,values in x,y and z matrices are NOT modified!!      
   subroutine XtoU(x,y,z,transx,transy,transz)
#ifdef USEFFTW
      use mod_fftw3
#endif
      use mod_general, only: idebug
      real(DP)  :: x(:,:), y(:,:), z(:,:)
      real(DP)  :: transx(:,:), transy(:,:), transz(:,:)
      integer   :: iat, iw

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

if(idebug.eq.1)then
       write(*,*)'complex coefficients'
       write(*,*)(cx(iw),iw=1,nwalk)
       write(*,*)(cy(iw),iw=1,nwalk)
       write(*,*)(cz(iw),iw=1,nwalk)
endif
       transx(iat,1) = realpart(cx(1)) / nwalk
       transy(iat,1) = realpart(cy(1)) / nwalk
       transz(iat,1) = realpart(cz(1)) / nwalk
       transx(iat,nwalk) = realpart(cx((nwalk+2)/2)) / nwalk
       transy(iat,nwalk) = realpart(cy((nwalk+2)/2)) / nwalk
       transz(iat,nwalk) = realpart(cz((nwalk+2)/2)) / nwalk
       do iw=2,nwalk/2
       transx(iat,2*iw-2) = realpart(cx(iw)) / nwalk
       transx(iat,2*iw-1) = imagpart(cx(iw)) / nwalk
       transy(iat,2*iw-2) = realpart(cy(iw)) / nwalk
       transy(iat,2*iw-1) = imagpart(cy(iw)) / nwalk
       transz(iat,2*iw-2) = realpart(cz(iw)) / nwalk
       transz(iat,2*iw-1) = imagpart(cz(iw)) / nwalk
       enddo
      enddo
!       write(*,*)'original cartesian to normal modes'
!       call printf(x/ang,y/ang,z/ang)
!       write(*,*)'transformed coordinates to normal modes'
!       call printf(transx/ang,transy/ang,transz/ang)

#else
      write(*,*)'FATAL ERROR: The program was not compiled with FFTW libraries.'
      write(*,*)'Normal mode transformations cannot be performed.'
      call abinerror('XtoU')
#endif


   end subroutine XtoU

end module mod_transform
