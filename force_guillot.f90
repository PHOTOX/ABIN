!-----BASED ON:
!Quantum effects in simulated water by the Feynman&Hibbs approach,
!Bertrand Guillot and Yves Guissani,J. Chem. Phys. 108, 10162 (1998);  

module mod_guillot
   use mod_array_size, ONLY: ANG,AUTOKCAL
   private
   public :: force_guillot,inames_guillot

   CONTAINS
      subroutine force_guillot(x,y,z,fx,fy,fz,eclas)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: inames
      implicit real*8(a-h,o-z)
      real*8,intent(in)  :: x(:,:),y(:,:),z(:,:)
      real*8,intent(out) :: fx(:,:),fy(:,:),fz(:,:)
      real*8,intent(out) :: eclas
      real*8 :: temp1,r,fr
!     integer :: inames(npartmax)

      eclas=0.0d0


      do i=1,natom
!       if(names(i).eq.'H')then
!        inames(i)=1
!       else if(names(i).eq.'O')then
!        inames(i)=0
!       else if (names(i).eq.'CL')then
!        inames(i)=2
!       endif
       do k=1,nwalk
       fx(i,k)=0.0d0
       fy(i,k)=0.0d0
       fz(i,k)=0.0d0
       enddo
      enddo

!TODO: IS OMP really working here?
!disabled for safety
!!$OMP PARALLEL DO PRIVATE(temp1,r,frb,fr)
      do k=1,nwalk
       temp1=0.0d0
!       tempoo=0.0d0
!       tempoh=0.0d0
!       temphh=0.0d0
       do i=1,natom
        do j=i+1,natom
         fr=0.0d0
         r=(x(i,k)-x(j,k))**2+(y(i,k)-y(j,k))**2+(z(i,k)-z(j,k))**2
         r=dsqrt(r)
         if(inames(i).eq.0.and.inames(j).eq.0)then
          temp1=temp1+voo(r)
!          tempoo=tempoo+voo(r)
          fr=froo(r) 
       else if(inames(i).eq.0.and.inames(j).eq.1.or.inames(i).eq.1&
      .and.inames(j).eq.0)then
          temp1=temp1+voh(r)
          fr=froh(r) 
!          tempoh=tempoh+voh(r)
        else if(inames(i).eq.2.and.inames(j).eq.0.or.inames(i).eq.0&
       .and.inames(j).eq.2)then
          temp1=temp1+vocl(r)
          fr=frocl(r) 
      else if(inames(i).eq.2.and.inames(j).eq.1.or.inames(i).eq.1&
      .and.inames(j).eq.2)then
          temp1=temp1+vhcl(r)
          fr=frhcl(r) 
      else if(inames(i).eq.1.and.inames(j).eq.1)then
          temp1=temp1+vhh(r)
          fr=frhh(r) 
!          temphh=temphh+vhh(r)
      else 
         write(*,*)'Unrecognized atom pair in force guillot!'
         call abinerror('force_guillot')
      end if
         fx(i,k)=fx(i,k)+fr*(x(i,k)-x(j,k))/(r)
         fx(j,k)=fx(j,k)+fr*(x(j,k)-x(i,k))/(r)
         fy(i,k)=fy(i,k)+fr*(y(i,k)-y(j,k))/(r)
         fy(j,k)=fy(j,k)+fr*(y(j,k)-y(i,k))/(r)
         fz(i,k)=fz(i,k)+fr*(z(i,k)-z(j,k))/(r)
         fz(j,k)=fz(j,k)+fr*(z(j,k)-z(i,k))/(r)
       enddo
      enddo
!!OMP ATOMIC 
       eclas=eclas+temp1

      enddo
!!$OMP END PARALLEL DO
      eclas=eclas/nwalk
!      write(*,*)"CLAS"
!      write(*,*)eclas,tempoo,tempoh,temphh
   
      end
                                        
      real*8 function vocl(r)
      use mod_array_size
      implicit real*8(a-h,o-z)
      r=r/ang
      qo=0.66d0
      vocl=(AUTOKCAL*qo/ANG)/r
      sigma=3.6835677d0
      rred=sigma/r
      epsilon=0.1230367d0
      vocl=vocl+4*epsilon*(rred**12-rred**6)

      r=r*ANG
      vocl=vocl/AUTOKCAL
      return
      end

      real*8 function voclIS(r)
      implicit real*8(a-h,o-z)
      r=r/ANG
      voclIS=0.0d0
      sigma=3.6835677d0
      rred=sigma/r
      epsilon=0.1230367d0
      voclIS=voclIS+4*epsilon*(rred**12-rred**6)


      r=r*ANG
      voclIS=voclIS/AUTOKCAL
      return
      end

      real*8 function frocl(r)
      implicit real*8(a-h,o-z)

      r=r/ANG
      qo=0.66d0
      frocl=-(AUTOKCAL*qo/ANG)/(r**2)
      sigma=3.6835677d0
      rred=sigma/r
      epsilon=0.1230367d0
      frocl=frocl+0.5d0*(-12*rred**13+6*rred**7)/sigma

      r=r*ANG
      frocl=frocl/AUTOKCAL
      frocl=-frocl/ANG
      return
      end

      real*8 function vhcl(r)
      implicit real*8(a-h,o-z)
      real*8 :: vocl
      vocl=0.0d0
      r=r/ANG
      qh=0.66d0/2.0d0
      vhcl=-(AUTOKCAL*qh/ANG)/r
      sigma=2.0000000d0
      rred=sigma/r
      epsilon=0.1230367d0
      vocl=vocl+4*epsilon*(rred**12-rred**6)
      r=r*ANG
      vhcl=vhcl/AUTOKCAL
      return
      end

      real*8 function vhclIS(r)
      implicit real*8(a-h,o-z)
      r=r/ANG
      vhclIS=0.0d0
      sigma=2.0000000
      rred=sigma/r
      epsilon=0.1230367
      vhclIS=vhclIS+4*epsilon*(rred**12-rred**6)
      r=r*ANG
      vhclIS=vhclIS/AUTOKCAL
      return
      end


      real*8 function frhcl(r)
      implicit real*8(a-h,o-z)
      r=r/ANG
      qh=0.66d0/2.0d0
      frhcl=(AUTOKCAL*qh/ANG)/(r**2)
      sigma=2.0000000d0
      rred=sigma/r
      epsilon=0.1230367d0
      frhcl=frhcl+0.5d0*(-12*rred**13+6*rred**7)/sigma
      r=r*ANG
      frhcl=frhcl/AUTOKCAL
      frhcl=-frhcl/ANG
      return
      end

                                                    
      real*8 function voo(r)
      implicit real*8(a-h,o-z)
      r=r/ANG

      voo=144.358d0/r
      sigma=3.74d0
      rred=sigma/r
      voo=voo+0.5d0*(rred**8-rred**6)

      r=r*ANG
      voo=voo/AUTOKCAL
!      write(*,*)'voo1',voo
      return
      end

      real*8 function froo(r)
      implicit real*8(a-h,o-z)
      r=r/ANG

      froo=-144.358d0/(r**2)
      sigma=3.74d0
      rred=sigma/r
      froo=froo+0.5d0*(-8*rred**9+6*rred**7)/sigma

      r=r*ANG
      froo=froo/AUTOKCAL
      froo=-froo/ANG
!      write(*,*)'froo',froo
      return
      end

      real*8 function voh(r)
      implicit real*8(a-h,o-z)
      r=r/ANG

!      do i=1,100
!      r=0.5d0*0.05d0*(i-1)

      c=72.269d0
      voh=-c/r
      expon=0.9d0*(r-2.2d0)
      expon1=-5.8d0*(r-1.07d0)
      voh=voh-4.3d0/(1+exp(expon))
      voh=voh+3.9d0*((exp(expon1)-1)**2-1)

      r=r*ANG
      voh=voh/AUTOKCAL
!      write(*,*)'voh',r,voh
!      enddo
!      stop
!      return
      end

      real*8 function froh(r)
      implicit real*8(a-h,o-z)
      r=r/ANG

      c=72.269d0
      froh=c/(r**2)
      expon=0.9d0*(r-2.2d0)
      expon1=-5.8d0*(r-1.07d0)
      froh=froh+(4.3d0/(1+exp(expon))**2)*0.9d0*exp(expon)
      froh=froh-3.9d0*2d0*5.8d0*(exp(expon1)-1)*exp(expon1)

      r=r*ANG
      froh=froh/AUTOKCAL
      froh=-froh/ANG
!      write(*,*)'froh',froh
      return
      end

      real*8 function vhh(r)
      implicit real*8(a-h,o-z)
      r=r/ANG

      c=36.1345d0
      vhh=c/r
      expon=3.1d0*(r-2.05d0)
      expon1=-6.0d0*(r-1.495d0)
      vhh=vhh+17/(1+exp(expon))
      vhh=vhh+13*((exp(expon1)-1)**2-1)

      r=r*ANG
      vhh=vhh/AUTOKCAL
!      write(*,*)'vhh',r,vhh
      return
      end

      real*8 function frhh(r)
      implicit real*8(a-h,o-z)
      r=r/ANG
!      write(*,*)'r in frhh',r

      c=36.1345d0
      frhh=-c/(r**2)
!      write(*,*)'frhh in frhh',frhh
      expon=3.1d0*(r-2.05)
      expon1=-6.0d0*(r-1.495d0)
      frhh=frhh-(17.0d0/(1+exp(expon))**2)*3.1d0*exp(expon)
!      write(*,*)'frhh in frhh',frhh
      frhh=frhh-13.0d0*2.0d0*6.0d0*(exp(expon1)-1)*exp(expon1)
!      write(*,*)'frhh in frhh',frhh

      r=r*ANG
      frhh=frhh/AUTOKCAL
      frhh=-frhh/ANG
!      write(*,*)'frhh',frhh
      return
   end
!---potentially useful for guillot and other empirical force fields, because string
!---comparison is very cpu demanding!!
!--TODO:move to guillot
   subroutine inames_guillot()
      use mod_general, ONLY:natom
      use mod_system, ONLY:inames,names
      implicit none
      integer :: i
      do i=1,natom
       if(names(i).eq.'H')then
        inames(i)=1
       else if(names(i).eq.'O')then
        inames(i)=0
       else if (names(i).eq.'CL')then
        inames(i)=2
       endif
      enddo
   end subroutine

END MODULE mod_GUILLOT



