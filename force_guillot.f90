!-----BASED ON:
!Quantum effects in simulated water by the Feynman&Hibbs approach,
!Bertrand Guillot and Yves Guissani,J. Chem. Phys. 108, 10162 (1998);  

      subroutine force_guillot(x,y,z,fx,fy,fz,eclas)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: inames
      implicit real*8(a-h,o-z)
      real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(out) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8,intent(out) :: eclas
      real*8 :: temp1,r,fr
!      integer :: inames(npartmax)

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
         r=(x(i,k)-x(j,k))**2+(y(i,k)-y(j,k))**2+(z(i,k)-z(j,k))**2
         r=dsqrt(r)
         if(inames(i).eq.0.and.inames(j).eq.0)then
          temp1=temp1+voo(r)
!          tempoo=tempoo+voo(r)
          fr=froo(r) 
         endif
         if(inames(i).eq.0.and.inames(j).eq.1.or.inames(i).eq.1&
      .and.inames(j).eq.0)then
          temp1=temp1+voh(r)
          fr=froh(r) 
!          tempoh=tempoh+voh(r)
         endif
         if(inames(i).eq.2.and.inames(j).eq.0.or.inames(i).eq.0&
       .and.inames(j).eq.2)then
          temp1=temp1+vocl(r)
          fr=frocl(r) 
         endif
         if(inames(i).eq.2.and.inames(j).eq.1.or.inames(i).eq.1&
      .and.inames(j).eq.2)then
          temp1=temp1+vhcl(r)
          fr=frhcl(r) 
         endif

         if(inames(i).eq.1.and.inames(j).eq.1)then
          temp1=temp1+vhh(r)
          fr=frhh(r) 
!          temphh=temphh+vhh(r)
         endif
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
      vocl=(627.14*qo/1.89)/r
      sigma=3.6835677
      rred=sigma/r
      epsilon=0.1230367
      vocl=vocl+4*epsilon*(rred**12-rred**6)

      r=r*1.89
      vocl=vocl/627.14
      return
      end

      real*8 function voclIS(r)
      implicit real*8(a-h,o-z)
      r=r/1.89
      voclIS=0.0d0
      sigma=3.6835677
      rred=sigma/r
      epsilon=0.1230367
      voclIS=voclIS+4*epsilon*(rred**12-rred**6)


      r=r*1.89
      voclIS=voclIS/627.14
      return
      end

      real*8 function frocl(r)
      implicit real*8(a-h,o-z)

      r=r/1.89
      qo=0.66d0
      frocl=-(627.14*qo/1.89)/(r**2)
      sigma=3.6835677
      rred=sigma/r
      epsilon=0.1230367
      frocl=frocl+0.5d0*(-12*rred**13+6*rred**7)/sigma

      r=r*1.89
      frocl=frocl/627.14
      frocl=-frocl/1.89d0
      return
      end

      real*8 function vhcl(r)
      implicit real*8(a-h,o-z)
      vocl=0.0d0
      r=r/1.89
      qh=0.66d0/2.0d0
      vhcl=-(627.14*qh/1.89)/r
      sigma=2.0000000
      rred=sigma/r
      epsilon=0.1230367
      vocl=vocl+4*epsilon*(rred**12-rred**6)
      r=r*1.89
      vhcl=vhcl/627.14
      return
      end

      real*8 function vhclIS(r)
      implicit real*8(a-h,o-z)
      r=r/1.89
      vhclIS=0.0d0
      sigma=2.0000000
      rred=sigma/r
      epsilon=0.1230367
      vhclIS=vhclIS+4*epsilon*(rred**12-rred**6)
      r=r*1.89
      vhclIS=vhclIS/627.14
      return
      end


      real*8 function frhcl(r)
      implicit real*8(a-h,o-z)
      r=r/1.89
      qh=0.66d0/2.0d0
      frhcl=(627.14*qh/1.89)/(r**2)
      sigma=2.0000000
      rred=sigma/r
      epsilon=0.1230367
      frhcl=frhcl+0.5d0*(-12*rred**13+6*rred**7)/sigma
      r=r*1.89
      frhcl=frhcl/627.14
      frhcl=-frhcl/1.89d0
      return
      end

                                                    
      real*8 function voo(r)
      implicit real*8(a-h,o-z)
      r=r/1.89

      voo=144.358/r
      sigma=3.74
      rred=sigma/r
      voo=voo+0.5d0*(rred**8-rred**6)

      r=r*1.89
      voo=voo/627.14
!      write(*,*)'voo1',voo
      return
      end

      real*8 function froo(r)
      implicit real*8(a-h,o-z)
      r=r/1.89

      froo=-144.358/(r**2)
      sigma=3.74
      rred=sigma/r
      froo=froo+0.5d0*(-8*rred**9+6*rred**7)/sigma

      r=r*1.89
      froo=froo/627.14
      froo=-froo/1.89d0
!      write(*,*)'froo',froo
      return
      end

      real*8 function voh(r)
      implicit real*8(a-h,o-z)
      r=r/1.89

!      do i=1,100
!      r=0.5*0.05*(i-1)

      c=72.269
      voh=-c/r
      expon=0.9*(r-2.2)
      expon1=-5.8*(r-1.07)
      voh=voh-4.3/(1+exp(expon))
      voh=voh+3.9*((exp(expon1)-1)**2-1)

      r=r*1.89
      voh=voh/627.14
!      write(*,*)'voh',r,voh
!      enddo
!      stop
!      return
      end

      real*8 function froh(r)
      implicit real*8(a-h,o-z)
      r=r/1.89

      c=72.269
      froh=c/(r**2)
      expon=0.9*(r-2.2)
      expon1=-5.8*(r-1.07)
      froh=froh+(4.3/(1+exp(expon))**2)*0.9*exp(expon)
      froh=froh-3.9*2*5.8*(exp(expon1)-1)*exp(expon1)

      r=r*1.89
      froh=froh/627.14
      froh=-froh/1.89
!      write(*,*)'froh',froh
      return
      end

      real*8 function vhh(r)
      implicit real*8(a-h,o-z)
      r=r/1.89

      c=36.1345
      vhh=c/r
      expon=3.1*(r-2.05)
      expon1=-6.0*(r-1.495)
      vhh=vhh+17/(1+exp(expon))
      vhh=vhh+13*((exp(expon1)-1)**2-1)

      r=r*1.89
      vhh=vhh/627.14
!      write(*,*)'vhh',r,vhh
      return
      end

      real*8 function frhh(r)
      implicit real*8(a-h,o-z)
      r=r/1.89
!      write(*,*)'r in frhh',r

      c=36.1345
      frhh=-c/(r**2)
!      write(*,*)'frhh in frhh',frhh
      expon=3.1*(r-2.05)
      expon1=-6.0*(r-1.495)
      frhh=frhh-(17.0d0/(1+exp(expon))**2)*3.1*exp(expon)
!      write(*,*)'frhh in frhh',frhh
      frhh=frhh-13.0*2*6.0*(exp(expon1)-1)*exp(expon1)
!      write(*,*)'frhh in frhh',frhh

      r=r*1.89
      frhh=frhh/627.14
      frhh=-frhh/1.89
!      write(*,*)'frhh',frhh
      return
      end



