      subroutine force_quantum(fx,fy,fz,x,y,z,amg,energy)
      use mod_array_size
      use mod_general
      use mod_nhc, ONLY: temp,inose
      implicit none
      real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8,intent(in) :: amg(npartmax,nwalkmax)
      real*8,intent(inout) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      real*8,intent(out) :: energy
      real*8 :: ak(npartmax,nwalkmax)
      real*8 :: equant
      integer :: iat,iw,i,j,kplus,kminus

      do iat=1,natom
       do iw=1,nwalk
        fx(iat,iw)=0.0d0
        fy(iat,iw)=0.0d0
        fz(iat,iw)=0.0d0
       enddo
      enddo

! TODO: ak parametry se nemuseji pocitat pokazde znova
! Setting the quantum force constants
! ak is defined is m*P/(beta^2*hbar^2)
        do iw=1,nwalk
         do i=1,natom
          ak(i,iw)=nwalk*amg(i,iw)*TEMP**2
         enddo
        enddo
! for PI+GLE we have different hamiltonian
       if(inose.eq.2)then
        do iw=1,nwalk
         do i=1,natom
          ak(i,iw)=nwalk*ak(i,iw)
         enddo
        enddo
       endif

      equant=0.0d0


! If the staging transformation is not used
      if(istage.eq.0)then
       do j=1,natom 
        do i=1,nwalk 
         kplus=i+1
         kminus=i-1
         if(i.eq.1)then
          kminus=nwalk
         endif
         if(i.eq.nwalk)then
          kplus=1
         endif
         fx(j,i)=(x(j,i)-x(j,kplus))
         fx(j,i)=fx(j,i)+(x(j,i)-x(j,kminus))
         fx(j,i)=-fx(j,i)*ak(j,i)
         fy(j,i)=(y(j,i)-y(j,kplus))
         fy(j,i)=fy(j,i)+(y(j,i)-y(j,kminus))
         fy(j,i)=-fy(j,i)*ak(j,i)
         fz(j,i)=(z(j,i)-z(j,kplus))
         fz(j,i)=fz(j,i)+(z(j,i)-z(j,kminus))
         fz(j,i)=-fz(j,i)*ak(j,i)
         equant=equant+0.5*ak(j,i)*(x(j,i)-x(j,kplus))**2
         equant=equant+0.5*ak(j,i)*(y(j,i)-y(j,kplus))**2
         equant=equant+0.5*ak(j,i)*(z(j,i)-z(j,kplus))**2
        enddo
       enddo
      endif

! If the staging transformation is used
       if(istage.eq.1.or.istage.eq.2)then
        do j=1,natom 
         do i=1,nwalk 
          fx(j,i)=-ak(j,i)*x(j,i)   
          fy(j,i)=-ak(j,i)*y(j,i)       
          fz(j,i)=-ak(j,i)*z(j,i) 
          equant=equant+0.5d0*ak(j,i)*x(j,i)**2
          equant=equant+0.5d0*ak(j,i)*y(j,i)**2
          equant=equant+0.5d0*ak(j,i)*z(j,i)**2
         enddo
        enddo
       endif 

      energy=equant
      
      end
                                        
