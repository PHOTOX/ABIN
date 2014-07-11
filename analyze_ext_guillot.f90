!-----BASED ON:
!Quantum effects in simulated water by the Feynman&Hibbs approach,
!Bertrand Guillot and Yves Guissani,J. Chem. Phys. 108, 10162 (1998);  

      module mod_analyze_ext
      use mod_array_size
      real*8  :: spec(nbinmax)
      real*8  :: emin=0,emax=0.35
      integer :: nbinen=400
      save
      contains
      subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
      use mod_array_size
      use mod_general
      use mod_system, ONLY: inames
      implicit real*8(a-h,o-z)
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
      real*8 amt(npartmax,nwalkmax)

      de=(emax-emin)/nbinen
      
      do k=1,nwalk
!       tempoo=0.0d0
!       tempoh=0.0d0
!       temphh=0.0d0
       tempGS=0.0d0
       tempIS=0.0d0
       do i=1,natom
        do j=i+1,natom
         r=(x(i,k)-x(j,k))**2+(y(i,k)-y(j,k))**2+(z(i,k)-z(j,k))**2
         r=sqrt(r)
!         if(inames(i).eq.0.and.inames(j).eq.0)then
!          tempoo=tempoo+voo(r)
!         endif
!         if(inames(i).eq.0.and.inames(j).eq.1.or.inames(i).eq.1&
!      .and.inames(j).eq.0)then
!          tempoh=tempoh+voh(r)
!         endif

!       inames=0 ....'O'
!       inames=1.... 'H'
!       inames=2 .... 'Cl'
! see modules.f90 inames_guillot() subroutine
         if(inames(i).eq.2.and.inames(j).eq.0.or.inames(i).eq.0&
       .and.inames(j).eq.2)then
          tempGS=tempGS+vocl(r)
          tempIS=tempIS+voclIS(r)
         endif
         if(inames(i).eq.2.and.inames(j).eq.1.or.inames(i).eq.1&
      .and.inames(j).eq.2)then
          tempGS=tempGS+vhcl(r)
          tempIS=tempIS+vhclIS(r)
         endif

!         if(inames(i).eq.1.and.inames(j).eq.1)then
!          temphh=temphh+vhh(r)
!         endif
       enddo
      enddo
      
       indexSPEC=(tempIS-tempGS)/de
       if(indexSPEC.gt.0.and.indexSPEC.lt.nbinen)then
!       write(*,*)'Energy difference',(tempIS-tempGS),indexSPEC
        spec(indexSPEC)=spec(indexSPEC)+1 
       else
        write(*,*)'Energy difference',(tempIS-tempGS),indexSPEC,de
        write(*,*)'Problem with indexSPEC'
        IndexSPEC=1
!       call abinerror('analyze_ext')
       endif
      enddo
      

      if(modulo(it,nwrite).eq.0)then
       open(200,file='spec.dat')
       do ibin=1,nbinen
        write(200,*)emin+(ibin-1)*de,spec(ibin)/it/nwalk
       enddo
       close(200)
      endif
   
      end subroutine analyze_ext



      end module analyze_ext
                                        

