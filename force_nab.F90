      module mod_nab
      use mod_const, only: DP
      implicit none
      real(DP) :: boxx,boxy,boxz
      real(DP) :: boxx2,boxy2,boxz2
      real(DP) :: alpha_pme=-1,kappa_pme=-1,cutoff=100.0d0
      real(DP) :: epsinf=3d33
      real(DP), allocatable :: charges(:)
#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4 
      integer, allocatable  :: natmol(:) !used for wrapping, independent of nmolt
#else
      integer  :: natmol(2000) !hardcoded for older gfortran
#endif
      integer  :: nmol=1 !used for wraping, independent of nmolt
      integer  :: ipbc=0, nsnb=1
      integer  :: ips=0 ! 1-both LJ and coul, 2-coul 3- LJ
      save
      contains

      subroutine wrap(x,y,z)
      use mod_general,only:nwalk
      real(DP) :: x(:,:),y(:,:),z(:,:)
      integer  :: i,iat,iat2,iw,iww,iwrap

      iwrap=0

      iat=1
      ! wrapping based on the 1st bead
      iww=1
      do i=1,nmol

       if (x(iat,iww).gt.boxx2)then
               iwrap=iwrap+1
               do iat2=iat,iat+natmol(i)-1
                do iw=1,nwalk
                x(iat2,iw)=x(iat2,iw)-boxx
                enddo
               enddo
       else if (x(iat,iww).lt.-boxx2)then
               iwrap=iwrap+1
               do iat2=iat,iat+natmol(i)-1
                do iw=1,nwalk
                x(iat2,iw)=x(iat2,iw)+boxx
                enddo
               enddo
       endif

       if (y(iat,iww).gt.boxy2)then
               iwrap=iwrap+1
               do iat2=iat,iat+natmol(i)-1
                do iw=1,nwalk
                y(iat2,iw)=y(iat2,iw)-boxy
                enddo
               enddo
       else if (y(iat,iww).lt.-boxy2)then
               iwrap=iwrap+1
               do iat2=iat,iat+natmol(i)-1
                do iw=1,nwalk
                y(iat2,iw)=y(iat2,iw)+boxy
                enddo
               enddo
       endif

       if (z(iat,iww).gt.boxz2)then
               iwrap=iwrap+1
               do iat2=iat,iat+natmol(i)-1
                do iw=1,nwalk
                z(iat2,iw)=z(iat2,iw)-boxz
                enddo
               enddo
       else if (z(iat,iww).lt.-boxz2)then
               iwrap=iwrap+1
               do iw=1,nwalk
                do iat2=iat,iat+natmol(i)-1
                z(iat2,iw)=z(iat2,iw)+boxz
                enddo
               enddo
       endif

       if(iwrap.gt.0)then
        write(*,*)'Wrapped molecule number:',i
        iwrap=0
       endif

       iat=iat+natmol(i)

      enddo

      end subroutine

      subroutine force_nab(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_const,       only: AUtoKCAL, ANG, AUtoKK
      use mod_general,     only: ihess, natom, ncalc, idebug, it
      use mod_estimators,  only: h
      use mod_harmon,      only: hess
      real(DP), intent(in)    :: x(:,:),y(:,:),z(:,:)
      real(DP), intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
      real(DP), intent(inout) :: eclas
      integer,  intent(in)    :: walkmax
      real(DP), parameter     :: fac=autokcal*ang
      real(DP), parameter     :: fac2=autokcal*ang*ang,fac3=autoKK*ang
      real(DP)                :: grad(size(x,1)*3),xyz(size(x,1)*3)
      real(DP)                :: dummy1(size(x,1)),dummy2(size(x,1)*3)
      character(len=size(x,1)):: dummy3
      real(DP)                :: energy,del
      integer :: iat,iat1,iat2,pom,iw
      integer :: idum1=1,idum2=1,idum3=1,idum4=1,idum5=1,idum6=1,idum7=1,idum8=1
      real(DP) :: mme, mme2
      integer  :: iter  !jak casto delame non-bonded list update?musi byt totozne s hodnotou v souboru sff_my.c

#ifdef NAB

      del = fac2 * walkmax
      iter = it   !dulezite pro update non-bond listu
      if (idebug.eq.1.and.it.ne.nsnb) iter=-1

!!    PRO PARALELNI IMPLEMENTACI JE TREBA  vztvorit wrapper pro update nblistu
!!$   if ((modulo(it,nsnb).eq.0.or.it.eq.1) )call nblistupdate()

!!$OMP PARALLEL DO PRIVATE(energy,pom,xyz,grad,idum4,dummy1,dummy2,dummy3) REDUCTION(+:eclas) !(asi neni bezpecne volat!gradhess paralelne)
      do iw = 1, walkmax
      idum4=1
      if( (modulo(it,nsnb).eq.0.or.it.eq.1).and.iw.gt.1 ) iter=-10  !update nblistu jen pro prvniho walkera

      !!$iter=-10  !u openmp,we update nblist before the loop
      
      pom=1
      do iat=1,natom
       xyz(pom)=x(iat,iw)/ang
       xyz(pom+1)=y(iat,iw)/ang
       xyz(pom+2)=z(iat,iw)/ang
       pom=pom+3
      enddo

      if(ihess.eq.1.and.modulo(it,ncalc).eq.0)then
!!$OMP CRITICAL (neco)
#ifdef NAB
       energy=mme2(xyz,grad,h,dummy1,dummy2,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,iter,dummy3)
#endif
!!$OMP END CRITICAL (neco)
       pom=1
       do iat1=1,natom*3
        do iat2=1,natom*3
        !if(iat2.le.natqm.and.iat1.le.natqm) cycle !uncomment this line if we have QM hessian for QM part
         hess(iat1,iat2,iw)=h(pom)/del
         pom=pom+1
        enddo
       enddo
!      TODO:z mme2 brat jen hessian...anebo zablokovat vypocet elstat. sil pro ipbc=1


      else    !only forces

!!$OMP CRITICAL (neco)
       energy = mme( xyz, grad, iter ); ! set iter=-1 for detailed energy info
!!$OMP END CRITICAL (neco)

      endif

      
      pom=1
      do iat=1,natom
       fx(iat,iw)=-grad(pom)/fac
       fy(iat,iw)=-grad(pom+1)/fac
       fz(iat,iw)=-grad(pom+2)/fac
       pom=pom+3
      enddo

!convert energies to au units
      energy = energy / autokcal

      eclas = eclas + energy
      enddo
!!$OMP END PARALLEL DO

      eclas = eclas / walkmax
#endif
   end subroutine force_nab

end module



