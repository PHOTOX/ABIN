      module mod_nab
      use mod_const, only: DP
      implicit none
      real(DP) :: boxx,boxy,boxz
      real(DP) :: boxx2,boxy2,boxz2
      real(DP) :: alpha_pme=-1,kappa_pme=-1,cutoff=100.0d0
      real(DP) :: epsinf=3d33
      real(DP), allocatable :: charges(:)
      integer, allocatable  :: natmol(:) !used for wraping, independent of nmolt
      integer  :: nmol=1 !used for wraping, independent of nmolt
      integer  :: ipbc=0,nsnb=1
      integer  :: ips=0 ! 1-both LJ and coul, 2-coul 3- LJ
      save
      contains

      subroutine wrap(x,y,z)
      use mod_general,only:nwalk
      real(DP) :: x(:,:),y(:,:),z(:,:)
      integer  :: i,iat,iat2,iw,iww,iwrap

      iwrap=0

      iat=1
      !delame wrapping na zaklade 1.walkera
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

      !if (iwrap.gt.0) write(*,*)'Number of wraps in this step:',iwrap

      end subroutine

      subroutine force_nab(x,y,z,fx,fy,fz,eclas)
      use mod_const,       only: AUtoKCAL, ANG, AUtoKK
      use mod_general,     only: ihess, natom, nwalk, ncalc, idebug, iqmmm, it
      use mod_estimators,  only: h
      use mod_harmon,      only: hess
      use mod_qmmm,        only: natqm
      real(DP), intent(in)    :: x(:,:),y(:,:),z(:,:)
      real(DP), intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
      real(DP), intent(inout) :: eclas
      real(DP), parameter     :: fac=autokcal*ang
      real(DP), parameter     :: fac2=autokcal*ang*ang,fac3=autoKK*ang
      real(DP)                :: grad(size(x,1)*3),xyz(size(x,1)*3)
      real(DP)                :: dummy1(size(x,1)),dummy2(size(x,1)*3)
      character(len=size(x,1)):: dummy3
      real(DP)                :: energy,del
      integer :: iat,iat1,iat2,pom,iw
      integer :: idum1=1,idum2=1,idum3=1,idum4=1,idum5=1,idum6=1,idum7=1,idum8=1
      !pro ewalda
      real(DP) ::  mme,mme2,en_ewald,mme_qmmm
      integer  :: iter  !jak casto delame non-bonded list update?musi byt totozne s hodnotou v souboru sff_my.c

      en_ewald=0.0d0


      del=fac2*nwalk
      iter=it   !dulezite pro update non-bond listu
      if (idebug.eq.1.and.it.ne.nsnb) iter=-1

!!    PRO PARALELNI IMPLEMENTACI JE TREBA  vztvorit wrapper pro update nblistu
!!$   if ((modulo(it,nsnb).eq.0.or.it.eq.1) )call nblistupdate()

!!$OMP PARALLEL DO PRIVATE(energy,pom,xyz,grad,idum4,dummy1,dummy2,dummy3) REDUCTION(+:eclas) !(asi neni bezpecne volat!gradhess paralelne)
      do iw=1,nwalk
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
       energy=mme2(xyz,grad,h,dummy1,dummy2,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,iter,dummy3)
!!$OMP END CRITICAL (neco)
       pom=1
       do iat1=1,natom*3
        do iat2=1,natom*3
        !if(iat2.le.natqm.and.iat1.le.natqm) cycle !uncomment this line if we have QM hessian for QM part
         hess(iat1,iat2,iw)=h(pom)/del
         pom=pom+1
        enddo
       enddo
!       TODO:z mme2 brat jen hessian...anebo zablokovat vypocet elstat. sil pro ipbc=1


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

! calculation of reciprocal part of Ewald sum
! grad se nemusi nulovat, prepisuje se uvnitr ewalda
      if(ipbc.eq.1)then
       call ewald(xyz,grad,charges,en_ewald,boxx,boxy,boxz,cutoff,alpha_pme,kappa_pme,epsinf,natom,ipbc)
       pom=1
       do iat=1,natom
        fx(iat,iw)=fx(iat,iw)+grad(pom)/fac3  
        fy(iat,iw)=fy(iat,iw)+grad(pom+1)/fac3
        fz(iat,iw)=fz(iat,iw)+grad(pom+2)/fac3
        pom=pom+3
       enddo
      endif

!     CORRECTION FOR QMMM, COMPUTING ONLY QM PART and SUBSTRACTING FROM THE WHOLE SYSTEM
!     currently does not work!!!
!     PBC NOT SUPPORTED!
      if(iqmmm.eq.2)then

        if(idebug.eq.1)then
                iter=-1
        else  
                iter=-5  !we dont update nblist for QM part!!!!
        endif

        energy=energy-mme_qmmm( xyz, grad, iter) !grad se nemusi nulovat

!       if(ihess.eq.1)then
!        energy=energy-mme2(xyz,grad,h,dummy1,dummy2,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,iter,dummy3)
!        pom=1
!        do iat1=1,natqm*3
!         do iat2=1,natqm*3
!           hess(iat1,iat2,iw)=hess(iat1,iat2,iw)-h(pom)/del
!           pom=pom+1
!         enddo
!        enddo
!      endif

       pom=1
       do iat=1,natqm
        fx(iat,iw)=fx(iat,iw)+grad(pom)/fac
        fy(iat,iw)=fy(iat,iw)+grad(pom+1)/fac
        fz(iat,iw)=fz(iat,iw)+grad(pom+2)/fac
        pom=pom+3
       enddo

      endif
!----END-OF-QMMM-----------------------       

!convert energies to au units
      en_ewald=en_ewald/autoKK ! K to a.u.
      energy=energy/autokcal

      eclas=eclas+energy+en_ewald
      enddo
!!$OMP END PARALLEL DO

      
!       open(100,file="hessian_nab.dat")
!       iw=1
!       do iat1=1,natom*3
!        do iat2=1,natom*3
!         write(100,*)hess(iat1,iat2,iw)
!    pom=pom+1
!   enddo
!  enddo
!  close(100)

      eclas=eclas/nwalk
   end subroutine force_nab

end module



