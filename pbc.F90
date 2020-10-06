! This module is a stub for a future
! Periodic Boundary Conditions functionality
! It is currently not used, not even compiled.
module PBC
      use mod_const, only: DP
      implicit none
      real(DP) :: boxx,boxy,boxz
      integer, allocatable  :: natmol(:) ! used for wrapping, independent of nmolt
      integer  :: nmol = 1 !used for wraping, independent of nmolt
      integer  :: ipbc = 0
      save

      contains

      subroutine wrap(x, y, z)
      use mod_general, only:nwalk
      real(DP) :: x(:,:), y(:,:), z(:,:)
      integer  :: i, iat, iat2, iw, iww, iwrap

      iwrap = 0

      iat = 1
      ! wrapping based on the 1st bead
      iww = 1
      do i = 1, nmol

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

       iat = iat + natmol(i)

      enddo

      end subroutine

end module PBC
