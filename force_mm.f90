
      subroutine force_LJCoul(x,y,z,fx,fy,fz,eclas)
      use mod_array_size
      use mod_general
      use mod_system, ONLY:inames 
      use mod_qmmm
      implicit none
      real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
      real*8 fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
      integer :: iat,iw,iat1,iat2,i1,i2
      real*8  :: eclas,r,fxLJ,fyLJ,fzLJ,fxC,fyC,fzC,kLJ,KC
      real*8  :: ri,ri3,dx,dy,dz

      do iw=1,nwalk
       do iat1=1,natom
        do iat2=iat1+1,natom
        if (iat2.le.natqm) cycle
        dx=x(iat1,iw)-x(iat2,iw)
        dy=y(iat1,iw)-y(iat2,iw)
        dz=z(iat1,iw)-z(iat2,iw)
        r=dx**2+dy**2+dz**2
        ri=1/r
        ri3=ri*ri*ri
        i1=inames(iat1)
        i2=inames(iat2)
        kLJ=ri3*(ri3*AIJ(i1,i2)-BIJ(i1,i2))*ri
        kC=q(i1)*q(i2)*sqrt(ri3)
        fx(iat1,iw)=fx(iat1,iw)+(kLJ+kC)*dx
        fx(iat2,iw)=fx(iat2,iw)-(kLJ+kC)*dx
        fy(iat1,iw)=fy(iat1,iw)+(kLJ+kC)*dy
        fy(iat2,iw)=fy(iat2,iw)-(kLJ+kC)*dy
        fz(iat1,iw)=fz(iat1,iw)+(kLJ+kC)*dz
        fz(iat2,iw)=fz(iat2,iw)-(kLJ+kC)*dz
        eclas=eclas+ri3*(ri3*AIJ(i1,i2)/12-BIJ(i1,i2)/6)/nwalk
        eclas=eclas+q(i1)*q(i2)/sqrt(r)/nwalk
        enddo
       enddo
      enddo


        end
