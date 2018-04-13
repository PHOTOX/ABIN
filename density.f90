
module mod_density
   use mod_const, only: DP, ANG
   use mod_array_size, only: ndistmax, nbinmax
   implicit none
   private
   public  :: ndist, nbin, nbin_ang, nang, ndih, xmin, xmax, shiftdih
   public  :: dist1, dist2, ang1, ang2, ang3, dih1, dih2, dih3, dih4
   public  :: dist_init, density, density_ang, density_dih, disterror

   integer :: ndist=0,nbin=1000, nbin_ang=180, disterror=0
   integer :: dist1(ndistmax),dist2(ndistmax)
   integer :: nang=0,ang1(ndistmax),ang2(ndistmax),ang3(ndistmax)
   integer :: ndih=0,dih1(ndistmax),dih2(ndistmax),dih3(ndistmax),dih4(ndistmax)
   real(DP)  :: dist(nbinmax,ndistmax),dist_ang(nbinmax,ndistmax),dist_dih(nbinmax,ndistmax) 
   real(DP)  :: xmin=0.5d0,xmax=5.0d0
   real(DP)  :: shiftdih=360.0d0 ! 0 for (-180,180), 360 for (0,360)
   save
   contains

   subroutine dist_init()
   dist=0.0d0
   dist_ang=0.0d0
   dist_dih=0.0d0
   end subroutine dist_init

   subroutine density(x,y,z)
   use mod_general,  only: it, nwalk, nwrite, pot,natom
   use mod_system,   only: dime
   use mod_utils,    only: get_distance, abinerror
   real(DP),intent(in)  :: x(:,:),y(:,:),z(:,:)
   real(DP)  :: r,anorm,dx,dbin
   integer :: idist,iw,ipom,ian,unit1

   unit1 = 600

   
   do idist=1,ndist
      do iw=1,nwalk
   
         r = get_distance(x, y, z, dist1(idist), dist2(idist), iw)

         if(dime.eq.1) r = x(1, iw)

         dbin = (xmax - xmin) / nbin

         ipom=ceiling(( (r/ang)-xmin )/dbin)
         if(ipom.gt.nbin.or.ipom.le.0)then
            write(*,*)'WARNING: Problems with distribution function.'
            write(*,*)'This may mean that your system is falling apart.'
            write(*,*)'Or maybe you should set xmin and xmax differently.'
            if(disterror.eq.1) call abinerror('density')
         else
               
            dist(ipom,idist) = dist(ipom,idist) + 1.0d0

         endif

      enddo
   enddo
   
   if(modulo(it,nwrite).eq.0)then
      open(unit1,file='dist.dat')
      do idist=1,ndist

         anorm=0.0d0
         dx=(xmax-xmin)/nbin
   
         do ian=1,nbin
            anorm=anorm+dist(ian,idist)
         enddo
     
         do ian=1,nbin
            write(unit1,*)ian*dx+xmin+dx/2,dist(ian,idist)/(anorm*dx)
         enddo
         write(unit1,*)
   
      enddo
   
      close(unit1)
   endif
   
   end subroutine density
                                     
   
   subroutine density_ang(x,y,z)
   use mod_general, only:it, nwalk, nwrite
   use mod_utils, only: abinerror
   real(DP) x(:,:),y(:,:),z(:,:)
   real(DP) anorm,dbin,angmin,angmax,alfa
   integer :: idist,iw,ipom,ian, iunit
   
   angmin=0.0d0
   angmax=180.0d0
   dbin=(angmax-angmin)/nbin_ang
   
   do idist=1,nang
      do iw=1,nwalk

      alfa=get_angle(x,y,z,ang1(idist),ang2(idist),ang3(idist),iw)
      ipom=ceiling( (alfa-angmin)/dbin )
   
      if(ipom.gt.nbin_ang.or.ipom.le.0)then
         write(*,*)'problems with angle distribution function'
         write(*,*)'For angle between atoms:',ang1(idist),ang2(idist),ang3(idist)
         write(*,*)'Value of ipom=',ipom
         call abinerror('density_ang')
      endif
   
      dist_ang(ipom,idist)=dist_ang(ipom,idist)+1.0d0
      enddo
   enddo
  
   iunit = 600 
   if(modulo(it,nwrite).eq.0)then
      open(iunit,file='angles.dat')
      do idist=1,nang
         anorm=0.0d0

         do ian=1,nbin_ang
            anorm=anorm+dist_ang(ian,idist)
         enddo
    
         do ian=1,nbin_ang
            write(iunit,*)ian*dbin+angmin,dist_ang(ian,idist)/(anorm*dbin)
         enddo
   
         write(iunit,*)
      enddo
      close(iunit)
   endif
   
   end subroutine density_ang
   



   subroutine density_dih(x,y,z)
   use mod_general, only:it, nwalk, nwrite
   use mod_utils, only: abinerror
   real(DP) x(:,:),y(:,:),z(:,:)
   real(DP) anorm,dbin,dihmin,dihmax,delta
   integer :: idist,iw,ipom,ian, iunit

   iunit = 600

   if (shiftdih.gt.0)then
     dihmin=0.0d0
     dihmax=360.0
   else
     dihmin=-180.0d0
     dihmax=180.0d0
   endif

   dbin=(dihmax-dihmin)/nbin_ang
   
   do idist=1,ndih
      do iw=1,nwalk
      delta=get_dihedral(x,y,z,dih1(idist),dih2(idist),dih3(idist),dih4(idist),iw)
      ipom=ceiling( (delta-dihmin)/dbin )
   
      if(ipom.gt.nbin_ang.or.ipom.le.0)then
         write(*,*)'problems with dihedral angle distribution function'
         write(*,*)'For dihedral between atoms:',dih1(idist),dih2(idist),dih3(idist),dih4(idist)
         write(*,*)'Value of ipom=',ipom
         call abinerror('density_ang')
      endif
      dist_dih(ipom,idist)=dist_dih(ipom,idist)+1.0d0
      enddo
   enddo
   
   if(modulo(it,nwrite).eq.0)then
      open(iunit,file='dihedrals.dat')
      do idist=1,ndih
         anorm=0.0d0

         do ian=1,nbin_ang
            anorm=anorm+dist_dih(ian,idist)
         enddo
    
         do ian=1,nbin_ang
            write(iunit,*)ian*dbin+dihmin,dist_dih(ian,idist)/(anorm*dbin)
         enddo
   
         write(iunit,*)
      enddo
      close(iunit)
   endif
   
   end subroutine  density_dih
   
   
   
   real(DP) function get_angle(x,y,z,at1,at2,at3,iw)
   use mod_const, only: PI
   real(DP) x(:,:),y(:,:),z(:,:)
   integer :: at1,at2,at3,iw
   real(DP)  :: vec1x,vec1y,vec1z
   real(DP)  :: vec2x,vec2y,vec2z
   
   vec1x=x(at1,iw)-x(at2,iw)
   vec1y=y(at1,iw)-y(at2,iw)
   vec1z=z(at1,iw)-z(at2,iw)
   vec2x=x(at3,iw)-x(at2,iw)
   vec2y=y(at3,iw)-y(at2,iw)
   vec2z=z(at3,iw)-z(at2,iw)
   get_angle=180/pi*acos((vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/ &
   (sqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2)))
   
   return 
   end function get_angle
   

   real(DP) function get_dihedral(x,y,z,at1,at2,at3,at4,iw)
   use mod_const, only: PI
   real(DP) x(:,:),y(:,:),z(:,:)
   integer :: at1,at2,at3,at4,iw
   real(DP)  :: vec1x,vec1y,vec1z
   real(DP)  :: vec2x,vec2y,vec2z
   real(DP)  :: vec3x,vec3y,vec3z,sign
   real(DP)  :: norm1x,norm1y,norm1z,norm2x,norm2y,norm2z
   
   vec1x=x(at1,iw)-x(at2,iw)
   vec1y=y(at1,iw)-y(at2,iw)
   vec1z=z(at1,iw)-z(at2,iw)
   vec2x=x(at3,iw)-x(at2,iw)
   vec2y=y(at3,iw)-y(at2,iw)
   vec2z=z(at3,iw)-z(at2,iw)
   vec3x=x(at4,iw)-x(at3,iw)
   vec3y=y(at4,iw)-y(at3,iw)
   vec3z=z(at4,iw)-z(at3,iw)
   
   norm1x=vec1y*vec2z-vec1z*vec2y
   norm1y=vec1z*vec2x-vec1x*vec2z
   norm1z=vec1x*vec2y-vec1y*vec2x
   norm2x=vec3y*vec2z-vec3z*vec2y
   norm2y=vec3z*vec2x-vec3x*vec2z
   norm2z=vec3x*vec2y-vec3y*vec2x

   sign = norm1x*vec3x+norm1y*vec3y+norm1z*vec3z
   get_dihedral = 180/pi*acos((norm1x*norm2x+norm1y*norm2y+norm1z*norm2z)/ &
   (sqrt(norm1x**2+norm1y**2+norm1z**2)*sqrt(norm2x**2+norm2y**2+norm2z**2)))
   
   if (sign.gt.0) get_dihedral = shiftdih-get_dihedral
   
   return
   end function get_dihedral
   
end module mod_density
   
