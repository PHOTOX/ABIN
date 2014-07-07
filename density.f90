subroutine density(x,y,z)
use mod_array_size
use mod_general
use mod_system
implicit none
real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8  :: r,anorm,dx,dbin
integer :: idist,iw,ipom,ian


do idist=1,ndist
   do iw=1,nwalk

   r=(x(dist1(idist),iw)-x(dist2(idist),iw))**2
   r=r+(y(dist1(idist),iw)-y(dist2(idist),iw))**2
   r=r+(z(dist1(idist),iw)-z(dist2(idist),iw))**2
   r=sqrt(r)

   if(dime.eq.1.and.pot.eq.'2dho') r=x(1,iw)

   dbin=(xmax-xmin)/nbin

   ipom=ceiling(( (r/ang)-xmin )/dbin)
   if(ipom.gt.nbin.or.ipom.le.0)then
      write(*,*)'problems with distribution function'
      stop 1
   endif

   dist(ipom,idist)=dist(ipom,idist)+1.0d0
   enddo
enddo

if(modulo(it,nwrite).eq.0)then
  open(128,file='dist.dat')
  do idist=1,ndist
  anorm=0.0d0
  dx=(xmax-xmin)/nbin

  do ian=1,nbin
   anorm=anorm+dist(ian,idist)
  enddo
  
  do ian=1,nbin
   write(128,*)ian*dx+xmin+dx/2,dist(ian,idist)/(anorm*dx)
  enddo
  write(128,*)

  enddo

  close(128)
endif

end subroutine density
                                  

subroutine density_ang(x,y,z)
use mod_array_size
use mod_general
use mod_system
implicit none
real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8 anorm,dbin,angmin,angmax,alfa,get_angle
integer :: idist,iw,ipom,ian

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
      stop 1
   endif

   dist_ang(ipom,idist)=dist_ang(ipom,idist)+1.0d0
   enddo
enddo

if(modulo(it,nwrite).eq.0)then
   open(10,file='angles.dat')
   do idist=1,nang
      anorm=0.0d0

      do ian=1,nbin_ang
         anorm=anorm+dist_ang(ian,idist)
      enddo
 
      do ian=1,nbin_ang
         write(10,*)ian*dbin+angmin,dist_ang(ian,idist)/(anorm*dbin)
      enddo

      write(10,*)
   enddo
   close(10)
endif

end




subroutine density_dih(x,y,z)
use mod_array_size
use mod_general
use mod_system
implicit none
real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
real*8 anorm,dbin,dihmin,dihmax,delta,get_dihedral
integer :: idist,iw,ipom,ian

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
      stop 1
   endif
   dist_dih(ipom,idist)=dist_dih(ipom,idist)+1.0d0
   enddo
enddo

if(modulo(it,nwrite).eq.0)then
   open(10,file='dihedrals.dat')
   do idist=1,ndih
      anorm=0.0d0

      do ian=1,nbin_ang
         anorm=anorm+dist_dih(ian,idist)
      enddo
 
      do ian=1,nbin_ang
         write(10,*)ian*dbin+dihmin,dist_dih(ian,idist)/(anorm*dbin)
      enddo

      write(10,*)
   enddo
   close(10)
endif

end



real*8 function get_angle(x,y,z,at1,at2,at3,iw)
use mod_array_size
implicit none
real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
integer :: at1,at2,at3,iw
real*8  :: vec1x,vec1y,vec1z
real*8  :: vec2x,vec2y,vec2z

vec1x=x(at1,iw)-x(at2,iw)
vec1y=y(at1,iw)-y(at2,iw)
vec1z=z(at1,iw)-z(at2,iw)
vec2x=x(at3,iw)-x(at2,iw)
vec2y=y(at3,iw)-y(at2,iw)
vec2z=z(at3,iw)-z(at2,iw)
get_angle=180/pi*acos((vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/ &
(sqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2)))

return 
end


real*8 function get_dihedral(x,y,z,at1,at2,at3,at4,iw)
use mod_array_size
use mod_system, only:shiftdih
implicit none
real*8 x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
integer :: at1,at2,at3,at4,iw
real*8  :: vec1x,vec1y,vec1z
real*8  :: vec2x,vec2y,vec2z
real*8  :: vec3x,vec3y,vec3z,sign
real*8  :: norm1x,norm1y,norm1z,norm2x,norm2y,norm2z

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
end

