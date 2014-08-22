!-Simple program for analysis of bonds, angles and dihedrals from xyz movies
!-Daniel Hollas     2014

module mod_analyze
implicit none
real*8, parameter  :: pi=3.14159265
real*8,allocatable :: x(:), y(:), z(:)
save
end module  


program analyze_movie
use mod_analyze
implicit none
real*8,allocatable  :: r(:), delta(:), alfa(:)
integer,allocatable :: bins_dist(:,:), bins_ang(:,:),bins_dih(:,:)
integer,allocatable :: dists(:,:),angles(:,:),dihs(:,:)
real*8  :: distmin=0.5d0,distmax=5.0d0,get_distance
real*8  :: dihmin=-180.0d0,dihmax=180.0d0,get_dihedral
real*8  :: angmin=0.0d0,angmax=180.0d0,get_angle
real*8  :: shiftdih=0.0d0
integer :: nbin_dist=1000,nbin_ang=180,nbin_dih=180
integer :: ndist=0,nang=0,ndih=0
real*8  :: anorm=0.0d0,dx
integer :: idist,ipom,ian,i
integer :: iat,natom,it,iost
character(len=100) :: chmovie
character(len=2),allocatable  :: names(:)

!- Default, can be change from command line
chmovie='movie.xyz'

!-Read optional command line options
call Get_cmdline(chmovie, nbin_dist, nbin_ang, nbin_dih, distmin, distmax, shiftdih)

!-Now read from stdin, what to analyze
write(*,*)'How many bonds?'
read(*,*,IOSTAT=iost)ndist
if (iost.ne.0) call PrintInputError()
if (ndist.gt.0)then
  allocate( dists(2,ndist) )
  allocate( r(ndist) )
  write(*,*)'Please, specify each bond by atom indices, one bond per line.'
  do i=1,ndist
    read(*,*,IOSTAT=iost)dists(1,i),dists(2,i)
    if (iost.ne.0) call PrintInputError()
  enddo
  allocate( bins_dist(nbin_dist, ndist) )
  bins_dist=0
end if

write(*,*)'How many angles?'
read(*,*,IOSTAT=iost)nang
if (iost.ne.0) call PrintInputError()

if (nang.gt.0)then
  allocate( angles(3,nang) )
  allocate( alfa(nang) )
  write(*,*)'Please, specify each angle by atom indices, one angle per line.'
  do i=1,nang
    read(*,*)angles(1,i),angles(2,i),angles(3,i)
  enddo
  allocate( bins_ang(nbin_ang, nang)    ) 
  bins_ang=0
end if

write(*,*)'How many dihedrals?'
read(*,*,IOSTAT=iost)ndih
if (iost.ne.0) call PrintInputError()

if (ndih.gt.0)then
  allocate( dihs(4,ndih) )
  allocate( delta(ndih) )
  write(*,*)'Please, specify each dihedral by atom indices, one dihedral per line.'
  do i=1,ndih
    read(*,*)dihs(1,i),dihs(2,i),dihs(3,i),dihs(4,i)
  enddo
  allocate( bins_dih(nbin_dih, ndih) )
  bins_dih=0
end if

open(100,file=chmovie, action = "read" , status = "old",iostat=iost)
if (iost.ne.0) call PrintInputError()
read(100, *, iostat=iost ) natom
if (iost.ne.0) call PrintInputError()
rewind(100)
!------END OF READING INPUT PARAMETERS------------------
  
allocate( x(natom) )
allocate( y(natom) )
allocate( z(natom) )
allocate( names(natom) )


dihmin=dihmin+shiftdih
dihmax=dihmax+shiftdih

if (ndist.gt.0 ) open(101,file="dist_time.dat")
if (nang.gt.0 )  open(102,file="ang_time.dat")
if (ndih.gt.0 )  open(103,file="dih_time.dat")
it=0
write(*,*)'Proccessing file ',chmovie
do 
it=it+1

read(100,*,iostat=iost)natom
if(iost.ne. 0) exit
read(100,*)
do iat=1,natom
 read(100,*)names(iat),x(iat),y(iat),z(iat)
enddo




dx=(distmax-distmin)/nbin_dist

do idist=1,ndist

 r(idist)=get_distance(dists(1,idist),dists(2,idist))

 ipom=ceiling( ( (r(idist)) - distmin )/dx )
 if(ipom.gt.nbin_dist.or.ipom.le.0)then
  write(*,*)'problems with distribution function'
  write(*,*)'For distance between atoms:',dists(1,idist),dists(2,idist)
  write(*,*)'Value of ipom=',ipom,'Geometry number=',it
  stop 1
 endif

 bins_dist(ipom,idist)=bins_dist(ipom,idist)+1
enddo
if (ndist.gt.0) write(101,*)it,(r(idist),idist=1,ndist)

dx=(angmax-angmin)/nbin_ang

do idist=1,nang

 alfa(idist)=get_angle(angles(1,idist),angles(2,idist),angles(3,idist))
 ipom=ceiling( (alfa(idist)-angmin)/dx )

 if(ipom.gt.nbin_ang.or.ipom.le.0)then
  write(*,*)'problems with angle distribution function'
  write(*,*)'For angle between atoms:',angles(1,idist),angles(2,idist),angles(3,idist)
  write(*,*)'Value of ipom=',ipom,'Geometry number=',it
  stop 1
 endif

 bins_ang(ipom,idist)=bins_ang(ipom,idist)+1
enddo
if ( nang.gt.0) write(102,*)it,(alfa(idist),idist=1,nang)

dx=(dihmax-dihmin)/nbin_dih

do idist=1,ndih
 delta(idist)=get_dihedral(dihs(1,idist),dihs(2,idist),dihs(3,idist),dihs(4,idist),shiftdih)
 ipom=ceiling((delta(idist)-dihmin)/dx)

 if(ipom.gt.nbin_dih.or.ipom.le.0)then
  write(*,*)'problems with dihedral angle distribution function'
  write(*,*)'For dihedral between atoms:',dihs(1,idist),dihs(2,idist),dihs(3,idist),dihs(4,idist)
  write(*,*)'Value of ipom=',ipom,'Geometry number=',it
  stop 1
 endif
 bins_dih(ipom,idist)=bins_dih(ipom,idist)+1
enddo

if( ndih.gt.0) write(103,*)it,(delta(idist),idist=1,ndih)

enddo

close(100)
if (ndist.gt.0 ) close(101)
if (nang.gt.0 )  close(102)
if (ndih.gt.0 )  close(103)

!--------PRINTING and Normalizing histograms---------------------
if (ndist.gt.0 )then
  open(128,file='dist_hist.dat')
  dx=(distmax-distmin)/nbin_dist
  do idist=1,ndist
    anorm=0.0d0

    do ian=1,nbin_dist
      anorm=anorm+bins_dist(ian,idist)
    enddo

    do ian=1,nbin_dist
      write(128,*)ian*dx+distmin+dx/2,bins_dist(ian,idist)/(anorm*dx)
    enddo
    write(128,*)

  enddo
  close(128)
endif

if (nang.gt.0 )then
  open(10,file='ang_hist.dat')
  dx=(angmax-angmin)/nbin_ang

  do idist=1,nang
    anorm=0.0d0

    do ian=1,nbin_ang
      anorm=anorm+bins_ang(ian,idist)
    enddo

    do ian=1,nbin_ang
      write(10,*)ian*dx+angmin+dx/2,bins_ang(ian,idist)/(anorm*dx)
    enddo

    write(10,*)
  enddo
  close(10)
endif


if (ndih.gt.0 )then
  open(11,file='dih_hist.dat')
  dx=(dihmax-dihmin)/nbin_dih

  do idist=1,ndih
    anorm=0.0d0

    do ian=1,nbin_dih
      anorm=anorm+bins_dih(ian,idist)
    enddo


    do ian=1,nbin_dih
      write(11,*)ian*dx+dihmin+dx/2,bins_dih(ian,idist)/(anorm*dx)
    enddo

    write(11,*)
  enddo
  close(11)

endif

if (ndist.gt.0)then
  deallocate( bins_dist )
  deallocate( dists )
  deallocate( r )
endif
if (nang.gt.0)then
  deallocate( bins_ang ) 
  deallocate( angles )
  deallocate( alfa )
endif
if (ndih.gt.0)then
  deallocate( bins_dih )
  deallocate( dihs )
  deallocate( delta )
endif
deallocate( x )
deallocate( y )
deallocate( z )

end





real*8 function get_distance( at1, at2)
use mod_analyze
implicit none
integer,intent(in) :: at1, at2
real*8  :: dx,dy,dz,temp

dx=x(at1)-x(at2)
dy=y(at1)-y(at2)
dz=z(at1)-z(at2)
temp=dx*dx+dy*dy+dz*dz
temp=dsqrt(temp)
get_distance=temp

return 
end


real*8 function get_angle(at1,at2,at3)
use mod_analyze
implicit none
integer,intent(in) :: at1,at2,at3
real*8  :: vec1x,vec1y,vec1z
real*8  :: vec2x,vec2y,vec2z

vec1x=x(at1)-x(at2)
vec1y=y(at1)-y(at2)
vec1z=z(at1)-z(at2)
vec2x=x(at3)-x(at2)
vec2y=y(at3)-y(at2)
vec2z=z(at3)-z(at2)
get_angle=180/pi*acos((vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/ &
(dsqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2)))

return 
end


real*8 function get_dihedral(at1, at2, at3, at4, shiftdih)
use mod_analyze
implicit none
integer,intent(in) :: at1,at2,at3,at4
real*8  :: shiftdih
real*8  :: vec1x,vec1y,vec1z
real*8  :: vec2x,vec2y,vec2z
real*8  :: vec3x,vec3y,vec3z,sign
real*8  :: norm1x,norm1y,norm1z,norm2x,norm2y,norm2z

vec1x=x(at1)-x(at2)
vec1y=y(at1)-y(at2)
vec1z=z(at1)-z(at2)
vec2x=x(at3)-x(at2)
vec2y=y(at3)-y(at2)
vec2z=z(at3)-z(at2)
vec3x=x(at4)-x(at3)
vec3y=y(at4)-y(at3)
vec3z=z(at4)-z(at3)

norm1x=vec1y*vec2z-vec1z*vec2y
norm1y=vec1z*vec2x-vec1x*vec2z
norm1z=vec1x*vec2y-vec1y*vec2x
norm2x=vec3y*vec2z-vec3z*vec2y
norm2y=vec3z*vec2x-vec3x*vec2z
norm2z=vec3x*vec2y-vec3y*vec2x

sign = norm1x*vec3x+norm1y*vec3y+norm1z*vec3z
get_dihedral = 180/pi*acos((norm1x*norm2x+norm1y*norm2y+norm1z*norm2z)/ &
(sqrt(norm1x**2+norm1y**2+norm1z**2)*sqrt(norm2x**2+norm2y**2+norm2z**2)))

if (sign.gt.0) get_dihedral = shiftdih - get_dihedral

return
end

subroutine Get_cmdline(chmovie, nbin_dist, nbin_ang, nbin_dih, distmin, distmax, shiftdih)
implicit none
real*8,intent(out)  :: distmin,distmax
real*8,intent(out)  :: shiftdih
integer,intent(out) :: nbin_dist,nbin_ang,nbin_dih
character(len=100),intent(out)   :: chmovie
character(len=100)   :: arg
integer             :: i

i=0
do while (i < command_argument_count())
  i=i+1
  call get_command_argument(i, arg)
  
  select case (arg)
  case ('-h', '--help')
    call PrintHelp()
    stop
  case ('-m')
    i=i+1
    call get_command_argument(i, arg)
    write(*,*)arg
    !-format specifier is needed here in case of slashes
    read(arg,'(A)')chmovie  
  case ('-distbin')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nbin_dist
  case ('-angbin')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nbin_ang
  case ('-dihbin')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nbin_dih
  case ('-dmin')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)distmin
  case ('-dmax')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)distmax
  case ('-shiftdih')
    shiftdih=360.0d0
  case default
    write(*,*)'Invalid command line argument!'
    call PrintHelp()
    stop 1
  end select

end do

end subroutine Get_cmdline

subroutine PrintHelp()
implicit none
    print '(a)', 'Program for bond, angle and dihedral analysis of xyz movies.'
    print '(a)', 'USAGE: ./analyze_movie [OPTIONS]'
    print '(a)', ''
    print '(a)', 'The bond and angle specifications are read from standard input.'
    print '(a)', ''
    print '(a)', 'Without cmdline options, the program uses sensible default values.'
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h, --help        print usage information and exit'
    print '(a)', '  -m <movie_file>   default: movie.xyz'
    print '(a)', '  -distbin <nbins>  number of bins for distance analysis, default: 1000'
    print '(a)', '  -angbin  <nbins>  number of bins for distance analysis, default: 180'
    print '(a)', '  -dihbin  <nbins>  number of bins for dihedral analysis, default: 180'
    print '(a)', '  -dmax    <number> maximum distance in Angstroms, default: 5.0'
    print '(a)', '  -dmin    <number> minimum distance in Angstroms, default: 0.5'
    print '(a)', '  -shiftdih         Calculate dihedral from 0-360 degrees rather then from -180 to 180'
end subroutine PrintHelp

subroutine PrintInputError()
  write(*,*)'Error during reading input. Exiting...'
  stop 1
end subroutine PrintInputError


