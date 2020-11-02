! Simple program for analysis of bonds, angles and dihedrals
! from XYZ trajectories.
! Writes out time evolution as well as histograms.

! Author: Daniel Hollas

module mod_analyze
implicit none
real*8, parameter  :: pi = 3.14159265d0
real*8, allocatable :: x(:), y(:), z(:)
save

CONTAINS

! The following function also handles the PBC conditions
real*8 function get_distance( at1, at2, boxx, boxy, boxz)
   implicit none
   integer,intent(in) :: at1, at2
   real*8, intent(in) :: boxx, boxy, boxz
   real*8  :: dx, dy, dz, temp

   if (at1.eq.at2)then
      dx = x(at1)
      dy = y(at1)
      dz = z(at1)
   else
      dx = x(at1) - x(at2)
      dy = y(at1) - y(at2)
      dz = z(at1) - z(at2)
   end if

   if (boxx.gt.0)then
      dx = dx - boxx * nint(dx/boxx)
      dy = dy - boxy * nint(dy/boxy)
      dz = dz - boxz * nint(dz/boxz)
   end if
   
   temp = dx*dx + dy*dy + dz*dz
   temp = dsqrt(temp)
   get_distance = temp
   return 
end function

real*8 function get_angle(at1, at2, at3, boxx, boxy, boxz)
   implicit none
   integer,intent(in) :: at1, at2, at3
   real*8, optional, intent(in) :: boxx, boxy, boxz
   real*8  :: vec1x, vec1y, vec1z
   real*8  :: vec2x, vec2y, vec2z
   
   vec1x = x(at1) - x(at2)
   vec1y = y(at1) - y(at2)
   vec1z = z(at1) - z(at2)
   vec2x = x(at3) - x(at2)
   vec2y = y(at3) - y(at2)
   vec2z = z(at3) - z(at2)
   
   if(boxx.gt.0)then
      vec1x = vec1x - boxx * nint(vec1x/boxx)
      vec1y = vec1y - boxy * nint(vec1y/boxy)
      vec1z = vec1z - boxz * nint(vec1z/boxz)
      vec2x = vec2x - boxx * nint(vec2x/boxx)
      vec2y = vec2y - boxy * nint(vec2y/boxy)
      vec2z = vec2z - boxz * nint(vec2z/boxz)
   end if

   get_angle=180/pi*acos((vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/ &
   (dsqrt(vec1x**2+vec1y**2+vec1z**2)*sqrt(vec2x**2+vec2y**2+vec2z**2)))

return 
end function get_angle

real*8 function get_dihedral(at1, at2, at3, at4, shiftdih, boxx, boxy, boxz)
   implicit none
   integer,intent(in) :: at1,at2,at3,at4
   real*8, optional, intent(in) :: boxx, boxy, boxz
   real*8  :: shiftdih
   real*8  :: vec1x,vec1y,vec1z
   real*8  :: vec2x,vec2y,vec2z
   real*8  :: vec3x,vec3y,vec3z,sign
   real*8  :: norm1x,norm1y,norm1z,norm2x,norm2y,norm2z
   
   vec1x = x(at1) - x(at2)
   vec1y = y(at1) - y(at2)
   vec1z = z(at1) - z(at2)
   vec2x = x(at3) - x(at2)
   vec2y = y(at3) - y(at2)
   vec2z = z(at3) - z(at2)
   vec3x = x(at4) - x(at3)
   vec3y = y(at4) - y(at3)
   vec3z = z(at4) - z(at3)

   if(boxx.gt.0)then
      vec1x = vec1x - boxx * nint(vec1x/boxx)
      vec1y = vec1y - boxy * nint(vec1y/boxy)
      vec1z = vec1z - boxz * nint(vec1z/boxz)
      vec2x = vec2x - boxx * nint(vec2x/boxx)
      vec2y = vec2y - boxy * nint(vec2y/boxy)
      vec2z = vec2z - boxz * nint(vec2z/boxz)
      vec3x = vec3x - boxx * nint(vec3x/boxx)
      vec3y = vec3y - boxy * nint(vec3y/boxy)
      vec3z = vec3z - boxz * nint(vec3z/boxz)
   end if
   
   norm1x = vec1y*vec2z - vec1z*vec2y
   norm1y = vec1z*vec2x - vec1x*vec2z
   norm1z = vec1x*vec2y - vec1y*vec2x
   norm2x = vec3y*vec2z - vec3z*vec2y
   norm2y = vec3z*vec2x - vec3x*vec2z
   norm2z = vec3x*vec2y - vec3y*vec2x
   
   sign = norm1x*vec3x+norm1y*vec3y+norm1z*vec3z
   get_dihedral = 180/pi*acos((norm1x*norm2x+norm1y*norm2y+norm1z*norm2z)/ &
   (sqrt(norm1x**2+norm1y**2+norm1z**2)*sqrt(norm2x**2+norm2y**2+norm2z**2)))
   
   if (sign.gt.0) get_dihedral = shiftdih - get_dihedral
   
   return
end function

end module  


program analyze_movie
use mod_analyze
implicit none
real*8,allocatable  :: r(:), delta(:), alfa(:)
real*8,allocatable  :: bins_dist(:,:), bins_ang(:,:),bins_dih(:,:)
integer,allocatable :: dists(:,:), angles(:,:), dihs(:,:)
integer,allocatable :: ignore(:)
real*8  :: distmin=0.5d0, distmax=5.0d0
real*8  :: dihmin=-180.0d0, dihmax=180.0d0
real*8  :: angmin=0.0d0, angmax=180.0d0
real*8  :: shiftdih=0.0d0
integer :: nbin_dist=1000, nbin_ang=180, nbin_dih=180
integer :: ndist=0, nang=0, ndih=0
real*8  :: anorm=0.0d0, dx
integer :: idist, ipom, ian, i
integer :: iat, iat2, natom, it, iost
integer :: nignore, rdfpairs, rdfparts1, rdfparts2
integer :: imod=0  !mode of action, =1 for RDF
character(len=100) :: chmovie
character(len=2),allocatable  :: names(:)
character(len=2)     :: rdfname1, rdfname2
logical              :: lecho=.true., lignore=.false.
real*8  :: boxx=-1, boxy=-1, boxz=-1, rho, dvol

!- Default, can be change from command line
chmovie = 'movie.xyz'

!-Read optional command line options
call Get_cmdline(chmovie, nbin_dist, nbin_ang, nbin_dih, distmin, distmax, shiftdih, lecho, &
               imod, boxx, boxy, boxz)

! Now read from stdin, what to analyze
select case (imod)
case (0)

if (lecho) write(*,*)'How many bonds?'
read(*, *, IOSTAT=iost) ndist
if (iost.ne.0) call PrintInputError()
if (ndist.gt.0)then
  allocate( dists(2,ndist) )
  allocate( r(ndist) )
  if (lecho) write(*,*)'Please, specify each bond by atom indices, one bond per line.'
  do i=1,ndist
    read(*,*, IOSTAT=iost)dists(1,i),dists(2,i)
    if (iost.ne.0) call PrintInputError()
  enddo
  allocate( bins_dist(nbin_dist, ndist) )
  bins_dist = 0.0d0
end if

if (lecho) write(*,*)'How many angles?'
read(*, *, IOSTAT=iost)nang
if (iost.ne.0) call PrintInputError()

if (nang.gt.0)then
  allocate( angles(3,nang) )
  allocate( alfa(nang) )
  if (lecho) write(*,*)'Please, specify each angle by atom indices, one angle per line.'
  do i=1,nang
    read(*, *)angles(1,i), angles(2,i), angles(3,i)
  enddo
  allocate( bins_ang(nbin_ang, nang)    ) 
  bins_ang = 0.0d0
end if

if (lecho) write(*,*)'How many dihedrals?'
read(*, '(I3)', IOSTAT = iost) ndih
if (iost.ne.0) call PrintInputError()

if (ndih.gt.0)then
  allocate( dihs(4,ndih) )
  allocate( delta(ndih) )
  if(lecho)  write(*,*)'Please, specify each dihedral by atom indices, one dihedral per line.'
  do i=1,ndih
    read(*,*)dihs(1,i),dihs(2,i),dihs(3,i),dihs(4,i)
  enddo
  allocate( bins_dih(nbin_dih, ndih) )
  bins_dih=0.0d0
end if

!! Input for RDF.
case (1)

   ndist=1
   allocate( dists(2,ndist) )
   allocate( r(ndist) )
   allocate( bins_dist(nbin_dist, ndist) )
   bins_dist=0.0d0
 
   if(lecho)  write(*,*)'Please, specify first atom element.'
   read(*,'(A2)', IOSTAT=iost)rdfname1
   if (iost.ne.0) call PrintInputError()
 
   if(lecho)  write(*,*)'Please, specify second atom element.'
   read(*,'(A2)')rdfname2
   if (iost.ne.0) call PrintInputError()
 
   if(lecho)  write(*,*)'How many atoms should I ignore?'
   read(*,*)nignore
   if (iost.ne.0) call PrintInputError()
 
   if(nignore.gt.0)then
      allocate( ignore(nignore) )
      if(lecho)  write(*,*)'Please, specify indices of atoms to ignore, one per line.'
      do i=1,nignore
         read(*,*)ignore(i)
         if (iost.ne.0) call PrintInputError()
      end do
   end if

case DEFAULT
   write(*,*)'FATAL ERROR: Invalid mode of action'
   call PrintHelp()
   stop 1

end select


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


dihmin = dihmin + shiftdih/2
dihmax = dihmax + shiftdih/2

if (ndist.gt.0.and.imod.ne.1 ) open(101,file="dist_time.dat")
if (nang.gt.0 )  open(102,file="ang_time.dat")
if (ndih.gt.0 )  open(103,file="dih_time.dat")


if (imod.eq.1)then
   ! determine the smallest edge of the box
   if (boxx.gt.0)then
      distmin=0.0d0
      if (boxx.lt.boxy)then
         distmax=boxx/2.0d0
      else
         distmax=boxy/2.0d0
      end if
      if (boxz/2.0d0.lt.distmax)then
         distmax=boxz/2.0d0
      end if
   end if
end if


if (lecho) write(*,*)'Proccessing file ',chmovie
it=0
do 

read(100,*,iostat=iost)natom
if(iost.ne. 0) exit
read(100,*)
do iat=1,natom
 read(100,*)names(iat),x(iat),y(iat),z(iat)
enddo

it=it+1




! RDF part
if (imod.eq.1)then
   dx = (distmax-distmin)/nbin_dist
   rdfpairs=0
   rdfparts1=0
   rdfparts2=0
   do iat=1,natom
      lignore=.false.
      ! first determine, whether to skip the current atom
      do i=1,nignore
         if (ignore(i).eq.iat) lignore=.true.
      end do
      if((names(iat).ne.rdfname1.and.names(iat).ne.rdfname2).or.lignore) cycle

      if(names(iat).eq.rdfname1)then
         rdfparts1=rdfparts1+1
      end if
      if (names(iat).eq.rdfname2)then
         rdfparts2=rdfparts2+1
      end if

      do iat2=iat+1,natom
         lignore=.false.
         ! first determine, whether to skip the current atom
         do i=1,nignore
            if (ignore(i).eq.iat2) lignore=.true.
         end do
         if((names(iat2).ne.rdfname1.and.names(iat2).ne.rdfname2).or.lignore) cycle
         !folowing line need to be if doing different elements
         if(names(iat).eq.names(iat2).and.rdfname1.ne.rdfname2)     cycle
         rdfpairs=rdfpairs+1
       
         r(1)=get_distance(iat, iat2, boxx, boxy, boxz)
       
         ! only half a box
         ! this line is here for non-cubic boxes
         if (r(1).gt.distmax.and.boxx.gt.0) cycle
            
         ipom=ceiling( ( (r(1)) - distmin )/dx )
         if(ipom.gt.nbin_dist.or.ipom.le.0)then
            write(*,*)'problems with distribution function'
            write(*,'(A,I3,I3)')'For distance between atoms:',iat, iat2
            write(*,*)'Value of ipom=',ipom,'Geometry number=',it
            stop 1
         endif
       
         bins_dist(ipom,1)=bins_dist(ipom,1)+2.0d0

      end do
   end do

end if

if (imod.eq.0)then

   dx = (distmax-distmin)/nbin_dist

   do idist = 1, ndist
   
      r(idist) = get_distance(dists(1,idist), dists(2,idist), boxx, boxy, boxz)
   
      ! TODO: Put histograming at the end
      ! so we can calculate optimal bin width based on Scott's rule
      ! https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
      ! and can calculate max and min data
      ipom = ceiling( ( (r(idist)) - distmin )/dx )
      if(ipom.gt.nbin_dist.or.ipom.le.0)then
         write(*,*)'Problems with bond distribution function'
         write(*,*)'For distance between atoms:',dists(1,idist), dists(2,idist)
         write(*,*)'Distance=',r(idist),'Geometry number=',it
         write(*,*)'Boxx, boxy, boxz:',boxx, boxy, boxz
         stop 1
      endif
 
      bins_dist(ipom,idist)=bins_dist(ipom,idist)+1.0d0
   end do

   if (ndist.gt.0) write(101,*)it,(r(idist),idist=1,ndist)
   
end if ! imod end if, the rest will be skipped anyway


do idist = 1, nang
   dx = (angmax-angmin) / nbin_ang

   alfa(idist) = get_angle(angles(1,idist), angles(2,idist), angles(3,idist), &
                        boxx, boxy, boxz)
   ipom = ceiling( (alfa(idist)-angmin)/dx )

   if(ipom.gt.nbin_ang.or.ipom.le.0)then
      write(*,*)'problems with angle distribution function'
      write(*,*)'For angle between atoms:',angles(1,idist),angles(2,idist),angles(3,idist)
      write(*,*)'Value of ipom=',ipom,'Geometry number=',it
      stop 1
   endif

   bins_ang(ipom,idist) = bins_ang(ipom,idist) + 1.0d0

enddo

if ( nang.gt.0) write(102,*)it,(alfa(idist),idist=1,nang)

do idist=1,ndih

   dx=(dihmax-dihmin)/nbin_dih
   delta(idist) = get_dihedral(dihs(1,idist), dihs(2,idist), dihs(3,idist), dihs(4,idist), &
                              shiftdih, boxx, boxy, boxz)
   ipom=ceiling((delta(idist)-dihmin)/dx)

   if(ipom.gt.nbin_dih.or.ipom.le.0)then
      write(*,*)'problems with dihedral angle distribution function'
      write(*,*)'For dihedral between atoms:',dihs(1,idist),dihs(2,idist),dihs(3,idist),dihs(4,idist)
      write(*,*)'Value of ipom=',ipom,'Geometry number=',it
      stop 1
   endif
   bins_dih(ipom,idist)=bins_dih(ipom,idist)+1.0d0
enddo

if( ndih.gt.0) write(103,*)it,(delta(idist),idist=1,ndih)

enddo

! Closing the trajectory file
close(100)
if (ndist.gt.0 ) close(101)
if (nang.gt.0 )  close(102)
if (ndih.gt.0 )  close(103)

write(*,*)'Succesfully processed ',it,' geometries.'

!--------PRINTING and Normalizing histograms---------------------
if (ndist.gt.0 )then
   open(128,file='dist_hist.dat')
   dx = (distmax-distmin)/nbin_dist
   do idist = 1,ndist
      anorm = 0.0d0

      do ian=1, nbin_dist
         anorm = anorm + bins_dist(ian,idist)
      enddo

      if (imod.eq.1)then
         write(*,*)'Number of particles type 1 in RDF:',rdfparts1
         write(*,*)'Number of particles type 2 in RDF:',rdfparts2
      end if

      if (imod.eq.0)then
         anorm = anorm * dx
      end if

      do ian = 1, nbin_dist
         r(1) = distmin+dx*(ian+0.5)
         if (imod.eq.1)then
            dvol = ((ian+1)**3-ian**3)*dx**3
            rho = 1.0d0
            if(boxx.gt.0) rho = rdfparts2/boxx/boxy/boxz
            anorm = 4.0d0/3.0d0*pi*dvol*rho
            anorm = anorm*it*rdfparts1
            if(rdfname1.ne.rdfname2)then
               anorm=anorm*2.0d0
            end if
         end if

         write(128,*) r(1), bins_dist(ian, idist) / anorm

      enddo

    write(128,*)

  enddo
  close(128)
endif

if (nang.gt.0)then
   open(10,file='ang_hist.dat')
   dx=(angmax-angmin)/nbin_ang

   do idist=1,nang
      anorm=0.0d0

      do ian=1,nbin_ang
         anorm = anorm + bins_ang(ian,idist)
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

end program

   
subroutine Get_cmdline(chmovie, nbin_dist, nbin_ang, nbin_dih, distmin, distmax, shiftdih, lecho, &
            imod, boxx, boxy, boxz)
implicit none
real*8,intent(inout) :: distmin, distmax, boxx, boxy, boxz
real*8,intent(inout) :: shiftdih
integer,intent(inout) :: nbin_dist,nbin_ang,nbin_dih
integer, intent(inout) :: imod
character(len=100),intent(out) :: chmovie
logical, intent(inout) :: lecho
character(len=100) :: arg
integer :: i, ngeom

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
  case ('-optimal_nbin')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)ngeom
    call print_bin_number_estimators(ngeom)
    stop 0
  case ('-box')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)boxx
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)boxy
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)boxz
  case ('-shiftdih')
    shiftdih=360.0d0
  case ('-noecho')
    lecho=.false.
  case ('-rdf')
    imod=1
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
    print '(a)', '  -noecho           Suppres questions (for non-interactive input). Useful for scripts.'
    print '(a)', '  -rdf              Calculate radial distribution function.'
    print '(a)', '  -box <size_x size_y size_z>  Dimensions of periodic boundary box.'
end subroutine PrintHelp

subroutine PrintInputError()
  write(*,*)'Error during reading input. Exiting...'
  stop 1
end subroutine PrintInputError

! Various estimates to number of bins
! See https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
function square_root_bin_number(n) result(k)
   integer, intent(in) :: n
   integer :: k
   k = ceiling(sqrt(REAL(n)))
end function square_root_bin_number

function sturges_bin_number(n) result(k)
   integer, intent(in) :: n
   integer :: k
   k = ceiling(log(REAL(n))/log(2.0d0)) + 1
end function sturges_bin_number

function rice_bin_number(n) result(k)
   integer, intent(in) :: n
   integer :: k
   k = ceiling(2*(REAL(n)**(1/3.d0)))
end function rice_bin_number

function scott_bin_width(n, sigma) result(h)
   integer, intent(in) :: n
   real*8, intent(in) :: sigma
   real*8 :: h
   h = 3.5d0 * sigma / n**(1/3.0d0)
end function scott_bin_width

! TODO Freedman-Diaconis choice


subroutine print_bin_number_estimators(ngeom)
   integer, intent(in) :: ngeom
   integer :: square_root_bin_number, sturges_bin_number, rice_bin_number

   write(*,'(A)')'Printing various estimators for optimal number of bins'
   write(*,'(A)')'Square root estimate'
   write(*,'(A4,I5)')'k = ', square_root_bin_number(ngeom)
   write(*,'(A)')'Sturges estimate'
   write(*,'(A4,I5)')'k = ', sturges_bin_number(ngeom)
   write(*,'(A)')'Rice estimate'
   write(*,'(A4,I5)')'k = ', rice_bin_number(ngeom)
end subroutine print_bin_number_estimators

