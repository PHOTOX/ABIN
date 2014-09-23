!# cut_sphere           D. Hollas, 2014
!- Program for cutting solvent molecules around the solute
!  from an xyz file.
!- It assumes, that the solute atoms are the first atoms in the file.
!- The coordinates are read from standard input.
!- See help in subroutine PrintHelp().

!  There are two modes of action:
!  i) Pick closest nmol solvent molecules.
!  ii) Pick all solvent molecules within specified radius. 

!- There are 3 ways how to determine the distance of solvent to solute.
!  1. Closest distance between ANY solvent atom and ANY solute atom.
!     This is the default if we want specific number of solvent molecules.
!  2. Distance is determined relative to the geometrical center of the solute.
!     This is the default if we cut specific radius. (option -com)
!  3. Distance is determined relative to one specific atom of the solute.
!     (option -i idx)

program cut_sphere
   implicit none
   integer, parameter :: maxatom=10000
   character(len=2)   :: at(maxatom)
   real*8    :: x(maxatom), y(maxatom), z(maxatom)
   real*8    :: rmin(maxatom),r(maxatom,maxatom), rmolmin(maxatom)
   integer   :: index(maxatom),index_rev(maxatom)
   integer   :: natmol=3, nmol=-1, nsolute=-1, atidx=-1
   integer   :: nmoltotal, natom
   real*8    :: rad=-1
   logical   :: lcom=.false., ltrans=.false.
   integer   :: i, j, iat, imol, idx
   real*8    :: xt, yt, zt


   call Get_cmdline(natmol, nmol, nsolute, rad, lcom, atidx, ltrans)

   read(*,*)natom
   write(*,*)'Total number of atoms:', natom
   if(natom.gt.maxatom)then
      write(*,*)'ERROR: number of atoms is greater than maximum.'
      write(*,*)'Adjust parameter maxatom and recompile'
      stop 1
   end if

   read(*,*)
   do i=1,natom
      read(*,*)at(i),x(i),y(i),z(i)
   enddo

   if (lcom.or.ltrans)then
      xt=0.0
      yt=0.0
      zt=0.0
      do i=1,nsolute
         xt=xt+x(i)
         yt=yt+y(i)
         zt=zt+z(i)
      enddo
      xt=xt/nsolute
      yt=yt/nsolute
      zt=zt/nsolute
      write(*,*)'Coordinates of the geometrical center.'
      write(*,*)xt, yt, zt
   end if

   if (ltrans)then
      write(*,*)'Translating the molecule to the origin.'
      do iat=1,natom
         x(iat)=x(iat)-xt
         y(iat)=y(iat)-yt
         z(iat)=z(iat)-zt
      end do
      xt=0.0d0; yt=0.0d0; zt=0.0d0
   end if

   if(atidx.gt.0)then
      xt=x(atidx)
      yt=y(atidx)
      zt=z(atidx)
      write(*,*)'Coordinates of the atom around which we cut:'
      write(*,*)xt, yt, zt
   end if
 

! determine the closest distance of every solvent atoms to solute
! based on all atoms in the solute
if(.not.lcom.and.atidx.lt.0)then
   do i=nsolute+1,natom
      rmin(i)=10000d0
      do j=1,nsolute
         r(i,j)=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2 
         r(i,j)=dsqrt(r(i,j))
         if(r(i,j).lt.rmin(i))then
            rmin(i)=r(i,j)
         endif
      enddo
   enddo
end if

! determine the distances to the geometrical center
! or to atom with index atidx
   if(lcom.or.atidx.gt.0)then
      do i=nsolute+1,natom
         rmin(i)=(x(i)-xt)**2+(y(i)-yt)**2+(z(i)-zt)**2 
         rmin(i)=dsqrt(rmin(i))
      enddo
   end if

! total number of molecules
   nmoltotal=(natom-nsolute)/natmol

!  Now we have to determine, which atom in each molecule is closest
!  i.e. determine the distance of each solvent molecule from solute
   do imol=1,nmoltotal
      rmolmin(imol)=10000d0
      idx=nsolute+1+(imol-1)*natmol
      do iat=idx,idx+natmol-1
         if(rmin(iat).lt.rmolmin(imol))then
            rmolmin(imol)=rmin(iat)
         endif
      enddo
   enddo
   
! determine, how many molecules are in the radius
   if (rad.gt.0)then
      nmol=0
      do imol=1,nmoltotal
         if(rmolmin(imol).lt.rad)then
            nmol=nmol+1
         end if
      end do
      write(*,*)'Number of solvent molecules within radius ',rad
      write(*,*)'is', nmol
   end if

! Here's where the index ordering happens.
! ###############################################
   do i=1,nmoltotal
      index(i)=1
      do j=1,nmoltotal

         if(i.eq.j) cycle

!        now we have to determine, which molecule is closer
         if(rmolmin(i).gt.rmolmin(j))then
          index(i)=index(i)+1
         endif

      enddo
      index_rev(index(i))=i
   enddo
       

!  WRITE THE RESULTS----------------------
   open(150,file='cut_qm.xyz')
   write(150,*)natmol*nmol+nsolute
   write(150,*)
   ! first write solute
   do idx=1,nsolute
      write(150,*)at(idx),x(idx),y(idx),z(idx)
   end do
   ! now write ordered solvent molecules
   do imol=1,nmol
      do iat=1,natmol
         idx=nsolute+(index_rev(imol)-1)*natmol+iat
         write(150,*)at(idx),x(idx),y(idx),z(idx)
      enddo
   enddo

   close(150)

!  Now print the solvent that we left out.
   open(150,file='cut_mm.xyz')
   write(150,*)natom-natmol*nmol-nsolute
   write(150,*)
   do imol=nmol+1,nmoltotal
      do iat=1,natmol
         idx=nsolute+(index_rev(imol)-1)*natmol+iat
         write(150,*)at(idx),x(idx),y(idx),z(idx)
      enddo
   enddo

   close(150)

end

subroutine Get_cmdline(natmol, nmol, nsolute, rad, lcom, atidx, ltrans)
implicit none
real*8,intent(inout)    :: rad
integer, intent(inout)  :: natmol, nmol, nsolute, atidx
logical, intent(inout)  :: lcom, ltrans
character(len=100)      :: arg
integer                 :: i

i=0
do while (i < command_argument_count())
  i=i+1
  call get_command_argument(i, arg)
  
  select case (arg)
  case ('-h', '--help')
    call PrintHelp()
    stop
  case ('-u')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nsolute
  case ('-v')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)nmol
    if (nmol.le.0)then
       call PrintHelp()
       write(*,*)'ERROR: Number of molecules must be a positive integer.'
       call PrintInputError()
    endif
  case ('-va')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)natmol
  case ('-r')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)rad
    if (rad.le.0)then
       call PrintHelp()
       write(*,*)'ERROR: Radius must be positive.'
       call PrintInputError()
    endif
  case ('-i')
    i=i+1
    call get_command_argument(i, arg)
    read(arg,*)atidx
    if (atidx.le.0)then
       call PrintHelp()
       write(*,*)'ERROR: Atom index must be positive.'
       call PrintInputError()
    endif
  case ('-com')
    lcom=.true.
  case ('-trans')
    ltrans=.true.
  case default
    call PrintHelp()
    write(*,*)'Invalid command line argument!'
    call PrintInputError()
  end select

end do

! INPUT SANITY CHECK
  if(rad.gt.0.and.nmol.gt.0)then
     write(*,*)'ERROR: Conflicting options -v and -r.'
     write(*,*)'You may cut EITHER constant number of solvent molecules,'
     write(*,*)'OR cut sphere of solvent of fixed radius.'
     call PrintHelp()
     call PrintInputError()
  end if

  if(nmol.lt.0.and.rad.lt.0)then
     call PrintHelp()
     write(*,*)'ERROR: Missing option -v or -r.'
     call PrintInputError()
  end if

  if(nsolute.lt.0)then
     call PrintHelp()
     write(*,*)'ERROR: Number of atoms in solute not specified.'
     call PrintInputError()
  end if

  if(natmol.lt.0)then
     call PrintHelp()
     write(*,*)'ERROR: Number of atoms in a solvent molecule must be positive.'
     call PrintInputError()
  end if

  if(rad.gt.0.and.atidx.le.0)then
     lcom=.true.
  end if


end subroutine Get_cmdline

subroutine PrintHelp()
implicit none
    print '(a)', ''
    print '(a)', 'Program for cutting solvent molecules around a given solute from a larger solvated system.'
    print '(a)', 'You may cut either specific number of closest solvent molecules'
    print '(a)', 'or all solvent mocules within given radius.'
    print '(a)', ''
    print '(a)', 'USAGE: ./cut_sphere [OPTIONS] < input.xyz'
    print '(a)', ''
    print '(a)', 'cmdline options:'
    print '(a)', ''
    print '(a)', '  -h or --help     Print this information and exit.'
    print '(a)', '  -u  <integer>    Number of atoms in the solute.'
    print '(a)', '  -va <integer>    Number of atoms in one solvent molecule (default=3).'
    print '(a)', '  -v  <integer>    Number solvent molecules to cut.'
    print '(a)', '  -r  <radius>     Radius of the sphere to cut.'
    print '(a)', '  -i  <atidx>      Cut around the atom with index atidx.'
    print '(a)', '  -com             Cut around the geometrical center of molecule.'
    print '(a)', '                   This is the default if -r is specified.'
    print '(a)', '  -trans           Translate the molecule so that the &
                  coordinates of geometrical center are (0, 0, 0).'
end subroutine PrintHelp

subroutine PrintInputError()
  write(*,*)'Error during reading input. Exiting...'
  stop 1
end subroutine PrintInputError


