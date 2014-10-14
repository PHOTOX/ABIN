!--- File with modules                                created by Daniel Hollas,9.2.2012

!!!!!! We are now using modules to initialize some variables and also
!!!!!! for passing global variables to different subroutines.
!------------------------------------------------------------------------------------



!-----It also contains some physical constants.
      ! Exact values for ANG,AUTOKCAL and AUTODEBYE,me,AUTOM,AUTOFS,AMU: 
      ! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
      ! and using the thermochemical calorie (1 cal = 4.184 J):'
      module mod_const
      implicit none
      INTEGER, PARAMETER  :: DP = KIND(1.0d0)
      real(DP), parameter :: AMU=1822.888484264545d0,ANG=1.889726132873d0
      real(DP), parameter :: AUTOFS=0.02418884326505d0,PI=3.14159265358979323846d0
      real(DP), parameter :: AUTOK=3.1577464d5,ME=9.10938215d-31 !electron mass
      real(DP), parameter :: AUTOM=5.2917720859d-11 ! atomic length
      real(DP), parameter :: AMBTOAU=0.84d0/15.3067320d0 !charges to  MACSIMUS in init
      real(DP), parameter :: AUTOKCAL=6.2750946943d2,AUTOKK=3.1577322d5,AUTOEV=27.21138386d0
      real(DP), parameter :: AUTODEBYE =  2.54174623d0
      real(DP), parameter :: KBinAU = 0.9999999748284666d0  !boltzmann units in au,not used
      save
      end module  

!-----mod_array_size contains various array limits. Modify here if you need larger arrays (for big systems)
      module mod_array_size
      use mod_const, only: DP 
      implicit none
      integer,parameter   :: MAXCHAIN=10, MAXTYPES = 10
      integer,parameter   :: NBINMAX=2000, NDISTMAX=30
      integer,parameter   :: NSTMAX=20, NTRAJMAX=1
      save
      end module  

!------------------------------------------------------
!------------------------------------------------------

!--General simulation parameters      
module mod_general
   use mod_const, only: DP
   implicit none
   integer :: it=0,ipimd=0,istage=0,nwalk=1,ihess=0,md=1  
   character(len=10) :: pot='none'
   integer :: imini=0,nstep=1,nabin=50
   integer :: nwrite=1,nwritex=1,ncalc=1,nrest=1,nwritev=0
   integer :: icv=0,irest=0,anal_ext=0,idebug=0
   integer :: irandom=156873,natom=0,pid
   integer :: iqmmm=0
   integer :: parrespa=0 
   save
end module
      
!------------------------------------------------------
!------------------------------------------------------

!-----Some information about simulated system, especially for distributions and shake
      module mod_system
      use mod_const, only: DP
      implicit none
      real(DP),allocatable  :: am(:)
      character(len=2), allocatable :: names(:)
      integer,allocatable  :: inames(:)
      integer              :: dime=3,f=3 !dimenze systemu a pocet zakonu zachovani 
      integer              :: conatom=0
!---distributions (distance,angle,dihedral)
      save
      CONTAINS
      subroutine mass_init(masses, massnames)
      use mod_const,    only: AMU
      use mod_general,  only: natom
      real(DP) :: masses(:)
      character(len=2) :: massnames(:)
      integer  :: i
      ! Accurate values for H1 and H2 taken from: 
      ! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
      ! Other atomic weights taken from Handbook of Chemistry and Physics, 2013
      ! Original citation: Wieser, M. E., et al., Pure Appl. Chem. 85, 1047, 2013
      allocate( am(natom) )
      am=-1.0d0
      print *,''
      print *,'-----------------ATOMIC MASSES----------------------------'
      do i=1,natom
       if(names(i).eq.'H')then
        am(i)=1.008d0
       else if(names(i).eq.'H1')then
        am(i)=1.00782503207d0 
       else if(names(i).eq.'H2'.or.names(i).eq.'D')then
        am(i)=2.01410177785d0
       else if(names(i).eq.'O')then
        am(i)=15.999d0
       else if (names(i).eq.'S')then
        am(i)=32.06d0
       else if (names(i).eq.'SE')then
        am(i)=78.971d0
       else if (names(i).eq.'TE')then
        am(i)=127.60d0
       else if (names(i).eq.'N')then
        am(i)=14.007d0
       else if (names(i).eq.'P')then
        am(i)=30.973761998d0
       else if (names(i).eq.'AS')then
        am(i)=74.921595d0
       else if (names(i).eq.'SB')then
        am(i)=121.760d0
       else if (names(i).eq.'BI')then
        am(i)=208.98040d0
       else if (names(i).eq.'F')then
        am(i)=18.998403163d0
       else if (names(i).eq.'CL')then
        am(i)=35.45d0
       else if (names(i).eq.'BR')then
        am(i)=79.904d0
       else if (names(i).eq.'I')then
        am(i)=126.90447d0
       else if (names(i).eq.'LI')then
        am(i)=6.94d0
       else if (names(i).eq.'NA')then
        am(i)=22.98976928d0
       else if (names(i).eq.'K')then
        am(i)=39.0983d0
       else if (names(i).eq.'BE')then
        am(i)=9.0121831d0
       else if (names(i).eq.'MG')then
        am(i)=24.305d0
       else if (names(i).eq.'CA')then
        am(i)=40.078d0
       else if (names(i).eq.'B')then
        am(i)=10.81d0
       else if (names(i).eq.'AL')then
        am(i)=26.9815385d0
       else if (names(i).eq.'C')then
        am(i)=12.0110d0
       else if (names(i).eq.'SI')then
        am(i)=28.085d0
       else if (names(i).eq.'GE')then
        am(i)=72.630d0
       else if (names(i).eq.'SN')then
        am(i)=118.710d0
       else if (names(i).eq.'PB')then
        am(i)=207.2d0
       else if (names(i).eq.'HE')then
        am(i)=4.002602d0
       else if (names(i).eq.'NE')then
        am(i)=20.1797d0
       else if (names(i).eq.'AR')then
        am(i)=39.948d0
       else if (names(i).eq.'KR')then
        am(i)=83.798d0
       else if (names(i).eq.'XE')then
        am(i)=131.293d0
       else if (names(i).eq.'FE')then
        am(i)=55.845d0
       else if (names(i).eq.'TI')then
        am(i)=47.867d0
       else if (names(i).eq.'V')then
        am(i)=50.9415d0
       else if (names(i).eq.'CR')then
        am(i)=51.9961d0
       else if (names(i).eq.'MN')then
        am(i)=54.938044d0
       else if (names(i).eq.'CO')then
        am(i)=58.933194d0
       else if (names(i).eq.'NI')then
        am(i)=58.6934d0
       else if (names(i).eq.'CU')then
        am(i)=63.546d0
       else if (names(i).eq.'ZN')then
        am(i)=65.38d0
       else if (names(i).eq.'AG')then
        am(i)=107.8682d0
       else if (names(i).eq.'AU')then
        am(i)=196.966569d0
       else if (names(i).eq.'PT')then
        am(i)=195.084d0
       else if (names(i).eq.'CD')then
        am(i)=112.414d0
       else if (names(i).eq.'HG')then
        am(i)=200.592d0
       else if (names(i).eq.'U')then
        am(i)=238.02891d0
       else if (names(i).eq.'TL')then
        am(i)=204.38d0
       else if (names(i).eq.'BA')then
        am(i)=137.327d0
       else if (names(i).eq.'CE')then
        am(i)=140.116d0
       else if (names(i).eq.'CS')then
        am(i)=132.90545196d0
       else if (names(i).eq.'DY')then
        am(i)=162.500d0
       else if (names(i).eq.'ER')then
        am(i)=167.259d0
       else if (names(i).eq.'EU')then
        am(i)=151.964d0
       else if (names(i).eq.'GD')then
        am(i)=157.25d0
       else if (names(i).eq.'GA')then
        am(i)=69.723d0
       else if (names(i).eq.'HF')then
        am(i)=178.49d0
       else if (names(i).eq.'HO')then
        am(i)=164.93033d0
       else if (names(i).eq.'IN')then
        am(i)=114.818d0
       else if (names(i).eq.'IR')then
        am(i)=192.217d0
       else if (names(i).eq.'LA')then
        am(i)=138.90547d0
       else if (names(i).eq.'LU')then
        am(i)=174.9668d0
       else if (names(i).eq.'MO')then
        am(i)=95.95d0
       else if (names(i).eq.'ND')then
        am(i)=144.242d0
       else if (names(i).eq.'NB')then
        am(i)=92.90637d0
       else if (names(i).eq.'OS')then
        am(i)=190.23d0
       else if (names(i).eq.'PD')then
        am(i)=106.42d0
       else if (names(i).eq.'PR')then
        am(i)=140.90766d0
       else if (names(i).eq.'PA')then
        am(i)=231.03588d0
       else if (names(i).eq.'RE')then
        am(i)=186.207d0
       else if (names(i).eq.'RH')then
        am(i)=102.90550d0
       else if (names(i).eq.'RB')then
        am(i)=85.4678d0
       else if (names(i).eq.'RU')then
        am(i)=101.07d0
       else if (names(i).eq.'SM')then
        am(i)=150.36d0
       else if (names(i).eq.'SC')then
        am(i)=44.955908d0
       else if (names(i).eq.'SR')then
        am(i)=87.62d0
       else if (names(i).eq.'TA')then
        am(i)=180.94788d0
       else if (names(i).eq.'TB')then
        am(i)=158.92535d0
       else if (names(i).eq.'TH')then
        am(i)=232.0377d0
       else if (names(i).eq.'TM')then
        am(i)=168.93422d0
       else if (names(i).eq.'W')then
        am(i)=183.84d0
       else if (names(i).eq.'YB')then
        am(i)=173.054d0
       else if (names(i).eq.'Y')then
        am(i)=88.90584d0
       else if (names(i).eq.'ZR')then
        am(i)=91.224d0
       else 
          print *,'Atom name ',names(i),' was not found in the library.'
          print *,'I hope you specified the mass in namelist "system"'
       endif
      enddo

      call init_usermass()

      do i=1,natom
        if(am(i).le.0)then
           write(*,*)'ERROR: Some masses were not specified. Exiting...'
           stop 1
        end if
      end do

      if (natom.le.100)then
         print *,'These are the atom names I found:'
         do i=1,natom
            write(*,'(A2, A1)',advance='no')names(i),' '
         end do
         write(*,*)
         print *,'The corresponding relative atomic masses are:'
         write(*,*)(am(i),i=1,natom)
      end if
      print *,'----------------------------------------------------------'
      print *,''

      am=am*AMU

      contains

         subroutine init_usermass()
         integer :: iat, j
            
         do i=1,size(massnames)
            do j=i+1,size(massnames)
               if (massnames(i).eq.massnames(j))then
                  if((masses(i).gt.0.or.masses(j).gt.0))then
                     write(*,*)'ERROR: ambiguous user input for masses.'
                     write(*,*)'Please, take a hard look at the input arrays "masses" and "massnames"'
                     stop 1
                  end if
               end if
             end do
         end do

         do iat=1,natom
            do j=1,size(massnames)
            if (names(iat).eq.massnames(j))then

               if (masses(j)>0)then
                  write(*,*)'Defining new atom ',names(iat),' with mass=',masses(j)
                  am(iat)=masses(j)
               else
                  write(*,*)'Mass cannot be negative. Please, fix arrays "masses" or "mass_names" in your input.'
               end if

            end if
            end do
         end do

      end subroutine init_usermass
   end subroutine

   subroutine constrainP (px,py,pz)
      use mod_general, only: nwalk
      real(DP),intent(inout)  :: px(:,:),py(:,:),pz(:,:)
      integer               :: iw,iat
      do iw=1,nwalk
         do iat=1,conatom
            px(iat,iw)=0.0d0
            py(iat,iw)=0.0d0
            pz(iat,iw)=0.0d0
         enddo
      enddo
   end subroutine constrainP

end module 
!------------------------------------------------------
!------------------------------------------------------
