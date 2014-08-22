! D. Hollas, P. Slavicek 2014

! Program for generation of ABIN restart files from WIGNER sampling
! It should be called within script create_wignertrajectories.sh
! It reads some vital data from file "initial".

! It cooperates with script wigner_sampling.sh
! Expects standard input from FMSTRAJS/TRAJ.$i
! Atom masses are determined from the library.
! Atom names must be provided in a separate file atoms.dat

! Call as: ./make_restart < FMSTRAJS/TRAJ.xx

   program MakeRestart
      implicit none
      integer,parameter :: npartmax=1000
      real*8,parameter    :: mu=1823.d0,ang=1.8897d0
      real*8  :: x(npartmax),p(npartmax),y(npartmax),amass(npartmax)
      character*2 :: atom(npartmax)
      integer :: i,natom,nstate,istate
      integer :: index_x,index_y,index_z
      real*8  :: time

      open(100,file='initial',status='old', action="read")
      read(100,*)istate,nstate,natom
      close(100)

      
      ! Read Wigner data, expecting atomic units
      read(*,*)time,(x(i),i=1,3*natom),(p(i),i=1,3*natom)
 
       
      do i=1,natom*3
        y(i)=x(i)/ang !Bohr to angstroms
      enddo
! INITIALIZING MASS
      open(103,file='atoms.dat',status='old', action="read")
      do i=1,natom
       read(103,*)atom(i)
      enddo
      close(103)
      call mass_init(natom,atom,npartmax,amass)
      do i=1,natom
       amass(i)=amass(i)*mu
      enddo

      open(201,file='restart.xyz', action="write")
      open(202,file='mini.dat', action="write")

      write(201,*)0
      write(201,*)' Cartesian Coordinates [au]'
      do i=1,natom
      index_x=(i-1)*3+1
      index_y=(i-1)*3+2
      index_z=(i-1)*3+3
      write(201,20)x(index_x),x(index_y),x(index_z)
      enddo
      write(201,*)' Cartesian Velocities [au]'
      do i=1,natom
      index_x=(i-1)*3+1
      index_y=(i-1)*3+2
      index_z=(i-1)*3+3
      write(201,20)p(index_x)/amass(i),p(index_y)/amass(i),p(index_z)/amass(i)
      enddo
      if(nstate.gt.1)then
       write(201,*)' Coefficients for SH'
       write(201,*)istate
       do i=1,nstate
        if (i.eq.istate) then
         write(201,*)'1.000000000000000E+000  0.000000000000000E+000'
        else
         write(201,*)'0.000000000000000E+000  0.000000000000000E+000'
        endif
       enddo
      endif
      write(201,*)' Cumulative averages of various estimators'
      write(201,*)'0.000000000000000E+000'
      write(201,*)'0.000000000000000E+000  0.000000000000000E+000'
      write(201,*)'0.000000000000000E+000'
20    format(3F15.10)

      write(202,*)natom
      write(202,*)' Cartesian Coordinates from Wigner sampling'
      do i=1,natom
      index_x=(i-1)*3+1
      index_y=(i-1)*3+2
      index_z=(i-1)*3+3
      write(202,30)atom(i),x(index_x)/ang,x(index_y)/ang,x(index_z)/ang
      enddo

30    format(A,3F15.10)
      close(201)
      close(202)
   end

   subroutine mass_init(natom,names,npartmax,am)
      implicit none
      character*2 :: names(npartmax)
      real*8  :: am(npartmax) 
      integer :: i,natom,npartmax
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
        write(*,*)'Unknown atom.',names(i),'Exiting..'
        stop 1
       endif
      enddo
   end subroutine
