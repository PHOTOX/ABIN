subroutine force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
   use mod_const, only: DP, ANG
   use mod_files, only: MAXUNITS
   use mod_general, only: ihess, ipimd, iqmmm, it, iremd, my_rank
   use mod_system, only: names
   use mod_harmon, only: hess
   use mod_sh_integ, only: nstate
   use mod_sh, only: tocalc, en_array, istate
   use mod_lz, only: nstate_lz, tocalc_lz, en_array_lz, istate_lz, nsinglet_lz, ntriplet_lz
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror, toupper
   use mod_io, only: read_forces
   use mod_interfaces, only: oniom
   implicit none
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: eclas
   integer, intent(in) :: walkmax
   character(len=*), intent(in) :: chpot
   real(DP) :: temp1
   character(len=100) :: chsystem
   character(len=30) :: chgeom, chforce, chhess, fgeom
   logical :: file_exists
   integer :: iat, iw, iat1, iat2, itest !,nthreads=1, ithread
   integer :: ist1, iost, ISTATUS
   integer :: system
!!$ integer :: omp_get_max_threads,OMP_get_thread_num

!!$ nthreads=omp_get_max_threads()

   eclas = 0.0D0

!  Format for geom.dat; needed, so that Molpro can read it
   fgeom = '(A2,3E25.17E2)'

!$OMP PARALLEL ! REDUCTION(+:eclas) alternativa k atomic
!$OMP DO PRIVATE(temp1,chsystem,chgeom,chforce,chhess,itest,file_exists,iost)

   do iw = 1, walkmax

!!$   ithread=OMP_get_thread_num()
      write (chgeom, '(A,I3.3)') 'geom.dat.', iw
      write (chforce, '(A,I3.3)') 'engrad.dat.', iw
      write (chhess, '(A,I3.3)') 'hessian.dat.', iw
      if (iremd == 1) then
         write (chgeom, '(A,I2.2)') trim(chgeom)//'.', my_rank
         write (chforce, '(A,I2.2)') trim(chforce)//'.', my_rank
         write (chhess, '(A,I2.2)') trim(chhess)//'.', my_rank
      end if

!     Delete the last geometry
      open (unit=MAXUNITS + iw, file=chgeom)
      close (unit=MAXUNITS + iw, status='delete')
!     WRITING GEOMETRY IN ANGSTROMS
      open (unit=MAXUNITS + iw, file=chgeom, action='write', access='SEQUENTIAL')
      do iat = 1, natqm
         write (MAXUNITS + iw, fgeom, iostat=iost) names(iat), &
                                                & x(iat, iw) / ANG, &
                                                & y(iat, iw) / ANG, &
                                                & z(iat, iw) / ANG
      end do
      close (unit=MAXUNITS + iw)

!     Surface hopping or Ehrenfest
      if (ipimd == 2 .or. ipimd == 4) then

         open (unit=MAXUNITS + iw + 2 * walkmax, file='state.dat')
         write (MAXUNITS + iw + 2 * walkmax, '(I2)') nstate

!        Diagonal of tocalc holds info about needed forces
!        tocalc(x,x)= 1 -> compute forces for electronic state X
!        totalc for Ehrenfest muset be set just for required states according to c coef. TO-DO in ehrenfest enrehfest_forces
         do ist1 = 1, nstate
            write (MAXUNITS + iw + 2 * walkmax, '(I1,A1)', advance='no') tocalc(ist1, ist1), ' '
         end do
         close (MAXUNITS + iw + 2 * walkmax)
      end if

!     Landau-Zener
      if (ipimd == 5) then

         open (unit=MAXUNITS + iw + 2 * walkmax, file='state.dat')
         write (MAXUNITS + iw + 2 * walkmax, '(I2)') nstate_lz !How many el. states

         !do ist1=1,nstate_lz
         !   write(MAXUNITS+iw+2*walkmax,'(I1,A1)',advance='no')tocalc_lz(ist1),' '
         !end do
         !First we have singlets, then triplets
         do ist1 = 1, nstate_lz
            if (tocalc_lz(ist1) == 1) write (MAXUNITS + iw + 2 * walkmax, '(I2,A1)') ist1, ' ' !Number of gradient state
         end do
         write (MAXUNITS + iw + 2 * walkmax, '(I2,A1)') nsinglet_lz, ' ' !Number of singlets
         write (MAXUNITS + iw + 2 * walkmax, '(I2,A1)') ntriplet_lz, ' ' !Number of triplets
         close (MAXUNITS + iw + 2 * walkmax)
      end if

!     HERE we decide which program we use to obtain gradients and energies
!     e.g. ./G09/r.g09
      chsystem = './'//trim(toupper(chpot))//'/r.'//chpot

      inquire (FILE=chsystem, EXIST=file_exists)
      if (.not. file_exists) then
         write (*, *) 'File ', chsystem
         write (*, *) 'does not exist! Exiting...'
         call abinerror('force_abin')
      end if

      ! Passing arguments to bash script
      ! First argument is time step
      ! Second argument is the bead index, neccessary for parallel calculations
      write (chsystem, '(A40,I13,I4.3)') chsystem, it, iw

      if (iremd == 1) write (chsystem, '(A,I2.2)') trim(chsystem)//'.', my_rank

!     for SH, pass the 4th parameter: precision of forces as 10^(-force_accu1)
!     TODO: This should not be hard-coded
      if (ipimd == 2 .or. ipimd == 4 .or. ipimd == 5) then
         write (chsystem, '(A60,I3,A12)') chsystem, 7, ' < state.dat'
      end if

      !-----MAKE THE CALL----------!
      ISTATUS = system(chsystem)

      ! Exit status 0 turns to 0
      ! For some reason, exit status 1 turns to 256
      ! However, this one we get by default from bash, don't know why...
      ! see this thread for explanation:
      ! http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/MAXUNITS07-01/msg00085.html
      ! If the bash script wants to notify ABIN, it can use e.g. exit 2
      if (ISTATUS /= 0 .and. ISTATUS /= 256) then
         write (*, '(A)') 'ERROR during the execution of the ab initio external program.'
         write (*, '(A)') 'Please inspect the output files in&
         & folder '//trim(toupper(chpot))//"/"
         call abinerror('force_abin')
      end if

!     make sure that the file exist and flush the disc buffer
      itest = 0
      inquire (FILE=chforce, EXIST=file_exists)
      do while (.not. file_exists .and. itest < 10)
         write (*, *) 'WARNING:File ', chforce, ' does not exist. Waiting..'
         ISTATUS = system('sync') !mel by zajistit flush diskoveho bufferu
         inquire (FILE=chforce, EXIST=file_exists)
         itest = itest + 1
      end do

      open (unit=MAXUNITS + iw, file=chforce, status='old', ACTION='READ', IOSTAT=iost)
      if (iost /= 0) then
         write (*, *) 'Fatal problem when trying to open the file ', chforce
         call abinerror('force_abin')
      end if

      ! READING ENERGY from engrad.dat
      read (MAXUNITS + iw, *, IOSTAT=iost) temp1
      if (iost /= 0) then
         write (*, *) 'ERROR: Could not read energy from file ', chforce
         write (*, *) 'Fortran ERROR = ', iost
         write (*, *) 'This usually means, that the ab initio program failed to converge.'
         write (*, *) 'See the appropriate output files from the external program in folder '
         write (*, *) trim(toupper(chpot))
         call abinerror('force_abin')
      end if
!$OMP ATOMIC
      eclas = eclas + temp1
! SH
      ! TODO: Have each state in different file?
      if (ipimd == 2 .or. ipimd == 4) then
         en_array(1, iw) = temp1
         do ist1 = 2, nstate
            read (MAXUNITS + iw, *) en_array(ist1, iw)
         end do
         ! TODO-EH: eclas must be correctly overwritten later
         eclas = en_array(istate(iw), iw)
      end if
! LZ
      if (ipimd == 5) then
         !Move old energies by 1
         en_array_lz(:, 3) = en_array_lz(:, 2)
         en_array_lz(:, 2) = en_array_lz(:, 1)
         !Store the new one
         en_array_lz(1, 1) = temp1
         do ist1 = 2, nstate_lz
            read (MAXUNITS + iw, *, IOSTAT=iost) en_array_lz(ist1, 1)
            if (iost /= 0) then
               write (*, *) 'ERROR: Could not read excited state energy from file ', chforce
               write (*, *) 'Fortran ERROR = ', iost
               call abinerror('force_abin')
            end if
         end do
         eclas = en_array_lz(istate_lz, 1)
      end if

!     TODO-EH: Read additional forces probably somewhere here, use second index (iw) for different states
!     Actually, we should make a routine read_engrad() and make it general for both EH and SH
!     always read energies, read forces based on tocalc

      if (ipimd == 2 .or. ipimd == 4) then
         iost = read_forces(fx, fy, fz, natqm, tocalc(istate(iw), istate(iw)), MAXUNITS + iw)
      else if (ipimd == 5) then
         iost = read_forces(fx, fy, fz, natqm, tocalc_lz(istate_lz), MAXUNITS + iw) !Save only the computed state
      else
!     reading energy gradients from engrad.dat
         iost = read_forces(fx, fy, fz, natqm, iw, MAXUNITS + iw)
      end if
      if (iost /= 0) then
         write (*, *) 'ERROR: Could not read gradients from file ', chforce
         write (*, *) 'Fortran ERROR = ', iost
         write (*, *) 'This usually means, that the ab initio program failed.'
         write (*, *) 'See the appropriate output files from external program in folder ' &
            //trim(toupper(chpot))//"/."
         call abinerror('force_abin')
      end if

      ! READING of HESSIAN
      if (ihess == 1) then

         inquire (FILE=chhess, EXIST=file_exists)
         do while (.not. file_exists .and. itest < 10)
            write (*, *) 'WARNING:File ', chhess, ' does not exist. Waiting..'
            ISTATUS = system('sync') !mel by zajistit flush diskoveho bufferu
            inquire (FILE=chhess, EXIST=file_exists)
            itest = itest + 1
         end do

         open (unit=MAXUNITS + iw + walkmax, file=chhess, status='old', ACTION='READ')

         do iat2 = 1, natqm * 3
            do iat1 = 1, natqm * 3, 3
               read (MAXUNITS + iw + walkmax, *) hess(iat1, iat2, iw), &
                                               & hess(iat1 + 1, iat2, iw), &
                                               & hess(iat1 + 2, iat2, iw)
               hess(iat1, iat2, iw) = hess(iat1, iat2, iw) / walkmax
               hess(iat1 + 1, iat2, iw) = hess(iat1 + 1, iat2, iw) / walkmax
               hess(iat1 + 2, iat2, iw) = hess(iat1 + 2, iat2, iw) / walkmax
            end do
         end do

      end if

      close (unit=MAXUNITS + iw, status='delete')
      if (ihess == 1) then
         close (unit=MAXUNITS + iw + walkmax, status='delete')
      end if

      if (iqmmm == 1) call oniom(x, y, z, fx, fy, fz, eclas, iw)

   end do
!$OMP END DO
!$OMP END PARALLEL

end

subroutine oniom(x, y, z, fx, fy, fz, eclas, iw)
   use mod_const, only: DP, ANG
   use mod_files, only: MAXUNITS
   use mod_general, only: natom, it
   use mod_system, only: names
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror
   implicit none
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(inout) :: eclas
   integer, intent(in) :: iw
   character(len=100) :: chsystem
   character(len=20) :: chgeom, chforce, fgeom
   real(DP) :: temp1, tempx, tempy, tempz
   logical :: file_exists
   integer :: iat, iost, itest

   write (chgeom, '(A,I3.3)') 'geom_mm.dat.', iw
   write (chforce, '(A,I3.3)') 'engrad_mm.dat.', iw
   write (chsystem, '(A)') './MM/r.mm '

   fgeom = '(A2,3E25.17E2)'

   inquire (FILE=chsystem, EXIST=file_exists)
   if (.not. file_exists) then
      write (*, *) 'File ', chsystem
      write (*, *) 'does not exist! Exiting...'
      call abinerror('oniom')
   end if

   write (chsystem, '(A20,I13,I4.3)') chsystem, it, iw

   ! WRITING GEOMETRY of the whole system
   open (unit=MAXUNITS + iw, file=chgeom, action='write')
   do iat = 1, natom
      write (MAXUNITS + iw, fgeom) names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
   end do
   close (unit=MAXUNITS + iw)

   call system(chsystem)

   ! make sure that the file exist and flush the disc buffer
   itest = 0
   inquire (FILE=chforce, EXIST=file_exists)
   do while (.not. file_exists .and. itest < 10)
      write (*, *) 'WARNING:File ', chforce, ' does not exist. Waiting..'
      call system('sync') !mel by zajistit flush diskoveho bufferu
      inquire (FILE=chforce, EXIST=file_exists)
      itest = itest + 1
   end do

   open (unit=MAXUNITS + iw, file=chforce, status='old', ACTION='READ')

   ! READING ENERGY from engrad_mm.dat
   read (MAXUNITS + iw, *, IOSTAT=iost) temp1
   if (iost /= 0) then
      write (*, *) 'Fatal problem with reading energy from file ', chforce
      write (*, *) 'This usually means, that the a program failed.'
      write (*, *) 'See the appropriate output files in folder MM/.'
      call abinerror('oniom')
   end if

!$OMP ATOMIC
   eclas = eclas + temp1

   ! READING energy gradients from engrad.dat
   do iat = 1, natom
      read (MAXUNITS + iw, *, IOSTAT=iost) tempx, tempy, tempz
      if (iost /= 0) then
         write (*, '(2A)') 'Fatal problem with reading gradients from file ', chforce
         write (*, *) 'This usually means, that the ab initio program failed.'
         write (*, *) 'See the appropriate output files in folder MM/.'
         call abinerror('oniom')
      end if
      ! Conversion from gradients to forces
      fx(iat, iw) = fx(iat, iw) - tempx
      fy(iat, iw) = fy(iat, iw) - tempy
      fz(iat, iw) = fz(iat, iw) - tempz
   end do

   close (unit=MAXUNITS + iw, status='delete')

!-----------------MM, only model QM part--------------------------

   ! WRITING GEOMETRY of the QM part
   open (unit=MAXUNITS + iw, file=chgeom)
   do iat = 1, natqm
      write (MAXUNITS + iw, fgeom) names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
   end do
   close (unit=MAXUNITS + iw)

   call system(chsystem)

   ! make sure that the file exist and flush the disc buffer
   itest = 0
   inquire (FILE=chforce, EXIST=file_exists)
   do while (.not. file_exists .and. itest < 10)
      write (*, *) 'WARNING:File ', chforce, ' does not exist. Waiting..'
      call system('sync') !mel by zajistit flush diskoveho bufferu
      inquire (FILE=chforce, EXIST=file_exists)
      itest = itest + 1
   end do

   open (unit=MAXUNITS + iw, file=chforce, status='old', ACTION='READ')

   ! READING ENERGY from engrad_mm.dat
   read (MAXUNITS + iw, *, IOSTAT=iost) temp1
   if (iost /= 0) then
      write (*, *) 'Fatal problem with reading energy from file ', chforce
      write (*, *) 'This usually means, that the external program failed.'
      write (*, *) 'See the appropriate output files in folder MM/.'
      call abinerror('oniom')
   end if

!$OMP ATOMIC
   eclas = eclas - temp1

   ! READING gradients from engrad_mm.dat
   do iat = 1, natqm
      read (MAXUNITS + iw, *, IOSTAT=iost) tempx, tempy, tempz
      if (iost /= 0) then
         write (*, '(2A)') 'Fatal problem with reading gradients from file ', chforce
         write (*, *) 'This usually means, that the external program failed.'
         write (*, *) 'See the appropriate output files in folder MM/.'
         call abinerror('force_abin')
      end if
      ! Conversion to forces
      fx(iat, iw) = fx(iat, iw) + tempx
      fy(iat, iw) = fy(iat, iw) + tempy
      fz(iat, iw) = fz(iat, iw) + tempz
   end do

   close (unit=MAXUNITS + iw, status='delete')

end subroutine oniom
