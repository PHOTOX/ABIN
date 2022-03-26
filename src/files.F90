! Module for permanent file handling

! Note that we're not trying to explicitly handle I/O errors here.
! If open() or write fails, we crash ungracefully.
module mod_files
   implicit none
   public
   private :: CHFILES, MAXFILENAME
   ! TODO: Make MAXUNITS private, currently used in force_abin.F90
   ! private :: MAXUNITS
   ! Defines maximum number of units available for permanently opened files
   integer, parameter :: MAXUNITS = 50, MAXFILENAME = 50
   character(len=MAXFILENAME) :: CHFILES(MAXUNITS)

   ! UNIT 1 is reserved for CP2K!!!
   integer, parameter :: UMOVIE = 10, UVELOC = 2, UFORCE = 15
   integer, parameter :: UENERGY = 3, UTEMPER = 4
   integer, parameter :: URADIUS = 11
   ! PIMD stuff
   integer, parameter :: UESTENERGY = 12, UCV = 13, UCVDCV = 14
   ! Surface hopping stuff
   integer, parameter :: UPOP = 20, UPROB = 21, UPES = 22
   integer, parameter :: UDOTPROD = 23, UNACME = 24, UWFCOEF = 25
   integer, parameter :: UPHASE = 26, UBKL = 27
   ! So far only for TeraChem
   integer, parameter :: UDOTPRODCI = 31, UCHARGES = 32
   integer, parameter :: UDIP = 33, UTDIP = 34
   ! Analysis output
   integer, parameter :: UDIST = 36, UANG = 37, UDIH = 38
   save

contains

   subroutine files_init(isbc, phase, ndist, nang, ndih)
      use mod_general
      use mod_system, only: names
      integer, intent(in) :: isbc, phase, ndist, nang, ndih
      character(len=10) :: chaccess
      integer :: i

      do i = 1, MAXUNITS
         chfiles(i) = ''
      end do

      chfiles(UVELOC) = 'velocities.xyz'
      chfiles(UFORCE) = 'forces.xyz'
      chfiles(UENERGY) = 'energies.dat'
      chfiles(UTEMPER) = 'temper.dat'

      ! radius for spherical boundary conditions
      chfiles(URADIUS) = 'radius.dat'

      ! Files for PIMD estimators
      chfiles(UESTENERGY) = 'est_energy.dat'
      chfiles(UCV) = 'cv.dat'
      ! file for advanced cv estimator
      chfiles(UCVDCV) = 'cv_dcv.dat'

      ! Files for Surface Hopping
      chfiles(UPOP) = 'pop.dat'
      chfiles(UPROB) = 'prob.dat'
      chfiles(UPES) = 'PES.dat'
      chfiles(UDOTPROD) = 'dotprod.dat'
      chfiles(UNACME) = 'nacm_all.dat'
      chfiles(UWFCOEF) = 'wfcoef.dat'
      chfiles(UPHASE) = 'phase.dat'
      chfiles(UBKL) = 'bkl.dat'
      ! Files for TeraChem SH interface
      chfiles(UCHARGES) = 'charges.dat'
      chfiles(UDIP) = 'dipoles.dat'
      chfiles(UTDIP) = 'trans_dipoles.dat'
      chfiles(UDOTPRODCI) = 'dotprodci.dat'

      ! Geometry analysis output
      chfiles(UDIST) = 'distances.dat'
      chfiles(UANG) = 'angles.dat'
      chfiles(UDIH) = 'dihedrals.dat'

      ! Here we ensure, that previous files are deleted
      if (irest == 0) then
         chaccess = 'SEQUENTIAL'
      else
         chaccess = 'APPEND'
      end if

      chfiles(UMOVIE) = 'movie.xyz'

      if (iremd == 1) then
         do i = 1, MAXUNITS
            write (chfiles(i), '(A,I2.2)') trim(chfiles(i))//'.', my_rank
         end do
      end if

!  OPEN trajectory file
!  Trajectory file is opened later in output function trajout
!  to prevent creating empty movie.xyz and then failing
!  We still open when using CP2K interface, since otherwise
!  strangely more MPI ranks can write to this file if not opened here...
!   if(pot.eq.'_cp2k_')then
!      open(UMOVIE,file=chfiles(UMOVIE),access=chaccess,action='write')
!   end if

      if (nwritev > 0) then
!      if(iremd.eq.1)then
!         write(chout,'(A,I2.2)')'vel.dat.',my_rank
!      else
!         chout='vel.dat'
!      end if
         open (UVELOC, file=chfiles(UVELOC), access=chaccess, action='write')
      end if

      if (nwritef > 0) open (UFORCE, file=chfiles(UFORCE), access=chaccess, action='write')

      if (ipimd /= 1) then
         open (UENERGY, file=chfiles(UENERGY), access=chaccess, action='write')
         write (UENERGY, *) '#        Time[fs] E-potential           E-kinetic     E-Total    E-Total-Avg'
      end if

      if (ipimd == 1) then
         open (UESTENERGY, file=chfiles(UESTENERGY), access=chaccess, action='write')
         write (UESTENERGY, *) '#     Time[fs] E-potential  E-primitive   E-virial  CumulAvg_prim  CumulAvg_vir'
      end if

      open (UTEMPER, file=chfiles(UTEMPER), access=chaccess, action='write')
      write (UTEMPER, *) '#      Time[fs] Temperature T-Average Conserved_quantity_of_thermostat'

      if (ipimd == 5) then
         open (UPES, file=chfiles(UPES), access=chaccess, action='write')
         write (UPES, *) '#    Time[fs] Potential energies (singlets, triplets)'
         open (UPOP, file=chfiles(UPOP), access=chaccess, action='write')
         write (UPOP, *) '#    Time[fs] CurrentState   Populations Sum-of-Populations'
      end if

      if (ipimd == 2) then
         open (UPOP, file=chfiles(UPOP), access=chaccess, action='write')
         write (UPOP, *) '#    Time[fs] CurrentState   Populations Sum-of-Populations'
         open (UPROB, file=chfiles(UPROB), access=chaccess, action='write')
         write (UPROB, *) '#    Time[fs] CurrentState   Probabilities'
         open (UPES, file=chfiles(UPES), access=chaccess, action='write')
         write (UPES, *) '#    Time[fs] Potential energies'
         open (UNACME, file=chfiles(UNACME), access=chaccess, action='write')
         open (UDOTPROD, file=chfiles(UDOTPROD), access=chaccess, action='write')
         write (UDOTPROD, *) '#    Time[fs] dotproduct(i,j) [i=1,nstate-1;j=i+1,nstate]'
         if (idebug > 1) then
            open (UBKL, file=chfiles(UBKL), access=chaccess, action='write')
            write (UBKL, *) '# Hopping probabilities - bkl(i) [i=1,nstate]'
            open (UWFCOEF, file=chfiles(UWFCOEF), access=chaccess, action='write', recl=250)
            write (UWFCOEF, *) '# WF coefficients c_real(i),i=1,nstate c_imag(i),i=1,nstate'
            if (phase == 1) then
               open (UPHASE, file=chfiles(UPHASE), access=chaccess, action='write')
               write (UPHASE, *) '# Lower triangular matrix of gamma (phase)  gamma(i,j) [i=1,nstate ;j=1,i-1]'
            end if
         end if

         if (pot == '_tera_') then
            open (UCHARGES, file=chfiles(UCHARGES), access=chaccess, action='write')
            write (UCHARGES, *) '# Atomic charges from current electronic state'
            write (UCHARGES, *) '# Time  state ', (names(i), i=1, natom)
            open (UDOTPRODCI, file=chfiles(UDOTPRODCI), access=chaccess, action='write')
            write (UDOTPRODCI, *) '# Dot products between current and previous CI vectors.'
            write (UDOTPRODCI, *) '# Time  cidotprod1  cidotprod2 ... '
            open (UDIP, file=chfiles(UDIP), access=chaccess, action='write')
            write (UDIP, *) '# Time dip_tot.1 dip_tot.2 ... dip_x.1 dip_y.1 dip_z.1 dip_x.2 dip_y.2 dip_z.2.'
            open (UTDIP, file=chfiles(UTDIP), access=chaccess, action='write')
            write (UTDIP, *) '# Time  st  tdip_tot.1 tdip_tot.2 ... tdip_x.1 tdip_y.1 tdip_z.1 tdip_x.2 tdip_y.2 tdip_z.2.'
         end if
      end if

      if (ipimd /= 2 .and. pot == '_tera_') then
         open (UCHARGES, file=chfiles(UCHARGES), access=chaccess, action='write')
         write (UCHARGES, *) '# Atomic Mulliken charges from current electronic state'
         write (UCHARGES, *) '# Time_step Bead_index ', (names(i), i=1, natom)

         open (UDIP, file=chfiles(UDIP), access=chaccess, action='write')
         write (UDIP, *) '# Time Bead_index |D| Dx Dy Dz'
      end if

      if (isbc == 1) then
         open (URADIUS, file=chfiles(URADIUS), access=chaccess, action='write')
         write (URADIUS, *) '#TimeStep     Radius[ANG]   approximate density[kg.m^3]'
      end if

      if (icv == 1) then
         open (UCV, file=chfiles(UCV), access=chaccess, action='write')
         write (UCV, *) '#         Time[fs]  Cv-prim   Cv-vir  Cv_cumul_prim  Cv_cumul_vir'
         if (ihess == 1) then
            open (UCVDCV, file=chfiles(UCVDCV), access=chaccess, action='write')
            write (UCVDCV, *) '#         Time[fs]  Cv-DCV   Cv_cumul_DCV'
         end if
      end if

      ! Analysis
      if (ndist > 0) then
         open (UDIST, file=chfiles(UDIST), access=chaccess, action='write')
         write (UDIST, '(A)') "# Distances [Angstrom]"
      end if
      if (nang > 0) then
         open (UANG, file=chfiles(UANG), access=chaccess, action='write')
         write (UANG, '(A)') "# Angles [Degree]"
      end if
      if (ndih > 0) then
         open (UDIH, file=chfiles(UDIH), access=chaccess, action='write')
         write (UDIH, '(A)') "# Dihedral Angles [Degree]"
      end if

   end subroutine files_init

   subroutine close_files()
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
      integer :: i
      logical :: lopen

      do i = 2, MAXUNITS
         inquire (unit=i, opened=lopen)
         if (lopen .and. i /= ERROR_UNIT .and. i /= OUTPUT_UNIT) then
            close (i)
         end if
      end do
   end subroutine close_files

end module mod_files
