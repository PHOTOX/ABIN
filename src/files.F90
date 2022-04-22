! Module for permanent file handling

! Note that we're not trying to explicitly handle I/O errors here.
! If open() or write fails, we crash ungracefully.
module mod_files
   use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
   implicit none
   public
   ! Defines maximum number of units available for permanently opened files
   integer, parameter, private :: MAXUNITS = 50

   integer, parameter, private :: MAXFILENAME = 50
   character(len=MAXFILENAME), private :: CHFILES(MAXUNITS)

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
   ! Other 
   integer, parameter :: UERMD = 39

   ! Default standard output and standard error units
   ! These are NOT parameters, we change them in REMD and in unit tests.
   integer, protected :: stdout = OUTPUT_UNIT
   integer, protected :: stderr = ERROR_UNIT
   save

contains

   subroutine stdout_to_devnull()
      use mod_error, only: fatal_error
      integer :: u, iost
      open (newunit=u, iostat=iost, file='/dev/null', action='write')
      if (iost == 0) then
         stdout = u
      else
         call fatal_error(__FILE__, __LINE__, 'Could not open file /dev/null')
      end if
   end subroutine stdout_to_devnull

   subroutine reset_stdout()
      if (stdout /= OUTPUT_UNIT) then
         close (stdout)
         stdout = OUTPUT_UNIT
      end if
   end subroutine reset_stdout

   subroutine stderr_to_stdout()
      stderr = stdout
   end subroutine stderr_to_stdout

   subroutine reset_stderr()
      stderr = ERROR_UNIT
   end subroutine reset_stderr

   subroutine files_init(isbc, phase, ndist, nang, ndih)
      use mod_general, only: ipimd, irest, iremd, pot, &
                           & icv, ihess, idebug, nwritev, nwritef, en_restraint
      use mod_error, only: fatal_error
      use mod_mpi, only: get_mpi_rank
      integer, intent(in) :: isbc, phase, ndist, nang, ndih
      character(len=10) :: chaccess
      integer :: i

      ! In this code we assume ERROR_UNIT == 0 and OUTPUT_UNIT == 6
      ! Other values might conflict from hard-coded values defined above,
      ! so we just stop early in that case.
      if (ERROR_UNIT /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Non-standard stderr unit, compiler not supported')
      end if
      if (OUTPUT_UNIT /= 6) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Non-standard stdout unit, compiler not supported')
      end if

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

      !Energy restratint
      chfiles(UERMD) = 'en_restraint.dat'

      ! Here we ensure, that previous files are deleted
      if (irest == 0) then
         chaccess = 'SEQUENTIAL'
      else
         chaccess = 'APPEND'
      end if

      chfiles(UMOVIE) = 'movie.xyz'

      if (iremd == 1) then
         do i = 1, MAXUNITS
            write (chfiles(i), '(A,I2.2)') trim(chfiles(i))//'.', get_mpi_rank()
         end do
      end if

      ! Trajectory file is opened later in output function trajout
      ! to prevent creating empty movie.xyz and then failing

      open (UTEMPER, file=chfiles(UTEMPER), access=chaccess, action='write')

      if (nwritev > 0) then
         open (UVELOC, file=chfiles(UVELOC), access=chaccess, action='write')
      end if

      if (nwritef > 0) then
         open (UFORCE, file=chfiles(UFORCE), access=chaccess, action='write')
      end if

      if (ipimd == 1) then
         open (UESTENERGY, file=chfiles(UESTENERGY), access=chaccess, action='write')
      else
         open (UENERGY, file=chfiles(UENERGY), access=chaccess, action='write')
      end if

      if (ipimd == 5) then
         open (UPES, file=chfiles(UPES), access=chaccess, action='write')
         open (UPOP, file=chfiles(UPOP), access=chaccess, action='write')
      end if

      if (ipimd == 2) then
         open (UPOP, file=chfiles(UPOP), access=chaccess, action='write')
         open (UPROB, file=chfiles(UPROB), access=chaccess, action='write')
         open (UPES, file=chfiles(UPES), access=chaccess, action='write')
         open (UNACME, file=chfiles(UNACME), access=chaccess, action='write')
         open (UDOTPROD, file=chfiles(UDOTPROD), access=chaccess, action='write')

         if (idebug > 1) then
            open (UBKL, file=chfiles(UBKL), access=chaccess, action='write')
            open (UWFCOEF, file=chfiles(UWFCOEF), access=chaccess, action='write', recl=250)
            if (phase == 1) then
               open (UPHASE, file=chfiles(UPHASE), access=chaccess, action='write')
            end if
         end if

         if (pot == '_tera_') then
            open (UCHARGES, file=chfiles(UCHARGES), access=chaccess, action='write')
            open (UDOTPRODCI, file=chfiles(UDOTPRODCI), access=chaccess, action='write')
            open (UDIP, file=chfiles(UDIP), access=chaccess, action='write')
            open (UTDIP, file=chfiles(UTDIP), access=chaccess, action='write')
         end if
      end if

      if (ipimd /= 2 .and. pot == '_tera_') then
         open (UCHARGES, file=chfiles(UCHARGES), access=chaccess, action='write')
         open (UDIP, file=chfiles(UDIP), access=chaccess, action='write')
      end if

      if (isbc == 1) then
         open (URADIUS, file=chfiles(URADIUS), access=chaccess, action='write')
      end if

      if (icv == 1) then
         open (UCV, file=chfiles(UCV), access=chaccess, action='write')
         if (ihess == 1) then
            open (UCVDCV, file=chfiles(UCVDCV), access=chaccess, action='write')
         end if
      end if

      ! Analysis
      if (ndist > 0) then
         open (UDIST, file=chfiles(UDIST), access=chaccess, action='write')
      end if
      if (nang > 0) then
         open (UANG, file=chfiles(UANG), access=chaccess, action='write')
      end if
      if (ndih > 0) then
         open (UDIH, file=chfiles(UDIH), access=chaccess, action='write')
      end if

      !Energy restraint
      if (en_restraint > 0) then
         open (UERMD, file=chfiles(UERMD), access=chaccess, action='write')
      end if

      if (irest == 0) then
         call print_file_headers()
      end if

   end subroutine files_init

   subroutine print_file_headers()
      use mod_general, only: ipimd, natom
      use mod_system, only: names
      character(len=200) :: headers(MAXUNITS)
      integer :: i
      logical :: lopened

      headers = ''
      headers(UTEMPER) = '#      Time[fs] Temperature T-Average Conserved_quantity_of_thermostat'
      headers(UENERGY) = '#        Time[fs] E-potential           E-kinetic     E-Total    E-Total-Avg'
      ! PIMD estimators
      headers(UESTENERGY) = '#     Time[fs] E-potential  E-primitive   E-virial  CumulAvg_prim  CumulAvg_vir'
      headers(UCV) = '#         Time[fs]  Cv-prim   Cv-vir  Cv_cumul_prim  Cv_cumul_vir'
      headers(UCVDCV) = '#         Time[fs]  Cv-DCV   Cv_cumul_DCV'

      ! Geometrical analysis
      headers(UDIST) = "# Distances [Angstrom]"
      headers(UANG) = "# Angles [Degree]"
      headers(UDIH) = "# Dihedral Angles [Degree]"

      headers(URADIUS) = '#TimeStep     Radius[ANG]   approximate density[kg.m^3]'

      ! Surface Hopping and Landau-Zener
      if (ipimd == 5) then
         headers(UPES) = '#    Time[fs] Potential energies (singlets, triplets)'
      else
         headers(UPES) = '#    Time[fs] Potential energies'
      end if

      headers(UPOP) = '#    Time[fs] CurrentState   Populations Sum-of-Populations'
      headers(UPHASE) = '# Lower triangular matrix of gamma (phase)  gamma(i,j) [i=1,nstate ;j=1,i-1]'
      headers(UPROB) = '#    Time[fs] CurrentState   Probabilities'
      headers(UDOTPROD) = '#    Time[fs] dotproduct(i,j) [i=1,nstate-1;j=i+1,nstate]'
      headers(UBKL) = '# Hopping probabilities - bkl(i) [i=1,nstate]'
      headers(UWFCOEF) = '# WF coefficients c_real(i),i=1,nstate c_imag(i),i=1,nstate'

      if (ipimd == 2 .or. ipimd == 5) then
         headers(UDIP) = '# Time dip_tot.1 dip_tot.2 ... dip_x.1 dip_y.1 dip_z.1 dip_x.2 dip_y.2 dip_z.2.'
      else
         headers(UDIP) = '# Time Bead_index |D| Dx Dy Dz'
      end if
      headers(UCHARGES) = '# Atomic charges from current electronic state'
      headers(UDOTPRODCI) = '# Time  cidotprod1  cidotprod2 ... '
      headers(UTDIP) = '# Time state tdip_tot.1..N tdip_x.1 tdip_y.1 tdip_z.1 ...'

      headers(UERMD) = '# Step, energy difference ES-GS (eV), deltaE (Ha), deltaE_next (Ha), &
                       &lambda multiplier'

      do i = 2, MAXUNITS
         inquire (unit=i, opened=lopened)
         if (lopened .and. i /= ERROR_UNIT .and. i /= OUTPUT_UNIT) then
            if (len_trim(headers(i)) /= 0) then
               write (i, *) trim(headers(i))
            end if
         end if
      end do

      ! Not sure if this is a good idea if there are many atoms
      inquire (unit=UCHARGES, opened=lopened)
      if (lopened) then
         if (ipimd == 2 .or. ipimd == 5) then
            write (UCHARGES, *) '# Time  state ', (names(i), i=1, natom)
         else
            write (UCHARGES, *) '# Time_step Bead_index ', (names(i), i=1, natom)
         end if
      end if
   end subroutine print_file_headers

   subroutine close_files()
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
