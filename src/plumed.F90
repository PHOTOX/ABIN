! Interface to PLUMED, implemented according to
! https://www.plumed.org/doc-v2.7/developer-doc/html/_how_to_plumed_your_m_d.html

! Compilation with PLUMED is optional.
! PLUMED is statically linked with ABIN, see Makefile.
! To install PLUMED libraries, see 'dev_scripts/install_plumed.sh'
module mod_plumed
   use iso_c_binding, only: C_NULL_PTR
   use mod_const, only: DP
#ifndef USE_PLUMED
   use mod_utils, only: not_compiled_with
#endif
   implicit none
   private
   public :: iplumed, plumedfile
   public :: plumed_init, finalize_plumed
   public :: force_plumed
   ! This is only public for unit tests.
   ! Since it is a parameter, it's safe to make it public.
   public :: PLUMED_OUTPUT_FILE

   ! We officially support PLUMED versions >= 2.5.2,
   ! i.e. the last officially supported PLUMED as of 14.1.2021.
   ! We run CI tests for PLUMED 2.5.x, 2.6.x and 2.7.x
   ! MIN_API_VERSION = 4 also allows versions >= 2.3.0,
   ! but you're on your own if you're treading these obsolete waters!
   integer, parameter :: MIN_API_VERSION = 4
   character(len=*), parameter :: PLUMED_OUTPUT_FILE = 'plumed.out'

   integer :: iplumed = 0
   ! PLUMED input file, can be rewritten in ABIN input.
   character(len=40) :: plumedfile = 'plumed.in'
   save

contains

#ifdef USE_PLUMED
   subroutine plumed_init(silent_output)
      use mod_general, only: natom, irest, dt0, nrest
      use mod_const, only: ANG, AUTOFS, AMU, AVOGADRO, AUTOJ
      use mod_utils, only: abinerror, c_string
      !use mod_nhc, only: temp
      implicit none
      ! Do not print any output (used in unit tests, unittest/test_plumed.pf)
      logical, optional :: silent_output
      logical :: silent
      integer, parameter :: PLUMED_RESTART = 1
      ! Conversion from ABIN to PLUMED units
      real(DP), parameter :: PLUMED_ENERGY_UNIT = AUTOJ * AVOGADRO * 1.0D-3 ! Ha -> kJ/mol
      real(DP), parameter :: PLUMED_LENGTH_UNIT = 0.1D0 / ANG ! Bohr -> nm
      real(DP), parameter :: PLUMED_TIME_UNIT = AUTOFS * 0.001D0 ! a.u. -> ps
      real(DP), parameter :: PLUMED_MASS_UNIT = 1.0D0 / AMU ! a.u. -> amu
      !real(DP) :: plumed_kbt
      integer :: api_version
      character(len=32) :: chint

      ! TODO: Fix this! Currently, we're calling plumed_init()
      ! before the temperature is ready, so temp=0.0d0 here!
      ! Also, this conversion only works if temp is in Kelvins,
      ! but we use temp in atomic units in ABIN!
      ! Verify kbt using https://www.plumed.org/doc-v2.5/user-doc/html/kt.html
      !plumed_kbt = temp * 8.3144598D0 ! in kJ/mol

      silent = .false.
      if (present(silent_output)) then
         silent = silent_output
      end if

      if (.not. silent) then
         write (*, *) 'PLUMED is ON'
         write (*, *) 'PLUMED input file is '//trim(plumedfile)
      end if

      ! Initialize the main plumed object.
      ! This must be the very first call.
      call plumed_f_gcreate()

      call plumed_f_gcmd(c_string("getApiVersion"), api_version)
      if (.not. silent) then
         write (*, *) "PLUMED API VERSION: ", api_version
      end if

      if (api_version < MIN_API_VERSION) then
         write (*, *) 'ERROR: PLUMED version not supported'
         write (*, *) 'Please, install a newer PLUMED version and recompile ABIN'
         call abinerror('plumed_init')
      end if

      ! Set units for conversion between PLUMED and ABIN.
      call plumed_f_gcmd(c_string("setRealPrecision"), DP)
      call plumed_f_gcmd(c_string("setMDEnergyUnits"), PLUMED_ENERGY_UNIT)
      call plumed_f_gcmd(c_string("setMDLengthUnits"), PLUMED_LENGTH_UNIT)
      call plumed_f_gcmd(c_string("setMDTimeUnits"), PLUMED_TIME_UNIT)
      call plumed_f_gcmd(c_string("setMDMassUnits"), PLUMED_MASS_UNIT)
      ! TODO: Fix this, see comment above
      !call plumed_f_gcmd(c_string("setKbT"), plumed_kbt)
      ! This is just a label.
      call plumed_f_gcmd(c_string("setMDEngine"), c_string("ABIN"))
      ! Set PLUMED input file.
      call plumed_f_gcmd(c_string("setPlumedDat"), c_string(plumedfile))

      call plumed_f_gcmd(c_string("setNatoms"), natom)
      call plumed_f_gcmd(c_string("setTimestep"), dt0)
      call plumed_f_gcmd(c_string("setLogFile"), c_string(PLUMED_OUTPUT_FILE))

      if (irest == 1) then
         if (.not. silent) then
            write (*, *) "PLUMED RESTART ON"
         end if
         call plumed_f_gcmd(c_string("setRestart"), PLUMED_RESTART)
      end if

      ! This needs to be the last call
      call plumed_f_gcmd(c_string("init"), C_NULL_PTR)
      ! Make sure PLUMED flushes file buffers with a reasonable frequency
      ! in case we suddenly die. We set the frequency to every nrest steps
      ! to ensure we have consistent data for restarting the simulation.
      write (chint, *) nrest
      ! NOTE: 'readInputLine' cannot be called before 'init'
      ! Only available for API_VERSION > 3
      ! https://www.plumed.org/doc-v2.6/user-doc/html/_f_l_u_s_h.html
      call plumed_f_gcmd(c_string('readInputLine'), 'FLUSH STRIDE='//c_string(chint))
   end subroutine plumed_init

   subroutine finalize_plumed()
      integer, parameter :: PLUMED_G_INITIALIZED = 1
      integer :: initialized

      call plumed_f_ginitialized(initialized)

      if (initialized == PLUMED_G_INITIALIZED) then
         call plumed_f_gfinalize()
      end if
   end subroutine finalize_plumed

   subroutine force_plumed(x, y, z, fx, fy, fz, eclas)
      use mod_general, only: natom, it, nwalk
      use mod_system, only: am
      use mod_utils, only: c_string
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(in) :: eclas
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP) :: xx(size(x, 1)), yy(size(y, 1)), zz(size(z, 1))
      real(DP) :: fxx(size(fx, 1)), fyy(size(fy, 1)), fzz(size(fz, 1))
      ! Currently not implemented
      !real(DP) :: pbox(3,3), pcharges(natom)
      real(DP) :: plumvirial(3, 3)
      integer :: iat, iw

      ! Apply PLUMED only to centroid/single bead in PIMD ("post-quantization restraint")
      ! https://aip.scitation.org/doi/abs/10.1063/1.4986915
      ! TODO: We should probably check that we're using normal modes
      iw = 1

      do iat = 1, natom
         xx(iat) = x(iat, iw)
         yy(iat) = y(iat, iw)
         zz(iat) = z(iat, iw)
      end do

      ! TODO: Is this correct?
      ! Passing zero, i.e. no info about current forces acting on the system,
      ! might influence some CV usage
      fxx = 0.0D0
      fyy = 0.0D0
      fzz = 0.0D0

      ! It seems that this must be the first call here
      call plumed_f_gcmd(c_string("setStep"), it)

      ! Do not pass charges at all since ABIN cannot do this now.
      ! PLUMED should die if we require CVs based on charges.
      ! Same goes for PBC.
      ! TODO: Make ABIN read optional Mulliken charges and pass them here.
      !call plumed_f_gcmd("setCharges"//char(0),pcharges)
      !call plumed_f_gcmd("setBox"//char(0),pbox)
      ! However, virial is required.
      plumvirial = 0.0D0
      call plumed_f_gcmd(c_string("setVirial"), plumvirial)

      call plumed_f_gcmd(c_string("setPositionsX"), xx)
      call plumed_f_gcmd(c_string("setPositionsY"), yy)
      call plumed_f_gcmd(c_string("setPositionsZ"), zz)
      call plumed_f_gcmd(c_string("setMasses"), am)
      call plumed_f_gcmd(c_string("setEnergy"), eclas)
      call plumed_f_gcmd(c_string("setForcesX"), fxx)
      call plumed_f_gcmd(c_string("setForcesY"), fyy)
      call plumed_f_gcmd(c_string("setForcesZ"), fzz)

      call plumed_f_gcmd(c_string("calc"), C_NULL_PTR)

      ! TODO: Is this correct? Should these forces be scaled?
      ! Adding PLUMED forces to the current ones, nwalk scaling for PIMD
      do iat = 1, natom
         fx(iat, iw) = fx(iat, iw) + fxx(iat) * nwalk
         fy(iat, iw) = fy(iat, iw) + fyy(iat) * nwalk
         fz(iat, iw) = fz(iat, iw) + fzz(iat) * nwalk
      end do
   end subroutine force_plumed

! USE_PLUMED
#else
   ! Stubbed subroutines for compilation without PLUMED
   ! We use this approach to avoid needing to have '#ifdef USE_PLUMED'
   ! elsewhere in the codebase
   subroutine plumed_init()
      call not_compiled_with('PLUMED', 'plumed_init')
   end subroutine plumed_init

   subroutine finalize_plumed()
      ! This must be a no-op!
      ! This routine is called from finalize()
      ! Which is in turn called by abinerror()
      ! So if we called `not_compiled_with()` here,
      ! We'd enter an infinite loop!
   end subroutine finalize_plumed

   subroutine force_plumed(x, y, z, fx, fy, fz, eclas)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas

      ! Just to get rid of compiler warnings :-(
      fx = x; fy = y; fz = z; eclas = 0.0D0

      call not_compiled_with('PLUMED', 'force_plumed')
   end subroutine force_plumed
#endif

end module mod_plumed
