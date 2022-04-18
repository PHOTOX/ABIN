! This module contains both massive and global version
! of the Nosé-Hoover chain (NHC) thermostat.
module mod_nhc
   use mod_const, only: DP
   implicit none
   private
   public :: tau0, inose, nchain, temp
   public :: nrespnose, nyosh
   public :: readNHC
   public :: imasst, nmolt, natmolt, nshakemol
   public :: nhc_init, finalize_nhc
   public :: get_nhcham, get_nhc_temp
   public :: shiftNHC_yosh_mass, shiftNHC_yosh
   public :: nhc_restout, nhc_restin

   ! Temperature (read in Kelvins and converted to a.u. in init.F90)
   real(DP) :: temp = 0.0D0
   ! Thermostat relaxation time in picoseconds
   real(DP) :: tau0 = -1

   ! Main switch for different thermostats
   ! inose == 1 turns on Nose-Hoover
   integer :: inose = -1
   ! Number of Nosé=Hoover chains
   integer :: nchain = 4
   ! Parameters for the Suzuki-Yoshida RESPA integration
   integer :: nrespnose = 3, nyosh = 7
   ! switch for massive thermostatting
   ! imasst == 0 turns on global thermostat
   integer :: imasst = 1
   ! for global thermostat, we can choose to thermostat
   ! molecules individually.
   integer :: nmolt = 1
   integer, allocatable :: natmolt(:)
   integer, allocatable :: nshakemol(:)

   ! Whether to read NHC params from the restart file
   ! Useful e.g. if we want to restart NVE simulation
   ! and turn on the thermostat.
   logical :: readNHC = .true.

   ! Internal variables and arrays
   integer, parameter :: MAXCHAIN = 10
   ! Auxiliary thermostat momenta
   real(DP), allocatable :: pnhx(:, :, :), pnhy(:, :, :), pnhz(:, :, :)
   real(DP), allocatable :: xi_x(:, :, :), xi_y(:, :, :), xi_z(:, :, :)
   ! Suzuki-Yoshida weights
   real(DP), allocatable :: w(:)
   ! Thermostat masses for global thermostat
   real(DP), allocatable :: ms(:, :)
   ! Thermostat masses for massive thermostat
   real(DP), allocatable :: Qm(:)
   save
contains

   ! Calculate conserved quantity of the NHC thermostat
   real(DP) function get_nhcham(natom, nwalk) result(nhcham)
      use mod_system, only: dime
      integer, intent(in) :: natom, nwalk
      integer iat, inh, iw

      nhcham = 0.0D0
      if (imasst == 1) then

         do inh = 1, nchain
            do iw = 1, nwalk
               do iat = 1, natom
                  nhcham = nhcham + &
                         & pnhx(iat, iw, inh) * pnhx(iat, iw, inh) * 0.5 / Qm(iw) + &
                         & temp * xi_x(iat, iw, inh)
                  if (dime > 1) nhcham = nhcham + &
                                       & pnhy(iat, iw, inh) * pnhy(iat, iw, inh) * 0.5 / Qm(iw) + &
                                       & temp * xi_y(iat, iw, inh)
                  if (dime > 2) nhcham = nhcham + &
                                       & pnhz(iat, iw, inh) * pnhz(iat, iw, inh) * 0.5 / Qm(iw) + &
                                       & temp * xi_z(iat, iw, inh)
               end do
            end do
         end do

      else

         iw = 1 !TODO: az bude shake+pimd,tak je tohle treba vyresit
         do iat = 1, nmolt
            nhcham = nhcham + &
                     & 0.5D0 * pnhx(iat, iw, 1) * pnhx(iat, iw, 1) / ms(iat, 1) + &
                     & dime * natmolt(iat) * temp * xi_x(iat, iw, 1)
            do inh = 2, nchain
               nhcham = nhcham + &
                        pnhx(iat, iw, inh) * pnhx(iat, iw, inh) * 0.5 / ms(iat, inh) + &
                        temp * xi_x(iat, iw, inh)
            end do
         end do

      end if
   end function get_nhcham

   subroutine check_nhc_parameters()
      use mod_utils, only: int_positive, int_nonnegative, int_switch, real_nonnegative
      use mod_files, only: stdout, stderr
      use mod_general, only: natom, ipimd
      use mod_const, only: AUTOK
      use mod_error, only: fatal_error
      use mod_chars, only: chknow
      use mod_general, only: iknow
      integer :: imol, ipom, iat
      logical :: error

      call int_nonnegative(inose, 'inose')
      call int_positive(nchain, 'nchain')
      call int_positive(nrespnose, 'nrespnose')
      call int_positive(nmolt, 'nmolt')
      call int_switch(imasst, 'imasst')
      call real_nonnegative(temp, 'temp')

      error = .false.
      if (nchain > maxchain) then
         write (stderr, *) 'Maximum number of Nose-Hoover chains exceeded'
         error = .true.
      end if
      if (nrespnose < 3) then
         write (stderr, *) 'Variable nrespnose < 3! Assuming this is an error in input and exiting.'
         write (stderr, *) 'Such low value would probably not produce stable results.'
         write (stderr, *) chknow
         if (iknow /= 1) error = .true.
      end if
      if (nyosh == 1) then
         write (stderr, *) 'Use Suzuki-Yoshida scheme for Nose-Hoover thermostat (nyosh=3 or 7).'
         write (stderr, *) chknow
         if (iknow /= 1) error = .true.
      end if
      if (imasst == 0 .and. ipimd == 1) then
         write (stderr, *) 'PIMD simulations must use massive thermostat (imasst=1)!'
         error = .true.
      end if
      if (nmolt > natom) then
         write (stdout, *) 'Input error: nmolt > natom, which is not possible. Consult the manual.'
         error = .true.
      end if
      if (imasst == 0) then
         do imol = 1, nmolt
            if (natmolt(imol) <= 0) then
               write (stderr, *) 'Number of atoms in molecules not specified! Set array natmolt properly.'
               error = .true.
            end if
         end do
      end if
      if (imasst == 0) then
         ipom = 0
         do iat = 1, nmolt
            ipom = ipom + natmolt(iat)
         end do
         if (ipom /= natom) then
            write (stderr, *) 'Number of atoms in thermostated molecules (natmolt) does not match natom.'
            write (stderr, *) chknow
            if (iknow /= 1) error = .true.
         end if
      end if
      if (temp * AUTOK < 1 .and. inose > 0) then
         write (stderr, *) 'Temperature below 1 Kelvin. Are you sure?'
         write (stderr, *) chknow
         if (iknow /= 1) error = .true.
      end if

      if (error) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Invalid NHC thermostat parameter')
      end if
   end subroutine check_nhc_parameters

   subroutine nhc_init()
      use mod_files, only: stdout
      if (imasst == 1) then
         write (stdout, *) 'Initializing massive Nosé-Hoover Chain thermostat'
      else
         write (stdout, *) 'Initializing global Nosé-Hoover Chain thermostat'
      end if

      call check_nhc_parameters()

      call set_yoshida_weights(nyosh)

      ! Not sure if we should have this default value
      if (tau0 < 0) then
         write (stdout, *) 'WARNING: tau0 not set.'
         write (stdout, *) 'Using default value tau0=0.001 picoseconds'
         tau0 = 0.001D0
      end if
      call set_nhc_masses(tau0)

      call initialize_nhc_momenta(temp)
   end subroutine nhc_init

   subroutine set_nhc_masses(tau_ps)
      use mod_const, only: AUtoFS
      real(DP), intent(in) :: tau_ps
      real(DP) :: tau_au, nhc_mass

      tau_au = tau_ps / AUtoFS * 1000

      ! See M. E. Tuckerman, Statistical mechanics, p.190
      nhc_mass = temp * tau_au * tau_au

      if (imasst == 0) then
         call set_nhc_global_masses(nhc_mass)
      else if (imasst == 1) then
         call set_nhc_massive_masses(nhc_mass)
      end if
   end subroutine set_nhc_masses

   subroutine set_nhc_massive_masses(nhc_mass)
      use mod_const, only: PI
      use mod_general, only: ipimd, nwalk, inormalmodes
      real(DP), intent(in) :: nhc_mass
      real(DP) :: omega
      integer :: iw

      allocate (Qm(nwalk))
      Qm = 0.0D0
      do iw = 1, nwalk
         Qm(iw) = nhc_mass
      end do

      ! For PIMD the Nose-Hoover mass is set within the code as
      ! 1/(beta*omega_p^2)
      ! where
      ! omega_p=sqrt(P)/(beta*hbar)
      if (ipimd == 1 .and. inormalmodes /= 1) then

         ! TODO: Is this optimal both with/wo staging transform?
         Qm(1) = nhc_mass ! see tuckermann,stat.mech.
         do iw = 2, nwalk
            Qm(iw) = 1 / (temp * nwalk)
         end do

      else if (ipimd == 1 .and. inormalmodes == 1) then

         ! so far, NHC with normal modes does not work
         temp = temp * nwalk
         Qm(1) = nhc_mass * 4
         do iw = 2, nwalk
            omega = 2 * temp * sin((iw - 1) * PI / nwalk)
            Qm(iw) = 1 / temp / omega**2
         end do

      end if
   end subroutine

   subroutine set_nhc_global_masses(mass)
      use mod_system, only: dime
      real(DP), intent(in) :: mass
      integer :: imol, inh

      allocate (ms(nmolt, nchain))
      ms = 0.0D0

      do imol = 1, nmolt
         ms(imol, 1) = (dime * natmolt(imol) - nshakemol(imol)) * mass
         do inh = 2, nchain
            ms(imol, inh) = mass
         end do
      end do
   end subroutine

   subroutine initialize_nhc_momenta(temp)
      use mod_general, only: natom, nwalk
      use mod_random, only: gautrg
      real(DP), intent(in) :: temp
      real(DP), allocatable :: ran(:)
      integer :: inh, iw, iat, ipom, imol

      if (imasst == 1) then
         allocate (pnhx(natom, nwalk, nchain)); pnhx = 0.0D0
         allocate (pnhy(natom, nwalk, nchain)); pnhy = 0.0D0
         allocate (pnhz(natom, nwalk, nchain)); pnhz = 0.0D0
         allocate (xi_x(natom, nwalk, nchain)); xi_x = 0.0D0
         allocate (xi_y(natom, nwalk, nchain)); xi_y = 0.0D0
         allocate (xi_z(natom, nwalk, nchain)); xi_z = 0.0D0
      else
         allocate (pnhx(nmolt, nwalk, nchain)); pnhx = 0.0D0
         allocate (xi_x(nmolt, nwalk, nchain)); xi_x = 0.0D0
      end if

      allocate (ran(natom * 3))

      if (imasst == 1) then
         do inh = 1, nchain
            do iw = 1, nwalk
               call gautrg(ran, natom * 3)
               ipom = 1
               do iat = 1, natom
                  pnhx(iat, iw, inh) = ran(ipom) * dsqrt(temp * Qm(iw))
                  pnhy(iat, iw, inh) = ran(ipom + 1) * dsqrt(temp * Qm(iw))
                  pnhz(iat, iw, inh) = ran(ipom + 2) * dsqrt(temp * Qm(iw))
               end do
               ipom = ipom + 3
            end do
         end do

      else if (imasst == 0) then

         do inh = 1, nchain
            do iw = 1, nwalk
               ! +1 if nmolt=1, gautrg needs array at least of length=2
               call gautrg(ran, nmolt + 1)
               do imol = 1, nmolt
                  pnhx(imol, iw, inh) = ran(imol) * dsqrt(temp * ms(imol, inh))
               end do
            end do
         end do

      end if

      deallocate (ran)
   end subroutine initialize_nhc_momenta

   ! Set weights for Suzuki-Yoshida integrator,
   ! depending on the integration order nyosh
   subroutine set_yoshida_weights(nyosh)
      use mod_error, only: fatal_error
      integer, intent(in) :: nyosh

      allocate (w(nyosh))

      select case (nyosh)
      case (1)
         w(1) = 1
      case (3)
         w(1) = 1.0D0 / (2.0D0 - 2**(1.0D0 / 3.0D0))
         w(3) = w(1)
         w(2) = 1 - w(1) - w(3)
      case (7)
         w(1) = 0.784513610477560_DP
         w(7) = w(1)
         w(2) = 0.235573213359357_DP
         w(6) = w(2)
         w(3) = -1.17767998417887_DP
         w(5) = w(3)
         w(4) = 1 - w(1) - w(2) - w(3) - w(5) - w(6) - w(7)
      case default
         call fatal_error(__FILE__, __LINE__, &
            & 'Invalid nyosh parameter. Allowed values are 1, 3 or 7')
      end select
   end subroutine set_yoshida_weights

   ! Read NHC momenta from the restart file
   subroutine nhc_restin(u)
      use mod_general, only: natom, nwalk
      ! Restart file unit
      integer, intent(in) :: u
      integer :: iat, iw, inh

      if (imasst == 1) then
         do inh = 1, nchain
            do iw = 1, nwalk
               do iat = 1, natom
                  read (u, *) pnhx(iat, iw, inh), pnhy(iat, iw, inh), pnhz(iat, iw, inh)
               end do
            end do
         end do

      else

         do inh = 1, nchain
            do iw = 1, nwalk
               do iat = 1, nmolt
                  read (u, *) pnhx(iat, iw, inh)
               end do
            end do
         end do

      end if
   end subroutine nhc_restin

   ! Write NHC momenta to the restart file
   subroutine nhc_restout(u)
      use mod_general, only: natom, nwalk
      ! Restart file unit
      integer, intent(in) :: u
      integer :: iat, iw, inh
      if (imasst == 1) then

         do inh = 1, nchain
            do iw = 1, nwalk
               do iat = 1, natom
                  write (u, '(3ES25.16E3)') pnhx(iat, iw, inh), pnhy(iat, iw, inh), pnhz(iat, iw, inh)
               end do
            end do
         end do

      else

         do inh = 1, nchain
            do iw = 1, nwalk
               do iat = 1, nmolt
                  write (u, '(ES24.16E3)') pnhx(iat, iw, inh)
               end do
            end do
         end do
      end if
   end subroutine nhc_restout

   subroutine finalize_nhc()
      if (allocated(w)) deallocate (w)
      if (allocated(Qm)) deallocate (Qm)
      if (allocated(ms)) deallocate (ms)
      if (allocated(pnhx)) deallocate (pnhx)
      if (allocated(pnhy)) deallocate (pnhy)
      if (allocated(pnhz)) deallocate (pnhz)
      if (allocated(xi_x)) deallocate (xi_x)
      if (allocated(xi_y)) deallocate (xi_y)
      if (allocated(xi_z)) deallocate (xi_z)
   end subroutine finalize_nhc

   ! Calculate the temperature of thermostat auxiliary momenta
   ! Currently not in use, works only for massive thermostat
   real(DP) function get_nhc_temp(natom, nwalk) result(nhc_temp)
      integer, intent(in) :: natom, nwalk
      real(DP) :: ekin
      integer :: iw, iat

      ekin = 0.0D0
      do iw = 1, nwalk
         do iat = 1, natom
            ekin = 0.5D0 / Qm(iw) * (pnhx(iat, iw, 1)**2 + &
                                & pnhy(iat, iw, 1)**2 + &
                                & pnhz(iat, iw, 1)**2)
         end do
      end do
      nhc_temp = 2 * ekin / 3 / natom / nwalk
   end function get_nhc_temp

   ! Suzuki-Yoshida split-operator integrator for global NHC
   subroutine shiftNHC_yosh(px, py, pz, amt, dt)
      use mod_array_size
      use mod_general
      use mod_system, only: dime
      real(DP) :: px(:, :), py(:, :), pz(:, :)
      real(DP) :: amt(:, :), G(MAXCHAIN)
      real(DP) :: dt, ekin2, AA
      real(DP) :: wdt, wdt2, wdt4, pscale
      integer :: iw, iat, inh
      integer :: nf, iresp, iyosh
      integer :: iat1, iat2, sumat, imol

      iw = 1 !only for CMD or centroid variable
      sumat = 0
      do imol = 1, nmolt
         iat1 = sumat + 1 !first atom in molecule
         sumat = sumat + natmolt(imol)
         iat2 = sumat
         nf = dime * natmolt(imol) - nshakemol(imol) !degrees of freedom
         pscale = 1.0D0

         ekin2 = 0.0D0
         do iat = iat1, iat2
            ekin2 = ekin2 + (px(iat, iw)**2 + py(iat, iw)**2 + pz(iat, iw)**2) / amt(iat, iw)
         end do

         G(1) = ekin2 - nf * temp
         do inh = 2, nchain
            G(inh) = pnhx(imol, iw, inh - 1) * pnhx(imol, iw, inh - 1) / ms(imol, inh - 1) - temp
         end do

         do iresp = 1, nrespnose
            do iyosh = 1, nyosh
               wdt = w(iyosh) * dt / nrespnose
               wdt2 = wdt / 2
               wdt4 = wdt2 / 2

               pnhx(imol, iw, nchain) = pnhx(imol, iw, nchain) + G(nchain) * wdt2
               ! UPDATE THERMOSTAT VELOCITIES
               do inh = 1, nchain - 1
                  AA = dexp(-wdt4 * pnhx(imol, iw, nchain - inh + 1) / ms(imol, nchain - inh + 1))
                  pnhx(imol, iw, nchain - inh) = pnhx(imol, iw, nchain - inh) * AA * AA + &
                                               & wdt2 * G(nchain - inh) * AA
               end do

               ! UPDATE PARTICLE VELOCITIES
               AA = dexp(-wdt * pnhx(imol, iw, 1) / ms(imol, 1))
               pscale = pscale * AA
               ! UPDATE FORCES
               G(1) = pscale * pscale * ekin2 - nf * temp

               ! UPDATE THERMOSTAT POSITIONS
               do inh = 1, nchain
                  xi_x(imol, iw, inh) = xi_x(imol, iw, inh) + &
                                      & pnhx(imol, iw, inh) / ms(imol, inh) * wdt
               end do

               ! UPDATE THERMOSTAT VELOCITIES
               do inh = 1, nchain - 1
                  AA = dexp(-wdt4 * pnhx(imol, iw, inh + 1) / ms(imol, inh + 1))
                  pnhx(imol, iw, inh) = pnhx(imol, iw, inh) * AA * AA + wdt2 * G(inh) * AA
                  G(inh + 1) = pnhx(imol, iw, inh) * pnhx(imol, iw, inh) / ms(imol, inh) - temp
               end do

               pnhx(imol, iw, nchain) = pnhx(imol, iw, nchain) + G(nchain) * wdt2
               ! nyosh enddo
            end do
            ! iresp enddo
         end do

         ! Update particle velocities
         do iat = iat1, iat2
            px(iat, iw) = px(iat, iw) * pscale
            py(iat, iw) = py(iat, iw) * pscale
            pz(iat, iw) = pz(iat, iw) * pscale
         end do

         ! imol enddo
      end do

      return
   end subroutine shiftNHC_yosh

   ! Suzuki-Yoshida split-operator integrator for massive NHC
   subroutine shiftNHC_yosh_mass(px, py, pz, amt, dt)
      use mod_general
      use mod_shake, only: nshake
      real(DP) :: px(:, :), py(:, :), pz(:, :)
      real(DP) :: amt(:, :)
      real(DP) :: Gx(MAXCHAIN), Gy(MAXCHAIN), Gz(MAXCHAIN)
      real(DP) :: dt, AA
      real(DP) :: wdt, wdt2, wdt4
      integer :: iw, iat, inh, istart
      integer :: iresp, iyosh

      istart = 1 !will be different with normal modes and shake
      if (nshake > 0) istart = 2

      do iw = istart, nwalk
         do iat = 1, natom

            Gx(1) = px(iat, iw) * px(iat, iw) / amt(iat, iw) - temp
            Gy(1) = py(iat, iw) * py(iat, iw) / amt(iat, iw) - temp
            Gz(1) = pz(iat, iw) * pz(iat, iw) / amt(iat, iw) - temp

            do inh = 2, nchain
               Gx(inh) = pnhx(iat, iw, inh - 1)**2 / Qm(iw) - temp
               Gy(inh) = pnhy(iat, iw, inh - 1)**2 / Qm(iw) - temp
               Gz(inh) = pnhz(iat, iw, inh - 1)**2 / Qm(iw) - temp
            end do

            do iresp = 1, nrespnose
               do iyosh = 1, nyosh
                  wdt = w(iyosh) * dt / nrespnose
                  wdt2 = wdt / 2
                  wdt4 = wdt2 / 2

                  pnhx(iat, iw, nchain) = pnhx(iat, iw, nchain) + Gx(nchain) * wdt2
                  pnhy(iat, iw, nchain) = pnhy(iat, iw, nchain) + Gy(nchain) * wdt2
                  pnhz(iat, iw, nchain) = pnhz(iat, iw, nchain) + Gz(nchain) * wdt2

                  ! UPDATE THERMOSTAT VELOCITIES
                  do inh = 1, nchain - 1
                     AA = dexp(-wdt4 * pnhx(iat, iw, nchain - inh + 1) / Qm(iw))
                     pnhx(iat, iw, nchain - inh) = pnhx(iat, iw, nchain - inh) * AA * AA + &
                                                 & wdt2 * Gx(nchain - inh) * AA
                     AA = dexp(-wdt4 * pnhy(iat, iw, nchain - inh + 1) / Qm(iw))
                     pnhy(iat, iw, nchain - inh) = pnhy(iat, iw, nchain - inh) * AA * AA + &
                                                 & wdt2 * Gy(nchain - inh) * AA
                     AA = dexp(-wdt4 * pnhz(iat, iw, nchain - inh + 1) / Qm(iw))
                     pnhz(iat, iw, nchain - inh) = pnhz(iat, iw, nchain - inh) * AA * AA + &
                                                 & wdt2 * Gz(nchain - inh) * AA
                  end do

                  ! UPDATE PARTICLE VELOCITIES
                  AA = dexp(-wdt * pnhx(iat, iw, 1) / Qm(iw))
                  px(iat, iw) = px(iat, iw) * AA
                  AA = dexp(-wdt * pnhy(iat, iw, 1) / Qm(iw))
                  py(iat, iw) = py(iat, iw) * AA
                  AA = dexp(-wdt * pnhz(iat, iw, 1) / Qm(iw))
                  pz(iat, iw) = pz(iat, iw) * AA

                  ! UPDATE FORCES
                  Gx(1) = (px(iat, iw) * px(iat, iw) / amt(iat, iw) - temp)
                  Gy(1) = (py(iat, iw) * py(iat, iw) / amt(iat, iw) - temp)
                  Gz(1) = (pz(iat, iw) * pz(iat, iw) / amt(iat, iw) - temp)

                  ! UPDATE THERMOSTAT POSITIONS
                  do inh = 1, nchain
                     xi_x(iat, iw, inh) = xi_x(iat, iw, inh) + pnhx(iat, iw, inh) / Qm(iw) * wdt
                     xi_y(iat, iw, inh) = xi_y(iat, iw, inh) + pnhy(iat, iw, inh) / Qm(iw) * wdt
                     xi_z(iat, iw, inh) = xi_z(iat, iw, inh) + pnhz(iat, iw, inh) / Qm(iw) * wdt
                  end do

                  ! UPDATE THERMOSTAT VELOCITIES
                  do inh = 1, nchain - 1
                     AA = dexp(-wdt4 * pnhx(iat, iw, inh + 1) / Qm(iw))
                     pnhx(iat, iw, inh) = pnhx(iat, iw, inh) * AA * AA + wdt2 * Gx(inh) * AA
                     Gx(inh + 1) = pnhx(iat, iw, inh)**2 / Qm(iw) - temp
                     AA = dexp(-wdt4 * pnhy(iat, iw, inh + 1) / Qm(iw))
                     pnhy(iat, iw, inh) = pnhy(iat, iw, inh) * AA * AA + wdt2 * Gy(inh) * AA
                     Gy(inh + 1) = pnhy(iat, iw, inh)**2 / Qm(iw) - temp
                     AA = dexp(-wdt4 * pnhz(iat, iw, inh + 1) / Qm(iw))
                     pnhz(iat, iw, inh) = pnhz(iat, iw, inh) * AA * AA + wdt2 * Gz(inh) * AA
                     Gz(inh + 1) = pnhz(iat, iw, inh)**2 / Qm(iw) - temp
                  end do

                  pnhx(iat, iw, nchain) = pnhx(iat, iw, nchain) + Gx(nchain) * wdt2
                  pnhy(iat, iw, nchain) = pnhy(iat, iw, nchain) + Gy(nchain) * wdt2
                  pnhz(iat, iw, nchain) = pnhz(iat, iw, nchain) + Gz(nchain) * wdt2

               end do

            end do

            ! iat enddo
         end do
         ! iw enddo
      end do

      return
   end subroutine shiftNHC_yosh_mass

end module mod_nhc
