! This module contains both massive and global version of the Nosé-Hoover chain thermostat
! It also contains all necessary variables and initialization routines.
module mod_nhc
   use mod_const, only: DP
   use mod_array_size, only: maxchain
   implicit none
   public :: ams, tau0, nhcham, inose, nchain, temp, temp0
   public :: nrespnose, nyosh, scaleveloc, readNHC, initNHC
   public :: imasst, nmolt, natmolt, nshakemol

   real(DP) :: ams = -1, nhcham = 0.0D0, tau0 = -1, tau ! in picoseconds
   real(DP) :: temp0 = -1, temp = 0.0D0
   integer :: inose = -1, nchain = 4, initNHC = -1
   integer :: imasst = 1 ! switch for massive thermostatting
   integer :: nrespnose = 3, nyosh = 7, nmolt = 1
   integer :: scaleveloc, readNHC = 1
#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || __GNUC__ > 4
   integer, allocatable :: natmolt(:)
   integer, allocatable :: nshakemol(:)
#else
   integer, parameter :: MAXMOL = 3000
   integer :: natmolt(MAXMOL)
   integer :: nshakemol(MAXMOL)
#endif
   real(DP), allocatable :: pnhx(:, :, :), pnhy(:, :, :), pnhz(:, :, :)
   real(DP), allocatable :: xi_x(:, :, :), xi_y(:, :, :), xi_z(:, :, :)
   real(DP), allocatable :: w(:), ms(:, :), Qm(:)
   save
contains

   subroutine calc_nhcham()
      use mod_general, only: natom, nwalk
      use mod_system, only: dime
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
   end subroutine

   subroutine nhc_init()
      use mod_const, only: AMU, AUtoFS, PI
      use mod_general, only: my_rank, ipimd, nwalk, natom, inormalmodes
      use mod_system, only: dime
      use mod_random, only: gautrg
      implicit none
      real(DP), allocatable :: ran(:)
      real(DP) :: omega
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

      allocate (w(nyosh))
      allocate (Qm(nwalk))
      allocate (ms(nmolt, nchain))
      allocate (ran(natom * 3))

      ! SETTING THERMOSTAT MASSES
      ams = ams * amu
      if (ams < 0 .and. tau0 < 0) then
         if (my_rank == 0) then
            write (*, *) 'Warning. Ams and tau0 not set.'
            write (*, *) 'Using default value tau0=0.001'
         end if
         tau = 0.001D0 / AUtoFS * 1000
      else
         tau = tau0 / AUtoFS * 1000
      end if
      ! pokud neni ams v inputu, tak ji pro clasickou simulaci priradime take
      ! automaticky, viz tuckermann, Statistical mechanics p.190
      if (ams < 0) then
         ams = temp * tau * tau
      end if
      do iw = 1, nwalk
         Qm(iw) = ams
      end do
      do imol = 1, nmolt
         ms(imol, 1) = (dime * natmolt(imol) - nshakemol(imol)) * ams
         do inh = 2, nchain
            ms(imol, inh) = ams
         end do
      end do
      ! in pimd the Nose-Hoover mass is set within the code as
      ! 1/(beta*omega_p^2)
      ! where
      ! omega_p=sqrt(P)/(beta*hbar)
      ! DH DEBUG
      if (ipimd == 1 .and. inormalmodes /= 1) then

         do iw = 1, nwalk
            Qm(iw) = 1 / (TEMP * NWALK)
         end do
         if (tau0 > 0) Qm(1) = ams ! see tuckermann,stat.mech.

      else if (ipimd == 1 .and. inormalmodes == 1) then

         ! so far, NHC with normal modes does not work
         temp = temp * nwalk
         Qm(1) = temp * tau * tau * 4
         do iw = 2, nwalk
            omega = 2 * TEMP * sin((iw - 1) * PI / NWALK)
            Qm(iw) = 1 / TEMP / omega**2
         end do

      end if

      if (nmolt <= 50 .and. my_rank == 0) then
         write (*, *) 'Thermostat masses'
         if (imasst == 1) write (*, *) (Qm(iw), iw=1, nwalk)
         if (imasst == 0) then
            do imol = 1, nmolt
               write (*, *) (ms(imol, inh), inh=1, nchain)
            end do
         end if
      end if

      ! Initialize thermostat momenta
      if (initNHC == 1 .and. imasst == 1) then
         do inh = 1, nchain
            do iw = 1, nwalk
               call gautrg(ran, natom * 3)
               ipom = 1
               do iat = 1, natom
                  pnhx(iat, iw, inh) = ran(ipom) * sqrt(temp * Qm(iw))
                  pnhy(iat, iw, inh) = ran(ipom + 1) * sqrt(temp * Qm(iw))
                  pnhz(iat, iw, inh) = ran(ipom + 2) * sqrt(temp * Qm(iw))
               end do
               ipom = ipom + 3
            end do
         end do

      else if (initNHC == 1 .and. imasst == 0) then

         do inh = 1, nchain
            do iw = 1, nwalk
               ! +1 if nmolt=1, gautrg needs array at least of length=2
               call gautrg(ran, nmolt + 1)
               do imol = 1, nmolt
                  pnhx(imol, iw, inh) = ran(imol) * sqrt(temp * ms(imol, inh))
               end do
            end do
         end do

      end if

      ! SET SUZUKI-YOSHIDA WEIGHTS
      if (nyosh == 3) then
         w(1) = 1.0D0 / (2.0D0 - 2**(1.0D0 / 3.0D0))
         w(3) = w(1)
         w(2) = 1 - w(1) - w(3)
      else if (nyosh == 1) then
         w(1) = 1
      else if (nyosh == 7) then
         w(1) = 0.784513610477560_DP
         w(7) = w(1)
         w(2) = 0.235573213359357_DP
         w(6) = w(2)
         w(3) = -1.17767998417887_DP
         w(5) = w(3)
         w(4) = 1 - w(1) - w(2) - w(3) - w(5) - w(6) - w(7)
      end if

      ! open(100,file='temper_nhc')
      ! close(100,status='delete')

      deallocate (ran)
   end subroutine

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

   subroutine nhc_temp() !currently not in use
      use mod_general, only: nwalk, natom
      implicit none
      integer :: iw, iat
      real(DP) :: ekin_mom = 0.0D0, temp1 = 0.0D0
      do iw = 1, nwalk
         do iat = 1, natom
            temp1 = pnhx(iat, iw, 1)**2 + pnhy(iat, iw, 1)**2 + pnhz(iat, iw, 1)**2
            temp1 = 0.5 * temp1 / ams
            ekin_mom = ekin_mom + temp1
         end do
      end do
      open (100, file='temper_nhc.dat', access='append')
      write (100, *) 2 * ekin_mom / 3 / natom / nwalk
      close (100)
   end subroutine

   subroutine shiftNHC_yosh(px, py, pz, amt, dt)
      use mod_array_size
      use mod_general
      use mod_system, only: dime
      implicit none
      real(DP) :: px(:, :), py(:, :), pz(:, :)
      real(DP) :: amt(:, :), G(maxchain)
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

   subroutine shiftNHC_yosh_mass(px, py, pz, amt, dt)
      use mod_array_size, only: MAXCHAIN
      use mod_general
      use mod_shake, only: nshake
      implicit none
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
