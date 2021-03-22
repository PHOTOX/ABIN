! Harmonic Spherical Potential
! (Spherical boundary conditions - SBC)
module mod_sbc
   use mod_const, only: DP, AUtoM, ME, PI
   implicit none
   private
   public :: isbc, rb_sbc, kb_sbc, rho
   public :: force_sbc, sbc_init
   integer :: isbc = 0
   real(DP) :: rb_sbc = -1, kb_sbc = 0.02D0, rho = -1
   real(DP), parameter :: fact = autom**3 / me
   real(DP) :: mass_total = 0.0D0
   integer, parameter :: ibag = 0 !elastic bag
   save
contains

   subroutine sbc_init(x, y, z)
      use mod_const, only: ANG
      use mod_general, only: natom !,nwalk
      use mod_system, only: am, names
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: r, rmax, xcm, ycm, zcm
      integer :: iw, iat

      mass_total = 0.0D0

      do iat = 1, natom
         mass_total = mass_total + am(iat)
      end do

      rmax = 0.0D0

      rb_sbc = rb_sbc * ANG

      ! calculation of Center-Of-Mass
      iw = 1
      call calc_com(x, y, z, iw, xcm, ycm, zcm)

      write (*, *) 'Center of mass is at [A]:', xcm / ang, ycm / ang, zcm / ang

      ! calculation of the cluster radius, taken from the center of mass
      do iat = 1, natom
         r = (x(iat, iw) - xcm)**2 + (y(iat, iw) - ycm)**2 + (z(iat, iw) - zcm)**2
         r = sqrt(r)
         if (r > rmax .and. names(iat) /= 'H') then
            rmax = r
         end if
      end do

      ! SETTING SIZE OF THE CLUSTER
      write (*, *) 'Cluster radius is r= ', rmax / ang, '  Angstrom'
      if (rb_sbc <= 0 .and. rho <= 0) then
         write (*, *) 'Cluster radius for spherical boundary conditions not specified.'
         if (rho < 0) then
            write (*, *) 'Will be set equal the initial radius of the system.'
            rb_sbc = rmax
         end if
      end if
      ! DETERMINATION of cluster size from given density
      if (rho > 0) then
         write (*, *) 'Calculating cluster radius from given densty.'
         rho = rho * fact !conversion from g/L to atomic units
         rb_sbc = mass_total / rho * 3 / 4 / PI
         rb_sbc = rb_sbc**(1 / 3.)
      end if

      if (rmax > rb_sbc) then
         write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
         write (*, *) 'WARNING: Cluster radius is bigger than rb_sbc.'
         write (*, *) 'Is this really what you want?'
         write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
      end if

      write (*, *) 'rb_sbc[A]=', rb_sbc / ang

   end subroutine

   subroutine force_sbc(x, y, z, fx, fy, fz)
      use mod_const, only: ANG
      use mod_general, only: natom, nwalk, it, nwrite
      use mod_system, only: names
      use mod_files, only: URADIUS
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP) :: r, frb, rmax, xcm, ycm, zcm
      integer :: iat, iw

      rmax = 0.0D0

      iw = 1
      ! Center-Of-Mass
      call calc_com(x, y, z, iw, xcm, ycm, zcm)

      ! write(*,*)'COM: ',xcm,ycm,zcm

      ! kb_sbc and rb_sbc must be specified in input.in
      do iw = 1, nwalk
         do iat = 1, natom
            r = (x(iat, iw) - xcm)**2 + (y(iat, iw) - ycm)**2 + (z(iat, iw) - zcm)**2
            r = sqrt(r)
            if (r > rb_sbc .and. names(iat) /= 'H') then
               frb = -kb_sbc * (r - rb_sbc)
               fx(iat, iw) = fx(iat, iw) + frb * (x(iat, iw) - xcm) / (r)
               fy(iat, iw) = fy(iat, iw) + frb * (y(iat, iw) - ycm) / (r)
               fz(iat, iw) = fz(iat, iw) + frb * (z(iat, iw) - zcm) / (r)
            end if

            if (r > rmax .and. names(iat) /= 'H') then
               rmax = r
            end if

         end do
      end do

      if (modulo(it, nwrite) == 0) then
         write (URADIUS, '(I12,2F15.3)') it, rmax / ang, &
                                       &  mass_total / (4.0D0 / 3.0D0 * pi * rmax**3) / fact
      end if

      return
   end subroutine

   ! TODO: This could be a general purpose routine,
   ! move to utils
   subroutine calc_com(x, y, z, iw, xcm, ycm, zcm)
      use mod_general, only: natom
      use mod_system, only: am
      real(DP), intent(out) :: xcm, ycm, zcm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw
      integer :: iat

      xcm = 0.0D0
      ycm = 0.0D0
      zcm = 0.0D0

      do iat = 1, natom
         xcm = xcm + x(iat, iw) * am(iat)
         ycm = ycm + y(iat, iw) * am(iat)
         zcm = zcm + z(iat, iw) * am(iat)
      end do
      xcm = xcm / mass_total
      ycm = ycm / mass_total
      zcm = zcm / mass_total

   end subroutine calc_com

end module mod_sbc
