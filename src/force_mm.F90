! Simple subroutine for calculation of
! Coulomb na Lennard-Jones forces and energies

! Currently, this module does not really support QMMM,
! Instead we can only do pure MM with LJ and Coulomb potential
module mod_force_mm
   use mod_const, only: DP
   use mod_array_size, only: MAXTYPES
   use mod_utils, only: abinerror
   use mod_qmmm, only: natmm, natqm
   implicit none
   private
   public :: inames_init, ABr_init, force_LJ_Coulomb
   ! TODO: These should all be private
   public :: attypes, q, rmin, eps
   character(len=2) :: attypes(MAXTYPES)
   ! L-J Combination rules
   character(len=10), parameter :: LJcomb = 'LB' !no other option for now
   real(DP) :: q(MAXTYPES), rmin(MAXTYPES), eps(MAXTYPES)
   real(DP), allocatable :: AIJ(:, :), BIJ(:, :)
   save
contains
   subroutine inames_init()
      use mod_general, only: natom
      use mod_system, only: inames, names
      integer :: i, iat2, pom

      do i = 1, natom
         pom = 0
         do iat2 = 1, natom
            if (names(i) == attypes(iat2)) then
               inames(i) = iat2
               pom = 1
               exit
            end if
         end do
         if (pom == 0) then
            write (*, *) 'Atom name does not have atom type for qmmm parameters. Exiting....'
            call abinerror('inames_init')
         end if
      end do

   end subroutine inames_init

   subroutine ABr_init()
      use mod_const, only: ANG
      use mod_general, only: natom
      use mod_system, only: inames
      integer :: iat1, iat2, i1, i2
      real(DP) :: epsij, rij
      allocate (AIJ(natom, natom))
      allocate (BIJ(natom, natom))
      do iat1 = 1, natom
         do iat2 = 1, natom
            i1 = inames(iat1)
            i2 = inames(iat2)
            if (LJcomb == 'LB') then
               rij = 0.5 * (rmin(i1) + rmin(i2)) * ang
               epsij = sqrt(eps(i1) * eps(i2))
            end if
            BIJ(i1, i2) = 2 * 6 * epsij * rij**6
            AIJ(i1, i2) = 12 * epsij * rij**12
         end do
      end do
   end subroutine

   subroutine force_LJ_Coulomb(x, y, z, fx, fy, fz, eclas)
      use mod_array_size
      use mod_general
      use mod_system, only: inames
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer :: iw, iat1, iat2, i1, i2
      real(DP) :: r, kLJ, kC
      real(DP) :: ri, ri3, dx, dy, dz

      do iw = 1, nwalk
         do iat1 = 1, natom
            do iat2 = iat1 + 1, natom
               if (iat2 <= natqm) cycle
               dx = x(iat1, iw) - x(iat2, iw)
               dy = y(iat1, iw) - y(iat2, iw)
               dz = z(iat1, iw) - z(iat2, iw)
               r = dx**2 + dy**2 + dz**2
               ri = 1 / r
               ri3 = ri * ri * ri
               i1 = inames(iat1)
               i2 = inames(iat2)
               kLJ = ri3 * (ri3 * AIJ(i1, i2) - BIJ(i1, i2)) * ri
               kC = q(i1) * q(i2) * sqrt(ri3)
               fx(iat1, iw) = fx(iat1, iw) + (kLJ + kC) * dx
               fx(iat2, iw) = fx(iat2, iw) - (kLJ + kC) * dx
               fy(iat1, iw) = fy(iat1, iw) + (kLJ + kC) * dy
               fy(iat2, iw) = fy(iat2, iw) - (kLJ + kC) * dy
               fz(iat1, iw) = fz(iat1, iw) + (kLJ + kC) * dz
               fz(iat2, iw) = fz(iat2, iw) - (kLJ + kC) * dz
               eclas = eclas + ri3 * (ri3 * AIJ(i1, i2) / 12 - BIJ(i1, i2) / 6) / nwalk
               eclas = eclas + q(i1) * q(i2) / sqrt(r) / nwalk
            end do
         end do
      end do

   end subroutine force_LJ_Coulomb

end module mod_force_mm
