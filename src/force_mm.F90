! A toy code for calculation of
! Coulomb and Lennard-Jones forces and energies

! Currently, this module does not support QMMM,
! Instead we can only do pure MM with LJ and Coulomb potential
! This is invoked via pot='_mm_'
module mod_force_mm
   use mod_const, only: DP
   use mod_files, only: stdout, stderr
   implicit none
   private
   public :: initialize_mm, force_LJ_Coulomb
   public :: attypes, q, rmin, eps

   ! Mapping of atomic types to atoms in the system
   integer, allocatable :: inames(:)

   ! User-specified atomic types read from the input file
   character(len=2), allocatable :: attypes(:)
   ! User-specified Lennard-Jones parameters
   real(DP), allocatable :: rmin(:), eps(:)
   ! User-specified Coulomb charges
   real(DP), allocatable :: q(:)

   ! L-J Combination rules
   character(len=10), parameter :: LJcomb = 'LB' !no other option for now
   ! Internal L-J parameters
   real(DP), allocatable :: AIJ(:, :), BIJ(:, :)
   save
contains

   subroutine initialize_mm(natom)
      use mod_utils, only: normalize_atom_name
      integer, intent(in) :: natom
      real(DP) :: q_total
      integer :: iat

      allocate (inames(natom))

      do iat = 1, size(attypes)
         if (attypes(iat) == '') exit
         attypes(iat) = normalize_atom_name(attypes(iat))
      end do

      call inames_init(natom)
      call ABr_init(natom)

      q_total = calc_total_charge(q, natom)
      write (stdout, '(A,F6.3,A)') 'Total system charge = ', q_total, ' a.u.'
      if (q_total /= 0.0D0) then
         write (stderr, *) 'WARNING: total charge is not zero!'
      end if
   end subroutine

   subroutine inames_init(natom)
      use mod_error, only: fatal_error
      use mod_system, only: names
      integer, intent(in) :: natom
      integer :: iat, iat2
      logical :: assigned

      do iat = 1, natom
         assigned = .false.
         do iat2 = 1, size(attypes)
            if (names(iat) == attypes(iat2)) then
               inames(iat) = iat2
               assigned = .true.
               exit
            end if
         end do
         if (.not. assigned) then
            call fatal_error(__FILE__, __LINE__,&
               & 'Atom name '//trim(names(iat))//' does not have corresponding atom type in MM parameters')
         end if
      end do
   end subroutine inames_init

   subroutine ABr_init(natom)
      use mod_const, only: ANG
      integer, intent(in) :: natom
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
               epsij = dsqrt(eps(i1) * eps(i2))
            end if
            BIJ(i1, i2) = 2 * 6 * epsij * rij**6
            AIJ(i1, i2) = 12 * epsij * rij**12
         end do
      end do
   end subroutine

   real(DP) function calc_total_charge(q, natom) result(charge)
      real(DP), intent(in) :: q(:)
      integer, intent(in) :: natom
      integer :: iat

      charge = 0.0D0
      do iat = 1, natom
         charge = charge + q(inames(iat))
      end do
   end function calc_total_charge

   subroutine force_LJ_Coulomb(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: natom
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw, iat1, iat2, i1, i2
      real(DP) :: r, kLJ, kC
      real(DP) :: ri, ri3, dx, dy, dz

      do iw = 1, walkmax
         do iat1 = 1, natom
            do iat2 = iat1 + 1, natom
               if (iat2 <= natqm) then
                  cycle
               end if
               dx = x(iat1, iw) - x(iat2, iw)
               dy = y(iat1, iw) - y(iat2, iw)
               dz = z(iat1, iw) - z(iat2, iw)
               r = dx**2 + dy**2 + dz**2
               ri = 1 / r
               ri3 = ri * ri * ri
               i1 = inames(iat1)
               i2 = inames(iat2)
               kLJ = ri3 * (ri3 * AIJ(i1, i2) - BIJ(i1, i2)) * ri
               kC = q(i1) * q(i2) * dsqrt(ri3)
               fx(iat1, iw) = fx(iat1, iw) + (kLJ + kC) * dx
               fx(iat2, iw) = fx(iat2, iw) - (kLJ + kC) * dx
               fy(iat1, iw) = fy(iat1, iw) + (kLJ + kC) * dy
               fy(iat2, iw) = fy(iat2, iw) - (kLJ + kC) * dy
               fz(iat1, iw) = fz(iat1, iw) + (kLJ + kC) * dz
               fz(iat2, iw) = fz(iat2, iw) - (kLJ + kC) * dz
               eclas = eclas + ri3 * (ri3 * AIJ(i1, i2) / 12 - BIJ(i1, i2) / 6) / walkmax
               eclas = eclas + q(i1) * q(i2) / dsqrt(r) / walkmax
            end do
         end do
      end do

   end subroutine force_LJ_Coulomb

end module mod_force_mm
