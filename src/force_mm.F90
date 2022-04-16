! A toy code for Coulomb and Lennard-Jones forces and energies

! Currently, this module does not support QMMM,
! only pure MM with LJ and Coulomb potentials (no other bonding terms).
! This is invoked via pot='_mm_'
module mod_force_mm
   use mod_const, only: DP
   use mod_files, only: stdout, stderr
   use mod_error, only: fatal_error
   implicit none
   private
   public :: initialize_mm, force_mm

   ! Lennard-Jones combination rules
   ! Only Lorentz-Berthelot rules are currently implemented
   ! https://en.wikipedia.org/wiki/Combining_rules#Lorentz-Berthelot_rules
   character(len=10), parameter :: LJcomb = 'LB'

   ! Mapping of atomic types to atoms in the system
   integer, allocatable :: inames(:)

   ! User-specified Coulomb charges
   real(DP), allocatable :: q(:)

   ! Internal two-particle L-J parameters
   real(DP), allocatable :: Aij(:, :), Bij(:, :)
   save
contains

   subroutine initialize_mm(natom, names, atom_types, charges, LJ_rmin, LJ_eps)
      integer, intent(in) :: natom
      character(len=2), intent(in) :: names(:), atom_types(:)
      real(DP), intent(in) :: LJ_rmin(:), LJ_eps(:)
      real(DP), intent(in) :: charges(:)
      character(len=2), allocatable :: normalized_atom_types(:)
      integer :: num_types

      write (stdout, '(A)') ''
      write (stdout, '(A)') 'Initializing MM potential'

      normalized_atom_types = normalize_atom_types(atom_types)
      num_types = size(normalized_atom_types)
      write (stdout, '(A,I0,A)') 'Found ', num_types, ' user-defined atom types'

      call initialize_inames(names, normalized_atom_types, natom)

      call initialize_charges(charges, num_types, natom)

      call initialize_LJ(LJ_rmin, LJ_eps, num_types, natom)

      call flush(stdout)
      call flush(stderr)
      write (stdout, '(A)') ''
   end subroutine initialize_mm

   integer function count_atom_types(atom_types) result(num_types)
      character(len=2), intent(in) :: atom_types(:)
      integer :: iat

      num_types = 0
      do iat = 1, size(atom_types)
         if (atom_types(iat) == '') cycle
         num_types = num_types + 1
      end do
   end function count_atom_types

   function normalize_atom_types(atom_types) result(attypes)
      use mod_utils, only: normalize_atom_name
      character(len=2), intent(in) :: atom_types(:)
      character(len=2), allocatable :: attypes(:)
      integer :: num_types
      integer :: iat

      num_types = count_atom_types(atom_types)
      allocate(attypes(num_types))
      attypes = ''

      do iat = 1, size(atom_types)
         if (atom_types(iat) == '') cycle
         attypes(iat) = normalize_atom_name(atom_types(iat))
      end do
   end function normalize_atom_types

   ! Here we create the mapping from atom index
   ! to corresponding atom type, stored in array inames.
   ! Example: Water molecule
   ! names = (/ 'O', 'H', 'H' /)
   ! attypes = (/ 'O', 'H' /)
   ! inames = (/ 1, 2, 2 /)
   ! q = (/ -0.84, 0.42 /)
   subroutine initialize_inames(names, attypes, natom)
      character(len=2), intent(in) :: names(:)
      character(len=2), intent(in) :: attypes(:)
      integer, intent(in) :: natom
      integer :: iat, iat2
      logical :: assigned

      allocate (inames(natom))
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
   end subroutine initialize_inames

   subroutine initialize_charges(charges, num_types, natom)
      real(DP), intent(in) :: charges(:)
      integer, intent(in) :: num_types, natom
      real(DP) :: q_total

      allocate(q(num_types))
      q(:) = charges(:num_types)

      q_total = calc_total_charge(q, natom)
      write (stdout, '(A,F6.3,A)') 'Total system charge = ', q_total, ' a.u.'
      if (q_total /= 0.0D0) then
         write (stderr, *) 'WARNING: total charge is not zero!'
      end if
   end subroutine initialize_charges

   real(DP) function calc_total_charge(q, natom) result(charge)
      real(DP), intent(in) :: q(:)
      integer, intent(in) :: natom
      integer :: iat

      charge = 0.0D0
      do iat = 1, natom
         charge = charge + q(inames(iat))
      end do
   end function calc_total_charge

   subroutine initialize_LJ(rmin, eps, num_types, natom)
      use mod_const, only: ANG
      real(DP), intent(in) :: rmin(:), eps(:)
      integer, intent(in) :: num_types, natom
      integer :: iat1, iat2, i, j
      real(DP) :: epsij, rij

      ! Aij and Bij are two-particle combined L-J parameters
      ! https://en.wikipedia.org/wiki/Combining_rules
      allocate (Aij(num_types, num_types))
      allocate (Bij(num_types, num_types))
      do iat1 = 1, natom
         do iat2 = 1, natom
            i = inames(iat1)
            j = inames(iat2)
            if (LJcomb == 'LB') then
               ! WARNING: We expect rmin in angstroms!
               rij = 0.5 * (rmin(i) + rmin(j)) * ANG
               epsij = dsqrt(eps(i) * eps(j))
            end if
            Bij(i, j) = 2 * 6 * epsij * rij**6
            Aij(i, j) = 12 * epsij * rij**12
         end do
      end do
   end subroutine initialize_LJ

   subroutine force_mm(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: natom
      use mod_qmmm, only: natqm
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      integer :: iw, iat1, iat2, i, j
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
               i = inames(iat1)
               j = inames(iat2)
               kLJ = ri3 * (ri3 * Aij(i, j) - Bij(i, j)) * ri
               kC = q(i) * q(j) * dsqrt(ri3)
               fx(iat1, iw) = fx(iat1, iw) + (kLJ + kC) * dx
               fx(iat2, iw) = fx(iat2, iw) - (kLJ + kC) * dx
               fy(iat1, iw) = fy(iat1, iw) + (kLJ + kC) * dy
               fy(iat2, iw) = fy(iat2, iw) - (kLJ + kC) * dy
               fz(iat1, iw) = fz(iat1, iw) + (kLJ + kC) * dz
               fz(iat2, iw) = fz(iat2, iw) - (kLJ + kC) * dz
               eclas = eclas + ri3 * (ri3 * Aij(i, j) / 12 - Bij(i, j) / 6) / walkmax
               eclas = eclas + q(i) * q(j) / dsqrt(r) / walkmax
            end do
         end do
      end do

      ! TODO: Divide by nwalkmax at the end

   end subroutine force_mm

end module mod_force_mm
