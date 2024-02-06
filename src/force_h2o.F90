! Interface to analytical H2O (single water molecule) potentials.
! This is invoked via pot='_h2o_', and the potential is selected via
! h2opot='schwenke' in namelist General in ABIN input file.
module mod_force_h2o
   use mod_const, only: DP
   use mod_files, only: stdout, stderr
   use mod_error, only: fatal_error
   use mod_interfaces, only: h2o_pot_schwenke
   implicit none
   private
   public :: initialize_h2o_pot, force_h2o, h2opot
   character(len=20) :: h2opot = 'schwenke'
   save
contains

   ! Perform some basic validation
   subroutine initialize_h2o_pot(natom, atom_names)
      integer, intent(in) :: natom
      character(len=2), intent(in) :: atom_names(natom)

      if (natom /= 3) then
         call fatal_error(__FILE__, __LINE__, "This is not a water molecule!")
      end if

      if (atom_names(1) /= 'O' .or. atom_names(2) /= 'H' .or. atom_names(3) /= 'H') then
         write (stderr, *) 'ERROR: Bad element type.'
         write (stderr, *) 'Water atoms must be ordered as "O H H"'
         call fatal_error(__FILE__, __LINE__, "This is not a water molecule!")
      end if
      ! TODO: Some PES initialization could be performed here
   end subroutine initialize_h2o_pot

   ! Compute internal coordinates from cartesians,
   ! OH distances in Bohrs, HOH angle in radians.
   subroutine get_internal_coords(x, y, z, iw, rOH1, rOH2, aHOH_rad)
      use mod_const, only: PI
      use mod_utils, only: get_distance, get_angle
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: iw
      real(DP), intent(out) :: rOH1, rOH2, aHOH_rad
      real(DP) :: aHOH_deg

      rOH1 = get_distance(x, y, z, 1, 2, iw)
      rOH2 = get_distance(x, y, z, 1, 3, iw)
      aHOH_deg = get_angle(x, y, z, 2, 1, 3, iw)
      aHOH_rad = aHOH_deg * PI / 180.0D0
   end subroutine get_internal_coords

   ! This is just a wrapper function that selects the H2O potential
   ! based on user input. For now only h2opot="schwenke" is supported
   subroutine force_h2o(x, y, z, fx, fy, fz, eclas, natom, nbeads)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: natom, nbeads

      ! TODO: Use function pointers to select the potential and pass it down
      if (h2opot == 'schwenke') then
         call force_h2o_schwenke(x, y, z, fx, fy, fz, eclas, natom, nbeads)
      else
         call fatal_error(__FILE__, __LINE__, 'Potential '//trim(h2opot)//' not implemented')
      end if
   end subroutine force_h2o

   subroutine force_h2o_schwenke(x, y, z, fx, fy, fz, Eclas, natom, nbeads)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: natom, nbeads
      ! Internal water coordinates
      real(DP) :: rOH1, rOH2, aHOH_rad
      real(DP) :: rij(nbeads, 3)
      real(DP) :: Epot(nbeads)
      integer :: iw

      ! The H2O potentials are evaluated in internal coordinates, but ABIN works in cartesians
      do iw = 1, nbeads
         call get_internal_coords(x, y, z, iw, rOH1, rOH2, aHOH_rad)
         rij(iw, 1) = rOH1
         rij(iw, 2) = rOH2
         rij(iw, 3) = aHOH_rad
      end do

      ! Calculated the potential for all PI beads at once
      ! The potential is implemented in h2o_schwenke.f
      call h2o_pot_schwenke(rij, Epot, nbeads)

      ! For Path Integrals, the final energy of the PI necklace
      ! is an average over all beads.
      ! TODO: Forces need to be appropriately scaled as well!
      do iw = 1, nbeads
         Eclas = Eclas + Epot(iw)
      end do
      Eclas = Eclas / nbeads

      call numerical_forces(x, y, z, fx, fy, fz, Epot, natom, nbeads)

   end subroutine force_h2o_schwenke

   ! TODO: Implement numerical forces generally for all potentials
   ! For now, they can be implemented here and hardcoded for a specific H2O potential

   subroutine numerical_forces(x, y, z, Epot, fx, fy, fz, natom, nbeads)
      real(DP) :: x(natom, nbeads)
      real(DP) :: y(natom, nbeads)
      real(DP) :: z(natom, nbeads)
      real(DP), intent(inout) :: fx(natom, nbeads)
      real(DP), intent(inout) :: fy(natom, nbeads)
      real(DP), intent(inout) :: fz(natom, nbeads)
      integer, intent(in) :: natom, nbeads
      
      ! This is the energy for the currrent geometry that has already been calculated
      real(DP), intent(in) :: Epot(nbeads)

      ! Internal water coordinates
      real(DP) :: rOH1, rOH2, aHOH_rad
      real(DP) :: rij(nbeads, 3)
      
      real(DP) :: Eclas_orig, Eclas_plus, Eclas_minus
      real(DP) :: delta = 0.00005_DP
      integer :: i, j, k, iw

      ! Calculate forces numerically using central differences
      do j = 1, nbeads
         ! Save the original energy
         Eclas_orig = Epot(j)
         do i = 1, natom
            do k = 1, 3 ! x, y, z
               ! Reset the energy
               Eclas_plus = 0.0_DP
               Eclas_minus = 0.0_DP

               ! Move the atom forwards
               select case (k)
                  case (1)
                     x(j, i) = x(j, i) + delta
                  case (2)
                     y(j, i) = y(j, i) + delta
                  case (3)
                     z(j, i) = z(j, i) + delta
               end select

               ! Calculate the energy for the forward perturbed geometry
               do iw = 1, nbeads
                  call get_internal_coords(x, y, z, iw, rOH1, rOH2, aHOH_rad)
                  rij(iw, 1) = rOH1
                  rij(iw, 2) = rOH2
                  rij(iw, 3) = aHOH_rad
               end do

               call h2o_pot_schwenke(rij, Epot, nbeads)

               do iw = 1, nbeads
                  Eclas_plus = Eclas_plus + Epot(iw)
               end do
               Eclas_plus = Eclas_plus / nbeads

               ! Move the atom backwards
               select case (k)
                  case (1)
                     x(j, i) = x(j, i) - 2.0 * delta
                  case (2)
                     y(j, i) = y(j, i) - 2.0 * delta
                  case (3)
                     z(j, i) = z(j, i) - 2.0 * delta
               end select
               
               ! Calculate the energy for the backward perturbed geometry
               do iw = 1, nbeads
                  call get_internal_coords(x, y, z, iw, rOH1, rOH2, aHOH_rad)
                  rij(iw, 1) = rOH1
                  rij(iw, 2) = rOH2
                  rij(iw, 3) = aHOH_rad
               end do

               call h2o_pot_schwenke(rij, Epot, nbeads)

               do iw = 1, nbeads
                  Eclas_plus = Eclas_plus + Epot(iw)
               end do
               Eclas_plus = Eclas_plus / nbeads

               ! Calculate the numerical force
               select case (k)
                  case (1)
                     fx(j, i) = (Eclas_plus - Eclas_minus) / (2.0 * delta)
                  case (2)
                     fy(j, i) = (Eclas_plus - Eclas_minus) / (2.0 * delta)
                  case (3)
                     fz(j, i) = (Eclas_plus - Eclas_minus) / (2.0 * delta)
               end select

               ! Restore coordinates
               select case (k)
                  case (1)
                     x(j, i) = x(j, i) + delta
                  case (2)
                     y(j, i) = y(j, i) + delta
                  case (3)
                     z(j, i) = z(j, i) + delta
               end select

            end do
         end do
      end do
   end subroutine numerical_forces

end module mod_force_h2o