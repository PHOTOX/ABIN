! This is invoked via pot='_h2o_'
module mod_h2o_pes
   use mod_const, only: DP
   use mod_error, only: fatal_error
   implicit none
   private
   public :: force_h2o_schwenke
   save
contains

   subroutine force_h2o_schwenke(x, y, z, fx, fy, fz, eclas, natom, walkmax)
      use mod_utils, only: get_distance, get_angle
      use mod_const, only: PI
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: natom, walkmax
      ! Internal water coordinates
      real(DP) :: rOH1, rOH2, aHOH_deg, aHOH_rad
      real(DP) :: rij(walkmax, 3)
      real(DP) :: Epot(walkmax)
      integer :: iw

      ! TODO: Check atom names
      if (natom /= 3) then
         call fatal_error(__FILE__, __LINE__, "This is not a water molecule!")
      end if

      ! Compute internal coordinates from cartesians
      ! OH distances in bohrs, HOH angle in radians
      do iw = 1, walkmax
         rOH1 = get_distance(x, y, z, 1, 2, iw)
         rOH2 = get_distance(x, y, z, 1, 3, iw)
         aHOH_deg = get_angle(x, y, z, 2, 1, 3, iw)
         aHOH_rad = aHOH_deg * PI / 180.0D0
         rij(iw, 1) = rOH1
         rij(iw, 2) = rOH2
         rij(iw, 3) = aHOH_rad
      end do

      call pes_h2o_schwenke(rij, Epot, walkmax)

      ! TODO: Implement numerical forces
      fx = 0.0D0; fy = 0.0D0; fz = 0.0D0

      do iw = 1, walkmax
         eclas = eclas + Epot(iw)
      end do
      eclas = eclas / walkmax
   end subroutine force_h2o_schwenke

end module mod_h2o_pes
