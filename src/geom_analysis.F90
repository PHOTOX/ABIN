module mod_analyze_geometry
   use mod_const, only: DP, ANG, AUtoFS
   use mod_array_size, only: NDISTMAX
   use mod_general, only: sim_time, nwalk
   implicit none
   private
   public :: ndist, nang, ndih, shiftdih
   public :: dist1, dist2, ang1, ang2, ang3, dih1, dih2, dih3, dih4
   public :: print_distances, print_angles, print_dihedrals

   ! TODO: Make all these allocatable, instead of using NDISTMAX
   integer :: ndist = 0, dist1(NDISTMAX), dist2(NDISTMAX)
   integer :: nang = 0, ang1(NDISTMAX), ang2(NDISTMAX), ang3(NDISTMAX)
   integer :: ndih = 0, dih1(NDISTMAX), dih2(NDISTMAX), dih3(NDISTMAX), dih4(NDISTMAX)
   real(DP) :: shiftdih = 360.0D0 ! 0 for (-180,180), 360 for (0,360)
   save

contains

   subroutine print_distances(x, y, z)
      use mod_files, only: UDIST
      use mod_general, only: natom
      use mod_utils, only: get_distance
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: r(NDISTMAX)
      integer :: idist, iw
      character(len=50) :: frmt

      write (frmt, '(A7,I3,A6)') '(F15.2,', ndist, 'F12.7)'

      do iw = 1, nwalk

         do idist = 1, ndist
            ! If we have only one atom, take the distance from the origin,
            ! i.e. [0, 0, 0]
            if (natom == 1) then
               r(idist) = dsqrt(x(1, iw)**2 + y(1, iw)**2 + z(1, iw)**2)
            else
               r(idist) = get_distance(x, y, z, dist1(idist), dist2(idist), iw)
            end if
            r(idist) = r(idist) / ANG
         end do

         write (UDIST, frmt) sim_time * AUtoFS, (r(idist), idist=1, ndist)

      end do
   end subroutine print_distances

   subroutine print_angles(x, y, z)
      use mod_utils, only: get_angle
      use mod_files, only: UANG
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: alfa(NDISTMAX)
      integer :: idist, iw
      character(len=50) :: frmt

      write (frmt, '(A7,I3,A6)') '(F15.2,', nang, 'F12.7)'

      do iw = 1, nwalk

         do idist = 1, nang
            alfa(idist) = get_angle(x, y, z, ang1(idist), ang2(idist), ang3(idist), iw)
         end do

         write (UANG, frmt) sim_time * AUtoFS, (alfa(idist), idist=1, nang)

      end do
   end subroutine print_angles

   subroutine print_dihedrals(x, y, z)
      use mod_utils, only: get_dihedral
      use mod_files, only: UDIH
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP) :: delta(NDISTMAX)
      integer :: idist, iw
      character(len=50) :: frmt

      write (frmt, '(A7,I3,A6)') '(F15.2,', ndih, 'F12.7)'

      do iw = 1, nwalk

         do idist = 1, ndih
            delta(idist) = get_dihedral(x, y, z, &
                                        dih1(idist), dih2(idist), dih3(idist), dih4(idist), &
                                        iw, shiftdih)
         end do

         write (UDIH, frmt) sim_time * AUtoFS, (delta(idist), idist=1, ndih)

      end do
   end subroutine print_dihedrals

end module mod_analyze_geometry
