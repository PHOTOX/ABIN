! Fortran interfaces created by hand
! for functions outside of modules
! TODO: Everything should be in modules!
! Some functions are currently outside of modules due to
! circular dependencies.
module mod_interfaces
   use, intrinsic :: iso_c_binding, only: C_INT, C_INT32_T
   use mod_const, only: DP
   implicit none
   public
   interface

      subroutine print_compile_info()
         implicit none
      end subroutine print_compile_info

      subroutine finish(error_code)
         implicit none
         integer, intent(in) :: error_code
      end subroutine finish

      subroutine force_clas(fx, fy, fz, x, y, z, eclas, chpot)
         import :: DP
         implicit none
         real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(out) :: eclas
         character(len=*), intent(in) :: chpot
      end subroutine force_clas

      subroutine force_quantum(fx, fy, fz, x, y, z, amg, energy)
         import :: DP
         implicit none
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(in) :: amg(:, :)
         real(DP), intent(out) :: energy
      end subroutine force_quantum

      subroutine force_wrapper(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
         import :: DP
         implicit none
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(out) :: eclas
         integer, intent(in) :: walkmax
         character(len=*), intent(in) :: chpot
      end subroutine force_wrapper

      subroutine omp_set_num_threads(nthreads)
         implicit none
         integer, intent(in) :: nthreads
      end subroutine omp_set_num_threads

      ! TODO: This interface currently doesn't work, probably because
      ! of how we pass the 2D arrays...
      !subroutine force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot) !bind(C)
      !   use, intrinsic :: iso_c_binding
      !   implicit none
      !   real(C_DOUBLE), intent(in) :: x(:, :), y(:, :), z(:, :)
      !   real(C_DOUBLE), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      !   real(C_DOUBLE), intent(inout) :: eclas
      !   integer(C_INT), intent(in) :: natom, walkmax, watpot
      !end subroutine

      ! https://cyber.dabamos.de/programming/modernfortran/sleep.html
      ! int usleep(useconds_t useconds)
      function usleep(useconds) bind(c, name='usleep')
         import :: C_INT, C_INT32_T
         implicit none
         integer(kind=C_INT32_T), intent(in), value :: useconds
         integer(kind=C_INT) :: usleep
      end function usleep

      ! Computes potential energy of a water molecule
      ! using Schwenke potential, see h2o_schwenke.f
      subroutine h2o_pot_schwenke(rij, v, n)
         import :: DP
         implicit none
         integer, intent(in) :: n
         real(DP), intent(in) :: rij(n, 3)
         real(DP), intent(out) :: v(n)
      end subroutine h2o_pot_schwenke

      subroutine h2o_pot_cvrqd(V, rOH1, rOH2, aHOH, mH, mO)
         import :: DP
         implicit none
         real(DP), intent(out) :: V
         real(DP), intent(in) :: rOH1, rOH2, aHOH, mH, mO
      end subroutine h2o_pot_cvrqd

   end interface

end module mod_interfaces
