! Fortran interfaces created by hand
! for functions outside of modules
! TODO: Everything should be in modules!
! Some functions are currently outside of modules due to
! circular dependencies.
module mod_interfaces
   use, intrinsic :: iso_c_binding, only: C_INT, C_INT32_T
   use mod_const, only: DP
   public
   interface

      subroutine print_compile_info()
      end subroutine print_compile_info

      subroutine finish(error_code)
         integer, intent(in) :: error_code
      end subroutine finish

      subroutine force_clas(fx, fy, fz, x, y, z, eclas, chpot)
         import :: DP
         real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(out) :: eclas
         character(len=*), intent(in) :: chpot
      end subroutine force_clas

      subroutine force_quantum(fx, fy, fz, x, y, z, amg, energy)
         import :: DP
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(in) :: amg(:, :)
         real(DP), intent(out) :: energy
      end subroutine force_quantum

      subroutine force_wrapper(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
         import :: DP
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(out) :: eclas
         integer, intent(in) :: walkmax
         character(len=*), intent(in) :: chpot
      end subroutine force_wrapper

      subroutine omp_set_num_threads(nthreads)
         integer, intent(in) :: nthreads
      end subroutine omp_set_num_threads

      ! TODO: This interface currently doesn't work, probably because
      ! of how we pass the 2D arrays...
      !subroutine force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot) !bind(C)
      !   use, intrinsic :: iso_c_binding
      !   real(C_DOUBLE), intent(in) :: x(:, :), y(:, :), z(:, :)
      !   real(C_DOUBLE), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      !   real(C_DOUBLE), intent(inout) :: eclas
      !   integer(C_INT), intent(in) :: natom, walkmax, watpot
      !end subroutine

      ! https://cyber.dabamos.de/programming/modernfortran/sleep.html
      ! int usleep(useconds_t useconds)
      function usleep(useconds) bind(c, name='usleep')
         import :: C_INT, C_INT32_T
         integer(kind=C_INT32_T), value :: useconds
         integer(kind=C_INT) :: usleep
      end function usleep

   end interface

end module mod_interfaces
