! Fortran interfaces created by hand
! for functions outside of modules
module mod_interfaces
   use mod_const, only: DP
   public
   interface

      subroutine init(dt)
         import :: DP
         real(DP), intent(out) :: dt
      end subroutine init

      subroutine finish(error_code)
         integer, intent(in) :: error_code
      end subroutine finish

      subroutine print_compile_info()

      end subroutine print_compile_info

      ! has to be here or in its own module, as it depends on mod_sh
      ! and mod_sh depends on force_clas
      subroutine force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
         import :: DP
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(out) :: eclas
         integer, intent(in) :: walkmax
         character(len=*), intent(in) :: chpot
      end subroutine force_abin

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

      subroutine propagate_nm(x, y, z, px, py, pz, amg, dt)
         import :: DP
         real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
         real(DP), intent(in) :: amg(:, :), dt
      end subroutine propagate_nm

      subroutine oniom(x, y, z, fx, fy, fz, eclas, iw)
         import :: DP
         real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
         real(DP), intent(inout) :: eclas
         integer, intent(in) :: iw
      end subroutine oniom

      subroutine omp_set_num_threads(nthreads)
         integer, intent(in) :: nthreads
      end subroutine omp_set_num_threads

      subroutine print_runtime_info()
      end subroutine

      ! TODO: This interface currently doesn't work, probably because
      ! of how we pass the 2D arrays...
      !subroutine force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot) !bind(C)
      !   use, intrinsic :: iso_c_binding
      !   real(C_DOUBLE), intent(in) :: x(:, :), y(:, :), z(:, :)
      !   real(C_DOUBLE), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      !   real(C_DOUBLE), intent(inout) :: eclas
      !   integer(C_INT), intent(in) :: natom, walkmax, watpot
      !end subroutine

   end interface

end module mod_interfaces
