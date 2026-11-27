! Surface-hopping adapter for the state-independent NaI model.
module mod_potentials_sh
   use mod_const, only: DP
   use mod_potentials_nai, only: nai_init, evaluate_nai
   implicit none
   private
   public :: nai_init, force_nai
contains
   subroutine force_nai(x, y, z, fx, fy, fz, eclas)
      use mod_sh, only: en_array, istate, nacx, nacy, nacz
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :), eclas
      real(DP) :: state_fx(size(x, 1), size(x, 2), 2)
      real(DP) :: state_fy(size(y, 1), size(y, 2), 2)
      real(DP) :: state_fz(size(z, 1), size(z, 2), 2)

      call evaluate_nai(x, y, z, en_array, state_fx, state_fy, state_fz, nacx, nacy, nacz)
      eclas = en_array(istate)
      fx = state_fx(:, :, istate)
      fy = state_fy(:, :, istate)
      fz = state_fz(:, :, istate)
   end subroutine force_nai
end module mod_potentials_sh
