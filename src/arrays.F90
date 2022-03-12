! Module for all dynamic arrays.
! Contains routines for allocation and deallocation.
module mod_arrays
   use mod_const, only: DP
   implicit none

   ! TODO: Stop this madness, put everything into Traj_Data derived type
   ! and then create type(Traj_Data) :: trajectory. trajectory_prev
   ! (to implement adaptive time step for SH)

   ! will require some refactor (i.e. everywhere where we import mod_arrays)
   real(DP), allocatable :: x(:, :), y(:, :), z(:, :)
   real(DP), allocatable :: x_old(:, :), y_old(:, :), z_old(:, :)
   real(DP), allocatable :: amt(:, :), amg(:, :)
   ! "classical" forces
   real(DP), allocatable :: fxc(:, :), fyc(:, :), fzc(:, :)
   real(DP), allocatable :: fxc_old(:, :), fyc_old(:, :), fzc_old(:, :)
   ! "quantum" forces for PIMD
   real(DP), allocatable :: fxq(:, :), fyq(:, :), fzq(:, :)
   real(DP), allocatable :: fxq_old(:, :), fyq_old(:, :), fzq_old(:, :)
   real(DP), allocatable :: vx(:, :), vy(:, :), vz(:, :)
   real(DP), allocatable :: vx_old(:, :), vy_old(:, :), vz_old(:, :)
   real(DP), allocatable :: px(:, :), py(:, :), pz(:, :)
   ! helper arrays used througout the code
   real(DP), allocatable :: transx(:, :), transy(:, :), transz(:, :)
   real(DP), allocatable :: transfxc(:, :), transfyc(:, :), transfzc(:, :)
   real(DP), allocatable :: transxv(:, :), transyv(:, :), transzv(:, :)
   ! This holds the difference between reference and full potential
   ! for multiple time-stepping and ring contraction ala O. Marsalek
   real(DP), allocatable :: fxc_diff(:, :), fyc_diff(:, :), fzc_diff(:, :)
   real(DP) :: epot_diff = 0.0D0 ! not really array, but for now let's keep it here
   save
contains

   subroutine allocate_arrays(natomalloc, nwalkalloc)
      use mod_general, only: pot_ref
      integer, intent(in) :: nwalkalloc
      integer, intent(in) :: natomalloc

      allocate (x(natomalloc, nwalkalloc))
      allocate (y(natomalloc, nwalkalloc))
      allocate (z(natomalloc, nwalkalloc))
      allocate (transx(natomalloc, nwalkalloc))
      allocate (transy(natomalloc, nwalkalloc))
      allocate (transz(natomalloc, nwalkalloc))
      allocate (x_old(natomalloc, nwalkalloc))
      allocate (y_old(natomalloc, nwalkalloc))
      allocate (z_old(natomalloc, nwalkalloc))

      x = 0.0D0
      y = x
      z = x
      x_old = 0.0D0
      y_old = x
      z_old = x
      transx = x
      transy = transx
      transz = transx

      allocate (amt(natomalloc, nwalkalloc))
      allocate (amg(natomalloc, nwalkalloc))
      amt = -1; amg = amt

      allocate (px(natomalloc, nwalkalloc))
      allocate (py(natomalloc, nwalkalloc))
      allocate (pz(natomalloc, nwalkalloc))
      allocate (vx(natomalloc, nwalkalloc))
      allocate (vy(natomalloc, nwalkalloc))
      allocate (vz(natomalloc, nwalkalloc))
      allocate (vx_old(natomalloc, nwalkalloc))
      allocate (vy_old(natomalloc, nwalkalloc))
      allocate (vz_old(natomalloc, nwalkalloc))
      allocate (transxv(natomalloc, nwalkalloc))
      allocate (transyv(natomalloc, nwalkalloc))
      allocate (transzv(natomalloc, nwalkalloc))

      px = 0.0D0; py = px; pz = px
      vx = px; vy = vx; vz = vx
      vx_old = vx; vy_old = vx_old; vz_old = vx_old
      transxv = vx; transyv = transxv; transzv = transxv

      allocate (fxc(natomalloc, nwalkalloc))
      allocate (fyc(natomalloc, nwalkalloc))
      allocate (fzc(natomalloc, nwalkalloc))
      allocate (fxc_old(natomalloc, nwalkalloc))
      allocate (fyc_old(natomalloc, nwalkalloc))
      allocate (fzc_old(natomalloc, nwalkalloc))
      allocate (fxq(natomalloc, nwalkalloc))
      allocate (fyq(natomalloc, nwalkalloc))
      allocate (fzq(natomalloc, nwalkalloc))
      allocate (fxq_old(natomalloc, nwalkalloc))
      allocate (fyq_old(natomalloc, nwalkalloc))
      allocate (fzq_old(natomalloc, nwalkalloc))
      allocate (transfxc(natomalloc, nwalkalloc))
      allocate (transfyc(natomalloc, nwalkalloc))
      allocate (transfzc(natomalloc, nwalkalloc))

      fxc = 0.0D0; fyc = fxc; fzc = fxc
      fxq = fxc; fyq = fxc; fzq = fxc
      transfxc = fxc; transfyc = fxc; transfzc = fxc
      if (pot_ref /= '_none_') then
         allocate (fxc_diff(natomalloc, nwalkalloc))
         allocate (fyc_diff(natomalloc, nwalkalloc))
         allocate (fzc_diff(natomalloc, nwalkalloc))
         fxc_diff = 0.0D0; fyc_diff = 0.0D0; fzc_diff = 0.0D0
      end if

   end subroutine allocate_arrays

   subroutine deallocate_arrays
      ! If x is allocated, all of them should be.
      ! If not, something is horribly wrong anyway
      if (allocated(x)) then
         deallocate (x); deallocate (y); deallocate (z); 
         deallocate (x_old); deallocate (y_old); deallocate (z_old); 
         deallocate (amt); deallocate (amg)
         deallocate (px); deallocate (py); deallocate (pz); 
         deallocate (vx); deallocate (vy); deallocate (vz); 
         deallocate (vx_old); deallocate (vy_old); deallocate (vz_old); 
         deallocate (transx); deallocate (transy); deallocate (transz); 
         deallocate (transxv); deallocate (transyv); deallocate (transzv); 
         deallocate (transfxc); deallocate (transfyc); deallocate (transfzc); 
         deallocate (fxc); deallocate (fyc); deallocate (fzc); 
         deallocate (fxq); deallocate (fyq); deallocate (fzq); 
         deallocate (fxc_old); deallocate (fyc_old); deallocate (fzc_old); 
         deallocate (fxq_old); deallocate (fyq_old); deallocate (fzq_old); 
      end if
      if (allocated(fxc_diff)) then
         deallocate (fxc_diff); deallocate (fyc_diff); deallocate (fzc_diff); 
      end if
   end subroutine deallocate_arrays

end module mod_arrays
