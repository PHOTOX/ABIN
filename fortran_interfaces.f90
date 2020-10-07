! Fortran interfaces created by hand
! for functions outside of modules
module mod_interfaces
   use mod_const, only: DP
   INTERFACE

   subroutine init(dt)
   IMPORT :: DP
   real(DP),intent(out) :: dt
   end subroutine init

   subroutine print_compile_info()

   end subroutine print_compile_info

   ! has to be here or in its own module, as it depends on mod_sh
   ! and mod_sh depends on force_clas
   subroutine force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
   IMPORT :: DP
   real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out)   :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(out)   :: eclas
   integer,intent(in)     :: walkmax
   character(len=*),intent(in) :: chpot
   end subroutine force_abin

   subroutine force_clas(fx, fy, fz, x, y, z, eclas, chpot)
   IMPORT :: DP
   real(DP),intent(inout) :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(out)   :: eclas
   character(len=*),intent(in) :: chpot
   end subroutine force_clas

   subroutine force_quantum(fx, fy, fz, x, y, z, amg, energy)
   IMPORT :: DP
   real(DP),intent(in)  :: x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout) :: fx(:,:), fy(:,:), fz(:,:)
   real(DP),intent(in)  :: amg(:,:)
   real(DP),intent(out) :: energy
   end subroutine force_quantum

   subroutine force_wrapper(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
   IMPORT :: DP
   real(DP),intent(in)    ::  x(:,:),  y(:,:),  z(:,:)
   real(DP),intent(inout) :: fx(:,:), fy(:,:), fz(:,:)
   real(DP),intent(out)   :: eclas
   integer, intent(in)    :: walkmax
   character(len=*),intent(in) :: chpot
   end subroutine force_wrapper

   subroutine propagate_nm(x, y, z, px, py, pz, amg, dt)
   IMPORT :: DP
   real(DP), intent(inout) :: x(:,:), y(:,:), z(:,:)
   real(DP), intent(inout) :: px(:,:), py(:,:), pz(:,:)
   real(DP), intent(in)    ::  amg(:,:), dt
   end subroutine propagate_nm

   subroutine oniom(x, y, z, fx, fy, fz, eclas, iw)
   IMPORT :: DP
   real(DP),intent(in)      :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   :: eclas
   integer, intent(in)      :: iw
   end subroutine oniom

   END INTERFACE

end module mod_interfaces

