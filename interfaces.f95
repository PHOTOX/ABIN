
module mod_interfaces
   use mod_const, only: DP
   INTERFACE

   subroutine init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt)
   IMPORT :: DP
   real(DP),intent(out) :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out) :: vx(:,:),vy(:,:),vz(:,:)
   real(DP),intent(out) :: fxc(:,:),fyc(:,:),fzc(:,:)
   real(DP),intent(out) :: fxq(:,:),fyq(:,:),fzq(:,:)
   real(DP),intent(out) :: dt
   end subroutine init

   !has to be here or in its own module, as it depends on mod_sh
   ! and mod_sh depends on mod_forces
   subroutine force_abin(x,y,z,fx,fy,fz,eclas)
   IMPORT :: DP
   real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(out)   :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(out)   :: eclas
   end subroutine force_abin

   END INTERFACE

end module mod_interfaces

