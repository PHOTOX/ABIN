
module mod_interfaces
   use mod_const, only: DP
   INTERFACE

   subroutine init(dt)
   IMPORT :: DP
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

   subroutine oniom(x, y, z, fx, fy, fz, eclas, iw)
   IMPORT :: DP
   real(DP),intent(in)      :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   :: eclas
   integer, intent(in)      :: iw
   end subroutine oniom

   END INTERFACE

end module mod_interfaces

