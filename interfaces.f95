
module mod_interfaces
   use mod_array_size 
   INTERFACE
      !TODO: do module minimize
   subroutine minimize(x,y,z,fx,fy,fz,eclas)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(inout) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(inout) :: eclas
   end subroutine minimize

   subroutine trajout(x,y,z,it)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   integer,intent(in) :: it
   end subroutine trajout

   subroutine init(x,y,z,vx,vy,vz,fxc,fyc,fzc,fxq,fyq,fzq,dt)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(out) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   real*8,intent(out) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(out) :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
   real*8,intent(out) :: dt
   end subroutine init

   subroutine force_clas(fx,fy,fz,x,y,z,energy)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) ::  x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(inout) ::  fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(out)   ::  energy
   end subroutine force_clas

   subroutine force_quantum(fx,fy,fz,x,y,z,amg,energy)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(in) :: amg(npartmax,nwalkmax)
   real*8,intent(inout) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(out) :: energy
   end subroutine force_quantum

   subroutine force_abin(x,y,z,fx,fy,fz,eclas)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)    :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)   :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(out)   :: eclas
   end subroutine force_abin

   subroutine analysis(x,y,z,vx,vy,vz,fxc,fyc,fzc,amt,eclas,equant,dt)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(in) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(in) :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   real*8,intent(in) :: amt(npartmax,nwalkmax)
   real*8,intent(in) :: eclas,equant
   real*8 :: dt  
   end subroutine analysis

   subroutine restout(x,y,z,vx,vy,vz,it)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(in)  :: vx(npartmax,nwalkmax),vy(npartmax,nwalkmax),vz(npartmax,nwalkmax)
   integer,intent(in) :: it
   end subroutine restout

   subroutine finish(values1,values2)
   integer,dimension(8),intent(in)  :: values1
   integer,dimension(8),intent(out) :: values2
   end subroutine finish

   END INTERFACE

end module mod_interfaces

