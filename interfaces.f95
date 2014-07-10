
module mod_interfaces
   use mod_array_size 
   INTERFACE
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

   subroutine init_mass(amg,amt)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(out) :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
   end subroutine 

   subroutine QtoX(x,y,z,transx,transy,transz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine QtoX
   subroutine XtoQ(x,y,z,transx,transy,transz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine XtoQ

   subroutine XtoU(x,y,z,transx,transy,transz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine XtoU
   subroutine UtoX(x,y,z,transx,transy,transz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)    :: transx(npartmax,nwalkmax),transy(npartmax,nwalkmax),transz(npartmax,nwalkmax)
   end subroutine UtoX

   subroutine FXtoFQ(fxab,fyab,fzab,fx,fy,fz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)     :: fxab(npartmax,nwalkmax),fyab(npartmax,nwalkmax),fzab(npartmax,nwalkmax)
   real*8,intent(inout)  :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   end subroutine FXtoFQ
   subroutine FQtoFX(fx,fy,fz,transfx,transfy,transfz)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)    :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(inout) :: transfx(npartmax,nwalkmax),transfy(npartmax,nwalkmax),transfz(npartmax,nwalkmax)
   end subroutine FQtoFX

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

   subroutine respastep(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
           fxc,fyc,fzc,fxq,fyq,fzq)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(inout)  :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(inout)  :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
   real*8,intent(inout)  :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8,intent(in)     :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
   real*8,intent(in)     :: dt
   real*8,intent(inout)  :: eclas,equant
   end subroutine respastep

   subroutine respashake(x,y,z,px,py,pz,amt,amg,dt,equant,eclas, &
              fxc,fyc,fzc,fxq,fyq,fzq)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(inout) :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8,intent(in)    :: amg(npartmax,nwalkmax),amt(npartmax,nwalkmax)
   real*8,intent(in)    :: dt
   real*8,intent(inout) :: eclas,equant
   real*8,intent(inout) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(inout) :: fxq(npartmax,nwalkmax),fyq(npartmax,nwalkmax),fzq(npartmax,nwalkmax)
   end subroutine respashake

   subroutine verletstep(x,y,z,px,py,pz,amt,dt,eclas,fxc,fyc,fzc)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(inout) :: fxc(npartmax,nwalkmax),fyc(npartmax,nwalkmax),fzc(npartmax,nwalkmax)
   real*8,intent(inout) :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8,intent(in)    :: amt(npartmax,nwalkmax)
   real*8,intent(in)    :: dt
   real*8,intent(inout) :: eclas
   end subroutine verletstep

   subroutine force_abin(x,y,z,fx,fy,fz,eclas)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)    :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out)   :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(out)   :: eclas
   end subroutine force_abin

   subroutine force_guillot(x,y,z,fx,fy,fz,eclas)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(in)  :: x(npartmax,nwalkmax),y(npartmax,nwalkmax),z(npartmax,nwalkmax)
   real*8,intent(out) :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(out) :: eclas
   end subroutine force_guillot

   SUBROUTINE shiftP (px,py,pz,fx,fy,fz,dt)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8,intent(in)    :: fx(npartmax,nwalkmax),fy(npartmax,nwalkmax),fz(npartmax,nwalkmax)
   real*8,intent(in)    :: dt
   END SUBROUTINE shiftP
   SUBROUTINE shiftX (rx,ry,rz,px,py,pz,mass,dt)
   IMPORT :: npartmax,nwalkmax
   real*8,intent(inout) :: rx(npartmax,nwalkmax),ry(npartmax,nwalkmax),rz(npartmax,nwalkmax)
   real*8,intent(in)    :: px(npartmax,nwalkmax),py(npartmax,nwalkmax),pz(npartmax,nwalkmax)
   real*8,intent(in)    :: mass(npartmax,nwalkmax)
   real*8,intent(in)    :: dt
   end subroutine shiftX

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

