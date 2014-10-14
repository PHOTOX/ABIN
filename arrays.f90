!- Module for all dynamic arrays.
!- Contains routines for allocation and deallocation.
module mod_arrays
   use mod_const, only: DP
   implicit none
   real(DP), allocatable  :: x(:,:),y(:,:),z(:,:)
   real(DP), allocatable  :: amt(:,:),amg(:,:)
   real(DP), allocatable  :: fxc(:,:),fyc(:,:),fzc(:,:)
   real(DP), allocatable  :: fxq(:,:),fyq(:,:),fzq(:,:)
   real(DP), allocatable  :: vx(:,:),vy(:,:),vz(:,:)
   real(DP), allocatable  :: px(:,:),py(:,:),pz(:,:)
   real(DP), allocatable  :: transx(:,:),transy(:,:),transz(:,:)
   real(DP), allocatable  :: transfxc(:,:),transfyc(:,:),transfzc(:,:)
   real(DP), allocatable  :: transxv(:,:),transyv(:,:),transzv(:,:)
   real(DP), allocatable  :: vx_old(:,:),vy_old(:,:),vz_old(:,:)
   contains

   subroutine allocate_arrays(natomalloc, nwalkalloc)
      integer,intent(in)  :: nwalkalloc
      integer  :: natomalloc

      ! This is to avoid segfault for natom=1 and ndist>0
      if (natomalloc.eq.1)then
         allocate( x(natomalloc+1, nwalkalloc) )
         allocate( y(natomalloc+1, nwalkalloc) )
         allocate( z(natomalloc+1, nwalkalloc) )
         allocate( transx(natomalloc+1, nwalkalloc) )
         allocate( transy(natomalloc+1, nwalkalloc) )
         allocate( transz(natomalloc+1, nwalkalloc) )
      else
         allocate( x(natomalloc, nwalkalloc) )
         allocate( y(natomalloc, nwalkalloc) )
         allocate( z(natomalloc, nwalkalloc) )
         allocate( transx(natomalloc, nwalkalloc) )
         allocate( transy(natomalloc, nwalkalloc) )
         allocate( transz(natomalloc, nwalkalloc) )
      end if

      x=0.0d0
      !automatic allocation
      y=x
      z=x
      transx=x
      transy=transx
      transz=transx

      allocate( amt(natomalloc, nwalkalloc) )
      allocate( amg(natomalloc, nwalkalloc) )
      amt=0.0d0; amg=amt


      allocate( px(natomalloc, nwalkalloc) )
      allocate( py(natomalloc, nwalkalloc) )
      allocate( pz(natomalloc, nwalkalloc) )
      allocate( vx(natomalloc, nwalkalloc) )
      allocate( vy(natomalloc, nwalkalloc) )
      allocate( vz(natomalloc, nwalkalloc) )
      allocate( transxv(natomalloc, nwalkalloc) )
      allocate( transyv(natomalloc, nwalkalloc) )
      allocate( transzv(natomalloc, nwalkalloc) )

      px=0.0d0;   py = px; pz=px 
      vx=px;      vy=vx;   vz=vx
      vx_old=vx;  vy_old=vx_old; vz_old=vx_old
      transxv=vx ;transyv=transxv; transzv=transxv

      allocate( fxc(natomalloc, nwalkalloc) )
      allocate( fyc(natomalloc, nwalkalloc) )
      allocate( fzc(natomalloc, nwalkalloc) )
      allocate( fxq(natomalloc, nwalkalloc) )
      allocate( fyq(natomalloc, nwalkalloc) )
      allocate( fzq(natomalloc, nwalkalloc) )
      allocate( transfxc(natomalloc, nwalkalloc) )
      allocate( transfyc(natomalloc, nwalkalloc) )
      allocate( transfzc(natomalloc, nwalkalloc) )

      fxc=0.0d0;  fyc=fxc; fzc=fxc
      fxq=fxc;    fyq=fxc; fzq=fxc
      transfxc=fxc;  transfyc=fxc;  transfzc=fxc

   end subroutine allocate_arrays

   subroutine deallocate_arrays
      deallocate( x ); deallocate( y ); deallocate( z );
      deallocate( amt ); deallocate( amg )
      deallocate( px ); deallocate( py ); deallocate( pz );
      deallocate( vx ); deallocate( vy ); deallocate( vz );
      deallocate( vx_old ); deallocate( vy_old ); deallocate( vz_old );
      deallocate( transx ); deallocate( transy ); deallocate( transz );
      deallocate( transxv ); deallocate( transyv ); deallocate( transzv );
      deallocate( transfxc ); deallocate( transfyc ); deallocate( transfzc );
      deallocate( fxc ); deallocate( fyc ); deallocate( fzc );
      deallocate( fxq ); deallocate( fyq ); deallocate( fzq );
   end subroutine deallocate_arrays
       

end module mod_arrays
