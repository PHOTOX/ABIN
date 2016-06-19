
module mod_fftw3
   use, intrinsic :: iso_c_binding
   include 'fftw3.f90'
   type(C_PTR) :: plan_utox,plan_xtou
   real(C_DOUBLE), dimension(:), allocatable :: x_in, y_in, z_in
   complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: cx, cy, cz
   save
   contains

   subroutine fftw_init(nwalk)
      use mod_const, only: DP
      use mod_utils, only: abinerror
      integer :: nwalk

      if(DP.ne.C_DOUBLE)then
         write(*,*)'WARNING: Kind DP is not equal kind C_DOUBLE'
         write(*,*)'Precision might be lost during normal mode transform.'
         write(*,*)'Set iknow=1 if you want to proceed.'
         if (iknow.ne.1) call abinerror('fftw_init')
      end if

      allocate( x_in(nwalk+1) )
      allocate( y_in(nwalk+1) )
      allocate( z_in(nwalk+1) )

      allocate( cx(nwalk+1) )
      allocate( cy(nwalk+1) )
      allocate( cz(nwalk+1) )

      plan_xtou=fftw_plan_dft_r2c_1d(nwalk,x_in,cx, FFTW_MEASURE)
      plan_utox=fftw_plan_dft_c2r_1d(nwalk,cx,x_in, FFTW_MEASURE)

   end subroutine

   subroutine fftw_end()
      if(allocated(x_in))then
         deallocate( x_in )
         deallocate( y_in )
         deallocate( z_in )
         deallocate( cx )
         deallocate( cy )
         deallocate( cz )
         call fftw_destroy_plan(plan_xtou)
         call fftw_destroy_plan(plan_utox)
      endif
   end subroutine fftw_end

end module mod_fftw3

