! Interface to FFTW library for Fast Fourier Transform
! Currently utilized for normal mode transformation in Path Integral simulations
module mod_fftw3
   use, intrinsic :: iso_c_binding
#ifndef USE_FFTW
   use mod_utils, only: not_compiled_with
#endif
   private
   public :: fftw_normalmodes_init, fftw_normalmodes_finalize
   public :: dft_normalmode2cart, dft_cart2normalmode
   ! TODO: Is there a better way to do this in newer FFTW versions?
#ifdef USE_FFTW
   include 'fftw3.F90'
   type(C_PTR) :: plan_utox, plan_xtou
#endif
   save
contains

#ifdef USE_FFTW
   subroutine fftw_normalmodes_init(nwalk)
      use mod_const, only: DP
      use mod_utils, only: abinerror
      integer, intent(in) :: nwalk
      real(C_DOUBLE), dimension(:), allocatable :: x_tmp
      complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: cx_tmp

      if (DP /= C_DOUBLE) then
         write (*, *) 'WARNING: Kind DP is not equal kind C_DOUBLE'
         write (*, *) 'Precision might be lost during normal mode transform.'
         write (*, *) 'Set iknow=1 if you want to proceed.'
         if (iknow /= 1) call abinerror('fftw_init')
      end if

      allocate (x_tmp(nwalk + 1))
      allocate (cx_tmp(nwalk + 1))

      plan_xtou = fftw_plan_dft_r2c_1d(nwalk, x_tmp, cx_tmp, FFTW_MEASURE)
      plan_utox = fftw_plan_dft_c2r_1d(nwalk, cx_tmp, x_tmp, FFTW_MEASURE)

      deallocate (x_tmp)
      deallocate (cx_tmp)
   end subroutine

   ! Simple wrapper functions around the FFTW interface.
   subroutine dft_normalmode2cart(nm, cart)
      complex(C_DOUBLE_COMPLEX), dimension(:), intent(inout) :: nm
      real(C_DOUBLE), dimension(:), intent(out) :: cart

      call fftw_execute_dft_c2r(plan_utox, nm, cart)
   end subroutine dft_normalmode2cart

   subroutine dft_cart2normalmode(cart, nm)
      real(C_DOUBLE), dimension(:), intent(inout) :: cart
      complex(C_DOUBLE_COMPLEX), dimension(:), intent(out) :: nm

      call fftw_execute_dft_r2c(plan_xtou, cart, nm)
   end subroutine dft_cart2normalmode

   subroutine fftw_normalmodes_finalize()
      call fftw_destroy_plan(plan_xtou)
      call fftw_destroy_plan(plan_utox)
   end subroutine fftw_normalmodes_finalize

#else

   ! Dummy functions when ABIN is not compiled with FFTW
   subroutine fftw_normalmodes_init(nwalk)
      use mod_const, only: DP
      use mod_utils, only: abinerror
      integer :: nwalk
      nwalk = 0
      write (*, *) 'ERROR: Normal mode transformations cannot be performed.'
      call not_compiled_with('FFTW library', 'fftw_normalmodes_init')
   end subroutine fftw_normalmodes_init

   subroutine dft_normalmode2cart(nm, cart)
      use mod_utils, only: abinerror
      complex(C_DOUBLE_COMPLEX), dimension(:) :: nm
      real(C_DOUBLE), dimension(:) :: cart
      cart = 0.0D0
      nm = (0.0D0, 0.0D0)
      call not_compiled_with('FFTW library', 'dft_normalmode2cart')
   end subroutine dft_normalmode2cart

   subroutine dft_cart2normalmode(cart, nm)
      use mod_utils, only: abinerror
      complex(C_DOUBLE_COMPLEX), dimension(:) :: nm
      real(C_DOUBLE), dimension(:) :: cart
      cart = 0.0D0
      nm = (0.0D0, 0.0D0)
      call not_compiled_with('FFTW library', 'dft_cart2normalmode')
   end subroutine dft_cart2normalmode

   ! This must be a no-op and must not call abinerror()
   ! since it is itself called from abinerror().
   subroutine fftw_normalmodes_finalize()
   end subroutine fftw_normalmodes_finalize

#endif
end module mod_fftw3
