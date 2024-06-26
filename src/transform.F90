! Transformation routines for PIMD and PI+GLE.
! We currently support staging (istage=1)
! and normal mode (inormalmodes=1) coordinates
!
! Here's the code for the function names:
! X - cartesian coordinates
! Q - staging coordinates
! U - normal mode coordinates
module mod_transform
   use, intrinsic :: iso_c_binding
   use mod_const, only: DP
   use mod_general, only: natom, nwalk
   use mod_error, only: fatal_error
   implicit none
   private
   public :: UtoX, XtoU, QtoX, XtoQ, FXtoFQ, FQtoFX
   public :: initialize_pi_masses, initialize_pi_transforms
   public :: equant_nm
   public :: finalize_normalmodes
   real(C_DOUBLE), dimension(:), allocatable :: x_in, y_in, z_in
   complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: cx, cy, cz

contains

   ! Staging/normal mode transformations for Path Integrals
   ! Masses, velocities and positions are transformed here into a new set of u variables
   ! See Tuckermann's article in "Quantum Simulations of Complex Many Body Systems'.
   subroutine initialize_pi_transforms(x, y, z, vx, vy, vz)
      use mod_arrays, only: transx, transy, transz
      use mod_general, only: istage, inormalmodes
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: vx(:, :), vy(:, :), vz(:, :)

      ! Staging transform
      if (istage == 1) then
         call XtoQ(x, y, z, transx, transy, transz)
         x = transx
         y = transy
         z = transz
         call XtoQ(vx, vy, vz, transx, transy, transz)
         vx = transx
         vy = transy
         vz = transz
         ! Normal mode transform
      else if (inormalmodes > 0) then
         call initialize_normalmodes()

         call XtoU(x, y, z, transx, transy, transz)
         x = transx
         y = transy
         z = transz
         call XtoU(vx, vy, vz, transx, transy, transz)
         vx = transx
         vy = transy
         vz = transz
      else
         call fatal_error(__FILE__, __LINE__, &
            & 'Invalid call to initialize_pi_transforms')
      end if
   end subroutine initialize_pi_transforms

   ! Allocate auxiliary arrays and FFTW library
   subroutine initialize_normalmodes()
      use mod_general, only: nwalk
      use mod_fftw3, only: fftw_normalmodes_init

      ! Initialize FFTW library
      call fftw_normalmodes_init(nwalk)

      allocate (x_in(nwalk))
      allocate (y_in(nwalk))
      allocate (z_in(nwalk))

      allocate (cx(nwalk))
      allocate (cy(nwalk))
      allocate (cz(nwalk))
   end subroutine initialize_normalmodes

   subroutine finalize_normalmodes()
      use mod_fftw3, only: fftw_normalmodes_finalize
      if (allocated(x_in)) then
         deallocate (x_in)
         deallocate (y_in)
         deallocate (z_in)
         deallocate (cx)
         deallocate (cy)
         deallocate (cz)
      end if
      call fftw_normalmodes_finalize()
   end subroutine finalize_normalmodes

   ! This routine transforms staging coordinates to cartesian coordinates, which are stored
   ! in trans matrices, values in x, y and z matrices are NOT modified!!
   subroutine QtoX(x, y, z, transx, transy, transz)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: transx(:, :), transy(:, :), transz(:, :)
      integer :: iat, iw

      do iat = 1, natom

         transx(iat, 1) = x(iat, 1)
         transy(iat, 1) = y(iat, 1)
         transz(iat, 1) = z(iat, 1)

         transx(iat, nwalk) = x(iat, nwalk) + ((nwalk - 1.0D0) / nwalk) * &
                            & transx(iat, 1) + x(iat, 1) / nwalk
         transy(iat, nwalk) = y(iat, nwalk) + ((nwalk - 1.0D0) / nwalk) * &
                            & transy(iat, 1) + y(iat, 1) / nwalk
         transz(iat, nwalk) = z(iat, nwalk) + ((nwalk - 1.0D0) / nwalk) * &
                            & transz(iat, 1) + z(iat, 1) / nwalk

         do iw = nwalk - 1, 2, -1
            transx(iat, iw) = x(iat, iw) + ((iw - 1.0D0) / iw) * &
                            & transx(iat, iw + 1) + x(iat, 1) / iw
            transy(iat, iw) = y(iat, iw) + ((iw - 1.0D0) / iw) * &
                            & transy(iat, iw + 1) + y(iat, 1) / iw
            transz(iat, iw) = z(iat, iw) + ((iw - 1.0D0) / iw) * &
                            & transz(iat, iw + 1) + z(iat, 1) / iw
         end do

      end do
   end subroutine QtoX

   ! This routine transforms cartesian to staging coordinates, which are stored
   ! in trans matrices, values in x,y and z matrices are NOT modified!
   subroutine XtoQ(x, y, z, transx, transy, transz)
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: transx(:, :), transy(:, :), transz(:, :)
      integer :: iat, iw

      do iat = 1, natom

         transx(iat, 1) = x(iat, 1)
         transy(iat, 1) = y(iat, 1)
         transz(iat, 1) = z(iat, 1)

         do iw = 2, nwalk - 1
            transx(iat, iw) = ((iw - 1.0D0) * x(iat, iw + 1) + x(iat, 1)) / iw
            transx(iat, iw) = x(iat, iw) - transx(iat, iw)
            transy(iat, iw) = ((iw - 1.0D0) * y(iat, iw + 1) + y(iat, 1)) / iw
            transy(iat, iw) = y(iat, iw) - transy(iat, iw)
            transz(iat, iw) = ((iw - 1.0D0) * z(iat, iw + 1) + z(iat, 1)) / iw
            transz(iat, iw) = z(iat, iw) - transz(iat, iw)
         end do

         transx(iat, nwalk) = x(iat, nwalk) - ((nwalk - 1.0D0) * x(iat, 1) + x(iat, 1)) / nwalk
         transy(iat, nwalk) = y(iat, nwalk) - ((nwalk - 1.0D0) * y(iat, 1) + y(iat, 1)) / nwalk
         transz(iat, nwalk) = z(iat, nwalk) - ((nwalk - 1.0D0) * z(iat, 1) + z(iat, 1)) / nwalk

      end do

   end subroutine XtoQ

   subroutine initialize_pi_masses(am, amg, amt)
      use mod_const, only: PI
      use mod_general, only: istage, inormalmodes, idebug
      ! Physical atomic masses (already in atomic units)
      real(DP), intent(in) :: am(:)
      ! For PIMD, we may transform physical masses
      ! into something else to improve sampling.
      ! Each bead may have different mass, hence the 2D array.
      ! This array is used throughout the code even for non-PI MD.
      real(DP), intent(out) :: amt(:, :)
      ! These are used for the propagation of PI bead necklaces,
      ! see force_quantum()
      real(DP), intent(out) :: amg(:, :)
      real(DP) :: lambda(size(amg, 2) + 1) !amg(natom,nwalk) - why +1?
      integer :: iat, iw, nmodes, odd

      do iw = 1, nwalk
         do iat = 1, natom
            amg(iat, iw) = am(iat)
            amt(iat, iw) = am(iat)
         end do
      end do

      ! Transforming masses for PIMD with staging coordinates..
      ! There are two sets of masses associated with the transformed coordinates:
      ! amg according to eq. 90 in Tuckermann and amt (which has non-zero first mass)
      ! amt is associated with the kinetic energy, amg with the potential energy.
      if (istage == 1) then
         do iat = 1, natom
            amg(iat, 1) = 0.0D0
            amt(iat, 1) = am(iat)
            do iw = 2, nwalk
               amg(iat, iw) = (iw / (iw - 1.0D0)) * am(iat)
               amt(iat, iw) = (iw / (iw - 1.0D0)) * am(iat)
            end do
         end do
      end if

      ! Transforming masses for PIMD with normal modes
      ! Mass rescaling for Centroid MD, not tested yet
      ! WARNING: PIMD with normal modes currently does not work,
      ! so there might be a bug in this code.
      ! NOTE: PI+GLE and PIGLET (inormalmodes==1) require physical masses!
      if (inormalmodes == 2) then
         do iat = 1, natom
            lambda(1) = 0
            do iw = 2, nwalk / 2
               lambda(2 * iw - 1) = 2 * nwalk * (1 - cos(2 * PI * (iw - 1) / nwalk))
               lambda(2 * iw - 2) = lambda(2 * iw - 1)
            end do
            nmodes = nwalk / 2
            odd = nwalk - 2 * nmodes ! 0 if even, 1 if odd
            if (odd == 0) then
               lambda(nwalk) = 4 * nwalk
            end if
            do iw = 2, nwalk
               amg(iat, iw) = am(iat) * lambda(iw) !potential
               amt(iat, iw) = am(iat) * lambda(iw) !kinetic
            end do
            amg(iat, 1) = 0.0D0
            amt(iat, 1) = am(iat)
            if (idebug > 0) then
               write (*, *) 'lambda', (lambda(iw), iw=1, nwalk)
               write (*, *) 'amg', (amg(iat, iw), iw=1, nwalk)
               write (*, *) 'amt', (amt(iat, iw), iw=1, nwalk)
               ! write(*,*)'amg2',(amg(iat,iw)*2*sin(pi*(iw-1)/nwalk),iw=1,nwalk)
            end if
         end do
      end if

   end subroutine initialize_pi_masses

   ! This routine transforms staging forces to cartesian forces, which are stored
   ! in trans matrices,values in fx,fy and fz matrices are NOT modified!!
   ! used in estimators.f90
   subroutine FQtoFX(fx, fy, fz, transfx, transfy, transfz)
      real(DP), intent(in) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: transfx(:, :), transfy(:, :), transfz(:, :)
      real(DP) :: sumx, sumy, sumz
      integer :: iat, iw

      do iat = 1, natom

         sumx = 0.0D0
         sumy = 0.0D0
         sumz = 0.0D0

         do iw = nwalk, 2, -1
            transfx(iat, iw) = fx(iat, iw) - (iw - 2.0D0) / (iw - 1.0D0) * fx(iat, iw - 1)
            transfy(iat, iw) = fy(iat, iw) - (iw - 2.0D0) / (iw - 1.0D0) * fy(iat, iw - 1)
            transfz(iat, iw) = fz(iat, iw) - (iw - 2.0D0) / (iw - 1.0D0) * fz(iat, iw - 1)

            sumx = sumx + transfx(iat, iw)
            sumy = sumy + transfy(iat, iw)
            sumz = sumz + transfz(iat, iw)
         end do

         transfx(iat, 1) = fx(iat, 1) - sumx
         transfy(iat, 1) = fy(iat, 1) - sumy
         transfz(iat, 1) = fz(iat, 1) - sumz

      end do
   end subroutine FQtoFX

   ! This routine transforms cartesian forces to staging forces, which are stored
   ! in fx matrices, values in fxab,fyab and fzab matrices are NOT modified!
   ! Called from force_clas()
   ! The classical force calculated in cartesian coordinates is transformed
   ! into staging coordinates.
   ! CAUTION, the formula 94 in Quantum simulations (Tuckermann) is incorrect
   ! dphi/dui on the output should be divided by P (nwalk) both for the i=1 and the others
   ! futhermore, the dphi/dx(i-1) should actually be dphi/du(i-1)
   ! See Tuckerman's lecture notes. The force is therefore transformed
   ! recursively. sumux-z contains the (i-1) force used for this purpose.
   subroutine FXtoFQ(fxab, fyab, fzab, fx, fy, fz)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(in) :: fxab(:, :), fyab(:, :), fzab(:, :)
      real(DP) :: sumux(nwalk), sumuy(nwalk), sumuz(nwalk)
      integer :: iat, iw

      do iat = 1, natom

         sumux = 0.0D0
         sumuy = 0.0D0
         sumuz = 0.0D0

         do iw = 1, nwalk
            fx(iat, 1) = fx(iat, 1) + fxab(iat, iw) / nwalk
            fy(iat, 1) = fy(iat, 1) + fyab(iat, iw) / nwalk
            fz(iat, 1) = fz(iat, 1) + fzab(iat, iw) / nwalk
            sumux(1) = sumux(1) + fxab(iat, iw)
            sumuy(1) = sumuy(1) + fyab(iat, iw)
            sumuz(1) = sumuz(1) + fzab(iat, iw)
         end do

         do iw = 2, nwalk
            sumux(iw) = sumux(iw) + fxab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumux(iw - 1)
            sumuy(iw) = sumuy(iw) + fyab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumuy(iw - 1)
            sumuz(iw) = sumuz(iw) + fzab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumuz(iw - 1)

            fx(iat, iw) = fx(iat, iw) + ( &
                        & fxab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumux(iw - 1) &
                        & ) / nwalk
            fy(iat, iw) = fy(iat, iw) + ( &
                        & fyab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumuy(iw - 1) &
                        & ) / nwalk
            fz(iat, iw) = fz(iat, iw) + ( &
                        & fzab(iat, iw) + (iw - 2.0D0) / (iw - 1.0D0) * sumuz(iw - 1) &
                        & ) / nwalk
         end do
      end do
   end subroutine FXtoFQ

   ! This routine transforms normal mode coordinates to cartesian coordinates,
   ! which are stored in trans matrices.
   ! Original arrays x, y, z are NOT modified.
   ! Masses are also modified.
   subroutine UtoX(x, y, z, transx, transy, transz)
      use mod_fftw3, only: dft_normalmode2cart
      use mod_general, only: nwalk, inormalmodes, idebug
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: transx(:, :), transy(:, :), transz(:, :)
      integer :: iat, iw
      real(DP) :: dnwalk, fac, equant
      integer :: nmodes, odd

      dnwalk = 1.0D0
      fac = 1.0D0
      if (inormalmodes == 1) then
         dnwalk = dsqrt(1.0D0 * nwalk)
         fac = dsqrt(2.0D0)
      end if

      nmodes = nwalk / 2
      odd = nwalk - 2 * nmodes ! 0 if even, 1 if odd

      if (idebug > 0) then
         call equant_nm(x, y, z, equant)
      end if

      do iat = 1, natom

         ! TODO: complex() is a GNU intrinsic, we need to figure out how to do this
         ! in a portable manner to support intel compilers.
#if __GNUC__ == 0
         call fatal_error(__FILE__, __LINE__, &
            & 'Normal mode transform not supported for non-GNU compilers')
#else
         cx(1) = complex(x(iat, 1), 0)
         cy(1) = complex(y(iat, 1), 0)
         cz(1) = complex(z(iat, 1), 0)
         do iw = 2, nmodes + odd
            cx(iw) = complex(x(iat, iw), x(iat, nwalk + 2 - iw)) / fac
            cy(iw) = complex(y(iat, iw), y(iat, nwalk + 2 - iw)) / fac
            cz(iw) = complex(z(iat, iw), z(iat, nwalk + 2 - iw)) / fac
         end do
         if (odd /= 1) then
            cx(nmodes + 1) = complex(x(iat, nmodes + 1), 0)
            cy(nmodes + 1) = complex(y(iat, nmodes + 1), 0)
            cz(nmodes + 1) = complex(z(iat, nmodes + 1), 0)
         end if
#endif

         call dft_normalmode2cart(cx, x_in)
         call dft_normalmode2cart(cy, y_in)
         call dft_normalmode2cart(cz, z_in)

         do iw = 1, nwalk
            transx(iat, iw) = x_in(iw) / dnwalk
            transy(iat, iw) = y_in(iw) / dnwalk
            transz(iat, iw) = z_in(iw) / dnwalk
         end do

      end do

      if (idebug > 0) then
         call equant_cart(transx, transy, transz, equant)
      end if
   end subroutine UtoX

   ! This routine transforms cartesian to normal coordinates, which are stored
   ! in trans matrices, values in x,y and z matrices are NOT modified!!
   ! see https://github.com/i-pi/i-pi/blob/2a09a6d652b1ffe5f485c4c078c1085db6fcf63a/ipi/utils/nmtransform.py
   subroutine XtoU(x, y, z, transx, transy, transz)
      use mod_fftw3, only: dft_cart2normalmode
      use mod_general, only: inormalmodes
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: transx(:, :), transy(:, :), transz(:, :)
      integer :: iat, iw
      real(DP) :: dnwalk, fac
      integer :: nmodes, odd

      dnwalk = 1.0D0
      fac = 1.0D0
      if (inormalmodes == 1) then
         dnwalk = dsqrt(1.0D0 * nwalk)
         fac = dsqrt(2.0D0)
      end if

      nmodes = nwalk / 2
      odd = nwalk - 2 * nmodes ! 0 if even, 1 if odd

      do iat = 1, natom
         do iw = 1, nwalk
            x_in(iw) = x(iat, iw)
            y_in(iw) = y(iat, iw)
            z_in(iw) = z(iat, iw)
         end do

         call dft_cart2normalmode(x_in, cx)
         call dft_cart2normalmode(y_in, cy)
         call dft_cart2normalmode(z_in, cz)

         ! TODO: realpart() and imagpart() are GNU intrinsics, we need to
         ! figure out a portable way to do this to support intel compilers.
#if __GNUC__ == 0
         call fatal_error(__FILE__, __LINE__, &
            & 'Normal mode transform not supported for non-GNU compilers')
#else
         transx(iat, 1) = realpart(cx(1))
         transy(iat, 1) = realpart(cy(1))
         transz(iat, 1) = realpart(cz(1))
         do iw = 2, nmodes + odd
            transx(iat, iw) = realpart(cx(iw)) * fac
            transx(iat, nwalk + 2 - iw) = imagpart(cx(iw)) * fac
            transy(iat, iw) = realpart(cy(iw)) * fac
            transy(iat, nwalk + 2 - iw) = imagpart(cy(iw)) * fac
            transz(iat, iw) = realpart(cz(iw)) * fac
            transz(iat, nwalk + 2 - iw) = imagpart(cz(iw)) * fac
         end do

         if (odd /= 1) then
            transx(iat, nmodes + 1) = realpart(cx((nwalk + 2) / 2))
            transy(iat, nmodes + 1) = realpart(cy((nwalk + 2) / 2))
            transz(iat, nmodes + 1) = realpart(cz((nwalk + 2) / 2))
         end if
#endif

      end do

      transx = transx / dnwalk
      transy = transy / dnwalk
      transz = transz / dnwalk

      ! write(*,*)'original cartesian to normal modes'
      ! call print_xyz_arrays(x/ang,y/ang,z/ang, natom, nwalk)
      ! write(*,*)'transformed coordinates to normal modes'
      ! call print_xyz_arrays(transx/ang,transy/ang,transz/ang, natom, nwalk)
   end subroutine XtoU

   subroutine equant_cart(x, y, z, equant)
      use mod_general, only: nwalk, natom, idebug, inormalmodes
      use mod_system, only: am
      use mod_nhc, only: temp
      use mod_utils, only: print_xyz_arrays
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: equant
      real(DP) :: equantx, equanty, equantz, omega_n
      integer :: iat, iw

      omega_n = NWALK * TEMP
      if (inormalmodes == 2) then
         omega_n = omega_n * dsqrt(NWALK * 1.D0)
      end if
      equantx = 0.0D0
      equanty = 0.0D0
      equantz = 0.0D0

      do iat = 1, natom
         do iw = 1, nwalk - 1
            equantx = equantx + 0.5D0 * am(iat) * omega_n**2 * (x(iat, iw) - x(iat, iw + 1))**2
            equanty = equanty + 0.5D0 * am(iat) * omega_n**2 * (y(iat, iw) - y(iat, iw + 1))**2
            equantz = equantz + 0.5D0 * am(iat) * omega_n**2 * (z(iat, iw) - z(iat, iw + 1))**2
         end do
         equantx = equantx + 0.5D0 * am(iat) * omega_n**2 * (x(iat, nwalk) - x(iat, 1))**2
         equanty = equanty + 0.5D0 * am(iat) * omega_n**2 * (y(iat, nwalk) - y(iat, 1))**2
         equantz = equantz + 0.5D0 * am(iat) * omega_n**2 * (z(iat, nwalk) - z(iat, 1))**2
         ! write(*,*) "Quantum energy per atom in cartesian coordinates"
         ! write(*,*) equantx, equanty, equantz, equantx + equanty + equantz
      end do

      equant = equantx + equanty + equantz

      if (idebug > 0) then
         write (*, *) "Quantum energy in cartesian coordinates"
         write (*, *) equantx, equanty, equantz, equant
         write (*, *) "Cartesian coordinates"
         call print_xyz_arrays(x, y, z, natom, nwalk)
      end if
   end subroutine equant_cart

   subroutine equant_nm(x, y, z, equant)
      use mod_const, only: PI
      use mod_general, only: nwalk, natom, idebug, inormalmodes
      use mod_system, only: am
      use mod_nhc, only: temp
      use mod_utils, only: print_xyz_arrays
      use mod_arrays, only: amg
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: equant
      real(DP) :: equantx, equanty, equantz, omega_n, omega(size(x, 2))
      integer :: iat, iw, k

      omega_n = NWALK * TEMP
      ! Tuckerman Hamiltonian
      if (inormalmodes == 2) then
         omega_n = omega_n * dsqrt(NWALK * 1.D0)
      end if
      equantx = 0.0D0
      equanty = 0.0D0
      equantz = 0.0D0

      do iw = 1, nwalk
         k = iw - 1
         omega(iw) = 2 * omega_n * sin(k * PI / nwalk)
      end do

      do iat = 1, natom
         do iw = 1, nwalk
            if (inormalmodes == 1) then
               equantx = equantx + 0.5D0 * am(iat) * omega(iw)**2 * (x(iat, iw))**2
               equanty = equanty + 0.5D0 * am(iat) * omega(iw)**2 * (y(iat, iw))**2
               equantz = equantz + 0.5D0 * am(iat) * omega(iw)**2 * (z(iat, iw))**2
            else
               equantx = equantx + 0.5D0 * amg(iat, iw) * omega_n**2 * (x(iat, iw))**2
               equanty = equanty + 0.5D0 * amg(iat, iw) * omega_n**2 * (y(iat, iw))**2
               equantz = equantz + 0.5D0 * amg(iat, iw) * omega_n**2 * (z(iat, iw))**2
            end if
         end do
!      write(*,*)"Quantum energy in normal modes coordinates"
!      write(*,*)equantx, equanty, equantz, equantx+equanty+equantz
      end do

      equant = equantx + equanty + equantz

      if (idebug > 0) then
         write (*, *) 'omegas', (omega(iw), iw=1, nwalk)
         write (*, *) "Quantum energy in normal modes coordinates"
         write (*, *) equantx, equanty, equantz, equant
         write (*, *) "Normal mode coordinates"
         call print_xyz_arrays(x, y, z, natom, nwalk)
      end if

   end subroutine equant_nm

end module mod_transform
