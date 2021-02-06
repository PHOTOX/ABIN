! A wrapper routine for getting forces and energies
! Ab-initio programs are called from force_abin routine
subroutine force_clas(fx, fy, fz, x, y, z, energy, chpot)
   use mod_const, only: DP
   use mod_general, only: natom, nwalk, istage, inormalmodes, iqmmm, &
                          pot, pot_ref, idebug
   use mod_force_mm, only: force_LJ_Coulomb
   use mod_sbc, only: force_sbc, isbc !,ibag
   use mod_system, only: conatom
   use mod_nhc, only: inose
   use mod_transform
   use mod_interfaces, only: force_wrapper
   use mod_plumed, only: iplumed, force_plumed
   implicit none
   real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: energy
   character(len=*), intent(in) :: chpot
   real(DP) :: transx(size(x, 1), size(x, 2))
   real(DP) :: transy(size(y, 1), size(y, 2))
   real(DP) :: transz(size(z, 1), size(z, 2))
   real(DP) :: fxab(size(fx, 1), size(fx, 2))
   real(DP) :: fyab(size(fy, 1), size(fy, 2))
   real(DP) :: fzab(size(fz, 1), size(fz, 2))
   integer :: iat, iw
   real(DP) :: eclas

   do iw = 1, nwalk
      do iat = 1, natom
         fx(iat, iw) = 0.0D0
         fy(iat, iw) = 0.0D0
         fz(iat, iw) = 0.0D0
         fxab(iat, iw) = 0.0D0
         fyab(iat, iw) = 0.0D0
         fzab(iat, iw) = 0.0D0
      end do
   end do

   eclas = 0.0D0

   ! Back stage transformation, Cartesian coordinates are kept in trans
   ! matrices (even if staging is OFF!)
   ! TODO: Rename transx to cart_x, cart_y, cart_z
   if (istage == 1) then
      call QtoX(x, y, z, transx, transy, transz)
   else if (inormalmodes > 0) then
      if (idebug > 0) write (*, *) 'Transforming coordinates back to cartesian'
      call UtoX(x, y, z, transx, transy, transz)
   else
      do iat = 1, natom
         do iw = 1, nwalk
            transx(iat, iw) = x(iat, iw)
            transy(iat, iw) = y(iat, iw)
            transz(iat, iw) = z(iat, iw)
         end do
      end do
   end if

   ! TODO: wraping molecules back to the PBC box
   ! The original code was used only for NAB force field potential,
   ! but perhaps it could be useful elsewhere
   ! (e.g. if we start doing QM/MM with PBC.
   ! Don't forget that if we do this here, we would need to propagate
   ! this change back to x, y, z (e.g. normal modes or staging in PIMD)
   ! if (ipbc.eq.1) call wrap(transx, transy, transz)

   ! LET'S GET FORCES!
   ! The ab initio interface is still deeper in force_abin()
   call force_wrapper(transx, transy, transz, fxab, fyab, fzab, eclas, chpot, nwalk)

   ! TODO: Change all one-line ifs to multiline to fix code coverage.

   ! TODO: It would probably make sense to put the following additional forces
   ! and QM/MM inside the force_wrapper as well.
   ! Spherical harmonic potential
   if (isbc == 1) call force_sbc(transx, transy, transz, fxab, fyab, fzab)

   ! QMMM SECTION
   ! ONIOM method (iqmmm=1) is called in force_abin
   ! The following are not working at the moment
   if (iqmmm == 3) call force_LJ_Coulomb(transx, transy, transz, fxab, fyab, fzab, eclas)

   if (iplumed == 1) call force_plumed(transx, transy, transz, fxab, fyab, fzab, eclas)

   ! For reference potential and ring-polymer contraction
   if (pot_ref /= 'none' .and. chpot == pot) then
      ! fxab now holds the full potential,
      ! but we need the difference force on the output
      fx = fxab; fy = fyab; fz = fzab
      fxab = 0.0D0; fyab = 0.0D0; fzab = 0.0D0
      energy = eclas
      eclas = 0.0D0

      call force_wrapper(transx, transy, transz, fxab, fyab, fzab, eclas, pot_ref, nwalk)
      if (isbc == 1) call force_sbc(transx, transy, transz, fxab, fyab, fzab)
      if (iqmmm == 3) call force_LJ_Coulomb(transx, transy, transz, fxab, fyab, fzab, eclas)
      if (iplumed == 1) call force_plumed(transx, transy, transz, fxab, fyab, fzab, eclas)

      fxab = fx - fxab
      fyab = fy - fyab
      fzab = fz - fzab
      ! fxab now holds the difference force
      ! we return the difference forces, but full energy
!      eclas = energy - eclas
      eclas = energy
   end if

   ! TRANSFORMING FORCES FROM CARTESIAN TO STAGING or NORMAL MODE COORDS
   ! Forces are divided by nwalk inside the FXtoFQ routine
   if (istage == 1) then

      call FXtoFQ(fxab, fyab, fzab, fx, fy, fz)

   else if (inormalmodes > 0) then
      ! for PIGLET

      call XtoU(fxab, fyab, fzab, fx, fy, fz)

   else if (inose == 2) then
      ! for PI+GLE
      do iw = 1, nwalk
         do iat = 1, natom
            fx(iat, iw) = fxab(iat, iw)
            fy(iat, iw) = fyab(iat, iw)
            fz(iat, iw) = fzab(iat, iw)
         end do
      end do

   else

      ! no transformation, classical MD
      do iw = 1, nwalk
         do iat = 1, natom
            fx(iat, iw) = fxab(iat, iw) / nwalk
            fy(iat, iw) = fyab(iat, iw) / nwalk
            fz(iat, iw) = fzab(iat, iw) / nwalk
         end do
      end do

   end if

   ! Constraining atoms
   ! Warning, this kills energy conservation!
   if (conatom > 0) then
      do iw = 1, nwalk
         do iat = 1, conatom
            fx(iat, iw) = 0.0D0
            fy(iat, iw) = 0.0D0
            fz(iat, iw) = 0.0D0
         end do
      end do
   end if

   energy = eclas

end subroutine force_clas

subroutine force_wrapper(x, y, z, fx, fy, fz, e_pot, chpot, walkmax)
   use mod_const, only: DP
   use mod_interfaces, only: force_abin
   use mod_general, only: natom, ipimd
   use mod_water, only: watpot
   use mod_force_mm, only: force_LJ_Coulomb
   use mod_harmon, only: force_harmon, force_2dho, force_morse, force_doublewell
   use mod_splined_grid
   use mod_cp2k, only: force_cp2k
   use mod_terampi, only: force_tera
   use mod_terampi_sh, only: force_terash
   implicit none
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: e_pot
   character(len=*), intent(in) :: chpot
   integer, intent(in) :: walkmax
   real(DP) :: eclas

   eclas = 0.0D0
!  Here we decide which forces we want.
!  By default we call external program in force_abin routine
!  TODO: All keywords for internal potentials should begin and end by '_'
!  such as pot='_tera_'
   select case (chpot)
   case ("mm")
      call force_LJ_Coulomb(x, y, z, fx, fy, fz, eclas)
   case ("mmwater")
      call force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot)
   case ("splined_grid")
      ! Only 1D spline grid supported at the moment
      call force_splined_grid(x, fx, eclas)
   case ("harm")
      call force_harmon(x, y, z, fx, fy, fz, eclas)
   case ("2dho")
      call force_2dho(x, y, z, fx, fy, fz, eclas)
   case ("morse")
      call force_morse(x, y, z, fx, fy, fz, eclas)
   case ("doublewell")
      call force_doublewell(x, y, fx, fy, eclas)
   case ("_cp2k_")
      call force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_tera_")
      if (ipimd == 2 .or. ipimd == 4 .or. ipimd == 5) then
         call force_terash(x, y, z, fx, fy, fz, eclas)
      else
         call force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
      end if
   case DEFAULT
      call force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
      eclas = eclas / walkmax
   end select

   e_pot = eclas
end subroutine force_wrapper

subroutine force_quantum(fx, fy, fz, x, y, z, amg, energy)
   use mod_const, only: DP
   use mod_array_size
   use mod_general, only: nwalk, inormalmodes, istage, natom, idebug
   use mod_nhc, only: temp, inose
   implicit none
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(in) :: amg(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: energy
   real(DP) :: ak(size(x, 1), size(x, 2))
   real(DP) :: equant
   integer :: i, j, kplus, kminus

   fx = 0.0D0
   fy = 0.0D0
   fz = 0.0D0

!  TODO: we should not calculate ak params each step...
!  Setting the quantum force constants
!  ak is defined is m*P/(beta^2*hbar^2)
   do j = 1, nwalk
      do i = 1, natom
         ak(i, j) = nwalk * amg(i, j) * TEMP**2
      end do
   end do

!  for PI+GLE we have different hamiltonian
   if (inose == 2) then
      ak = nwalk * ak
   end if

!   if(inormalmodes.eq.2)then
!      ak = ak / dsqrt(nwalk*1.0d0)
!   end if
   ! Tuckerman normal modes Hamiltonian
   if (inormalmodes == 2) then
      ak = NWALK * TEMP**2 * amg / sqrt(nwalk * 1.0D0)
   end if

   equant = 0.0D0

!  If the staging transformation is not used
   if (istage == 0 .and. inormalmodes == 0) then
      do j = 1, natom
         do i = 1, nwalk
            kplus = i + 1
            kminus = i - 1
            if (i == 1) then
               kminus = nwalk
            end if
            if (i == nwalk) then
               kplus = 1
            end if
            fx(j, i) = (x(j, i) - x(j, kplus))
            fx(j, i) = fx(j, i) + (x(j, i) - x(j, kminus))
            fx(j, i) = -fx(j, i) * ak(j, i)
            fy(j, i) = (y(j, i) - y(j, kplus))
            fy(j, i) = fy(j, i) + (y(j, i) - y(j, kminus))
            fy(j, i) = -fy(j, i) * ak(j, i)
            fz(j, i) = (z(j, i) - z(j, kplus))
            fz(j, i) = fz(j, i) + (z(j, i) - z(j, kminus))
            fz(j, i) = -fz(j, i) * ak(j, i)
            equant = equant + 0.5 * ak(j, i) * (x(j, i) - x(j, kplus))**2
            equant = equant + 0.5 * ak(j, i) * (y(j, i) - y(j, kplus))**2
            equant = equant + 0.5 * ak(j, i) * (z(j, i) - z(j, kplus))**2
         end do
      end do
   end if

!  If the staging or normal mode transformation is used
   if (istage == 1 .or. inormalmodes > 0) then
      do j = 1, natom
         do i = 1, nwalk
            fx(j, i) = -ak(j, i) * x(j, i)
            fy(j, i) = -ak(j, i) * y(j, i)
            fz(j, i) = -ak(j, i) * z(j, i)
            equant = equant + 0.5D0 * ak(j, i) * x(j, i)**2
            equant = equant + 0.5D0 * ak(j, i) * y(j, i)**2
            equant = equant + 0.5D0 * ak(j, i) * z(j, i)**2
         end do
      end do
   end if

   energy = equant

   if (idebug == 1) then
      write (*, *) 'EQUANT', equant
   end if
end subroutine force_quantum

! Based on:
! Efficient stochastic thermostatting of path integral molecular dynamics
! J. Chem. Phys. 133, 124104 ?2010?
subroutine propagate_nm(x, y, z, px, py, pz, m, dt)
   use mod_const, only: DP, PI
   use mod_array_size
   use mod_general, only: nwalk, natom
   use mod_nhc, only: temp
   implicit none
   real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
   real(DP), intent(in) :: m(:, :), dt
   real(DP) :: omega(size(x, 2)), omega_n, om, omt, c, s
   real(DP) :: x_old(size(x, 1), size(x, 2))
   real(DP) :: y_old(size(x, 1), size(x, 2))
   real(DP) :: z_old(size(x, 1), size(x, 2))
   real(DP) :: px_old(size(px, 1), size(px, 2))
   real(DP) :: py_old(size(px, 1), size(px, 2))
   real(DP) :: pz_old(size(px, 1), size(px, 2))
   integer :: iat, iw, k

   ! expecting m = am = physical masses

   x_old = x
   y_old = y
   z_old = z
   px_old = px
   py_old = py
   pz_old = pz

   omega_n = TEMP * NWALK

   do iw = 2, nwalk
      k = iw - 1
      omega(iw) = 2 * omega_n * sin(k * PI / nwalk)
   end do

   ! First, propagate centroid
   iw = 1
   do iat = 1, natom
      X(iat, iw) = X(iat, iw) + dt * PX(iat, iw) / M(iat, iw)
      Y(iat, iw) = Y(iat, iw) + dt * PY(iat, iw) / M(iat, iw)
      Z(iat, iw) = Z(iat, iw) + dt * PZ(iat, iw) / M(iat, iw)
   end do

   ! eq 23 from J. Chem. Phys. 133, 124104 2010
   ! exact propagation of a free ring polymer in normal mode coordinates
   do iw = 2, nwalk
      om = omega(iw)
      omt = omega(iw) * dt
      c = cos(omt)
      s = sin(omt)
      do iat = 1, natom
         ! Propagate positions
         x(iat, iw) = x_old(iat, iw) * c &
                      + px_old(iat, iw) * s / m(iat, iw) / om
         y(iat, iw) = y_old(iat, iw) * c &
                      + py_old(iat, iw) * s / m(iat, iw) / om
         z(iat, iw) = z_old(iat, iw) * c &
                      + pz_old(iat, iw) * s / m(iat, iw) / om

         ! propagate momenta
         px(iat, iw) = px_old(iat, iw) * c &
                       - x_old(iat, iw) * s * m(iat, iw) * om
         py(iat, iw) = py_old(iat, iw) * c &
                       - y_old(iat, iw) * s * m(iat, iw) * om
         pz(iat, iw) = pz_old(iat, iw) * c &
                       - z_old(iat, iw) * s * m(iat, iw) * om

      end do

   end do

end subroutine propagate_nm
