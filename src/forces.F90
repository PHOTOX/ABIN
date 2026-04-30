! A wrapper routine for getting forces and energies
! Ab-initio programs are called from force_abin routine
subroutine force_clas(fx, fy, fz, x, y, z, energy, chpot)
   use mod_const, only: DP
   use mod_general, only: natom, nwalk, istage, inormalmodes, &
                          pot, pot_ref
   use mod_system, only: conatom
   use mod_nhc, only: inose
   use mod_transform
   use mod_interfaces, only: force_wrapper
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
   real(DP) :: eclas, eclas_ref

   fx = 0.0D0
   fy = 0.0D0
   fz = 0.0D0
   fxab = 0.0D0
   fyab = 0.0D0
   fzab = 0.0D0

   eclas = 0.0D0

   ! Back stage transformation, Cartesian coordinates are kept in trans
   ! matrices (even if staging is OFF!)
   if (istage == 1) then
      call QtoX(x, y, z, transx, transy, transz)
   else if (inormalmodes > 0) then
      call UtoX(x, y, z, transx, transy, transz)
   else
      transx = x
      transy = y
      transz = z
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

   ! For reference potential (ring-polymer contraction not yet implemented)
   if (pot_ref /= '_none_' .and. chpot == pot) then
      ! fxab now holds the full potential,
      ! but we need the difference force on the output
      fx = fxab; fy = fyab; fz = fzab
      fxab = 0.0D0; fyab = 0.0D0; fzab = 0.0D0
      eclas_ref = 0.0D0

      ! Calculate reference (cheap potential)
      call force_wrapper(transx, transy, transz, fxab, fyab, fzab, eclas_ref, pot_ref, nwalk)

      fxab = fx - fxab
      fyab = fy - fyab
      fzab = fz - fzab
      ! fxab now holds the difference force
      ! we return the difference forces, but full energy
      ! eclas = eclas - eclas_ref
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
   use mod_general, only: natom, ipimd
   use mod_water, only: watpot
   use mod_force_h2o, only: force_h2o
   use mod_force_mm, only: force_mm
   use mod_sbc, only: force_sbc, isbc
   use mod_plumed, only: iplumed, force_plumed
   use mod_potentials, only: force_harmonic_rotor, force_harmonic_oscillator, force_morse, force_doublewell
   use mod_potentials_sh, only: force_nai
   use mod_splined_grid, only: force_splined_grid
   use mod_cp2k, only: force_cp2k
   use mod_shell_interface, only: force_abin
   use mod_force_tcpb, only: force_tcpb
   use mod_force_tera, only: force_tera
   use mod_terampi_sh, only: force_terash
   use mod_force_mace, only: force_mace
   implicit none
   ! allow(external-procedure) ! fortitude linter
   external :: force_water
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: e_pot
   character(len=*), intent(in) :: chpot
   integer, intent(in) :: walkmax
   real(DP) :: eclas

   eclas = 0.0D0
   ! Here we decide which forces we want.
   ! By default we call an external program in force_abin routine
   select case (chpot)
   case ("_mm_")
      call force_mm(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_mmwater_")
      call force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot)
   case ("_h2o_")
      call force_h2o(x, y, z, fx, fy, fz, eclas, natom, walkmax)
   case ("_splined_grid_")
      ! Only 1D spline grid supported at the moment
      call force_splined_grid(x, fx, eclas, walkmax)
   case ("_harmonic_rotor_")
      call force_harmonic_rotor(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_harmonic_oscillator_")
      call force_harmonic_oscillator(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_morse_")
      call force_morse(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_doublewell_")
      call force_doublewell(x, y, fx, fy, eclas, walkmax)
   case ("_cp2k_")
      call force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_tcpb_")
      call force_tcpb(x, y, z, fx, fy, fz, eclas, walkmax)
   case ("_nai_")
      call force_nai(x, y, z, fx, fy, fz, eclas)
   case ("_tera_")
      if (ipimd == 2 .or. ipimd == 5) then
         call force_terash(x, y, z, fx, fy, fz, eclas)
      else
         call force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
      end if
   case ("_mace_")
      call force_mace(x, y, z, fx, fy, fz, eclas, walkmax)
   case DEFAULT
      call force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
   end select

   ! Spherical harmonic potential
   if (isbc == 1) then
      call force_sbc(x, y, z, fx, fy, fz, walkmax)
   end if

   ! External forces from PLUMED library
   if (iplumed == 1) then
      call force_plumed(x, y, z, fx, fy, fz, eclas)
   end if

   e_pot = eclas
end subroutine force_wrapper

subroutine force_quantum(fx, fy, fz, x, y, z, amg, quantum_energy)
   use mod_const, only: DP
   use mod_general, only: nwalk, inormalmodes, istage, natom
   use mod_nhc, only: temp, inose
   implicit none
   real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
   real(DP), intent(in) :: amg(:, :)
   real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
   real(DP), intent(out) :: quantum_energy
   real(DP) :: ak(natom, nwalk)
   real(DP) :: equant
   integer :: i, j, kplus, kminus

   fx = 0.0D0
   fy = 0.0D0
   fz = 0.0D0

   ! Setting the quantum force constants
   ! ak is defined as m*P/(beta^2*hbar^2)
   ak = 0.0D0
   do j = 1, nwalk
      do i = 1, natom
         ak(i, j) = nwalk * amg(i, j) * TEMP**2
      end do
   end do

   ! for PI+GLE we have different hamiltonian
   if (inose == 2) then
      ak = nwalk * ak
   end if

   ! Tuckerman normal modes Hamiltonian
   ! Not tested
   if (inormalmodes == 2) then
      ak = NWALK * TEMP**2 * amg / dsqrt(nwalk * 1.0D0)
   end if

   ! This is the energy coming from the Path Integral harmonic forces
   equant = 0.0D0

   ! If the staging transformation is not used
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
            equant = equant + 0.5D0 * ak(j, i) * (x(j, i) - x(j, kplus))**2
            equant = equant + 0.5D0 * ak(j, i) * (y(j, i) - y(j, kplus))**2
            equant = equant + 0.5D0 * ak(j, i) * (z(j, i) - z(j, kplus))**2
         end do
      end do
   end if

   ! If staging or normal mode transformation is used
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

   quantum_energy = equant
end subroutine force_quantum

