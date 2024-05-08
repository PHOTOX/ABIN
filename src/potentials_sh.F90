! Various analytical potentials for surface hopping:
! - NaI potential and couplings in adiabatic basis

! Created by Jiri Janos

! User-definer parameters are set input section 'system'
module mod_potentials_sh
   use mod_const, only: DP
   implicit none
   public
   private :: nai

   ! We use derived types to encapsulate potential parameters
   ! https://fortran-lang.org/learn/quickstart/derived_types

   ! Parameters for the NaI model as defined in the following articles
   ! https://doi.org/10.1063/1.4919780
   ! https://doi.org/10.1063/1.456377
   ! These articles contain diabatic basis while here we use adiabatic
   type :: nai_params
      ! Diabatic state VX parameters
      real(DP) :: a2 = 2760d0
      real(DP) :: b2 = 2.398d0
      real(DP) :: c2 = 11.3d0
      real(DP) :: lp = 0.408d0
      real(DP) :: lm = 6.431d0
      real(DP) :: rho = 0.3489d0
      real(DP) :: de0 = 2.075d0
      real(DP) :: e2 = 14.399613877582553d0 ! elementary charge in eV/angs units
      ! Diabatic state VA parameters
      real(DP) :: a1 = 0.813d0
      real(DP) :: beta1 = 4.08d0
      real(DP) :: r0 = 2.67d0
      ! Diabatic coupling VXA parameters
      real(DP) :: a12 = 0.055d0
      real(DP) :: beta12 = 0.6931d0
      real(DP) :: rx = 6.93d0
   end type
   type(nai_params) :: nai 
   save
contains

   ! NaI potential
   subroutine nai_init(natom, nwalk, ipimd, nstate)
      use mod_error, only: fatal_error
      use mod_utils, only: real_positive
      integer, intent(in) :: natom, nwalk, ipimd, nstate

      if (natom /= 2) then 
         call fatal_error(__FILE__, __LINE__, &
            & 'NaI potential is only for 2 particles.')
      end if

      if (nwalk /= 1) then
         call fatal_error(__FILE__, __LINE__, &
            & 'NaI model requires only 1 walker.')
      end if

      if (ipimd /= 2) then
         call fatal_error(__FILE__, __LINE__, &
            & 'NaI model works only with SH dynamics.')
      end if

      if (nstate /= 2) then
         call fatal_error(__FILE__, __LINE__, &
            & 'NaI model is implemented for 2 states only.')
      end if

   end subroutine nai_init

   subroutine force_nai(x, y, z, fx, fy, fz, eclas)
      use mod_sh, only: en_array, istate, nacx, nacy, nacz
      use mod_const, only: ANG, AUTOEV
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      real(DP) :: VX, VA, VXA, dVX, dVA, dVXA ! diabatic hamiltonian and its derivatives
      real(DP) :: E1, E2, d12, dE1, dE2 ! adiabatic hamiltonian and derivatives of energies
      real(DP) :: r, fr, dx, dy, dz 

      ! NOTE: The original potential is defined in diabatic basis. For ABIN's purposes, we need to transfer to adiabatic basis,
      ! which can be easily done for two states. However, the adiabatic formulas are very long and complicated to both code and
      ! read. Therefore, we calculate the diabatic quantities and their derivatives and then construct the adiabatic states and
      ! couplings from diabatic ones. It should be easier to read and also more efficient. Note that the diabatic potential is in eV
      ! and Angstrom so conversions are necessary.

      ! distance between atoms
      dx = x(2, 1) - x(1, 1)
      dy = y(2, 1) - y(1, 1)
      dz = z(2, 1) - z(1, 1)
      r = dx**2 + dy**2 + dz**2
      r = dsqrt(r)

      ! normalized distance
      dx = dx/r
      dy = dy/r
      dz = dz/r

      ! converting distance to Angstrom as the potential was done
      r = r/ANG

      ! calculating diabatic hamiltonian 
      VX = (nai%a2 + (nai%b2/r)**8)*dexp(-r/nai%rho) - nai%e2/r - nai%e2*(nai%lp + nai%lm)/2/r**4 - &
              nai%c2/r**6 - 2*nai%e2*nai%lm*nai%lp/r**7 + nai%de0
      VA = nai%a1*dexp(-nai%beta1*(r - nai%r0))
      VXA = nai%a12*dexp(-nai%beta12*(r - nai%rx)**2)
      ! calculating derivatives of diabatic hamiltonian
      dVX = -(nai%a2 + (nai%b2/r)**8)*dexp(-r/nai%rho)/nai%rho - 8.0d0*nai%b2**8/r**9*dexp(-r/nai%rho) + &
              14.0d0*nai%e2*nai%lm*nai%lp/r**8 + 6.0d0*nai%c2/r**7 + 2.0d0*nai%e2*(nai%lp + nai%lm)/r**5 + nai%e2/r**2
      dVA = -nai%beta1*VA
      dVXA = -2.0d0*nai%beta12*(r - nai%rx)*VXA
      ! converting them to a.u. as they are in eV or eV/Angstrom
      VX = VX/AUTOEV
      VA = VA/AUTOEV
      VXA = VXA/AUTOEV
      dVX = dVX/AUTOEV/ANG
      dVA = dVA/AUTOEV/ANG
      dVXA = dVXA/AUTOEV/ANG


      ! calculating derivatives of diabatic quantities
      ! adiabatic potentials
      E1 = (VA + VX)/2.0d0 - dsqrt((VA - VX)**2.0d0 + 4.0d0*VXA**2.0d0)/2.0d0
      E2 = (VA + VX)/2.0d0 + dsqrt((VA - VX)**2.0d0 + 4.0d0*VXA**2.0d0)/2.0d0
      ! nonadiabatic coupling vector in the reduced system
      d12 = -(VXA*(dVA - dVX) + (-VA + VX)*dVXA)/(VA**2.0d0 - 2.0d0*VA*VX + VX**2.0d0 + 4.0d0*VXA**2.0d0)
      d12 = d12/ANG ! one more conversion necessary
      ! derivatives of energies
      dE1 = (dVA + dVX)/2.0d0 - (2.0d0*(VA - VX)*(dVA - dVX) + 8.0d0*VXA*dVXA)/(4.0d0*dsqrt((VA - VX)**2.0d0 + 4.0d0*VXA**2.0d0))
      dE2 = (dVA + dVX)/2.0d0 + (2.0d0*(VA - VX)*(dVA - dVX) + 8.0d0*VXA*dVXA)/(4.0d0*dsqrt((VA - VX)**2.0d0 + 4.0d0*VXA**2.0d0))


      ! saving electronic energies for SH
      en_array(1) = E1
      en_array(2) = E2

      ! saving classical energy for Verlet
      eclas = en_array(istate)

      ! calculate force (-dE/dr) on the propagated state
      if (istate == 1) fr = -dE1
      if (istate == 2) fr = -dE2

      ! projecting the force on the cartesian coordinates
      fx(1, 1) = -fr * dx
      fx(2, 1) = -fx(1, 1)
      fy(1, 1) = -fr * dy
      fy(2, 1) = -fy(1, 1)
      fz(1, 1) = -fr * dz
      fz(2, 1) = -fz(1, 1)

      ! projecting coupling to the cartesian coordinates
      nacx(1, 1, 2) = -d12 * dx
      nacy(1, 1, 2) = -d12 * dy
      nacz(1, 1, 2) = -d12 * dz
      nacx(1, 2, 1) = -nacx(1, 1, 2)
      nacy(1, 2, 1) = -nacy(1, 1, 2)
      nacz(1, 2, 1) = -nacz(1, 1, 2)
      nacx(2, 1, 2) = d12 * dx
      nacy(2, 1, 2) = d12 * dy
      nacz(2, 1, 2) = d12 * dz
      nacx(2, 2, 1) = -nacx(2, 1, 2)
      nacy(2, 2, 1) = -nacy(2, 1, 2)
      nacz(2, 2, 1) = -nacz(2, 1, 2)

   end subroutine force_nai

end module mod_potentials_sh
