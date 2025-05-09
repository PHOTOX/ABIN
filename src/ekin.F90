! initial version                         P. Slavicek a kol., Mar 25, 2009
module mod_kinetic
   use mod_const, only: DP
   implicit none
   private
   public :: temperature, ekin_v, est_temp_cumul, entot_cumul
   real(DP) :: est_temp_cumul = 0.0D0, entot_cumul = 0.0D0
   save
contains
   subroutine temperature(px, py, pz, amt, eclas)
      use mod_const, only: AUtoFS, AUtoK
      use mod_general, only: it, sim_time, natom, nwalk, nwrite, inormalmodes, ipimd
      use mod_system, only: dime, f, conatom
      use mod_files, only: UTEMPER, UENERGY
      use mod_nhc, only: inose, get_nhcham
      use mod_gle, only: get_langham
      use mod_shake, only: nshake
      integer :: iw, iat
      real(DP), intent(in) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      real(DP), intent(in) :: eclas
      real(DP) :: est_temp, temp1, ekin_mom

      ekin_mom = 0.0D0
      do iw = 1, nwalk
         do iat = 1, natom
            temp1 = px(iat, iw)**2 + py(iat, iw)**2 + pz(iat, iw)**2
            temp1 = 0.5D0 * temp1 / amt(iat, iw)
            ekin_mom = ekin_mom + temp1
         end do
      end do

      est_temp = 2 * ekin_mom / (dime * natom * nwalk - nshake - f - dime * conatom)
      if (inormalmodes == 1 .or. inose == 2) est_temp = est_temp / nwalk
      est_temp_cumul = est_temp_cumul + est_temp

      if (ipimd == 0) then
         entot_cumul = entot_cumul + eclas + ekin_mom
      end if

      if (modulo(it, nwrite) == 0) then

         ! Temperature and conserved quantities of thermostats
         write (UTEMPER, '(F15.2,F10.2)', advance="no") sim_time * AUtoFS, est_temp * autok
         if (ipimd /= 2) write (UTEMPER, '(F10.2)', advance="no") est_temp_cumul * autok / it

         if (inose == 1) then
            write (UTEMPER, '(E20.10)', advance="no") get_nhcham(natom, nwalk) + ekin_mom + eclas
         else if (inose == 2 .or. inose == 3 .or. inose == 4) then
            write (UTEMPER, '(E20.10)', advance="no") get_langham() + ekin_mom + eclas
         end if

         write (UTEMPER, *)

         if (ipimd /= 1) then
            ! Printing to file energies.dat
            write (UENERGY, '(F15.2,3E20.10)', advance='no') sim_time * AUtoFS, &
                                                           & eclas, ekin_mom, eclas + ekin_mom
            if (ipimd == 0) write (UENERGY, '(E20.10)', advance="no") entot_cumul / it
            write (UENERGY, *)
         end if

      end if
   end subroutine temperature

   ! WARNING: We assume physical masses here!
   real(DP) function ekin_v(vx, vy, vz)
      use mod_general, only: nwalk, natom
      use mod_system, only: am
      real(DP), intent(in) :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: temp1, ekin_mom
      integer :: iw, iat

      ekin_mom = 0.0D0

      do iw = 1, nwalk
         do iat = 1, natom
            temp1 = vx(iat, iw)**2 + vy(iat, iw)**2 + vz(iat, iw)**2
            temp1 = 0.5D0 * temp1 * am(iat)
            ekin_mom = ekin_mom + temp1
         end do
      end do

      ekin_v = ekin_mom
   end function ekin_v

end module mod_kinetic
