! Module for PIMD estimators (energy, heat capacity)
module mod_estimators
   use mod_const, only: DP, AUtoFS
   use mod_files, only: UCV, UCVDCV, UESTENERGY
   implicit none
   real(DP) :: est_prim_cumul = 0.0D0, est_vir_cumul = 0.0D0
   real(DP) :: est_prim2_cumul = 0.0D0, est_prim_vir = 0.0D0, est_vir2_cumul = 0.0D0
   real(DP) :: cv_prim_cumul = 0.0D0, cv_vir_cumul = 0.0D0, cv_dcv_cumul = 0.0D0
   real(DP), allocatable :: cvhess_cumul(:)
   integer :: enmini = 100
   save
contains
   ! Expecting cartesian coordinates and forces!
   subroutine estimators(x, y, z, fxab, fyab, fzab, eclas)
      use mod_general, only: ipimd, pot, natom, nwalk, it, sim_time, &
                           & ncalc, nwrite, inormalmodes, imini, ihess, icv
      use mod_nhc, only: temp, inose
      use mod_system, only: am, dime
      use mod_potentials, only: hessian_harmonic_oscillator, hessian_morse, hessian_harmonic_rotor
      use mod_shake, only: nshake
      use mod_error, only: fatal_error
      use mod_arrays, only: hess
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      ! fxab array is classical force in cartesian coordinates
      real(DP), intent(in) :: fxab(:, :), fyab(:, :), fzab(:, :)
      real(DP), intent(in) :: eclas
      real(DP) :: xc(size(x, 1)), yc(size(x, 1)), zc(size(x, 1))
      real(DP) :: cvhess(size(x, 2)), dc(size(x, 1) * 3, size(x, 2))
      real(DP) :: est_vir, est_prim, cv_prim, cv_vir, cv_dcv
      integer :: iat, iw, ipom, iat1, iat2, nf
      real(DP) :: it2, itnc

      if (inormalmodes == 1 .and. inose == 1) then
         temp = temp / nwalk
      end if

      ! We calculate all quantities only every ncalc steps
      ! also we begin to accumulate energies only after first enmini steps to avoid
      ! large initial oscilations
      itnc = (it - enmini) / ncalc
      nf = dime * natom - nshake !degrees of freedom

      ! We begin to accumulate averages of heat capacities only after it > imini
      ! This auxiliary variable is for cumulative averaging if imini > 0
      it2 = (it - imini) / ncalc

      ! We cannot accumulate heat capacity without accumulating energy first
      if (enmini > imini .and. icv == 1) then
         call fatal_error(__FILE__, __LINE__, 'enmini > imini not allowed')
      end if

      ! CALCULATING PRIMITIVE ESTIMATORS
      est_prim = 0.0D0
      cv_prim = 0.0D0

      do iat = 1, natom
         do iw = 1, nwalk - 1
            est_prim = est_prim + am(iat) * (x(iat, iw) - x(iat, iw + 1))**2
            est_prim = est_prim + am(iat) * (y(iat, iw) - y(iat, iw + 1))**2
            est_prim = est_prim + am(iat) * (z(iat, iw) - z(iat, iw + 1))**2
         end do
         est_prim = est_prim + am(iat) * (x(iat, nwalk) - x(iat, 1))**2
         est_prim = est_prim + am(iat) * (y(iat, nwalk) - y(iat, 1))**2
         est_prim = est_prim + am(iat) * (z(iat, nwalk) - z(iat, 1))**2
      end do

      cv_prim = (nwalk * temp * temp**2) * est_prim

      est_prim = nf * nwalk * temp * 0.5D0 - 0.5D0 * nwalk * temp**2 * est_prim + eclas
      est_prim_cumul = est_prim_cumul + est_prim

      if (icv == 1 .and. itnc > 0) then
         est_prim2_cumul = est_prim2_cumul + est_prim * est_prim
         cv_prim = (1 / temp**2) * (est_prim2_cumul / itnc - (est_prim_cumul / itnc)**2 + &
                                    0.5D0 * nf * nwalk * temp**2 - cv_prim)
         cv_prim_cumul = cv_prim_cumul + cv_prim
      end if
      ! END OF PRIMITIVE ESTIMATORS

      ! CALCULATING VIRIAL ESTIMATORS
      est_vir = 0.0D0
      cv_vir = 0.0D0

      ! Calculating centroids
      do iat = 1, natom
         xc(iat) = 0.0D0
         yc(iat) = 0.0D0
         zc(iat) = 0.0D0
         do iw = 1, nwalk
            xc(iat) = xc(iat) + x(iat, iw)
            yc(iat) = yc(iat) + y(iat, iw)
            zc(iat) = zc(iat) + z(iat, iw)
         end do
         xc(iat) = xc(iat) / nwalk
         yc(iat) = yc(iat) / nwalk
         zc(iat) = zc(iat) / nwalk
      end do

      do iat = 1, natom
         do iw = 1, nwalk
            est_vir = est_vir - (x(iat, iw) - xc(iat)) * fxab(iat, iw)
            est_vir = est_vir - (y(iat, iw) - yc(iat)) * fyab(iat, iw)
            est_vir = est_vir - (z(iat, iw) - zc(iat)) * fzab(iat, iw)
         end do
      end do

      if (inose == 2 .or. inormalmodes == 1) then
         est_vir = est_vir / nwalk
      end if

      est_vir = 0.5D0 * est_vir + nf * temp * 0.5D0 + eclas
      est_vir_cumul = est_vir_cumul + est_vir

      if (icv == 1 .and. itnc > 0) then

         est_vir2_cumul = est_vir2_cumul + est_vir * est_vir
         est_prim_vir = est_prim_vir + est_prim * est_vir
         cv_vir = est_prim_vir / itnc - (est_vir_cumul / itnc)**2
         cv_vir = cv_vir / temp**2 + 0.5D0 * nf
         cv_vir_cumul = cv_vir_cumul + cv_vir

         ! DOUBLE centroid virial Cv estimator
         ! external hessian is read in force_abin because of parallelization and through mod_estimators
         if (ihess == 1) then

            if (pot == '_harmonic_rotor_') then
               call hessian_harmonic_rotor(x, y, z, nwalk, hess)
            else if (pot == '_morse_') then
               call hessian_morse(x, y, z, nwalk, hess)
            else if (pot == '_harmonic_oscillator_') then
               call hessian_harmonic_oscillator(nwalk, hess)
            end if

            cv_dcv = 0.0D0
            do iw = 1, nwalk
               cvhess(iw) = 0.0D0
               ipom = 0
               do iat = 1, natom * 3, +3
                  ipom = ipom + 1
                  dc(iat, iw) = x(ipom, iw) - xc(ipom)
                  dc(iat + 1, iw) = y(ipom, iw) - yc(ipom)
                  dc(iat + 2, iw) = z(ipom, iw) - zc(ipom)
               end do
            end do
            ! summation over hessian elements etc.
            do iw = 1, nwalk

               do iat = 1, natom
                  cvhess(iw) = cvhess(iw) - (x(iat, iw) - xc(iat)) * fxab(iat, iw) * 1.5D0
                  cvhess(iw) = cvhess(iw) - (y(iat, iw) - yc(iat)) * fyab(iat, iw) * 1.5D0
                  cvhess(iw) = cvhess(iw) - (z(iat, iw) - zc(iat)) * fzab(iat, iw) * 1.5D0
               end do
               ! PIGLE --- different hamiltonian, fxc not divided by nwalk
               if (inose == 2 .or. inormalmodes == 1) then
                  cvhess(iw) = cvhess(iw) / nwalk
               end if

               do iat1 = 1, natom * 3
                  do iat2 = 1, natom * 3
                     cvhess(iw) = cvhess(iw) + &
                                & 0.5D0 * dc(iat1, iw) * dc(iat2, iw) * hess(iat1, iat2, iw)
                  end do
               end do
               cvhess_cumul(iw) = cvhess_cumul(iw) + cvhess(iw)
               cv_dcv = cv_dcv + cvhess_cumul(iw) / itnc

               ! summation enddo
            end do

            cv_dcv = est_vir2_cumul / itnc - (est_vir_cumul / itnc)**2 + &
                     0.5D0 * nf * temp**2 - 0.5D0 * temp * cv_dcv
            cv_dcv = cv_dcv / temp**2
            cv_dcv_cumul = cv_dcv_cumul + cv_dcv

            ! ihess endif
         end if
         ! icv endif
      end if

      ! END OF VIRIAL ESTIMATORS CALCULATIONS

      if (imini >= it) then
         cv_prim_cumul = 0.0D0
         cv_vir_cumul = 0.0D0
         cv_dcv_cumul = 0.0D0
         ! dirty hack to avoid NaN when it.eq.imin
         it2 = 1
      end if

      if (enmini >= it) then
         cv_prim_cumul = 0.0D0
         cv_vir_cumul = 0.0D0
         cv_dcv_cumul = 0.0D0
         est_prim_cumul = 0.0D0
         est_vir_cumul = 0.0D0
         est_prim2_cumul = 0.0D0
         est_prim_vir = 0.0D0
         est_vir2_cumul = 0.0D0
         if (ihess == 1) cvhess_cumul = 0.0D0
         ! dirty hack to avoid NaN when it.eq.enmin
         itnc = 1
      end if

      if (modulo(it, nwrite) == 0) then

         if (icv == 1 .and. enmini < it) then
            write (UCV, '(F15.2,4E20.10)') sim_time * AUtoFS, &
                                         & cv_prim, cv_vir, cv_prim_cumul / it2, cv_vir_cumul / it2

            if (ihess == 1) then
               write (UCVDCV, '(F15.2,2E20.10)') sim_time * AUtoFS, cv_dcv, cv_dcv_cumul / it2
            end if

            ! icv endif
         end if

         if (ipimd == 1) then
            write (UESTENERGY, '(F15.2,5E20.10)') sim_time * AUtoFS, eclas, &
                                                & est_prim, est_vir, &
                                                & est_prim_cumul / itnc, est_vir_cumul / itnc
         end if

      end if

      if (inormalmodes == 1 .and. inose == 1) then
         temp = temp * nwalk
      end if

   end subroutine estimators

end module mod_estimators
