! Too simple steepest descent gradient minimization
! gamm = coefficient, by which we multiply forces
! gamm_thr = threshold for gamm, which ends the minimization
! when we overshoot and energy gets higher, we divide gamm by 2 and try again
! This is extremely naive and stupid implementation,
! DON'T USE THIS!
module mod_minimize
   use mod_const, only: DP
   private
   public :: minimize, gamm, gammthr
   real(DP) :: gamm = 20.D0, gammthr = 1D-10 !minthr=1e-15
   save

contains

   subroutine minimize(x, y, z, fx, fy, fz, eclas)
      use mod_const, only: ANG
      use mod_general, only: natom, nwrite, nwritex, nstep, pot
      use mod_system, only: names, conatom
      use mod_analysis, only: trajout
      use mod_interfaces, only: force_clas
      implicit none
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      real(DP) :: x_new(size(x, 1), size(x, 2)), y_new(size(x, 1), &
                                                       size(x, 2)), z_new(size(x, 1), size(x, 2))
      real(DP) :: fx_new(size(x, 1), size(x, 2)), fy_new(size(x, 1), &
                                                         size(x, 2)), fz_new(size(x, 1), size(x, 2))
      real(DP) :: eclas_new
      integer :: iat, iw, iter

      iw = 1
      open (100, file='minimize.dat')
      write (100, *) '#Iteration    gamma   eclas  deltaE '
      call force_clas(fx, fy, fz, x, y, z, eclas, pot)
      write (100, '(A5,F10.4,1E20.8)') '#  0  ', gamm, eclas

      do iat = 1, conatom
         fx(iat, iw) = 0.0D0
         fy(iat, iw) = 0.0D0
         fz(iat, iw) = 0.0D0
      end do

      do iter = 1, nstep

         do iat = 1, natom
            x_new(iat, iw) = x(iat, iw) + gamm * fx(iat, iw)
            y_new(iat, iw) = y(iat, iw) + gamm * fy(iat, iw)
            z_new(iat, iw) = z(iat, iw) + gamm * fz(iat, iw)
         end do

         call force_clas(fx_new, fy_new, fz_new, x_new, y_new, z_new, eclas_new, pot)

         if (eclas_new < eclas) then
            do iat = conatom + 1, natom
               x(iat, iw) = x_new(iat, iw)
               y(iat, iw) = y_new(iat, iw)
               z(iat, iw) = z_new(iat, iw)
               fx(iat, iw) = fx_new(iat, iw)
               fy(iat, iw) = fy_new(iat, iw)
               fz(iat, iw) = fz_new(iat, iw)
            end do

            write (100, '(I8,F10.4,2E20.8)') iter, gamm, eclas_new, eclas_new - eclas
            if (modulo(iter, nwrite) == 0) then
               write (*, '(I8,F10.4,2E20.8)') iter, gamm, eclas_new, eclas_new - eclas
            end if
            eclas = eclas_new

            if (modulo(iter, nwritex) == 0) then
               call trajout(x, y, z, iter)
            end if

         else

            ! OK, this is really stupid
            nstep = nstep + 1
            gamm = gamm / 2
            if (gamm < gammthr) then
               write (*, *) '#Gamma smaller than ', gammthr, 'Optimization stopped...'
               write (100, *) '#Gamma smaller than ', gammthr, 'Optimization stopped...'
               exit
            end if
            write (*, *) '#Energy rose.Reducing gamma to:', gamm
            write (100, *) '#Energy rose.Reducing gamma to:', gamm

         end if

         if (iter >= 100000) then
            write (*, *) '#Minimization did 100000 steps. Exiting...'
            exit
         end if

      end do
!---------------------------------------------------------------------------

      close (100)
      open (100, file='geom.mini.xyz')
      write (100, *) natom
      write (100, *) 'Step:', iter
      do iat = 1, natom
         write (100, *) names(iat), x(iat, iw) / ang, y(iat, iw) / ang, z(iat, iw) / ang
      end do
      close (100)

   end subroutine minimize

end module mod_minimize
