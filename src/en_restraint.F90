module mod_en_restraint
   use mod_const, only: DP, AUTOEV
   use mod_files, only: stderr, stdout
   use mod_error, only: fatal_error
   use mod_general, only: en_restraint
   implicit none
   public
   real(DP) :: en_diff, en_kk
   real(DP), allocatable :: fxr(:, :), fyr(:, :), fzr(:, :)
   character(len=200) :: restrain_pot = 'none'

contains

! ENERGY RESTRAINT subroutine
! Fixing the excitation energy of molecule - adapts forces based on knowledge
! of ground and excited state gradients

! Using:
! A) Lagrange multipliers (default) - using energy gradient as approximation for next step energy change
! Paper: On the Importance of Initial Conditions for Excited-State Dynamics, 10.1039/C8FD00088C

! B) Quadratic potential around target value (en_kk must be set)

   subroutine en_rest_init(natom)
      integer, intent(in) :: natom

      if (en_restraint == 1) then
         write (stdout, *) 'Energy restraint is ON(1): Using method of Lagrange multipliers.'
      else if (en_restraint == 2 .and. en_kk >= 0) then
         write (stdout, *) 'Energy restraint is ON(2): Using quadratic potential restraint.'
      else
         call fatal_error(__FILE__, __LINE__, &
            & 'en_restraint must be either 0, 1 (Lagrange multipliers) or 2 (umbrella, define en_kk)')
      end if

      allocate (fxr(natom, 2)) ! two states, not beads.
      allocate (fyr(natom, 2))
      allocate (fzr(natom, 2))
   end subroutine en_rest_init

   subroutine energy_restraint(x, y, z, px, py, pz, eclas)
      use mod_general, only: natom, nwalk, dt0, it, idebug 
      use mod_system, only: am
      use mod_sh, only: en_array
      use mod_terampi_sh, only: force_terash
      use mod_files, only: UERMD
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      ! TODO: Eclas is not modified in this routine, but probably should be
      real(DP), intent(inout) :: eclas
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), dimension(natom) :: fxgs, fygs, fzgs, fxes, fyes, fzes
      real(DP) :: eclasexc, eclasground, excE, deltaE, lambda, lsum, deltaEnext, convercrit, deltaD
      integer :: ios, iat, iat2, iw
      character(len=30) :: formt, chforce_ground, chforce_exc

      do iw = 1, nwalk

         if (restrain_pot == '_tera_') then
            call force_terash(x, y, z, fxr, fyr, fzr, eclasground)
            do iat = 1, natom
               fxgs(iat) = fxr(iat, 1)
               fygs(iat) = fyr(iat, 1)
               fzgs(iat) = fzr(iat, 1)
               fxes(iat) = fxr(iat, 2)
               fyes(iat) = fyr(iat, 2)
               fzes(iat) = fzr(iat, 2)
            end do
            eclasground = en_array(1)
            eclasexc = en_array(2)
         else

            ! Should be done as a separate call eventually..

            write (chforce_ground, '(A,I3.3)') 'engrad.ground.dat.', iw
            write (chforce_exc, '(A,I3.3)') 'engrad.exc.dat.', iw

            !-----READING energy of groud state (engrad.ground.dat)
            open (901, file=chforce_ground, status='OLD', iostat=ios, action='read')
            if (ios /= 0) then
               call fatal_error(__FILE__, __LINE__, 'Could not open file '//chforce_ground)
            end if
            read (901, *) eclasground

            ! Reading gradient of ground state
            do iat = 1, natom
               read (901, *, IOSTAT=ios) fxgs(iat), fygs(iat), fzgs(iat)
               if (ios /= 0) then
                  call fatal_error(__FILE__, __LINE__, &
                     & 'Fatal problem with reading gradients from file engrad.ground.dat')
               end if
               ! Conversion to forces
               fxgs(iat) = -fxgs(iat)
               fygs(iat) = -fygs(iat)
               fzgs(iat) = -fzgs(iat)
            end do
            close (901)

            ! Reading energy of excited state (engrad.exc.dat)
            open (901, file=chforce_exc, status='OLD', iostat=ios, action='read')
            if (ios /= 0) then
               call fatal_error(__FILE__, __LINE__, 'Could not open file '//chforce_exc)
            end if
            read (901, *) eclasexc

            ! Reading gradient of excited state
            do iat = 1, natom
               read (901, *, IOSTAT=ios) fxes(iat), fyes(iat), fzes(iat)
               if (ios /= 0) then
                  call fatal_error(__FILE__, __LINE__, &
                                   'Could not read gradients from file '//chforce_exc)
               end if
               ! Conversion to forces
               fxes(iat) = -fxes(iat)
               fyes(iat) = -fyes(iat)
               fzes(iat) = -fzes(iat)
            end do
            close (901)

            ! restraint_pot endif
         end if

         ! Energy difference
         excE = eclasexc - eclasground
         deltaE = (excE - en_diff)

         if (en_restraint == 1) then
            !======== A) Computing lagrange multiplier lambda =========

            !Iterative procedure
            convercrit = 100
            do while (convercrit > 0.00001)
               ! deltaE(t+dt) prediction - improves the accuracy
               deltaEnext = 0
               do iat = 1, natom
                  deltaD = px(iat, iw) * dt0 / am(iat)
                  deltaEnext = deltaEnext - px(iat, iw) * dt0 / am(iat) * (fxes(iat) - fxgs(iat))
                  deltaEnext = deltaEnext - py(iat, iw) * dt0 / am(iat) * (fyes(iat) - fygs(iat))
                  deltaEnext = deltaEnext - pz(iat, iw) * dt0 / am(iat) * (fzes(iat) - fzgs(iat))
               end do

               ! lambda computation

               lsum = 0
               do iat2 = 1, natom
                  lsum = lsum + 1 / (am(iat2)) * (fxes(iat2) - fxgs(iat2))**2
                  lsum = lsum + 1 / (am(iat2)) * (fyes(iat2) - fygs(iat2))**2
                  lsum = lsum + 1 / (am(iat2)) * (fzes(iat2) - fzgs(iat2))**2
               end do

               deltaEnext = deltaE + deltaEnext
               lambda = deltaEnext / (lsum * dt0 * dt0)
               ! For debugging or computing excitation energy on the fly uncomment:
               ! lambda = 0
               ! Applying new forces
               ! p = p0 + f * dt

               do iat = 1, natom
                  px(iat, iw) = px(iat, iw) + dt0 * lambda * (fxes(iat) - fxgs(iat))
                  py(iat, iw) = py(iat, iw) + dt0 * lambda * (fyes(iat) - fygs(iat))
                  pz(iat, iw) = pz(iat, iw) + dt0 * lambda * (fzes(iat) - fzgs(iat))
               end do

               convercrit = deltaEnext
               if (idebug > 1) then
                   write (*, *) 'deltaEnext', deltaEnext
               end if

            end do

            !Output to en_restraint.dat
            write (formt, '(A30)') '(I8,F16.8,E20.10,E20.10,F16.8)' 
            write (UERMD, fmt=formt) it, excE * AUTOEV, deltaE, deltaEnext, lambda

         else if (en_restraint == 2) then
            !======= B) Quadratic restraint =============

            ! Applying new forces

            do iat = 1, natom
               px(iat, iw) = px(iat, iw) + dt0 * en_kk * deltaE * (fxes(iat) - fxgs(iat))
               py(iat, iw) = py(iat, iw) + dt0 * en_kk * deltaE * (fyes(iat) - fygs(iat))
               pz(iat, iw) = pz(iat, iw) + dt0 * en_kk * deltaE * (fzes(iat) - fzgs(iat))
            end do

            ! TODO:
            eclas = eclas ! + quadratic_restraint_energy

            !Output to en_restraint.dat
            write (formt, '(A30)') '(I8,F16.8,E20.10,E20.10,F16.8)'
            write (UERMD, fmt=formt) it, excE * AUTOEV, deltaE, 0.0, 0.0

         end if

      end do

   end subroutine energy_restraint

   subroutine en_rest_finalize()
      deallocate (fxr, fyr, fzr)
   end subroutine en_rest_finalize

end module mod_en_restraint
