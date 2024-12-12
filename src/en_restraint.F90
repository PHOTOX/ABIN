! Driver routines for Energy Restraint Molecular Dynamics (ERMD)
! J.Suchan, https://doi.org/10.1039/C8FD00088C
!
! Fixing the excitation energy of a molecule during dynamics
! by adapting forces based on knowledge of ground and excited state gradients.
! Originally developed to sample initial conditions excitable by CW laser.

! Usage:
! A) Lagrange multipliers (default, en_restraint=1) - using energy gradient as approximation for next step energy change
! B) Quadratic potential around target value (en_restraint=2, force constant en_kk must be set)

module mod_en_restraint
   use mod_const, only: DP, AUTOEV
   use mod_files, only: stderr, stdout
   use mod_error, only: fatal_error
   use mod_general, only: en_restraint
   implicit none
   private
   public :: en_rest_init, energy_restraint, en_rest_finalize
   public :: en_diff, en_kk, restrain_pot

   real(DP) :: en_diff, en_kk
   real(DP), allocatable :: fxr(:, :), fyr(:, :), fzr(:, :)
   character(len=200) :: restrain_pot = 'none'

contains

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
      use mod_shell_interface_private, only: read_energy, read_forces, open_engrad_file
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      ! TODO: Eclas is not modified in this routine, but probably should be
      real(DP), intent(inout) :: eclas
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), dimension(natom, 1) :: fxgs, fygs, fzgs, fxes, fyes, fzes
      real(DP) :: eclasexc, eclasground, excE, deltaE, lambda, lsum, deltaEnext, convercrit, deltaD
      integer :: iat, iat2, iw, ugs, ues
      character(len=30) :: chforce_ground, chforce_exc
      logical :: abort

      abort = .false.

      do iw = 1, nwalk
         if (restrain_pot == '_tera_') then
            call force_terash(x, y, z, fxr, fyr, fzr, eclasground)
            do iat = 1, natom
               fxgs(iat, 1) = fxr(iat, 1)
               fygs(iat, 1) = fyr(iat, 1)
               fzgs(iat, 1) = fzr(iat, 1)
               fxes(iat, 1) = fxr(iat, 2)
               fyes(iat, 1) = fyr(iat, 2)
               fzes(iat, 1) = fzr(iat, 2)
            end do
            eclasground = en_array(1)
            eclasexc = en_array(2)
         else
            write (chforce_ground, '(A,I3.3)') 'engrad.ground.dat.', iw
            write (chforce_exc, '(A,I3.3)') 'engrad.exc.dat.', iw

            ugs = open_engrad_file(chforce_ground, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not open file '//trim(chforce_ground))

            ues = open_engrad_file(chforce_exc, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not open file '//trim(chforce_exc))

            eclasground = read_energy(ugs, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not read energy from file '//trim(chforce_ground))

            call read_forces(fxgs, fygs, fzgs, natom, 1, ugs, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not read gradients from file '//trim(chforce_ground))

            eclasexc = read_energy(ues, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not read energy from file '//trim(chforce_exc))

            call read_forces(fxes, fyes, fzes, natom, 1, ues, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, 'Could not read gradients from file '//trim(chforce_exc))

            close (ugs, status='delete')
            close (ues, status='delete')
         end if

         ! Energy difference
         excE = eclasexc - eclasground
         deltaE = (excE - en_diff)

         if (en_restraint == 1) then
            !======== A) Computing lagrange multiplier lambda =========

            !Iterative procedure
            convercrit = 100
            do while (convercrit > 0.00001D0)
               ! deltaE(t+dt) prediction - improves the accuracy
               deltaEnext = 0
               do iat = 1, natom
                  deltaD = px(iat, iw) * dt0 / am(iat)
                  deltaEnext = deltaEnext - px(iat, iw) * dt0 / am(iat) * (fxes(iat, 1) - fxgs(iat, 1))
                  deltaEnext = deltaEnext - py(iat, iw) * dt0 / am(iat) * (fyes(iat, 1) - fygs(iat, 1))
                  deltaEnext = deltaEnext - pz(iat, iw) * dt0 / am(iat) * (fzes(iat, 1) - fzgs(iat, 1))
               end do

               ! lambda computation

               lsum = 0
               do iat2 = 1, natom
                  lsum = lsum + 1 / (am(iat2)) * (fxes(iat2, 1) - fxgs(iat2, 1))**2
                  lsum = lsum + 1 / (am(iat2)) * (fyes(iat2, 1) - fygs(iat2, 1))**2
                  lsum = lsum + 1 / (am(iat2)) * (fzes(iat2, 1) - fzgs(iat2, 1))**2
               end do

               deltaEnext = deltaE + deltaEnext
               lambda = deltaEnext / (lsum * dt0 * dt0)
               ! For debugging or computing excitation energy on the fly uncomment:
               ! lambda = 0
               ! Applying new forces
               ! p = p0 + f * dt

               do iat = 1, natom
                  px(iat, iw) = px(iat, iw) + dt0 * lambda * (fxes(iat, 1) - fxgs(iat, 1))
                  py(iat, iw) = py(iat, iw) + dt0 * lambda * (fyes(iat, 1) - fygs(iat, 1))
                  pz(iat, iw) = pz(iat, iw) + dt0 * lambda * (fzes(iat, 1) - fzgs(iat, 1))
               end do

               convercrit = deltaEnext
               if (idebug > 1) then
                  write (*, *) 'deltaEnext', deltaEnext
               end if

            end do

            !Output to en_restraint.dat
            write (UERMD, '(I8,F16.8,E20.10,E20.10,F16.8)') it, excE * AUTOEV, deltaE, deltaEnext, lambda

         else if (en_restraint == 2) then
            !======= B) Quadratic restraint =============

            ! Applying new forces

            do iat = 1, natom
               px(iat, iw) = px(iat, iw) + dt0 * en_kk * deltaE * (fxes(iat, 1) - fxgs(iat, 1))
               py(iat, iw) = py(iat, iw) + dt0 * en_kk * deltaE * (fyes(iat, 1) - fygs(iat, 1))
               pz(iat, iw) = pz(iat, iw) + dt0 * en_kk * deltaE * (fzes(iat, 1) - fzgs(iat, 1))
            end do

            ! TODO:
            eclas = eclas ! + quadratic_restraint_energy

            !Output to en_restraint.dat
            write (UERMD, '(I8,F16.8,E20.10,E20.10,F16.8)') it, excE * AUTOEV, deltaE, 0.0D0, 0.0D0

         end if

      end do

   end subroutine energy_restraint

   subroutine en_rest_finalize()
      if (allocated(fxr)) deallocate (fxr, fyr, fzr)
   end subroutine en_rest_finalize

end module mod_en_restraint
