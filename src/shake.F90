module mod_shake
   use mod_const, only: DP
   use mod_utils, only: abinerror
   implicit none
   private
   public :: shake_init, shake_tol, nshake, shake, ishake1, ishake2
   real(DP), allocatable :: dshake(:)
#if ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6 ) || ( __GNUC__ > 4 )
   ! TODO: merge ishake1 and ishake2
   integer, allocatable :: ishake1(:), ishake2(:)
#else
   integer, parameter :: MAXSHAKE = 6000
   integer :: ishake1(maxshake), ishake2(maxshake)
#endif
   real(DP) :: shake_tol = 0.001D0
   integer :: nshake = 0
   save
contains

   subroutine find_hbonds(x, y, z, num_shake)
      use mod_const, only: ANG
      use mod_general, only: natom
      use mod_system, only: names
      use mod_utils, only: get_distance, count_atoms_by_name
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(out) :: num_shake
      real(DP) :: r, hbond_len_thr
      integer :: iat1, iat2, iw, ish

      ! DH fixed threshold for now, we should make it a parameter
      ! Or determine bonds in a better way
      hbond_len_thr = 2.0 * ANG
      ! Here we assume that each H atom has a single bond
      num_shake = count_atoms_by_name(names, 'H', natom)
      ! TODO: Not sure about this, but shake does not
      ! work with PIMD anyway at this point
      iw = 1

      ish = 0
      do iat1 = 1, natom

         if (trim(names(iat1)) /= 'H') cycle

         do iat2 = 1, natom

            ! TODO: Is there a nicer way to do it?
            if (iat1 == iat2) cycle

            r = get_distance(x, y, z, iat1, iat2, iw)
            if (r <= hbond_len_thr) then
               ish = ish + 1
               ishake1(ish) = iat1
               ishake2(ish) = iat2
               ! TODO: Maybe we should continue and check
               ! whether we find more than one bond?
               exit
            end if

            ! We expect each H to be bonded!
            if (iat2 == natom) then
               write (*, *) 'ERROR: Could not find bond for hydrogen atom! index = ', iat1
               call abinerror('find_bonds')
            end if
         end do

      end do

      if (ish /= nshake) then
         write (*, *) 'ERROR: Whoops, something wrong when determining H bonds!'
         call abinerror('find_hbonds')
      end if

      write (*, *) 'Found the following bonds to H atoms'
      do iat1 = 1, num_shake
         write (*, '(A2, A2, 2I6)') names(ishake1(iat1)), names(ishake1(iat1)), &
                                  & ishake1(iat1), ishake2(iat1)
      end do

   end subroutine find_hbonds

   subroutine shake_init(x, y, z)
      real(DP) x(:, :), y(:, :), z(:, :)
      real(DP) xi, yi, zi, xj, yj, zj
      integer :: ixshake, i, j

      ! TODO: This is a temporary swith for SHAKEing H atoms
      if (nshake == -1) then
         call find_hbonds(x, y, z, nshake)
      end if

      allocate (dshake(nshake))
      do ixshake = 1, nshake
         i = ishake1(ixshake)
         j = ishake2(ixshake)
         xi = x(i, 1)
         yi = y(i, 1)
         zi = z(i, 1)
         xj = x(j, 1)
         yj = y(j, 1)
         zj = z(j, 1)
         dshake(ixshake) = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
      end do
   end subroutine shake_init

   subroutine shake(x, y, z, px, py, pz, amt, iq, iv)
      use mod_general, only: natom, nwalk, inormalmodes
      use mod_system, only: am
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: px(:, :), py(:, :), pz(:, :)
      real(DP), intent(in) :: amt(:, :)
      integer, intent(in) :: iq, iv
      real(DP), allocatable, save :: vx(:, :), vy(:, :), vz(:, :)
      real(DP) :: mi, mj, mij, agama
      integer :: iat, iw, i, j, ixshake, iiter, maxcycle, itest
      real(DP) :: dijiter2, xij, yij, zij, rij2
      real(DP) :: xdotij, ydotij, zdotij, dot
      real(DP) :: xiiter, yiiter, ziiter, xjiter, yjiter, zjiter

      if (.not. allocated(vx)) then
         allocate (vx(natom, nwalk))
         allocate (vy(natom, nwalk))
         allocate (vz(natom, nwalk))
      end if

      if (iv == 1) then
         do iw = 1, nwalk
            do iat = 1, natom
               vx(iat, iw) = px(iat, iw) / amt(iat, iw)
               vy(iat, iw) = py(iat, iw) / amt(iat, iw)
               vz(iat, iw) = pz(iat, iw) / amt(iat, iw)
            end do
         end do
      end if

      maxcycle = 1000
      do iw = 1, nwalk

         if (iw == 2 .and. inormalmodes > 1) exit !shake only centroid variable i.e. first normal mode

! SHAKE algorithm implemented, iteratively solved until convergence
! criterion agama < shake_tol are met for all bonds.
! We stop program if we do not converge after 1000 steps.
! Velocities recalculated in a primitive way by method of finite differences.

         if (iq == 1) then

            do IIter = 1, maxcycle

               ! variable which controls convergence,0=converged, 1=not converged
               itest = 0

               do ixshake = 1, NShake
                  i = IShake1(ixshake)
                  j = IShake2(ixshake)
                  mi = am(i)
                  mj = am(j)
                  mij = ((1 / mi)) + (1 / (mj))
                  xiiter = x(i, iw)
                  yiiter = y(i, iw)
                  ziiter = z(i, iw)
                  xjiter = x(j, iw)
                  yjiter = y(j, iw)
                  zjiter = z(j, iw)

                  dijiter2 = (xiiter - xjiter)**2 + (yiiter - yjiter)**2 + (ziiter - zjiter)**2
                  agama = dijiter2 - dshake(ixshake)
                  agama = agama / (4 * dijiter2 * mij)

                  if (abs(agama) > shake_tol) then
                     itest = 1
                     ! Rapaport, p.275, here with masses(small modification)
                     x(i, iw) = xiiter - agama * (xiiter - xjiter) / mi
                     y(i, iw) = yiiter - agama * (yiiter - yjiter) / mi
                     z(i, iw) = ziiter - agama * (ziiter - zjiter) / mi
                     x(j, iw) = xjiter + agama * (xiiter - xjiter) / mj
                     y(j, iw) = yjiter + agama * (yiiter - yjiter) / mj
                     z(j, iw) = zjiter + agama * (ziiter - zjiter) / mj
                  end if

                  ! ixshake loop
               end do

               if (itest == 0) then
                  exit
               end if
               ! IITER loop
            end do

            if (iiter >= maxcycle) then
               write (*, *) 'Error: Shake not converged after', maxcycle, 'iterations. Exiting...'
               call abinerror('shake')
            end if

         end if

         ! velocity update
         if (iv == 1) then
            do IIter = 1, maxcycle

               ! variable which controls convergence, 0=converged, 1=not converged
               itest = 0

               do ixshake = 1, NShake
                  i = IShake1(ixshake)
                  j = IShake2(ixshake)
                  mi = am(i)
                  mj = am(j)
                  mij = ((1 / mi)) + (1 / (mj))

                  xij = x(j, iw) - x(i, iw)
                  yij = y(j, iw) - y(i, iw)
                  zij = z(j, iw) - z(i, iw)
                  xdotij = vx(j, iw) - vx(i, iw)
                  ydotij = vy(j, iw) - vy(i, iw)
                  zdotij = vz(j, iw) - vz(i, iw)

                  dot = xij * xdotij + yij * ydotij + zij * zdotij
                  rij2 = xij**2 + yij**2 + zij**2
                  agama = -dot / (2 * rij2 * mij)

                  if (abs(agama) > shake_tol) then
                     itest = 1
                     vx(i, iw) = vx(i, iw) - agama * xij / mi
                     vy(i, iw) = vy(i, iw) - agama * yij / mi
                     vz(i, iw) = vz(i, iw) - agama * zij / mi
                     vx(j, iw) = vx(j, iw) + agama * xij / mj
                     vy(j, iw) = vy(j, iw) + agama * yij / mj
                     vz(j, iw) = vz(j, iw) + agama * zij / mj
                  end if

                  ! ixshake loop
               end do

               if (itest == 0) then
                  exit
               end if
               ! IITER loop
            end do
            if (iiter >= maxcycle) then
               write (*, *) 'Velocity shake not converged after', maxcycle, 'iterations. Exiting...'
               call abinerror('shake')
            end if
         end if

      end do

      if (iv == 1) then
         do iw = 1, nwalk
            do iat = 1, natom
               px(iat, iw) = vx(iat, iw) * amt(iat, iw)
               py(iat, iw) = vy(iat, iw) * amt(iat, iw)
               pz(iat, iw) = vz(iat, iw) * amt(iat, iw)
            end do
         end do
      end if

      return
   end subroutine shake

end module mod_shake
