module mod_vinit
   use mod_const, only: DP
   implicit none
   private
   public :: vinit, scalevelocities, remove_rotations, remove_comvel
   contains
!------------------------------------------------------------------------
!    
! INITIAL VELOCITY DISTRIBUTION                
! =============================
!
! This subroutine returns a velocity vector VX/VY/VZ for an 
! atom of given MASS chosen randomly from a Gaussian velocity 
! distribution corresponding to a given temperature TEMP. The
! mean of this distribution is zero, the variance is given by
!
!      2     k_B*T
! sigma   = -------
!              m
!
! Original version is code f24 from the book by 
! M. P. Allen & D. J. Tildesley:
! "Computer simulation of Liquids", Appendix G
!
! TEMP      desired temperature (atomic units: E_h / k_B)
! MASS      mass of atom (atomic units: m_e)
! VX,VY,VZ  velocity vector (output)
! NATOM     number of atoms
! NOUT      FORTRAN output channel
!
! Adapted (real(DP), no angular velocities, mass)   B. Schmidt, Apr 6, 1995
!
!------------------------------------------------------------------------
SUBROUTINE vinit(TEMP, MASS, vx, vy, vz)
   USE mod_general,  ONLY: natom, pot, nwalk, my_rank
   USE mod_random,   ONLY: gautrg
   real(DP),intent(out)    :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: mass(:)
   real(DP)    :: rans(3*size(mass))
   real(DP)    :: TEMP,SIGMA
   integer     :: iw,iat,pom

   do iw=1,nwalk

      call gautrg(rans,natom*3,0,6)  

      pom=1
      do iat=1,natom
! Variance of distribution
         sigma = SQRT ( TEMP/MASS(iat) )
 
! Velocity vector of i-th atom
         vx(iat,iw) = sigma * rans(pom) 
         vy(iat,iw) = sigma * rans(pom+1) 
         vz(iat,iw) = sigma * rans(pom+2)
         pom=pom+3

         if(abs(vx(iat,iw)).gt.1d-2.or.abs(vy(iat,iw)).gt.1d-2.or.abs(vz(iat,iw)).gt.1d-2)then
            write(*,*)'WARNING: The initial velocity of atom',iat,'is very large.'
            write(*,*)'Your system might blow up. Maybe try different random seed.'
            stop 1
         end if

      end do

      if(pot.eq.'2dho')then
         vx(1,iw)=0.0d0
         vy(1,iw)=0.0d0
         vz(1,iw)=0.0d0
      endif

   end do

END subroutine vinit

 
! Scaling of velocities to correct temperature after COM velocity removal           
! or if restarting to different temperature
subroutine ScaleVelocities(vx,vy,vz)
   use mod_const,   only: autok
   use mod_general, only: natom, nwalk, my_rank
   use mod_system,  only: dime, f, conatom
   use mod_nhc,     only: scaleveloc, temp
   use mod_kinetic, only: ekin_v
   use mod_shake,   only: nshake
   real(DP),intent(out)    :: vx(:,:), vy(:,:), vz(:,:)
   real(DP)  :: ekin_mom, temp_mom, scal

   ekin_mom=ekin_v(vx, vy, vz)

   if(ekin_mom.gt.0.1d-10)then
      temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
   else
      temp_mom=0.0d0
   end if

   if (my_rank.eq.0) write(*,*)'Initial temperature (K):',temp_mom*autok

! TODO: pro normal modes nemusi nutne fungovat!
   if(scaleveloc.eq.1.and.temp_mom.gt.0.1e-10)then

      if (my_rank.eq.0) write(*,*) 'Scaling velocities to correct temperature.'
      scal = sqrt(temp/temp_mom)
      vx = vx * scal 
      vy = vy * scal
      vz = vz * scal

      ekin_mom=ekin_v(vx, vy, vz)
      temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
      if (my_rank.eq.0) write(*,*)'Temperature after scaling (K):',temp_mom*autok
   end if

end subroutine ScaleVelocities

subroutine remove_comvel(vx, vy, vz, mass, lremove)
   use mod_general,     only: natom, nwalk, my_rank
   real(DP),intent(inout)  :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: mass(:)
   ! This variable decides whether we actually remove the momenta
   logical, intent(in)     :: lremove
   real(DP)                :: vcmx, vcmy, vcmz, tm
   integer                 :: iat, iw


   do iw=1,nwalk

!     Velocity vector of center of mass
      vcmx = 0.d0
      vcmy = 0.d0
      vcmz = 0.d0
      tm   = 0.d0
      do iat=1,natom
         vcmx = vcmx + vx(iat,iw)*mass(iat)
         vcmy = vcmy + vy(iat,iw)*mass(iat)
         vcmz = vcmz + vz(iat,iw)*mass(iat)
         tm = tm + mass(iat)
      end do
      vcmx = vcmx / tm
      vcmy = vcmy / tm
      vcmz = vcmz / tm

!     Shift velocities such that momentum of center of mass is zero
      if(lremove)then
         if(my_rank.eq.0) write(*,*)'Removing center of mass velocity.'
         do iat=1,natom
            vx(iat,iw) = vx(iat,iw) - vcmx
            vy(iat,iw) = vy(iat,iw) - vcmy
            vz(iat,iw) = vz(iat,iw) - vcmz
         end do
      end if

   end do

END SUBROUTINE remove_comvel


subroutine remove_rotations(x, y, z, vx, vy, vz, masses, lremove)
   use mod_system,      only: dime
   use mod_general,     only: natom, nwalk, my_rank, pot
   real(DP),intent(in)     :: x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout)  :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: masses(:)
   ! This variable decides whether we actually remove the momenta
   logical, intent(in)     :: lremove
   real(DP) :: Lx, Ly, Lz, L_tot
   real(DP) :: xx, xy, xz, yy, yz, zz
   real(DP) :: I(0:8), Iinv(0:8)
   real(DP) :: Omx, Omy, Omz, Om_tot
   integer  :: iat, iw

   if (dime.eq.1.or.natom.eq.1) return


   ! It would probably be more correct to calculate angular momentum
   ! for the whole bead necklace, same with COM velocity
   do iw=1,nwalk

      ! Angular momentum
      Lx = 0.0d0 
      Ly = 0.0d0 
      Lz = 0.0d0
 
      do iat=1, natom
         Lx = Lx + (y(iat,iw)*vz(iat,iw) - z(iat,iw)*vy(iat,iw)) * masses(iat)
         Ly = Ly + (z(iat,iw)*vx(iat,iw) - x(iat,iw)*vz(iat,iw)) * masses(iat)
         Lz = Lz + (x(iat,iw)*vy(iat,iw) - y(iat,iw)*vx(iat,iw)) * masses(iat)
      end do
 
      L_tot = dsqrt(Lx**2+Ly**2+Lz**2)
 
      write(*,*)"Angular momenta:", Lx, Ly, Lz, L_tot
 
      ! tensor of inertia 
 
      xx = 0.0d0
      xy = 0.0d0
      xz = 0.0d0
      yy = 0.0d0
      yz = 0.0d0
      zz = 0.0d0
 
      do iat=1, natom
         xx = xx + x(iat,iw)*x(iat,iw) * masses(iat)
         xy = xy + x(iat,iw)*y(iat,iw) * masses(iat)
         xz = xz + x(iat,iw)*z(iat,iw) * masses(iat)
         yy = yy + y(iat,iw)*y(iat,iw) * masses(iat)
         yz = yz + y(iat,iw)*z(iat,iw) * masses(iat)
         zz = zz + z(iat,iw)*z(iat,iw) * masses(iat)
      end do
 
      I(0) = (yy+zz)
      I(1) = -xy
      I(2) = -xz
      I(3) = -xy
      I(4) = (xx+zz)
      I(5) = -yz
      I(6) = -xz
      I(7) = -yz
      I(8) = (xx+yy)
 
      ! inverse of tensor 
      call mat_inv_3x3(I, Iinv)
 
      ! rotation: angular velocity
      Omx = 0.0d0
      Omy = 0.0d0
      Omz = 0.0d0
 
      Omx = Iinv(0)*Lx + Iinv(1)*Ly + Iinv(2)*Lz
      Omy = Iinv(3)*Lx + Iinv(4)*Ly + Iinv(5)*Lz
      Omz = Iinv(6)*Lx + Iinv(7)*Ly + Iinv(8)*Lz
      Om_tot = dsqrt(Omx**2+Omy**2+Omz**2)
 
      ! We probably need to diagonalize I to get kinetic energy...
      !write(*,*)"Angular velocity:", Omx, Omy, Omz, Om_tot
      !Ekin = 0.5 * (Omx**2
      !write(*,*)"Angular kinetic energy:", 0.5 * Om* Om_tot**2
      !write(*,*)"Angular kinetic energy:", 0.5 * * Om_tot**2
 
      if(lremove)then
         if(my_rank.eq.0) write(*,*)'Removing angular momentum.'
         do iat =1, natom
            vx(iat,iw) = vx(iat,iw) - Omy*z(iat,iw) + Omz*y(iat,iw)
            vy(iat,iw) = vy(iat,iw) - Omz*x(iat,iw) + Omx*z(iat,iw)
            vz(iat,iw) = vz(iat,iw) - Omx*y(iat,iw) + Omy*x(iat,iw)
         end do
      end if

   end do

   CONTAINS
   ! brute force inverse of 3x3 matrix
   SUBROUTINE mat_inv_3x3(a, ainv)
      real(DP), intent(in)    :: a(0:8)
      real(DP), intent(out)   :: ainv(0:8)
      real(DP) :: det, minor_det(0:8)
 
      minor_det(0) = a(4) * a(8) - a(5) * a(7)
      minor_det(1) = a(3) * a(8) - a(5) * a(6)
      minor_det(2) = a(3) * a(7) - a(4) * a(6)
      minor_det(3) = a(1) * a(8) - a(2) * a(7)
      minor_det(4) = a(0) * a(8) - a(2) * a(6)
      minor_det(5) = a(0) * a(7) - a(1) * a(6)
      minor_det(6) = a(1) * a(5) - a(2) * a(4)
      minor_det(7) = a(0) * a(5) - a(2) * a(3)
      minor_det(8) = a(0) * a(4) - a(1) * a(3)
 
      det = a(0)* minor_det(0) - a(1) * minor_det(1) + a(2) * minor_det(2)
 
      ainv(0) =  minor_det(0) / det
      ainv(1) = -minor_det(3) / det
      ainv(2) =  minor_det(6) / det
      ainv(3) = -minor_det(1) / det
      ainv(4) =  minor_det(4) / det
      ainv(5) = -minor_det(7) / det
      ainv(6) =  minor_det(2) / det
      ainv(7) = -minor_det(5) / det
      ainv(8) =  minor_det(8) / det
   END SUBROUTINE mat_inv_3x3

END SUBROUTINE REMOVE_ROTATIONS


end module mod_vinit
