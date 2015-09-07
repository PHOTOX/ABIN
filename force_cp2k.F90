module mod_cp2k
  use mod_const,    only: DP, ANG
  public
  private :: pos, force
  integer :: f_env_id, ierr
  real(DP), dimension(:), pointer :: pos, force 


#ifdef CP2K
CONTAINS

subroutine cp2k_init()
    use mod_general, only: natom
    CALL cp_init_cp2k(1,ierr)
    IF (ierr/=0) STOP "init_cp2k"
    CALL cp_create_fenv(f_env_id,"cp2k.inp","cp2k.out",ierr)
    IF (ierr/=0) STOP "create_force_env"

    ALLOCATE(pos(natom*3),force(natom*3),stat=ierr)
      IF (ierr/=0) STOP "Could not allocate memory"

end subroutine cp2k_init

subroutine cp2k_finalize()
!   deallocate( pos )
!   deallocate( force )
   CALL cp_finalize_cp2k(1, ierr)
   IF (ierr/=0) STOP "finalize_cp2k"
end subroutine cp2k_finalize

subroutine force_cp2k(x,y,z,fx,fy,fz,eclas)
      use mod_general
      use mod_system,   only: names
      use mod_qmmm,     only: natqm
      use mod_utils,    only: abinerror, printf
      use mod_interfaces, only: oniom
      implicit none
      real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
      real(DP),intent(out)   :: fx(:,:),fy(:,:),fz(:,:)
      real(DP),intent(out)   :: eclas
      real(DP)  :: e0
      character(len=100) :: chsystem
      character(len=20) :: chgeom,chforce,chhess,fgeom
      logical :: file_exists
      integer :: iat,iw,iat1,iat2,itest 
      integer :: ist1,ist2,iost, ISTATUS, ierr, ind

     eclas=0.0d0

     do iw=1,nwalk
!     call cp_get_pos(f_env_id, pos, natom*3,ierr) 

!     do iat=1,natom*3
!     write(*,*)pos(iat)
!     end do

       ind=1
       do iat=1,natom
        pos(ind)   = x(iat,iw) 
        pos(ind+1) = y(iat,iw) 
        pos(ind+2) = z(iat,iw) 
        ind=ind+3
       end do

     call cp_set_pos(f_env_id, pos, natom*3,ierr) 
     IF (ierr/=0) STOP "Error when passing coordinates to CP2K"

     call cp_calc_energy_force(f_env_id, 1, ierr)
     IF (ierr/=0) STOP "Error when calculating CP2K forces"

     call cp_get_energy(f_env_id, e0, ierr)
     IF (ierr/=0) STOP "Error when getting CP2K energy."

     call cp_get_force(f_env_id, force, natom*3, ierr)
       IF (ierr/=0) STOP "Error when getting CP2K forces."

     eclas = eclas + e0
!     write(*,*)e0

     ind = 1
     do iat = 1, natom
        fx(iat, iw) = force(ind)
        fy(iat, iw) = force(ind + 1)
        fz(iat, iw) = force(ind + 2)
        ind = ind + 3
     end do

!     call printf(fx, fy, fz)

     if (iqmmm.eq.1) call oniom(x, y, z, fx, fy, fz, eclas, iw)

    end do

end subroutine force_cp2k

#endif

end module mod_cp2k


