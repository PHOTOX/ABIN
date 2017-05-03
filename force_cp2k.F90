module mod_cp2k
   use iso_c_binding
   use mod_const,    only: DP
   implicit none

! We are actually connecting to C-interface, hence this mess
INTERFACE 
   SUBROUTINE CP2K_SET_POSITIONS(env_id, new_pos, sz)
      IMPORT :: C_DOUBLE
      INTEGER, VALUE ::  env_id, sz
      REAL(C_DOUBLE), dimension(1:sz),intent(in) :: new_pos
   END SUBROUTINE CP2K_SET_POSITIONS

   SUBROUTINE CP2K_CALC_ENERGY_FORCE(env_id)
      INTEGER, VALUE ::  env_id
   END SUBROUTINE CP2K_CALC_ENERGY_FORCE

   SUBROUTINE CP2K_DESTROY_FORCE_ENV(env_id)
      INTEGER, VALUE ::  env_id
   END SUBROUTINE CP2K_DESTROY_FORCE_ENV

   SUBROUTINE CP2K_GET_POTENTIAL_ENERGY(env_id, E)
      import :: C_DOUBLE
      INTEGER, VALUE ::  env_id
      REAL(C_DOUBLE), intent(out) :: E
   END SUBROUTINE CP2K_GET_POTENTIAL_ENERGY

   SUBROUTINE CP2K_GET_FORCES(env_id, f, sz)
      IMPORT :: C_DOUBLE
      INTEGER, VALUE ::  env_id, sz
      REAL(C_DOUBLE), dimension(1:sz), intent(out) :: f
   END SUBROUTINE CP2K_GET_FORCES

END INTERFACE

  public
  private :: pos, force, f_env_id
  integer :: f_env_id
  real(C_DOUBLE), dimension(:), pointer :: pos, force 
  save

CONTAINS

#ifdef CP2K
subroutine init_cp2k()
    use mod_general, only: natom, idebug
    use iso_c_binding, only: C_CHAR,c_null_char
#ifdef MPI
    include "mpif.h"
#endif
    integer :: ierr
    integer :: cp2k_mpicomm
    character(len=*, KIND=C_CHAR), parameter :: cp2k_input_file='cp2k.inp'//c_null_char
    character(len=*, KIND=C_CHAR), parameter :: cp2k_output_file='cp2k.out'//c_null_char

#ifdef MPI
    CALL cp2k_init()

    !call MPI_Comm_dup(MPI_COMM_WORLD, cp2k_mpicomm, ierr)
    !if (idebug.gt.0) write(*,*)'Created new communicator ', cp2k_mpicomm

    ! TODO: we might want to create different communicators for different beads
    !CALL cp2k_create_force_env_comm(f_env_id, cp2k_input_file, cp2k_output_file, cp2k_mpicomm)
    CALL cp2k_create_force_env(f_env_id, cp2k_input_file, cp2k_output_file)
#else

    CALL cp2k_init_without_mpi()

    CALL cp2k_create_force_env(f_env_id, cp2k_input_file, cp2k_output_file)

#endif

    write(*,*)'Created CP2K force environment. id=', f_env_id

    ALLOCATE(pos(natom*3),force(natom*3),stat=ierr)
    IF (ierr/=0) STOP "Could not allocate memory"

end subroutine init_cp2k

subroutine finalize_cp2k()
   deallocate( pos )
   deallocate( force )
   call cp2k_destroy_force_env(f_env_id)

#ifdef MPI
   CALL cp2k_finalize()
#else
   call cp2k_finalize_without_mpi()
#endif

end subroutine finalize_cp2k

! CP2K endif
#endif

subroutine force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general,  only: natom, iqmmm, idebug
      use mod_utils,    only: abinerror
      use mod_interfaces, only: oniom
      implicit none
      real(DP),intent(in)    :: x(:,:),  y(:,:),  z(:,:)
      real(DP),intent(out)   :: fx(:,:), fy(:,:), fz(:,:)
      real(DP),intent(out)   :: eclas
      integer, intent(in)    :: walkmax
      real(DP)  :: e0
      integer :: iat,iw, ind

     eclas=0.0d0

     do iw = 1, walkmax

       ind=1
       do iat=1,natom
        pos(ind)   = x(iat,iw) 
        pos(ind+1) = y(iat,iw) 
        pos(ind+2) = z(iat,iw) 
        ind=ind+3
       end do

#ifdef CP2K
     if (idebug.gt.0) write(*,*)'Setting positions into CP2K engine.', f_env_id
     call cp2k_set_positions(f_env_id, pos, natom*3) 

     if (idebug.gt.0) write(*,*)'CP2K engine now calculates forces and energies.'
     call cp2k_calc_energy_force(f_env_id)

     call cp2k_get_potential_energy(f_env_id, e0)

     call cp2k_get_forces(f_env_id, force, natom*3)
#endif

     eclas = eclas + e0

     ind = 1
     do iat = 1, natom
        fx(iat, iw) = force(ind)
        fy(iat, iw) = force(ind + 1)
        fz(iat, iw) = force(ind + 2)
        ind = ind + 3
     end do

     if (iqmmm.eq.1) call oniom(x, y, z, fx, fy, fz, eclas, iw)

    end do

end subroutine force_cp2k

end module mod_cp2k


