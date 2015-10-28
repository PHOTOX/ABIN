
module mod_remd
  use mod_const, only: DP
  use mod_utils, only: abinerror
  implicit none
  include "mpif.h"
  public
!  private :: 
  integer :: iremd=0, nreplica, nswap=-1
  real(DP) :: deltaT=-1
  save

CONTAINS
   subroutine remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
   real(DP),intent(inout)   ::  x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout)   ::  px(:,:), py(:,:), pz(:,:)
   real(DP),intent(inout)   ::  fxc(:,:), fyc(:,:), fzc(:,:)
   real(DP),intent(inout)   ::  eclas

   end subroutine remd_swap

   subroutine remd_init()
   use mod_general, only: pot, my_rank
   integer :: ierr
   if(nswap.lt.0)then
      call abinerror('remd_init')
   end if
   if(deltaT.lt.0)then
      call abinerror('remd_init')
   end if
   ! determine number of replicas via MPI call
   call MPI_Comm_size(MPI_COMM_WORLD, nreplica, ierr)
   if (my_rank.eq.0) write(*,*)'Number of REMD replicas: ', nreplica
   if(nreplica.le.1)then
      write(*,*)'You cannot do REMD with just one replica!'
      call abinerror('remd_init')
   end if
!   if(pot.eq.'_tera_'.or.pot.eq.'_cp2k_')then
!      call abinerror('remd_init')
!   end if
   end subroutine remd_init

end module mod_remd

