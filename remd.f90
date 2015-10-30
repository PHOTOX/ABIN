
module mod_remd
  use mod_const, only: DP
  use mod_utils, only: abinerror
  implicit none
  include "mpif.h"
  public
  private :: swap_replicas
  integer :: nreplica, nswap=-1
  integer :: MAX_REPLICA=50
  real(DP) :: deltaT=-1
  save

CONTAINS
   !TODO: make general subroutine MPI_ERROR in utils.f90 and check every ierr
   !TODO: reading multiple geometries for each replica
   !TODO: writing output from all replicas...mainly energies
   subroutine remd_swap(x, y, z, px, py, pz, fxc, fyc, fzc, eclas)
   use mod_general, only: my_rank 
   use mod_nhc,     only: temp
   use mod_random,  only: vranf
   real(DP),intent(inout)   ::  x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout)   ::  px(:,:), py(:,:), pz(:,:)
   real(DP),intent(inout)   ::  fxc(:,:), fyc(:,:), fzc(:,:)
   real(DP),intent(inout)   ::  eclas
   real(DP) :: eclas_new, prob, ran(MAX_REPLICA), rn
   integer :: source, ierr, src, tag_en=10, tag_swap=1
   integer             :: status(MPI_STATUS_SIZE), i
   logical             :: lswap

!   eclas = 1.0d0 / (my_rank+1)

   ! Broadcast array of random numbers from rank 0 
   if (my_rank.eq.0)then
      call vranf(ran, nreplica-1, 0, 6)
!      write(*,*)'rans', (ran(i),i=1,nreplica)
   end if
   call MPI_Bcast(ran, MAX_REPLICA, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   if (my_rank.ne.0) rn = ran(my_rank)

   ! Main Swapping loop, starting with the second replica and going up
   do source=1, nreplica-1

      if(my_rank.eq.source-1)then

         call MPI_Recv(eclas_new, 1, MPI_REAL8, source, tag_en, MPI_COMM_WORLD, status, ierr)
         call MPI_Send(eclas, 1, MPI_DOUBLE_PRECISION, source, tag_en, MPI_COMM_WORLD, ierr)
         call MPI_Recv(lswap, 1, MPI_LOGICAL, source, tag_swap, MPI_COMM_WORLD, status,  ierr)

      else if (my_rank.eq.source)then
         ! Replica with the higher temperature is the master here
         call MPI_Send(eclas, 1, MPI_REAL8, source-1, tag_en, MPI_COMM_WORLD, ierr)
         call MPI_Recv(eclas_new, 1, MPI_DOUBLE_PRECISION, source-1, tag_en, MPI_COMM_WORLD, status, ierr)
         prob = (eclas-eclas_new)*(1/temp-1/(temp-deltaT))
         if(rn.lt.prob.or.eclas.lt.eclas_new)then
            lswap=.True.
         else
            lswap=.False.
         end if
!         write(*,*)'Probability, ran. number, lswap', prob, rn, lswap
         call MPI_Send(lswap, 1, MPI_LOGICAL, source-1,tag_swap, MPI_COMM_WORLD, ierr)

      end if

      if (my_rank.eq.source.or.my_rank.eq.source-1)then
         if(lswap)then
            eclas = eclas_new
            call swap_replicas(x, y, z, px, py, pz, fxc, fyc, fzc, source-1, source)
         end if
      end if

   end do 


   end subroutine remd_swap

   subroutine swap_replicas(x, y, z, px, py, pz, fxc, fyc, fzc, rank1, rank2)
   use mod_general, only: my_rank 
   use mod_nhc,     only: temp
   real(DP),intent(inout)  ::  x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout)  ::  px(:,:), py(:,:), pz(:,:)
   real(DP),intent(inout)  ::  fxc(:,:), fyc(:,:), fzc(:,:)
   integer ,intent(in)     ::  rank1, rank2
   real(DP),allocatable    ::  x_new(:,:), y_new(:,:), z_new(:,:)
   real(DP),allocatable    ::  px_new(:,:), py_new(:,:), pz_new(:,:)
   real(DP),allocatable    ::  fxc_new(:,:), fyc_new(:,:), fzc_new(:,:)
   real(DP)                ::  scal
   integer                 ::  reqs(18), stats(18)
   integer                 ::  status(MPI_STATUS_SIZE)

   integer :: tag_x=11, tag_y=12, tag_z=13
   integer :: tag_px=114, tag_py=115, tag_pz=116
   integer :: tag_fx=17, tag_fy=18, tag_fz=19
   integer :: ierr, dest, src, size1, size2

   ! First, rank1 sends data and rank2 receives
   size1=size(x(:,1))
   size2=size(x(1,:))
   allocate( x_new(size1, size2)   )
   allocate( y_new(size1, size2)   )
   allocate( z_new(size1, size2)   )
   allocate( px_new(size1, size2)  )
   allocate( py_new(size1, size2)  )
   allocate( pz_new(size1, size2)  )
   allocate( fxc_new(size1, size2) )
   allocate( fyc_new(size1, size2) )
   allocate( fzc_new(size1, size2) )

!   write(*,*)'pz ', pz, my_rank
   if (my_rank.eq.rank1)then
      dest = rank2
      call MPI_Send(x, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(y, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(z, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(px, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(py, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(pz, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fxc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fyc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fzc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr )
      call MPI_Recv(x_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(y_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(z_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(px_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(py_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(pz_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fxc_new ,size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fyc_new ,size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fzc_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr )
   else if (my_rank.eq.rank2)then
      dest = rank1
      call MPI_Recv(x_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(y_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(z_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(px_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(py_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(pz_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fxc_new ,size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fyc_new ,size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr )
      call MPI_recv(fzc_new, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(x, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_x, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(y, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_y, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(z, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_z, MPI_COMM_WORLD, status,  ierr )
      call MPI_Send(px, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_px, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(py, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_py, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(pz, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_pz, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fxc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fx, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fyc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fy, MPI_COMM_WORLD, status, ierr )
      call MPI_Send(fzc, size1*size2, MPI_DOUBLE_PRECISION, dest, tag_fz, MPI_COMM_WORLD, status, ierr )
   end if

!   call MPI_WaitAll(18, reqs, stats, ierr)
!   if(ierr.ne.0)then
!      write(*,*)'Unspecified MPI problem during swapping with code and rank', &
!      ierr, my_rank
!      call abinerror('swap_replicas')
!   end if

!   Now all ranks have required data, so they can swap them
   x=x_new; y=y_new; z=z_new

   ! Momenta are scaled according to the new temperature
   ! p_new = p_old * sqrt(T_new/T_old)
   if (my_rank.eq.rank1)then
      scal = dsqrt(temp/(temp+deltaT))
   else if (my_rank.eq.rank2)then
      scal = dsqrt(temp/(temp-deltaT))
   end if
   px = px_new*scal;  py = py_new*scal;  pz = pz_new*scal

   fxc=fxc_new; fyc=fyc_new; fzc=fzc_new

!   write(*,*)'pz_new ', pz, my_rank

   deallocate( x_new ); deallocate( y_new ); deallocate ( z_new )
   deallocate( px_new ); deallocate( py_new );  deallocate( pz_new )
   deallocate( fxc_new ); deallocate( fyc_new ); deallocate ( fzc_new )

   if(my_rank.eq.rank1) write(*,*)'Succesfull swap between replicas ', rank1, rank2

   end subroutine swap_replicas

   subroutine remd_init()
   use mod_general, only: pot, my_rank, ipimd
   integer :: ierr

   ! First, do sanity check
   if(nswap.lt.0)then
      write(*,*)'ERROR: nswap must be a positive integer!'
      call abinerror('remd_init')
   end if

   if(deltaT.lt.0)then
      write(*,*)'ERROR: deltaT must be a positive real number!'
      call abinerror('remd_init')
   end if
   
   if(pot.eq.'_tera_'.or.pot.eq.'_cp2k_')then
      call abinerror('remd_init')
   end if

   if(ipimd.gt.0)then
      write(*,*)'REMD is currently implemented only for classical dynamics (ipimd.eq.0)!'
      call abinerror('remd_init')
   end if

   ! determine number of replicas via MPI call
   call MPI_Comm_size(MPI_COMM_WORLD, nreplica, ierr)
   if (my_rank.eq.0) write(*,*)'Number of REMD replicas: ', nreplica

   if(nreplica.le.1)then
      write(*,*)'You cannot do REMD with just one replica!'
      call abinerror('remd_init')
   end if

   end subroutine remd_init

end module mod_remd

