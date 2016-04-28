module mod_terampi
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD.
! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
! So far, I was not able to make it work with OpenMPI.
! (but now that we use file based tera_port, it should work as well)
!
! Currently supports:
! pure QM and ONIOM
! Adapted from Sander (Amber14)
!
! Original Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
! Modified by Dan Hollas hollas@vscht.cz
! Date: September 2014 - September 2015
! ----------------------------------------------------------------
   use mod_const, only: DP
   implicit none
   private
   public :: teraport, finalize_terachem, initialize_terachem, &
          connect_terachem, force_tera, newcomm, &
          qmcharges, qmcoords, dxyz_all ! for terash
   ! This does not work at this moment.
   character*50  ::  teraport = 'terachem_port'
   integer     ::  newcomm ! Communicator, initialized in mpi_init subroutine
!  DH WARNING, initial hack, we do not support TeraChem-based QM/MM yet
   integer, parameter     ::  natmm_tera=0
   character(len=2), allocatable :: names_qm(:)
   real(DP), allocatable :: qmcharges(:)  ! QM charges from population analysis
   real(DP), allocatable :: mmcharges(:)  ! QM charges from population analysis
   real(DP), allocatable :: qmcoords(:,:)
   real(DP), allocatable :: mmcoords(:,:) 
   real(DP), allocatable :: dxyz_all(:,:)
   save

contains

subroutine force_tera(x, y, z, fx, fy, fz, eclas)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, iqmmm, nwalk, DP
!   use mod_system, only: names
   use mod_qmmm, only: natqm
   use mod_utils, only: printf, abinerror
   use mod_interfaces, only: oniom
   include 'mpif.h'
   real(DP),intent(in)      ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   ::  eclas
!   integer, intent(in) :: qmtypes(natqm)    ! QM atom types (nuclear charge in au)
   real(DP) :: escf                 ! SCF energy
   real(DP) :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
   real(DP) :: dummy=1.0d0
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iw, iat

! DHnote: we cannot use Niklasson's propagator if nwalk > 1
! This is the responsibility of the user

! For parallel terachem servers, we would need to read multiple port.txt
! files and newcomm would be an array. That's not implemented at this moment.
!!$ nteraservers=1
!!$ call OMP_SET_NESTED(TRUE)
!!$ if (OMP_GET_NESTED().ne.TRUE.and.parallel_qmmm.eq.1)then
!!$   write(*,*)'Nested parallelism is not supported by this compiler.'
!!$   write(*,*)'Please set parallel_qmmm=0.'
!!$   call abinerror('force_tera')
!!$ end if

!!$ OMP PARALLEL DO PRIVATE(qmcoords, mmcoords, ierr, iat, dxyz_all, escf) & !qmcharges,mmcharges and dipmom are not used at this point 
!!$     IF(nwalk.gt.1.and.nteraservers.gt.1) NUM_THREADS(nteraservers)
   do iw=1,nwalk

!!$OMP PARALLEL SECTIONS IF(parallel_qmmm.eq.true) NUM_THREADS(2)
!!$OMP_SECTION

      do iat=1,natqm
         qmcoords(1,iat) = x(iat,iw)/ANG
         qmcoords(2,iat) = y(iat,iw)/ANG
         qmcoords(3,iat) = z(iat,iw)/ANG
      end do
      do iat=1,natmm_tera
         mmcoords(1,iat) = x(iat+natqm,iw)
         mmcoords(2,iat) = y(iat+natqm,iw)
         mmcoords(3,iat) = z(iat+natqm,iw)
      end do

   ! -----------------------------------------
   ! Begin sending data each step to terachem
   ! -----------------------------------------


   ! Send natqm and the type of each qmatom
   if ( idebug > 1 ) then
      write(6,'(/, a, i0)') 'Sending natqm = ', natqm
      call flush(6)
   end if
   call MPI_SSend( natqm, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

   if ( idebug > 1 ) then
      write(6,'(/,a)') 'Sending QM atom types: '
      write(*,*)(names_qm(iat), iat=1,natqm)
      call flush(6)
   end if
   call MPI_SSend( names_qm, 2*size(names_qm), MPI_CHARACTER, 0, 2, newcomm, ierr )

   ! Send QM coordinate array
   if ( idebug > 1 ) then
      write(6,'(a)') 'Sending QM coords: '
   end if
   if ( idebug > 1) then
      do iat=1, natqm
         write(6,*) 'Atom ',iat,': ',qmcoords(:,iat)
         call flush(6)
      end do 
   end if
   call MPI_SSend( qmcoords, size(qmcoords), MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

if(natmm_tera.gt.0)then
   ! Send natmm and the charge of each atom
   if ( idebug > 1 ) then
      write(6,'(a, i0)') 'Sending natmm = ', natmm_tera
      call flush(6)
   end if
   call MPI_Send( natmm_tera, 1, MPI_INTEGER, 0, 2, newcomm, ierr ) 

   if ( idebug > 1 ) then
      write(6,'(a)') 'Sending charges: '
   end if
   call MPI_Send( mmcharges, natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

   ! Send MM point charge coordinate array
   if ( idebug > 1 ) then
      write(6,'(a)') 'Sending charges coords: '
   end if

   call MPI_Send( mmcoords, 3*natmm_tera, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 
end if

   ! -----------------------------------
   ! Begin receiving data from terachem
   ! -----------------------------------

   ! DH TODO: we need to somehow make sure that we don't wait forever if
   ! terachem crashes
   ! At this moment, this is ensured at the bash script level

   ! Energy
   if ( idebug > 2 ) then
      write(6,'(a)') 'Waiting to receive scf energy from TeraChem...'
      call flush(6)
   end if
   call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   ! Checking for TAG=1, which means that SCF did not converge
   if (status(MPI_TAG).eq.1)then
      write(*,*)'GOT TAG 1 from TeraChem: SCF probably did not converge.'
      call abinerror('force_tera')
   end if
   if ( idebug > 1 ) then
      write(6,'(a,es15.6)') 'Received scf energy from server:', escf
      call flush(6)
   end if

   ! Charges (Mulliken or other)
   if ( idebug > 2 ) then
      write(6,'(a)') 'Waiting to receive charges...'
   end if
   call MPI_Recv( qmcharges(:), natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 2 ) then
      write(6,'(a)') 'Received the following charges from server:'
      do iat=1, natqm
         write(6,*) 'Atom ',iat, ': ', qmcharges(iat)
      end do
      call flush(6)
   end if

   ! Dipole moment
   if ( idebug > 2 ) then
      write(6,'(a)') 'Waiting to receive dipole moment...'
   end if
   ! QM dipole moment
   call MPI_Recv( dipmom(:,1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 1 ) then
      write(6,'(a,4es15.6)') 'Received QM  dipole moment from server:', dipmom(:,1)
      call flush(6)
   end if
   ! MM dipole moment, disabled for now
!   call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
!   if ( idebug > 1 ) then
!      write(6,'(a,4es15.6)') 'Received MM  dipole moment from server:', dipmom(:,2)
!      call flush(6)
!   end if
   ! TOT dipole moment
!   call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
!   if ( idebug > 1 ) then
!      write(6,'(a,4es15.6)') 'Received TOT dipole moment from server:', dipmom(:,3)
!      call flush(6)
!   end if
   
   ! QM gradients
   if ( idebug > 1 ) then
         write(*,'(a)') 'Waiting to receive gradients...'
   end if
   call MPI_Recv( dxyz_all, 3*(natqm+natmm_tera), MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 1 ) then
      write(6,'(a)') 'Received the following gradients from server:'
      do iat=1, natqm+natmm_tera
         write(6,*) 'Atom ',iat, ': ',dxyz_all(:,iat)
      end do
      call flush(6)
   end if
!!$OMP CRITICAL
   do iat=1,natqm+natmm_tera
      fx(iat,iw)=-dxyz_all(1,iat)
      fy(iat,iw)=-dxyz_all(2,iat)
      fz(iat,iw)=-dxyz_all(3,iat)
   end do
!!$OMP END CRITICAL

!!$OMP ATOMIC
   eclas = eclas + escf / nwalk

!!$OMP SECTION

!!$OMP SECTION
   ! ONIOM was not yet tested!!
   if (iqmmm.eq.1) call oniom(x, y, z, fx, fy, fz, eclas, iw)
!!$OMP SECTION
!!$OMP END PARALLEL SECTIONS
   ! nwalk end do
   end do
!!$OMP END PARALLEL DO


end subroutine force_tera



subroutine connect_terachem( )
   use mod_utils, only: abinerror
   include 'mpif.h'
   character(255)  :: port_name
   integer         :: ierr
!   real*8          :: timer
!   logical         :: done=.false.
!   character(len=50) :: server_name

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------
!     server_name = trim(teraport)  !//'.'//trim(id)
!     write(*,*)''
!     write(6,'(2a)') 'Looking up TeraChem server under name:', trim(server_name)
!     call flush(6)

!    timer = MPI_WTIME(ierr)
!    do while (done .eqv. .false.)

!    done = .true.

!      call MPI_LOOKUP_NAME(server_name, MPI_INFO_NULL, port_name, ierr)
!      if (ierr == MPI_SUCCESS) then
!        if ( idebug > 1 ) then
!          write(6,'(2a)') 'Found port: ', trim(port_name)
!          call flush(6)
!        end if
!        done=.true.

!      else
!         write(*,*)'Error in MPI_LOOKUP_NAME. Error code:', ierr
!      end if

!      if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
!              write(*,*)'Port"'//trim(server_name)//'" not found. Timed out after 60 seconds.'
!              call abinerror("connect_to_terachem")
!      end if

!    end do

    ! DH hack, since MPI_LOOKUP_NAME does not work:
    write(6,'(A)') 'Reading TeraChem port name from file port.txt...'
    open(500, file="port.txt", action="read")
    read(500, *)port_name
    close(500)
    write(6,'(2a)') 'Looking up TeraChem port under name:', trim(port_name)
    ! ----------------------------------------
    ! Establish new communicator via port name
    ! ----------------------------------------
    write(*,*)'Establishing connection...'
    call flush(6)
    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
    write(6,'(a,i0)') 'Established new communicator:', newcomm

  end subroutine connect_terachem

   subroutine initialize_terachem()
   use mod_qmmm,  only: natqm
   use mod_system,only: names
   use mod_utils, only: abinerror
   include 'mpif.h'
   integer :: ierr, iat

   if (natmm_tera.gt.0)then
      allocate(mmcharges(natmm_tera), mmcoords(3, natmm_tera))
   end if
   allocate( qmcoords(3, natqm), qmcharges(natqm))
   allocate( dxyz_all(3, natqm+natmm_tera) )

   write(*,*)'Sending initial number of QM atoms to TeraChem.'
   call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

   ! DH WARNING: FOR QMMM WE WILL NEED TO send only part of names()
   ! Is it safe to send the whole array? (it should be, it is passed by reference)
   allocate(names_qm(natqm))
   do iat=1,natqm
      names_qm(iat) = names(iat)
   end do
   write(*,*)'Sending initial QM atom names to TeraChem.'
   call MPI_Send( names_qm, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )

   end subroutine initialize_terachem


  subroutine finalize_terachem(error_code)
  include 'mpif.h'
  integer, intent(in) :: error_code
  integer :: request
  integer :: ierr
  integer :: empty=0

  if(allocated(names_qm)) deallocate( names_qm )
  write(*,*)'Shutting down TeraChem.'
  if (error_code.eq.0)then
     call MPI_Send( empty, 1, MPI_INTEGER, 0, 0, newcomm, ierr )
  else
     ! Not sure whether this can lead to deadlock when terachem sends and not
     ! receive. Maybe we should avoid this.
     call MPI_Send( empty, 1, MPI_INTEGER, 0, 13, newcomm, ierr )
     !call MPI_ISend( empty, 1, MPI_INTEGER, 0, 13, newcomm, request, ierr )
     ! for some reason, non-blocking ISend does not work :(
  end if
  if (ierr.ne.0)then
     write(*,*)'I got a MPI Error when I tried to shutdown TeraChem'
     write(*,*)'Please, verify manually that the TeraChem server was terminated.'
     write(*,*)'Error code was:', ierr
  end if
  ! Ugly, but portable for different compilers
!  call system("sleep 5")
  end subroutine finalize_terachem


end module mod_terampi

