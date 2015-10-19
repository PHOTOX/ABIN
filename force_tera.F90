module mod_terampi
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD.
! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
! So far, I was not able to make it work with OpenMPI.
! (but now that we use file based tera_port, it should work as well)
!
! Currently supports:
! pure QM
! Adapted from Sander (Amber14)
!
! Original Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
! Modified by Dan Hollas hollas@vscht.cz
! Date: September 2014
! ----------------------------------------------------------------
  implicit none
  private
  public :: teraport, tc_finalize, connect_to_terachem, force_tera
  ! This does not work at this moment.
  character*50  ::  teraport = 'terachem_port'
  integer     ::  newcomm ! Communicator, initialized in mpi_init subroutine
  save

contains

subroutine force_tera(x, y, z, fx, fy, fz, eclas)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, nwalk, DP
   use mod_system, only: names
   use mod_qmmm, only: natqm, natmm
   use mod_utils, only: printf
   include 'mpif.h'
   real(DP),intent(in)      ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   ::  eclas
!   integer, intent(in) :: qmtypes(natqm)    ! QM atom types (nuclear charge in au)
   real*8 :: escf                 ! SCF energy
   real*8 :: dxyzqm(3,natqm)   ! SCF QM force
   real*8              :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
   real*8              :: qmcharges(natqm)  ! QM charges from population analysis
   real*8              :: mmcharges(natmm)  ! QM charges from population analysis
   real*8               :: qmcoords(3,natqm) 
   real*8               :: mmcoords(3,natmm) 
!   real*8,   :: clcoords(4,nqmatoms)

   real*8              :: dxyz_all(3,natmm+natqm)
   integer             :: i, status(MPI_STATUS_SIZE)
   integer             :: ierr, iw, iat


   !DH WARNING, initial hack, we do not support QM/MM yet
   natmm=0

   ! DHnote: we cannot use Niklasson's propagator if nwalk > 1
   ! This is the responsibility of the user
   do iw=1,nwalk

      do iat=1,natqm
         qmcoords(1,iat)=x(iat,iw)/ANG
         qmcoords(2,iat)=y(iat,iw)/ANG
         qmcoords(3,iat)=z(iat,iw)/ANG
      end do
      do iat=1,natmm
         mmcoords(1,iat)=x(iat+natqm,iw)
         mmcoords(2,iat)=y(iat+natqm,iw)
         mmcoords(3,iat)=z(iat+natqm,iw)
      end do

   ! -----------------------------------------
   ! Begin sending data each step to terachem
   ! -----------------------------------------


   ! Send natqm and the type of each qmatom
   if ( idebug > 1 ) then
      write(6,'(/, a, i0)') 'Sending natqm = ', natqm
      call flush(6)
   end if
   call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

   if ( idebug > 1 ) then
      write(6,'(/,a)') 'Sending QM atom types: '
   end if
   call MPI_Send( names, 2*size(names), MPI_CHARACTER, 0, 2, newcomm, ierr )

   ! Send QM coordinate array
   if ( idebug > 1 ) then
      write(6,'(a)') 'Sending QM coords: '
   end if
   do i=1, natqm
      if ( idebug > 2 ) then
         write(6,*) 'Atom ',i,': ',qmcoords(:,i)
         call flush(6)
      end if
   end do 
   call MPI_Send( qmcoords, 3*natqm, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

   ! Send natmm and the charge of each atom
   if ( idebug > 2 ) then
      write(6,'(a, i0)') 'Sending natmm = ', natmm
      call flush(6)
   end if
   call MPI_Send( natmm, 1, MPI_INTEGER, 0, 2, newcomm, ierr ) 

   if ( idebug > 2 ) then
      write(6,'(a)') 'Sending charges: '
   end if
   call MPI_Send( mmcharges, natmm, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

   ! Send MM point charge coordinate array
   if ( idebug > 2 ) then
      write(6,'(a)') 'Sending charges coords: '
   end if

   call MPI_Send( mmcoords, 3*natmm, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

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
      do i=1, natqm
         write(6,*) 'Atom ',i, ': ', qmcharges(i)
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
   ! MM dipole moment
   call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 1 ) then
      write(6,'(a,4es15.6)') 'Received MM  dipole moment from server:', dipmom(:,2)
      call flush(6)
   end if
   ! TOT dipole moment
   call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 1 ) then
      write(6,'(a,4es15.6)') 'Received TOT dipole moment from server:', dipmom(:,3)
      call flush(6)
   end if
   
   ! QM gradients
   if ( idebug > 1 ) then
         write(*,'(a)') 'Waiting to receive gradients...'
   end if
   call MPI_Recv( dxyz_all, 3*(natqm+natmm), MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   if ( idebug > 1 ) then
      write(6,'(a)') 'Received the following gradients from server:'
      do i=1, natqm+natmm
         write(6,*) 'Atom ',i, ': ',dxyz_all(:,i)
      end do
      call flush(6)
   end if

   do iat=1,natqm+natmm
      fx(iat,iw)=-dxyz_all(1,iat)
      fy(iat,iw)=-dxyz_all(2,iat)
      fz(iat,iw)=-dxyz_all(3,iat)
   end do

   eclas = eclas + escf / nwalk
   ! nwalk end do
   end do

end subroutine force_tera

subroutine connect_to_terachem( )
   use mod_qmmm,  only: natqm
   use mod_system,only: names
   use mod_utils, only: abinerror
   include 'mpif.h'
!  character(len=3) , intent(in) :: iw

   character(255)  :: port_name
   real*8          :: timer
   integer         :: ierr
   logical         :: done=.false.
   character(len=50) :: server_name

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------
     server_name = trim(teraport)  !//'.'//trim(id)
     write(*,*)''
     write(6,'(2a)') 'Looking up TeraChem server under name:', trim(server_name)
     call flush(6)

    timer = MPI_WTIME(ierr)
    do while (done .eqv. .false.)

    ! DH hack, since MPI_LOOKUP_NAME does not work:
    open(500, file="port.txt", action="read")
    read(500, *)port_name
    close(500)
    done = .true.

!      call MPI_LOOKUP_NAME(server_name, MPI_INFO_NULL, port_name, ierr)
!      if (ierr == MPI_SUCCESS) then
!        if ( idebug > 1 ) then
!          write(6,'(2a)') 'Found port: ', trim(port_name)
!          call flush(6)
 !       end if
!        done=.true.

!      else
!         write(*,*)'Error in MPI_LOOKUP_NAME. Error code:', ierr
!      end if

      if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
              write(*,*)'Port"'//trim(server_name)//'" not found. Timed out after 60 seconds.'
              call abinerror("connect_to_terachem")
      end if

    end do

    ! ----------------------------------------
    ! Establish new communicator via port name
    ! ----------------------------------------
    write(*,*)'Establishing connection...'
    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
    write(6,'(a,i0)') 'Established new communicator:', newcomm


    write(*,'(/, a, i0)') 'Sending number of QM atoms...' 
    call flush(6)
    call MPI_Send( natqm, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

    write(*,'(/,a)') 'Sending QM atom types: '
    ! DH WARNING: FOR QMMM WE WILL NEED TO send only part of names()
    ! Is it safe to send the whole array?(it should be, it is passed by reference)
    call MPI_Send( names, 2*size(names), MPI_CHARACTER, 0, 2, newcomm, ierr )

  end subroutine connect_to_terachem


  subroutine tc_finalize()
   include 'mpif.h'
    integer :: ierr
    real*8  :: empty

    write(*,*)'Shutting down TeraChem.'
    call MPI_Send( empty, 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr )
  end subroutine tc_finalize


end module mod_terampi

