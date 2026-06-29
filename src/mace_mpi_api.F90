! MPI interface with MACE Python server
! Modeled after mod_terampi for TeraChem
module mod_mace_mpi
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
#ifdef USE_MPI
   use mpi
   use mod_mpi, only: handle_mpi_error, get_mpi_error_string
#endif
   implicit none
   private

   ! MPI tags for MACE protocol
   integer, parameter :: MACE_TAG_EXIT = 666
   integer, parameter :: MACE_TAG_DATA = 2
   integer, parameter :: MACE_TAG_ERROR = 13

   ! Port file for MACE MPI connection
   character(len=*), parameter :: MACE_PORT_FILE_NAME = 'mace_port.txt'

#ifdef USE_MPI
   integer :: mace_comm = MPI_COMM_NULL
   logical :: mace_communication_established = .false.
#endif

   public :: MACE_TAG_EXIT, MACE_TAG_ERROR, MACE_TAG_DATA
#ifdef USE_MPI
   public :: get_mace_communicator
   public :: send_mace_atom_types, send_mace_coordinates
#endif
   public :: initialize_mace_interface, initialize_mace_server, finalize_mace
   save

contains

#ifdef USE_MPI

   subroutine initialize_mace_interface()
      integer :: ierr

      mace_comm = MPI_COMM_NULL
      mace_communication_established = .false.

      ! Set MPI error handler to allow retries
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      call handle_mpi_error(ierr)

      call connect_mace_server()
   end subroutine initialize_mace_interface

   subroutine connect_mace_server()
      character(len=MPI_MAX_PORT_NAME) :: port_name
      integer :: ierr, newcomm

      call read_mpi_port_from_file(MACE_PORT_FILE_NAME, port_name)

      write (stdout, '(2A)') 'Found MACE port: ', trim(port_name)
      write (stdout, '(A)') 'Establishing connection to MACE server...'
      call flush (stdout)

      call MPI_Comm_connect(trim(port_name), MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
      call handle_mpi_error(ierr)
      write (stdout, '(A)') 'Connection to MACE server established!'

      mace_comm = newcomm
      mace_communication_established = .true.
   end subroutine connect_mace_server

   ! Read MPI port from a file.
   subroutine read_mpi_port_from_file(portfile, port_name)
      character(len=*), intent(in) :: portfile
      character(len=MPI_MAX_PORT_NAME), intent(out) :: port_name
      integer :: iunit, iost

      write (stdout, '(A)') 'Reading MPI port name from file '//portfile
      port_name = ''

      open (newunit=iunit, file=portfile, action="read", status="old", iostat=iost)
      read (iunit, '(A)', iostat=iost) port_name
      close (iunit, status='delete')
   end subroutine read_mpi_port_from_file

   integer function get_mace_communicator() result(comm)
      comm = mace_comm
   end function get_mace_communicator

   subroutine initialize_mace_server()
      use mod_qmmm, only: natqm
      use mod_system, only: names

      call send_natom(natqm, mace_comm)
      call send_mace_atom_types(names, natqm, mace_comm)
   end subroutine initialize_mace_server

   subroutine finalize_mace()
      integer :: ierr, empty

      ! Set error handler to return so we can handle errors gracefully
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)

      if (.not. mace_communication_established) return

      write (stdout, '(A)') 'Shutting down MACE server'

      call MPI_Send(empty, 0, MPI_INTEGER, 0, MACE_TAG_EXIT, mace_comm, ierr)

      if (ierr /= MPI_SUCCESS) then
         write (stderr, '(A)') 'MPI ERROR during shutdown of MACE server'
         write (stderr, '(A)') 'Verify manually that the MACE server was terminated.'
         write (stderr, *) get_mpi_error_string(ierr)
      end if

      call MPI_Comm_free(mace_comm, ierr)
      if (ierr /= MPI_SUCCESS) then
         write (stderr, *) get_mpi_error_string(ierr)
      end if
   end subroutine finalize_mace

   subroutine send_natom(num_atom, comm)
      use mod_general, only: idebug
      integer, intent(in) :: num_atom
      integer, intent(in) :: comm
      integer :: ierr

      if (idebug > 1) then
         write (stdout, '(A, I0)') 'MACE: Sending number of atoms = ', num_atom
         call flush (stdout)
      end if
      call MPI_Send(num_atom, 1, MPI_INTEGER, 0, MACE_TAG_DATA, comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_natom

   subroutine send_mace_atom_types(at_names, num_atom, comm)
      use mod_general, only: idebug
      character(len=2), intent(in) :: at_names(:)
      integer, intent(in) :: num_atom
      integer, intent(in) :: comm
      character(len=2*num_atom) :: buffer
      integer :: ierr, offset, iat

      buffer = ''
      offset = 1
      do iat = 1, num_atom
         write (buffer(offset:offset + 1), '(A2)') at_names(iat)
         offset = offset + 2
      end do

      if (idebug > 1) then
         write (stdout, '(A)') 'MACE: Sending atom types: '
         write (stdout, '(A)') trim(buffer)
         call flush (stdout)
      end if

      call MPI_Send(buffer, num_atom * 2, MPI_CHARACTER, 0, MACE_TAG_DATA, comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_mace_atom_types

   subroutine send_mace_coordinates(x, y, z, num_atom, iw, comm)
      use mod_general, only: idebug
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: num_atom
      integer, intent(in) :: iw
      integer, intent(in) :: comm
      real(DP), allocatable :: coords(:, :)
      integer :: ierr, iat

      allocate (coords(3, num_atom))

      ! Send coordinates in Bohr (server converts to Angstrom)
      do iat = 1, num_atom
         coords(1, iat) = x(iat, iw)
         coords(2, iat) = y(iat, iw)
         coords(3, iat) = z(iat, iw)
      end do

      if (idebug > 1) then
         write (stdout, '(A)') 'MACE: Sending coordinates [Bohr]: '
         do iat = 1, num_atom
            write (stdout, *) 'Atom ', iat, ': ', coords(:, iat)
            call flush (stdout)
         end do
      end if
      call MPI_Send(coords, num_atom * 3, MPI_DOUBLE_PRECISION, 0, MACE_TAG_DATA, comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine send_mace_coordinates

#else

   subroutine initialize_mace_interface()
      use mod_error, only: not_compiled_with
      call not_compiled_with('MPI')
   end subroutine initialize_mace_interface

   subroutine initialize_mace_server()
      use mod_error, only: not_compiled_with
      call not_compiled_with('MPI')
   end subroutine initialize_mace_server

   ! This must be a no-op, since it is called from finish()
   subroutine finalize_mace(error_code)
      integer, intent(in) :: error_code
      integer :: i
      i = error_code
   end subroutine finalize_mace

! USE_MPI
#endif

end module mod_mace_mpi
