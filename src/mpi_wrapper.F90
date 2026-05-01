! Wrapper functions for MPI calls.
!  - MPI initialization / finalization
!  - Common MPI function with proper error handling
!  - MPI error handling helpers
module mod_mpi
#ifdef USE_MPI
   use mpi
#endif
   implicit none
   private
   public :: initialize_mpi, finalize_mpi
   public :: mpi_barrier_wrapper
   public :: get_mpi_rank, get_mpi_size

#ifdef USE_MPI
   public :: check_recv_count, handle_mpi_error
   public :: get_mpi_error_string
#endif

contains

#ifdef USE_MPI
   subroutine initialize_mpi(pot, pot_ref, nteraservers)
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
      character(len=*), intent(in) :: pot, pot_ref
      integer, intent(in) :: nteraservers
      logical :: initialized
      integer :: ierr, ithread

      call MPI_Initialized(initialized, ierr)
      if (initialized) then
         return
      end if

      if ((pot == "_tera_" .or. pot_ref == "_tera_") .and. nteraservers > 1) then
         ! We will be calling TS servers concurently
         ! via OpenMP parallelization, hence we need MPI_Init_thread().
         ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node303.htm
         call MPI_Init_thread(MPI_THREAD_MULTIPLE, ithread, ierr)
         if (ierr /= 0) then
            write (ERROR_UNIT, *) 'Bad signal from MPI_Init_thread:', ierr
            stop 1
         end if
         if (ithread /= MPI_THREAD_MULTIPLE) then
            write (ERROR_UNIT, *) 'Provided safety level is not MPI_THREAD_MULTIPLE'
            write (ERROR_UNIT, '(A,I0,A,I0)') 'Requested ', MPI_THREAD_MULTIPLE, ' got: ', ithread
            stop 1
         end if
      else
         call MPI_Init(ierr)
         if (ierr /= 0) then
            write (ERROR_UNIT, *) 'Bad signal from MPI_INIT:', ierr
            stop 1
         end if
      end if
   end subroutine initialize_mpi

   subroutine finalize_mpi(error_code)
      integer, intent(in) :: error_code
      logical :: initialized, finalized
      integer :: ierr

      call MPI_Initialized(initialized, ierr)
      if (.not. initialized) then
         return
      end if

      call MPI_Finalized(finalized, ierr)

      if (.not. finalized) then
         ! NOTE: For some reason, MPI_Abort does not play nicely with code coverage.
         ! Therefore, for the case of just one MPI replica, e.g. for pot == _tera_,
         ! we only call MPI_Finalize() and terminate by `stop error_code`.
         ! For REMD, we need to ensure that all replicas are stopped,
         ! so MPI_Abort is safer.
         if (error_code /= 0 .and. get_mpi_size() > 1) then
            call MPI_Abort(MPI_COMM_WORLD, error_code, ierr)
         else
            call MPI_Finalize(ierr)
         end if
      end if
   end subroutine finalize_mpi

   function get_mpi_error_string(mpi_err) result(error_string)
      integer, intent(in) :: mpi_err
      character(len=MPI_MAX_ERROR_STRING) :: error_string
      integer :: result_len
      integer :: ierr

      call MPI_Error_string(mpi_err, error_string, result_len, ierr)
      if (ierr /= MPI_SUCCESS) then
         error_string = 'Unspecified MPI error'
      end if
   end function get_mpi_error_string

   subroutine handle_mpi_error(mpi_err)
      use mod_error, only: fatal_error
      integer, intent(in) :: mpi_err

      if (mpi_err /= MPI_SUCCESS) then
         call fatal_error(__FILE__, __LINE__, get_mpi_error_string(mpi_err))
      end if
   end subroutine handle_mpi_error

   subroutine check_recv_count(mpi_status, expected_count, datatype)
      use mod_error, only: fatal_error
      integer, intent(in) :: mpi_status(:)
      integer, intent(in) :: expected_count
      integer, intent(in) :: datatype ! e.g. MPI_INTEGER
      character(len=100) :: error_msg
      integer :: recv_count
      integer :: ierr

      call MPI_Get_count(mpi_status, datatype, recv_count, ierr)
      call handle_mpi_error(ierr)
      if (recv_count /= expected_count) then
         write (error_msg, '(A,I0,A,I0)') 'Unexpected MPI_Recv count.&
            & Received ', recv_count, ' bytes, expected ', expected_count
         call fatal_error(__FILE__, __LINE__, error_msg)
      end if
   end subroutine check_recv_count

   subroutine mpi_barrier_wrapper(communicator)
      integer, intent(in), optional :: communicator
      integer :: comm
      integer :: ierr

      comm = MPI_COMM_WORLD
      if (present(communicator)) then
         comm = communicator
      end if

      call MPI_Barrier(comm, ierr)
      call handle_mpi_error(ierr)
   end subroutine mpi_barrier_wrapper

   integer function get_mpi_size(communicator) result(mpi_size)
      integer, intent(in), optional :: communicator
      integer :: comm
      integer :: ierr

      comm = MPI_COMM_WORLD
      if (present(communicator)) then
         comm = communicator
      end if

      call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)
      call handle_mpi_error(ierr)
   end function get_mpi_size

   integer function get_mpi_rank(communicator) result(rank)
      integer, intent(in), optional :: communicator
      integer :: comm
      integer :: ierr

      comm = MPI_COMM_WORLD
      if (present(communicator)) then
         comm = communicator
      end if

      call MPI_Comm_rank(comm, rank, ierr)
      call handle_mpi_error(ierr)
   end function get_mpi_rank

#else

   subroutine initialize_mpi(pot, pot_ref, nteraservers)
      use mod_general, only: iremd
      use mod_error, only: not_compiled_with
      character(len=*), intent(in) :: pot, pot_ref
      integer, intent(in) :: nteraservers
      integer :: i
      i = nteraservers
      if (iremd == 1 .or. &
        & pot == '_tera_' .or. pot_ref == '_tera_' .or. &
        & pot == '_mace_' .or. &
        & pot == '_cp2k_' .or. pot_ref == '__cp2k__') then
         call not_compiled_with('MPI')
      end if
   end subroutine initialize_mpi

   subroutine finalize_mpi(error_code)
      integer, intent(in) :: error_code
      integer :: e
      e = error_code
   end subroutine finalize_mpi

   subroutine mpi_barrier_wrapper()
   end subroutine mpi_barrier_wrapper

   integer function get_mpi_rank() result(rank)
      rank = 0
   end function get_mpi_rank

   integer function get_mpi_size() result(mpi_size)
      mpi_size = 1
   end function get_mpi_size

#endif

end module mod_mpi
