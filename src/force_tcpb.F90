module mod_force_tcpb
! ----------------------------------------------------------------
! Interface for TeraChem based ab initio MD based on Protocol Buffers.
!
! Amber MPI interface that we've been using for ground state MD (force_tera.F90)
! was removed in 2021 from TeraChem.
! This module implements a replacement interface based on plain unix socket
! communication and data serialization through Google's Protocol buffers library.
!
! We don't need to deal with the details of the interface on ABIN side
! thanks to the tcpb-cpp library that provides a thin Fortran API. See
! https://github.com/mtzgroup/tcpb-cpp
!
! The simulations with this interface are triggered by pot='_tcpb_',
! see utils/r.tcpbabin for full setup.
!
! WARNING:
!  - QM/MM not implemented yet
!  - REMD is not supported as we cannot connect to multiple TC servers concurrently.
! ----------------------------------------------------------------
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
   use mod_utils, only: c_string
   implicit none
   private
   public :: initialize_tcpb, finalize_tcpb, force_tcpb
   ! TODO: Provide type bound function for validation
   ! https://fortran-lang.org/en/learn/quickstart/derived_types/
   type :: tcpb_params
      integer :: port = -1
      character(len=1024) :: hostname = ''
      character(len=1024) :: input_file = ''
      ! This make TC reuse WF from previous step.
      integer :: globaltreatment = 0
   end type
   type(tcpb_params) :: tcpb
   character(len=5), allocatable :: qmattypes(:)
   save

contains

   subroutine initialize_tcpb(natqm, at_names, tc_port, tc_host, tc_file)
      use mod_error, only: not_compiled_with
      use mod_utils, only: file_exists_or_exit
      integer, intent(in) :: natqm
      character(len=2), intent(in) :: at_names(:)
      integer, intent(in) :: tc_port
      character(len=*), intent(in) :: tc_host, tc_file
      integer, parameter :: natmm = 0
      integer, parameter :: MIN_PORT_NUMBER = 1025, MAX_PORT_NUMBER = 65536
      integer :: status
      integer :: i

#ifndef USE_TCPB
      call not_compiled_with("TCPB interface")
#endif
      if (tc_port < MIN_PORT_NUMBER .and. tc_port > 0) then
         call fatal_error(__FILE__, __LINE__, &
            & "Invalid TCPB port. Ports between 0 and 1024 are reserved by the system")
      else if (tc_port < 0) then
         call fatal_error(__FILE__, __LINE__, &
            & "Invalid TCPB port. Port cannot be negative and must be between 1025 and 65536")
      else if (tc_port > MAX_PORT_NUMBER) then
         call fatal_error(__FILE__, __LINE__, &
            & "TCPB port out of range. Port must be between 1025 and 65536")
      end if

      if (trim(tc_file) == "") then
         call fatal_error(__FILE__, __LINE__, &
            & "TCPB input file not provided")
      end if
      call file_exists_or_exit(tc_file)

      if (trim(tc_host) == "") then
         call fatal_error(__FILE__, __LINE__, &
            & "TCPB hostname not provided")
      end if

      tcpb = tcpb_params(port=tc_port, input_file=c_string(tc_file), hostname=c_string(tc_host))

      allocate (qmattypes(natqm))
      do i = 1, natqm
         qmattypes(i) = c_string(at_names(i))
      end do

      write (*, '(A,I0,A)') "Connecting to TeraChem TCPB server at " &
          & //trim(tc_host)//":", tcpb%port, " with input file "//trim(tc_file)
      status = -1

#ifdef USE_TCPB
      call tc_connect(tcpb%hostname, tcpb%port, status)
#endif

      if (status == 0) then
         write (stdout, *) "Successfully connected to TeraChem server."
      else if (status == 1) then
        call fatal_error(__FILE__, __LINE__, &
           & "Connection to TeraChem TCPB server failed! Is it running?")
      else if (status == 2) then
        call fatal_error(__FILE__, __LINE__, &
           & "Connection to TeraChem server succeed, but the "&
           &  //" server is not available!")
      else
        call fatal_error(__FILE__, __LINE__, &
           & "Could not connect to TCPB server. Is it running?")
      end if

      ! Setup TeraChem
      status = -1
#ifdef USE_TCPB
      call tc_setup(tcpb%input_file, qmattypes, natqm, status)
#endif
      status = 0
      if (status == 0) then
        write (*,*) "TeraChem setup completed with success."
      else if (status == 1) then
        call fatal_error(__FILE__, __LINE__, &
           & "TCPB: No options read from TeraChem input file or mismatch in the input options!")
      else if (status == 2) then
        call fatal_error(__FILE__, __LINE__, &
           & "TCPB: Failed to setup TeraChem.")
      else
        call fatal_error(__FILE__, __LINE__, &
           & "TCPB: Status on tc_setup function is not recognized!")
      end if
   end subroutine initialize_tcpb

   subroutine finalize_tcpb()
      ! TODO: I am not sure whether this actually does anything.
      ! TCPB server output is silent when this is called.
#ifdef USE_TCPB
      call tc_finalize()
#endif
   end subroutine finalize_tcpb

   subroutine force_tcpb(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: iqmmm
      use mod_qmmm, only: natqm
      use mod_shell_interface, only: oniom
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      real(DP), dimension(:), allocatable :: qmcoords, qmcharges, qmgrad
      ! For now we do not support QM/MM
      real(DP), dimension(1) :: mmcoords, mmcharges, mmgrad
      integer, parameter :: natmm = 0
      integer :: iw, iat, status
      logical :: abort

      allocate (qmcharges(natqm))
      allocate (qmcoords(3 * natqm))
      allocate (qmgrad(3* natqm))

      qmgrad = 0.0D0
      mmgrad = 0.0D0
      mmcoords = 0.0D0
      mmcharges= 0.0D0
      qmcharges= 0.0D0

      do iw = 1, walkmax

         do iat = 1, natqm
            qmcoords(3*iat - 2) = x(iat, iw)
            qmcoords(3*iat - 1) = y(iat, iw)
            qmcoords(3*iat) = z(iat, iw)
         end do

         status = -1
#ifdef USE_TCPB
         call tc_compute_energy_gradient(qmattypes, qmcoords, natqm, eclas, qmgrad, &
             & mmcoords, mmcharges, natmm, mmgrad, tcpb%globaltreatment, status)
#endif
         if (status == 0) then
            write (stdout, *) "TCPB Computed energy and gradient with success."
         else if (status == 1) then
            call fatal_error(__FILE__, __LINE__, &
               & "TCPB: Mismatch in the variables passed to tc_compute_energy_gradient()!")
         else if (status == 2) then
            call fatal_error(__FILE__, __LINE__, &
               & "TCPB: Problem to compute energy and gradient!")
         else
            call fatal_error(__FILE__, __LINE__, &
               & "TCPB: Status on tc_compute_energy_gradient function is not recognized!")
         end if

         do iat = 1, natqm
            fx(iat, iw) = -qmgrad(3*iat - 2)
            fy(iat, iw) = -qmgrad(3*iat - 1)
            fz(iat, iw) = -qmgrad(3*iat)
         end do

         ! ONIOM was not yet tested!!
         if (iqmmm == 1) then
            call oniom(x, y, z, fx, fy, fz, eclas, iw, abort)
            if (abort) call fatal_error(__FILE__, __LINE__, "ONIOM calculation failed")
         end if

      end do
   end subroutine force_tcpb

end module mod_force_tcpb
