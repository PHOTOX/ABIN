module mod_force_tcpb
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD.
! ----------------------------------------------------------------
   use mod_const, only: DP
   use mod_error, only: fatal_error
   use mod_files, only: stdout, stderr
   use mod_utils, only: c_string
   implicit none
   private
   public :: initialize_tcpb, finalize_tcpb, force_tcpb
   ! TODO: Pack these into a derived type
   integer :: port, status, globaltreatment
   character(len=80)  :: host
   character(len=256) :: tcfile
   character(len=5), allocatable :: qmattypes(:)
   save

contains

   subroutine initialize_tcpb(natqm, at_names)
      integer, intent(in) :: natqm
      character(len=2), intent(in) :: at_names(:)
      integer, parameter :: natmm = 0
      integer :: i

      ! TODO: These all need to be either cmdline arguments or input file params
      host = c_string("localhost")
      port = 12345
      tcfile = c_string("terachem.inp")
      ! TODO: Find out what this does
      globaltreatment = 0

      allocate (qmattypes(natqm))
      do i = 1, natqm
         qmattypes(i) = c_string(at_names(i))
      end do

      write (*, '(A,I0)') "Attempting to connect to TeraChem server using host " &
          & //trim(host)//" and port ", port
      status = -1

      call tc_connect(host, port, status)

      if (status == 0) then
         write (stdout, *) "Successfully connected to TeraChem server."
      else if (status == 1) then
        call fatal_error(__FILE__, __LINE__, &
           & "ERROR: Connection to TeraChem server failed!")
      else if (status == 2) then
        call fatal_error(__FILE__, __LINE__, &
           & "ERROR: Connection to TeraChem server succeed, but the "&
           &  //" server is not available!")
      else
        call fatal_error(__FILE__, __LINE__, &
           & "ERROR: Status on tc_connect function is not recognized!")
      end if

      ! Setup TeraChem
      status = -1
      call tc_setup(tcfile, qmattypes, natqm, status)
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
      call tc_finalize()
   end subroutine finalize_tcpb

   subroutine force_tcpb(x, y, z, fx, fy, fz, eclas, walkmax)
      use mod_general, only: iqmmm
      use mod_qmmm, only: natqm, natmm
      use mod_shell_interface, only: oniom
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: walkmax
      real(DP), dimension(:), allocatable :: qmcoords, qmcharges, qmgrad
      real(DP), dimension(1) :: mmcoords, mmcharges, mmgrad
      integer :: iw, iat, status
      logical :: abort

      ! TODO: Preallocate these in the init function
      allocate (qmcharges(natqm))
      allocate (qmcoords(3 * natqm))
      allocate (qmgrad(3* natqm))

      do iw = 1, walkmax

         do iat = 1, natqm
            qmcoords(3*iat - 2) = x(iat, iw)
            qmcoords(3*iat - 1) = y(iat, iw)
            qmcoords(3*iat) = z(iat, iw)
         end do

         status = -1
         call tc_compute_energy_gradient(qmattypes, qmcoords, natqm, eclas, qmgrad, &
             & mmcoords, mmcharges, natmm, mmgrad, globaltreatment, status)
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
