! Ehrenfest electronic-structure interfaces.
! Backends receive coordinates and return all adiabatic states in atomic units.
module mod_eh_backend
   use mod_const, only: DP, ANG
   use mod_error, only: fatal_error
   use mod_files, only: stderr
   use mod_shell_io, only: write_geom, get_shellscript, open_engrad_file, read_energy, read_forces
   implicit none
   private
   public :: eh_evaluate_backend
contains

   subroutine eh_evaluate_backend(x, y, z, energies, fx, fy, fz, nacx, nacy, nacz, pot, nac_accu1, nac_accu2)
      use mod_potentials_nai, only: evaluate_nai
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: energies(:), fx(:, :, :), fy(:, :, :), fz(:, :, :)
      real(DP), intent(out) :: nacx(:, :, :), nacy(:, :, :), nacz(:, :, :)
      character(len=*), intent(in) :: pot
      integer, intent(in) :: nac_accu1, nac_accu2

      select case (pot)
      case ('_nai_')
         call evaluate_nai(x, y, z, energies, fx, fy, fz, nacx, nacy, nacz)
      case default
         ! Add direct program interfaces here; shell backends use the contract below.
         call eh_shell_backend(x, y, z, energies, fx, fy, fz, nacx, nacy, nacz, pot, nac_accu1, nac_accu2)
      end select
   end subroutine eh_evaluate_backend

   subroutine eh_shell_backend(x, y, z, energies, fx, fy, fz, nacx, nacy, nacz, pot, nac_accu1, nac_accu2)
      use mod_system, only: names
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: energies(:), fx(:, :, :), fy(:, :, :), fz(:, :, :)
      real(DP), intent(out) :: nacx(:, :, :), nacy(:, :, :), nacz(:, :, :)
      character(len=*), intent(in) :: pot
      integer, intent(in) :: nac_accu1, nac_accu2
      character(:), allocatable :: script
      integer :: u, state, nstates, natoms, status, gradient_accuracy
      logical :: abort

      natoms = size(x, 1)
      nstates = size(energies)
      script = trim(get_shellscript(pot))
      call write_geom('geom.dat.001', names, x, y, z, natoms, 1)
      call write_eh_requests(nstates, nac_request=.false.)
      call discard_output('engrad.dat.001')
      gradient_accuracy = 7
      ! BAGEL-EH evaluates gradients and NACs together and caches the NAC output.
      if (pot == 'bagel-eh' .and. nstates > 1) gradient_accuracy = max(gradient_accuracy, nac_accu1)
      status = run_eh_script(script, gradient_accuracy)
      if (status /= 0) call fatal_error(__FILE__, __LINE__, 'Ehrenfest gradient interface failed: '//script)

      abort = .false.
      u = open_engrad_file('engrad.dat.001', abort)
      if (abort) call fatal_error(__FILE__, __LINE__, 'Missing Ehrenfest energies and gradients')
      do state = 1, nstates
         energies(state) = read_energy(u, abort)
         if (abort) exit
      end do
      if (.not. abort) then
         do state = 1, nstates
            call read_forces(fx(:, :, state), fy(:, :, state), fz(:, :, state), natoms, 1, u, abort)
            if (abort) exit
         end do
      end if
      close (u, status='delete')
      if (abort) call fatal_error(__FILE__, __LINE__, 'Invalid Ehrenfest energies or state gradients')

      nacx = 0.0_DP; nacy = 0.0_DP; nacz = 0.0_DP
      if (nstates == 1) return
      call write_eh_requests(nstates, nac_request=.true.)
      call discard_output('nacm.dat')
      status = run_eh_script(script//'.nacm', nac_accu1)
      if (status == 0) status = read_eh_nac(nacx, nacy, nacz)
      if (status /= 0 .and. nac_accu1 > nac_accu2) then
         write (stderr, '(A)') 'WARNING: Retrying Ehrenfest NAC calculation with decreased accuracy'
         call discard_output('nacm.dat')
         status = run_eh_script(script//'.nacm', nac_accu2)
         if (status == 0) status = read_eh_nac(nacx, nacy, nacz)
      end if
      if (status /= 0) call fatal_error(__FILE__, __LINE__, 'Invalid or missing Ehrenfest NAC vectors')
   end subroutine eh_shell_backend

   ! eh_state.dat is private to this interface. All states/pairs are requested,
   ! independently of populations, an active state or SH screening thresholds.
   subroutine write_eh_requests(nstates, nac_request)
      integer, intent(in) :: nstates
      logical, intent(in) :: nac_request
      integer :: u, count, i
      count = nstates
      if (nac_request) count = nstates * (nstates - 1) / 2
      open (newunit=u, file='eh_state.dat', status='replace', action='write')
      write (u, '(I0)') nstates
      write (u, '(*(I0,1X))') (1, i=1, count)
      close (u)
   end subroutine write_eh_requests

   integer function run_eh_script(script, accuracy) result(status)
      use mod_general, only: it
      character(len=*), intent(in) :: script
      integer, intent(in) :: accuracy
      character(len=64) :: arguments
      integer :: command_status
      logical :: exists
      inquire (file=script, exist=exists)
      if (.not. exists) call fatal_error(__FILE__, __LINE__, 'Ehrenfest interface does not exist: '//script)
      write (arguments, '(1X,I0,A,I0,A)') it, ' 001 ', accuracy, ' < eh_state.dat'
      call execute_command_line(script//trim(arguments), exitstat=status, cmdstat=command_status)
      if (command_status /= 0) status = command_status
   end function run_eh_script

   integer function read_eh_nac(nacx, nacy, nacz) result(status)
      real(DP), intent(out) :: nacx(:, :, :), nacy(:, :, :), nacz(:, :, :)
      integer :: u, i, j, atom
      nacx = 0.0_DP; nacy = 0.0_DP; nacz = 0.0_DP
      open (newunit=u, file='nacm.dat', status='old', action='read', iostat=status)
      if (status /= 0) return
      do i = 1, size(nacx, 2) - 1
         do j = i + 1, size(nacx, 3)
            do atom = 1, size(nacx, 1)
               read (u, *, iostat=status) nacx(atom, i, j), nacy(atom, i, j), nacz(atom, i, j)
               if (status /= 0) then
                  close (u, status='delete')
                  return
               end if
            end do
            nacx(:, j, i) = -nacx(:, i, j)
            nacy(:, j, i) = -nacy(:, i, j)
            nacz(:, j, i) = -nacz(:, i, j)
         end do
      end do
      close (u, status='delete')
   end function read_eh_nac

   ! Never accept an output file left behind by an earlier backend call.
   subroutine discard_output(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=filename, exist=exists)
      if (.not. exists) return
      open (newunit=u, file=filename, status='old')
      close (u, status='delete')
   end subroutine discard_output

end module mod_eh_backend
