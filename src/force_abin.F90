! This file defines the shell interface that we use to get forces
! and energies from external ab initio programs. Individual interfaces
! are implemented as BASH scripts in ../interfaces/
!
! The basic workflow is very simple:
! 1. ABIN writes current geometry into file geom.dat
!   in a XYZ format without the header.
! 2. ABIN launches the shell script POT/r.pot
! 3. The shellscript does a few things:
!    i) Takes the input geometry and prepares the input file
!    ii) Launches the QM program
!    iii) extracts energies and forces from the program output
!         into file engrad.dat
! 4. ABIN reads energies and forces from file engrad.dat
!
! NOTE: We append bead index to every file name so that we can
! call the interface in parallel in PIMD simulations.
! NOTE: Interface for Surface Hopping is a bit more complicated,
! see interfaces/MOLPRO-SH/r.molpro-sh
module mod_shell_interface_private
   use mod_const, only: DP, ANG
   use mod_error, only: fatal_error
   use mod_files, only: stderr
   implicit none
   ! Everything is public in this module for unit tests,
   ! but in the program one should only use mod_shell_interface,
   ! which re-exports public interface, at the end of this file.
   public

contains

   subroutine write_geom(fname, at_names, x, y, z, natom, iw)
      character(len=*), intent(in) :: fname
      character(len=2), dimension(:), intent(in) :: at_names
      real(DP), dimension(:, :), intent(in) :: x, y, z
      integer, intent(in) :: natom, iw
      ! Format for geom.dat; needed, so that Molpro can read it
      character(len=*), parameter :: fgeom = '(A2,3E25.17E2)'
      integer :: u, iat

      ! Delete the last geometry
      open (newunit=u, file=fname)
      close (u, status='delete')

      open (newunit=u, file=fname, status='new', action='write', access='sequential')
      do iat = 1, natom
         write (u, fgeom) at_names(iat), x(iat, iw) / ANG, y(iat, iw) / ANG, z(iat, iw) / ANG
      end do
      close (u)
   end subroutine write_geom

   subroutine write_sh_data(nstate, tocalc)
      integer, intent(in) :: nstate
      integer, dimension(:, :), intent(in) :: tocalc
      integer :: ist
      integer :: u

      open (newunit=u, file='state.dat')
      write (u, '(I0)') nstate
      ! Diagonal of tocalc holds info about needed forces
      ! tocalc(x,x) = 1 -> compute forces for electronic state X
      ! off-diagonal elements correspond to non-adiabatic couplings
      ! Here we only request forces, so only printing diagonal elements.
      ! NACMs are computed separately.
      do ist = 1, nstate
         write (u, '(I1,A1)', advance='no') tocalc(ist, ist), ' '
      end do
      close (u)
   end subroutine write_sh_data

   subroutine write_lz_data(nstate, nsinglet, ntriplet, tocalc)
      integer, intent(in) :: nstate
      integer, intent(in) :: nsinglet, ntriplet
      integer, dimension(:), intent(in) :: tocalc
      integer :: ist
      integer :: u

      if (nstate /= nsinglet + ntriplet) then
         call fatal_error(__FILE__, __LINE__, &
            & 'LZ: nstate /= nsinglet + ntriplet')
      end if

      open (newunit=u, file='state.dat')
      write (u, '(I0)') nstate

      ! Print for which state we need gradient
      ! First we have singlets, then triplets
      do ist = 1, nstate
         if (tocalc(ist) == 1) write (u, '(I0)') ist
      end do
      write (u, '(I0)') nsinglet
      write (u, '(I0)') ntriplet
      close (u)
   end subroutine write_lz_data

   function get_shellscript(potential) result(shellscript)
      use mod_utils, only: toupper
      character(len=*), intent(in) :: potential
      character(:), allocatable :: shellscript
      logical :: exists

      shellscript = './'//trim(toupper(potential))//'/r.'//potential
      inquire (file=shellscript, exist=exists)
      if (.not. exists) then
         call fatal_error(__FILE__, __LINE__, 'Shell executable '//shellscript//' does not exist')
      end if
   end function get_shellscript

   subroutine call_shell(shellscript, it, iw, ipimd)
      use mod_utils, only: append_rank
      character(:), allocatable, intent(in) :: shellscript
      ! Time step
      integer, intent(in) :: it
      ! Bead index (perhaps appended by REMD index)
      integer, intent(in) :: iw
      integer, intent(in) :: ipimd
      integer :: istatus
      character(len=300) :: call_cmd

      call_cmd = ''
      ! Passing arguments to bash script
      ! First argument is time step
      ! Second argument is the bead index, neccessary for parallel calculations
      write (call_cmd, '(A,I0,I4.3)') './'//trim(shellscript)//' ', it, iw
      call_cmd = append_rank(call_cmd)

      ! For SH, pass the 4th parameter: precision of forces as 10^(-force_accu1)
      ! TODO: This threshold should not be hard-coded
      if (ipimd == 2 .or. ipimd == 5) then
         write (call_cmd, '(A,I0,A)') trim(call_cmd)//' ', 7, ' < state.dat'
      end if

      ! Call the shell interface script
      call execute_command_line(trim(call_cmd), exitstat=istatus)

      ! TODO: Test whether this is true for execute_command_line
      ! NOTE: For some reason, shell exit status 1 turns into 256
      ! However, this one we also get by default from BASH, I don't know why.
      ! If the BASH script wants to exit with an error, it should use e.g. 'exit 2'
      if (istatus /= 0 .and. istatus /= 256) then
         write (stderr, *) 'Error during the execution of the external ab initio program.'
         write (stderr, *) 'Inspect the output files generated by '//trim(shellscript)
         call fatal_error(__FILE__, __LINE__, 'External forces error')
      end if
   end subroutine call_shell

   integer function open_engrad_file(fname) result(uengrad)
      character(len=*) :: fname
      character(len=300) :: errmsg
      logical :: exists
      integer :: iost

      inquire (file=fname, exist=exists)
      if (.not. exists) then
         write (stderr, *) 'WARNING: File '//trim(fname)//' does not exist. Waiting..'
         ! Should flush HDD buffer
         call execute_command_line('sync')
         call execute_command_line('sleep 0.5')
      end if

      open (newunit=uengrad, file=fname, status='old', action='read', iostat=iost, iomsg=errmsg)
      if (iost /= 0) then
         write (stderr, *) trim(errmsg)
         call fatal_error(__FILE__, __LINE__, 'Could not open file '//fname)
      end if
   end function open_engrad_file

   real(DP) function read_energy(engrad_unit) result(energy)
      integer, intent(in) :: engrad_unit
      character(len=300) :: errmsg, fname
      integer :: iost
      logical :: lopened

      ! Read electronic energy from engrad.dat
      read (engrad_unit, *, iostat=iost, iomsg=errmsg) energy
      if (iost /= 0) then
         write (stderr, *) trim(errmsg)
         inquire (unit=engrad_unit, opened=lopened, name=fname)
         call fatal_error(__FILE__, __LINE__, 'Could not read energy from file '//trim(fname))
      end if
   end function

   subroutine read_forces(fx, fy, fz, num_atom, iw, engrad_unit)
      use mod_files, only: stderr
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      integer, intent(in) :: iw, engrad_unit, num_atom
      character(len=300) :: errmsg, fname
      logical :: lopened
      integer :: iat, iost

      ! WARNING: The engrad file contains energy gradients, we need to convert to forces.
      do iat = 1, num_atom
         read (engrad_unit, *, iostat=iost) fx(iat, iw), fy(iat, iw), fz(iat, iw)
         if (iost /= 0) then
            write (stderr, *) trim(errmsg)
            inquire (unit=engrad_unit, opened=lopened, name=fname)
            call fatal_error(__FILE__, __LINE__, 'Could not read gradients from file '//trim(fname))
            return
         end if
         ! Convert gradients to forces
         fx(iat, iw) = -fx(iat, iw)
         fy(iat, iw) = -fy(iat, iw)
         fz(iat, iw) = -fz(iat, iw)
      end do
   end subroutine read_forces

   subroutine force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
      use mod_mpi, only: get_mpi_rank
      use mod_general, only: ipimd, iqmmm, it
      use mod_system, only: names
      use mod_sh_integ, only: nstate
      use mod_sh, only: tocalc, en_array, istate
      use mod_lz, only: nstate_lz, tocalc_lz, en_array_lz, istate_lz, nsinglet_lz, ntriplet_lz
      use mod_qmmm, only: natqm
      use mod_utils, only: toupper, append_rank
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(out) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(out) :: eclas
      integer, intent(in) :: walkmax
      character(len=*), intent(in) :: chpot
      real(DP) :: energy
      character(:), allocatable :: shellscript
      character(len=50) :: chgeom, chforce
      integer :: iw
      integer :: ist1
      integer :: engrad_unit

      eclas = 0.0D0

      ! Here we decide which program we use to obtain gradients and energies
      ! e.g. ./G09/r.g09
      shellscript = get_shellscript(chpot)

!$OMP PARALLEL
!$OMP DO PRIVATE(engrad_unit, energy, chgeom, chforce)
      do iw = 1, walkmax

         ! Write XYZ geometry for the external program
         write (chgeom, '(A,I3.3)') 'geom.dat.', iw
         chgeom = append_rank(chgeom)
         call write_geom(chgeom, names, x, y, z, natqm, iw)

         ! Number of states and current state for Surface Hopping
         if (ipimd == 2) then
            call write_sh_data(nstate, tocalc)
         end if

         ! Landau-Zener
         if (ipimd == 5) then
            call write_lz_data(nstate_lz, nsinglet_lz, ntriplet_lz, tocalc_lz)
         end if

         ! Call the external program
         call call_shell(shellscript, it, iw, ipimd)

         write (chforce, '(A,I3.3)') 'engrad.dat.', iw
         chforce = append_rank(chforce)
         engrad_unit = open_engrad_file(chforce)

         ! Read electronic energies from engrad.dat
         if (ipimd == 2) then

            do ist1 = 1, nstate
               en_array(ist1) = read_energy(engrad_unit)
            end do
            eclas = en_array(istate)

         else if (ipimd == 5) then

            ! Move old energies by 1
            en_array_lz(:, 3) = en_array_lz(:, 2)
            en_array_lz(:, 2) = en_array_lz(:, 1)

            ! Store the new on
            do ist1 = 1, nstate_lz
               en_array_lz(ist1, 1) = read_energy(engrad_unit)
            end do
            eclas = en_array_lz(istate_lz, 1)

         else
            energy = read_energy(engrad_unit)
!$OMP ATOMIC
            eclas = eclas + energy
         end if

         ! Read gradients from engrad.dat
         if (ipimd == 5) then
            ! Read only the computed state
            call read_forces(fx, fy, fz, natqm, tocalc_lz(istate_lz), engrad_unit)
         else
            call read_forces(fx, fy, fz, natqm, iw, engrad_unit)
         end if

         close (unit=engrad_unit, status='delete')

         if (iqmmm == 1) then
            call oniom(x, y, z, fx, fy, fz, eclas, iw)
         end if

      end do
!$OMP END DO
!$OMP END PARALLEL

      eclas = eclas / walkmax

   end subroutine force_abin

   subroutine oniom(x, y, z, fx, fy, fz, eclas, iw)
      use mod_general, only: natom, it, ipimd
      use mod_utils, only: append_rank
      use mod_system, only: names
      use mod_qmmm, only: natqm
      use mod_sh_integ, only: nstate
      use mod_sh, only: en_array, istate
      implicit none
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: iw
      real(DP), dimension(natom, 1) :: fx_tmp, fy_tmp, fz_tmp
      character(:), allocatable :: shellscript
      character(len=100) :: chgeom, chforce
      real(DP) :: energy
      integer :: iat, ist
      integer :: engrad_unit

      write (chgeom, '(A,I3.3)') 'geom_mm.dat.', iw
      write (chforce, '(A,I3.3)') 'engrad_mm.dat.', iw
      chgeom = append_rank(chgeom)
      chforce = append_rank(chforce)

      shellscript = get_shellscript('mm')

      ! Write geometry of the whole system
      call write_geom(chgeom, names, x, y, z, natom, iw)

      call call_shell(shellscript, it, iw, ipimd)

      engrad_unit = open_engrad_file(chforce)

      if (ipimd == 2) then
         ! TODO: Surface Hopping with ONIOM not tested!
         do ist = 1, nstate
            en_array(ist) = en_array(ist) + read_energy(engrad_unit)
         end do
         eclas = en_array(istate)
      else if (ipimd == 5) then
         call fatal_error(__FILE__, __LINE__, 'Landau-Zener with ONIOM not implemented')
      else
         energy = read_energy(engrad_unit)
!$OMP ATOMIC
         eclas = eclas + energy
      end if

      ! Read gradients from engrad_mm.dat
      call read_forces(fx_tmp, fy_tmp, fz_tmp, natom, 1, engrad_unit)
      do iat = 1, natom
         fx(iat, iw) = fx(iat, iw) + fx_tmp(iat, 1)
         fy(iat, iw) = fy(iat, iw) + fy_tmp(iat, 1)
         fz(iat, iw) = fz(iat, iw) + fz_tmp(iat, 1)
      end do

      close (engrad_unit, status='delete')

      ! MM, only QM part
      call write_geom(chgeom, names, x, y, z, natqm, iw)

      call call_shell(shellscript, it, iw, ipimd)

      engrad_unit = open_engrad_file(chforce)

      ! Here we use the substractive QM/MM scheme,
      ! so we are substracting results from the QM-only part
      if (ipimd == 2) then
         do ist = 1, nstate
            en_array(ist) = en_array(ist) - read_energy(engrad_unit)
         end do
         eclas = en_array(istate)
      else if (ipimd == 5) then
         call fatal_error(__FILE__, __LINE__, 'Landau-Zener with ONIOM not implemented')
      else
         energy = read_energy(engrad_unit)
!$OMP ATOMIC
         eclas = eclas - energy
      end if

      ! Read gradients from engrad_mm.dat
      call read_forces(fx_tmp, fy_tmp, fz_tmp, natqm, 1, engrad_unit)
      do iat = 1, natqm
         ! NOTE: Substracting the MM forces for the QM part
         fx(iat, iw) = fx(iat, iw) - fx_tmp(iat, 1)
         fy(iat, iw) = fy(iat, iw) - fy_tmp(iat, 1)
         fz(iat, iw) = fz(iat, iw) - fz_tmp(iat, 1)
      end do

      close (engrad_unit, status='delete')

   end subroutine oniom

end module mod_shell_interface_private

! Re-export public interface
module mod_shell_interface
   use mod_shell_interface_private
   implicit none
   private
   public :: force_abin, oniom
end module mod_shell_interface
