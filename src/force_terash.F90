module mod_terampi_sh
   !----------------------------------------------------------------
   ! Interface for TeraChem based Surface Hopping.
   ! Based on the FMS interface from FMS90 (TerachemModule.f90)
   !
   ! Original Authors: Basile Curchod, J. Snyder and Ed Hohenstein
   !----------------------------------------------------------------
   use mod_const, only: DP
   use mod_files, only: stdout
   use mod_terampi, only: TC_TAG
   implicit none
   private
   public :: force_terash
   public :: init_terash, finalize_terash
   public :: write_wfn, read_wfn, move_new2old_terash
   ! These are public for unit tests
   public :: allocate_tc_arrays, set_blob_size, set_civec_size, set_nbf
   real(DP), allocatable :: CIvecs(:, :), MO(:, :), blob(:), NAC(:)
   real(DP), allocatable :: CIvecs_old(:, :), MO_old(:, :), blob_old(:)
   real(DP), allocatable :: SMatrix(:)
   integer :: civec, nbf, blobsize, oldWfn = 0
   save

contains

#ifdef USE_MPI
   subroutine force_terash(x, y, z, fx, fy, fz, eclas)
      use mod_terampi, only: get_tc_communicator
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer :: tc_comm

      ! For SH always use only one TC server.
      tc_comm = get_tc_communicator(1)

      call send_terash(x, y, z, tc_comm)

      call receive_terash(fx, fy, fz, eclas, tc_comm)
   end subroutine force_terash

   subroutine receive_terash(fx, fy, fz, eclas, tc_comm)
      use mod_terampi, only: wait_for_terachem
      use mod_const, only: DP, ANG
      use mod_error, only: fatal_error
      use mod_general, only: idebug, natom, en_restraint, ipimd
      use mod_mpi, only: handle_mpi_error, check_recv_count
      use mod_qmmm, only: natqm
      use mod_io, only: print_charges, print_dipoles, print_transdipoles
      use mod_sh_integ, only: nstate
      use mod_sh, only: check_CIVector, en_array, istate, nacx, nacy, nacz
      use mod_lz, only: en_array_lz
      use mpi
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      integer, intent(in) :: tc_comm
      real(DP) :: dip(nstate * 3), tdip((nstate - 1) * 3) ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
      real(DP) :: qmcharges(size(fx, 1))
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr, iat, iw, ist1, ist2, ipom, i

      iw = 1

      call wait_for_terachem(tc_comm)

      ! Receive energies from TC
      if (idebug > 0) write (*, '(a)') 'Receiving energies from TC.'
      call MPI_Recv(en_array, nstate, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, nstate, MPI_DOUBLE_PRECISION)

      eclas = en_array(istate)

      ! Landau-Zener arrays
      if (ipimd == 5) then
         !Move old energies by 1
         en_array_lz(:, 3) = en_array_lz(:, 2); 
         en_array_lz(:, 2) = en_array_lz(:, 1); 
         !Store the new one
         en_array_lz(:, 1) = en_array(:)
      end if

      if (idebug > 0) write (*, '(a)') 'Receiving transition dipoles from TC.'
      call MPI_Recv(TDip, (nstate - 1) * 3, &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, (nstate - 1) * 3, MPI_DOUBLE_PRECISION)

      ! TODO: these things should be printed in analysis.F90
      ! TODO: move charges and dipoles to array module and make them universal
      ! TODO: move TDIP to surface hopping module
      ! allow reading this stuff from other programs as well
      call print_transdipoles(TDip, istate, nstate - 1)

      if (idebug > 0) write (*, '(a)') 'Receiving dipole moments from TC.'
      call MPI_Recv(Dip, nstate * 3, &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, nstate * 3, MPI_DOUBLE_PRECISION)

      call print_dipoles(Dip, iw, nstate)

      if (idebug > 0) write (*, '(a)') 'Receiving atomic charges from TC.'
      call MPI_Recv(qmcharges, natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, natqm, MPI_DOUBLE_PRECISION)

      call print_charges(qmcharges, istate)

      if (idebug > 0) write (*, '(a)') 'Receiving MOs from TC.'
      call MPI_Recv(MO, nbf * nbf, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, nbf * nbf, MPI_DOUBLE_PRECISION)

      if (idebug > 0) write (*, '(a)') 'Receiving CI vectors from TC.'
      call MPI_Recv(CIvecs, nstate * civec, &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, nstate * civec, MPI_DOUBLE_PRECISION)

      if (idebug > 0) write (*, *) "Receiving wavefunction overlap."
      call MPI_Recv(SMatrix, nstate * nstate, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr); 
      call handle_mpi_error(ierr)
      call check_recv_count(status, nstate * nstate, MPI_DOUBLE_PRECISION)

      ! Should change the following according to what is done in TeraChem
      if (oldwfn /= 0) then
         i = Check_CIVector(CIvecs, CIvecs_old, civec, nstate)
      end if

      CIVecs_old = Civecs

      if (idebug > 0) write (*, '(a)') 'Receiving blob.'
      call MPI_Recv(blob, blobsize, &
                    MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)
      call check_recv_count(status, blobsize, MPI_DOUBLE_PRECISION)

      ! TODO: Extract all this to a function.
      if (idebug > 0) write (*, '(A)') 'Receiving gradients and NACME.'
      do ist1 = 1, nstate
         do ist2 = ist1, nstate

            if (idebug > 0) write (*, '(A,i0,1X,i0)') 'Receiving derivatives between states ', ist1, ist2

            ! NOTE: We do not filter here based on tocalc because TC always sends the whole
            ! derivative matrix, including zero elements, see 'terachem/fms.cpp:'
            ! Is TC sending zero arrays for NAC that we did not want it to compute???
            call MPI_Recv(NAC, 3 * natom, MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG, tc_comm, status, ierr)
            call handle_mpi_error(ierr)
            call check_recv_count(status, 3 * natom, MPI_DOUBLE_PRECISION)

            if (idebug > 0) write (*, *) (NAC(i), i=1, 3 * natom)

            ipom = 1
            ! Gradients
            if (ist1 == ist2 .and. istate == ist1) then
               do iat = 1, natom
                  fx(iat, iw) = -NAC(ipom)
                  fy(iat, iw) = -NAC(ipom + 1)
                  fz(iat, iw) = -NAC(ipom + 2)
                  ipom = ipom + 3
               end do
            else if (ist1 == ist2) then
               ! DH2Jirka: here we will read excited state forces..
               ! perhaps we can use the iw index for the excited state force e.g.
               ! (this assumes, that the initial state is ground state)
               if (en_restraint >= 1) then
                  if (ist1 > 2) then
                     call fatal_error(__FILE__, __LINE__, &
                        & ': Energy restraint not implemented for more than 2 states!')
                  end if
                  do iat = 1, natom
                     fx(iat, 2) = -NAC(ipom)
                     fy(iat, 2) = -NAC(ipom + 1)
                     fz(iat, 2) = -NAC(ipom + 2)
                     ipom = ipom + 3
                  end do
               else
                  cycle
               end if
            else
               ! NACME
               do iat = 1, natom
                  nacx(iat, ist1, ist2) = NAC(ipom)
                  nacy(iat, ist1, ist2) = NAC(ipom + 1)
                  nacz(iat, ist1, ist2) = NAC(ipom + 2)
                  nacx(iat, ist2, ist1) = -nacx(iat, ist1, ist2)
                  nacy(iat, ist2, ist1) = -nacy(iat, ist1, ist2)
                  nacz(iat, ist2, ist1) = -nacz(iat, ist1, ist2)
                  ipom = ipom + 3
               end do
            end if

         end do
      end do

      oldWfn = 1

   end subroutine receive_terash

   subroutine send_terash(x, y, z, tc_comm)
      use mod_terampi, only: send_coordinates
      use mod_const, only: DP, ANG, AUTOFS
      use mod_mpi, only: handle_mpi_error
      use mod_general, only: natom, idebug, sim_time, en_restraint
      use mod_qmmm, only: natqm
      use mod_sh_integ, only: nstate
      use mod_sh, only: istate, tocalc
      use mpi
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer, intent(in) :: tc_comm
      real(DP) :: bufdoubles(100)
      real(DP) :: vels(3, size(x, 1))
      integer :: ierr, iw, i, ist1, ist2
      integer :: bufints(nstate * (nstate - 1) / 2 + nstate + 12)
      integer, parameter :: FMSInit = 0

      iw = 1

      ! Send ESinit
      bufints(1) = FMSinit
      bufints(2) = natom
      bufints(3) = 1 ! doCoup
      bufints(4) = 0 ! TrajID=0 for SH
      bufints(5) = 0 ! T_FMS%CentID(1)
      bufints(6) = 0 ! T_FMS%CentID(2)
      bufints(7) = istate - 1 ! T_FMS%StateID ! currently not used in fms.cpp
      bufints(8) = oldWfn ! does ABIN have info about WF?
      bufints(9) = istate - 1 ! iCalcState-1 ! TC Target State
      bufints(10) = istate - 1 ! jCalcState-1
      bufints(11) = 0 ! first_call, not used
      bufints(12) = 0 ! FMSRestart, not used

      if (idebug > 0) write (*, *) 'Sending surface hopping configuration.'
      call MPI_Send(bufints, 12, MPI_INTEGER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      ! The following bit is not in FMS code
      ! let ABIN decide which derivatives should TC compute
      i = 1
      do ist1 = 1, nstate
         do ist2 = ist1, nstate
            if (ist1 == ist2 .and. ist1 == istate) then
               bufints(i) = 1
            else if (ist1 == ist2) then
               ! DH hack for jirka
               ! this will work only if we compute only S0 and S1 states
               if (en_restraint >= 1) then
                  bufints(i) = 1
               else
                  bufints(i) = 0
               end if
            else
               bufints(i) = tocalc(ist1, ist2)
            end if
            i = i + 1
         end do
      end do

      if (idebug > 0) then
         write (*, *) 'Sending derivative matrix logic.'
         write (*, *) (bufints(i), i=1, nstate * (nstate - 1) / 2 + nstate)
      end if
      call MPI_SSend(bufints, nstate * (nstate - 1) / 2 + nstate, MPI_INTEGER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      ! Send Time
      bufdoubles(1) = sim_time
      call MPI_Send(bufdoubles, 1, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      call send_coordinates(x, y, z, natqm, iw, tc_comm, 'bohr')

      ! Send previous diabatic MOs
      if (idebug > 0) write (*, *) 'Sending previous orbitals.', nbf * nbf
      call MPI_Send(MO, nbf * nbf, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      ! Send previous CI vecs
      if (idebug > 0) write (*, *) 'Sending CI vector of size ', civec * nstate
      call MPI_Send(CIvecs, civec * nstate, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      if (idebug > 0) write (*, *) 'Sending blob.'
      call MPI_Send(blob, blobsize, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      if (idebug > 0) write (*, *) 'Sending velocities'
      ! Only needed for numerical NACME, so send 0 instead for now
      vels = 0.0D0
      call MPI_Send(vels, 3 * natom, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
      ! Imaginary velocities for FMS, not needed here, sending zeros...
      if (idebug > 0) write (*, *) 'Sending imaginary velocities'
      call MPI_Send(vels, 3 * natom, MPI_DOUBLE_PRECISION, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)

      if (idebug > 0) write (*, *) 'Succesfully sent all data to TeraChem-FMS'
   end subroutine send_terash

   subroutine init_terash(x, y, z)
      use mpi
      use mod_const, only: DP, ANG
      use mod_general, only: idebug, DP, natom
      use mod_system, only: names
      use mod_qmmm, only: natqm
      use mod_sh_integ, only: nstate
      use mod_mpi, only: handle_mpi_error
      use mod_terampi, only: get_tc_communicator, &
                             send_natom, &
                             send_atom_types_and_scrdir, &
                             send_coordinates
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr, iw, tc_comm
      integer, parameter :: FMS_INIT = 1
      logical, parameter :: send_scrdir = .false.
      integer :: bufints(3)
      ! QMMM currently not supported
      integer, parameter :: natmm_tera = 0

      ! use only one TC server !
      tc_comm = get_tc_communicator(1)

      iw = 1

      bufints(1) = FMS_INIT
      bufints(2) = natqm
      bufints(3) = natmm_tera
      call MPI_SSend(bufints, 3, MPI_INTEGER, 0, TC_TAG, tc_comm, ierr)
      call handle_mpi_error(ierr)
      if (idebug > 0) write (stdout, '(a)') 'Sent initial FMSinit.'

      call send_atom_types_and_scrdir(names, natqm, iw, tc_comm, send_scrdir)

      call send_coordinates(x, y, z, natqm, iw, tc_comm, 'bohr')

      ! Receive initial info TeraChem.
      ! Receive nbf, CI length and blob size
      ! No energies or gradients yet.
      call MPI_Recv(bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                    MPI_ANY_TAG, tc_comm, status, ierr)
      call handle_mpi_error(ierr)

      call set_civec_size(bufints(1))
      call set_nbf(bufints(2))
      call set_blob_size(bufints(3))

      write (stdout, *) 'size of CI vector, number of AOs, blob size:', bufints(1), bufints(2), bufints(3)
      call allocate_tc_arrays(nstate, natom)
   end subroutine init_terash
! USE_MPI
#endif

   subroutine set_civec_size(n)
      use mod_error, only: fatal_error
      integer, intent(in) :: n
      if (n <= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'invalid array size in set_civec_size')
         return
      end if
      civec = n
   end subroutine set_civec_size

   subroutine set_blob_size(n)
      use mod_error, only: fatal_error
      integer, intent(in) :: n
      if (n <= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'invalid array size in set_blob_size')
         return
      end if
      blobsize = n
   end subroutine set_blob_size

   subroutine set_nbf(n)
      use mod_error, only: fatal_error
      integer, intent(in) :: n
      if (n <= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'invalid array size in set_nbf')
         return
      end if
      nbf = n
   end subroutine set_nbf

   subroutine allocate_tc_arrays(nstate, natom)
      integer, intent(in) :: nstate
      integer, intent(in) :: natom
      allocate (MO(nbf, nbf))
      allocate (MO_old(nbf, nbf))
      allocate (CiVecs(civec, nstate))
      allocate (CiVecs_old(civec, nstate))
      allocate (NAC(natom * 3))
      allocate (blob(blobsize))
      allocate (blob_old(blobsize))
      allocate (SMatrix(nstate * nstate))
      MO = 0.0D0
      MO_old = 0.0D0
      CIVecs = 0.0D0
      CIVecs_old = 0.0D0
      blob = 0.0D0
      blob_old = 0.0D0
      SMatrix = 0.0D0
      NAC = 0.0D0
   end subroutine allocate_tc_arrays

   subroutine finalize_terash()
      if (allocated(MO)) then
         deallocate (MO, MO_old)
         deallocate (blob, blob_old)
         deallocate (CiVecs, CiVecs_old)
         deallocate (NAC, SMatrix)
      end if
   end subroutine finalize_terash

   subroutine write_wfn()
      use mod_general, only: it, sim_time, narchive
      use mod_sh_integ, only: nstate
      use mod_utils, only: rename_file, archive_file
      character(len=*), parameter :: fname = 'wfn.bin'
      integer :: uwfn

      call rename_file(fname, fname//'.old')

      open (newunit=uwfn, file=fname, action='write', status="new", access="sequential", form="unformatted")

      write (uwfn) it, sim_time
      write (uwfn) nbf
      write (uwfn) MO
      write (uwfn) civec, nstate
      write (uwfn) Civecs
      write (uwfn) blobsize
      write (uwfn) blob

      close (uwfn)

      if (modulo(it, narchive) == 0) then
         call archive_file('wfn.bin', it)
      end if
   end subroutine write_wfn

   subroutine read_wfn()
      use mod_general, only: iknow, it
      use mod_files, only: stderr
      use mod_chars, only: chknow
      use mod_error, only: fatal_error
      use mod_utils, only: archive_file
      use mod_sh_integ, only: nstate
      character(len=*), parameter :: fname = 'wfn.bin'
      logical :: file_exists
      integer :: temp, temp2, time_step
      integer :: uwfn
      real(DP) :: stime

      inquire (file=fname, exist=file_exists)
      if (.not. file_exists) then
         close (uwfn)
         write (stderr, *) 'Wavefunction restart file '//trim(fname)//' does not exist!'
         if (iknow /= 1) then
            write (stderr, *) chknow
            call fatal_error(__FILE__, __LINE__, &
               & 'missing restart file '//trim(fname))
         end if
         return
      end if

      open (newunit=uwfn, file=fname, action='read', status="old", access="sequential", form="unformatted")

      read (uwfn) time_step, stime
      read (uwfn) temp
      if (temp /= nbf) then
         close (uwfn)
         call fatal_error(__FILE__, __LINE__, &
            & 'Number of MOs in restart file '//trim(fname)//' is inconsistent')
         return
      end if

      read (uwfn) MO

      read (uwfn) temp, temp2
      if (temp /= civec) then
         close (uwfn)
         call fatal_error(__FILE__, __LINE__, &
            & 'Size of CI vectors in restart file '//trim(fname)//' is inconsistent')
         return
      end if
      if (temp2 /= nstate) then
         close (uwfn)
         call fatal_error(__FILE__, __LINE__, &
            & 'Number of states in restart file '//trim(fname)//' is inconsistent')
         return
      end if

      read (uwfn) CIVecs

      read (uwfn) temp
      if (temp /= blobsize) then
         close (uwfn)
         call fatal_error(__FILE__, __LINE__, &
            & 'Size of blob in restart file '//trim(fname)//' is inconsistent')
         return
      end if

      read (uwfn) blob

      close (uwfn)

      call move_new2old_terash()
      oldWFN = 1
      call archive_file('wfn.bin', it)
   end subroutine read_wfn

   subroutine move_new2old_terash
      MO_old = MO
      CIVecs_old = CIVecs
      blob_old = blob
   end subroutine move_new2old_terash

#ifndef USE_MPI
   subroutine init_terash(x, y, z)
      use mod_error, only: not_compiled_with
      real(DP), intent(inout) :: x(:, :), y(:, :), z(:, :)
      ! Assignments just to squash compiler warnings
      x = 0.0D0; y = 0.0D0; z = 0.0D0
      call not_compiled_with('MPI')
   end subroutine init_terash

   subroutine force_terash(x, y, z, fx, fy, fz, eclas)
      use mod_const, only: DP
      use mod_error, only: not_compiled_with
      real(DP), intent(in) :: x(:, :), y(:, :), z(:, :)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      real(DP), intent(inout) :: eclas
      ! Assignments just to squash compiler warnings
      fx = x; fy = y; fz = z; eclas = 0.0D0
      call not_compiled_with('MPI')
   end subroutine force_terash
#endif
end module mod_terampi_sh
