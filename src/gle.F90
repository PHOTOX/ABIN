! *********************************************************
! * Here we define the functions to be called for         *
! * colored-noise thermostatting. Incidentally, we also   *
! * write here a couple of functions for white-noise.     *
! *                                                       *
! * Code is licensed under GPLv3 [www.gnu.org]            *
! * please consider citing the relevant papers (listed    *
! * below) if you use GLE in your simulations.            *
! *                                                       *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************

! The original code was modified to allow for PI+GLE and PIGLET simulations
! in ABIN by Daniel Hollas, danekhollas at gmail dot com

! NOTE: In the ABIN code, use mod_gle, defined at the end of this
! file, which exports public interface to the GLE module.
! mod_gle_private should be used only in unit tests (test_gle.pf).
module mod_gle_private
   use mod_const, only: DP
   use mod_files, only: stdout
   implicit none
   public
   ! Relaxation time of the white noise Langevin thermostat (PILE)
   ! In ABIN input, it should be set in picoseconds
   real(DP) :: tau0_langevin = -1.0D0
   ! Special input flag for GLE tests
   logical :: gle_test
   ! Size of the GLE auxiliary arrays
   integer :: ns
   ! GLE auxiliary array
   real(DP), allocatable :: ps(:, :, :)
   ! GLE propagators
   real(DP), allocatable :: gS(:, :), gT(:, :)
   ! PIGLET propagators for path integral centroid
   real(DP), allocatable :: gS_centroid(:, :), gT_centroid(:, :)
   ! Temporary helper array
   real(DP), allocatable :: gp(:, :)
   ! Cumulative conserved quantity
   real(DP) :: langham = 0.0D0
   ! internal arrays for PILE thermostat
   real(DP), allocatable :: c1(:), c2(:)
   save

contains

   ! Getter for a conserved quantity
   ! Used both by white-noise and GLE thermostats.
   real(DP) function get_langham()
      get_langham = langham
   end function get_langham

   ! Initialize white-noise PILE thermostat,
   ! which can be used both for classical MD and PIMD.
   subroutine pile_init(dt, tau0)
      use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
      use mod_const, only: PI, AUtoFS
      use mod_general, only: nwalk, ipimd, inormalmodes
      use mod_nhc, only: temp
      use mod_error, only: fatal_error
      real(DP), intent(in) :: dt, tau0
      real(DP) :: gam, omega, tau
      integer :: iw

      ! Sanity checks
      if (ipimd == 1 .and. inormalmodes <= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'You must use normal mode coordinates with PILE thermostat.')
      end if

      if (tau0 <= 0.0D0) then
         write (ERROR_UNIT, *) 'ERROR: tau0_langevin for PILE thermostat was not set or was negative.'
         write (ERROR_UNIT, *) 'Set "tau0_langevin" in picoseconds in the ABIN input in section "nhcopt".'
         call fatal_error(__FILE__, __LINE__, 'invalid tau0_langevin')
      end if

      langham = 0.D0 ! sets to zero accumulator for langevin 'conserved' quantity

      ! c1 and c2 as defined in the paper, but
      ! c2 is further multiplied by sqrt(temp*nwalk)
      allocate (c1(nwalk), c2(nwalk))

      ! Centroid
      tau = tau0 / AUtoFS * 1000
      gam = 1 / tau
      c1(1) = exp(-dt * gam)
      c2(1) = dsqrt(1 - c1(1)**2) * dsqrt(temp * nwalk)

      ! Now normal modes of bead necklace
      do iw = 2, nwalk
         omega = 2 * temp * nwalk * sin((iw - 1) * PI / nwalk)
         gam = 2 * omega
         c1(iw) = exp(-dt * gam)
         c2(iw) = dsqrt(1 - c1(iw)**2) * dsqrt(temp * nwalk)
      end do
   end subroutine pile_init

   ! White-noise propagator. time-step has been set in pile_init
   subroutine pile_step(px, py, pz, m)
      use mod_general, only: natom, nwalk
      use mod_utils, only: ekin_p
      use mod_random, only: gautrg
      real(DP), intent(inout) :: px(:, :)
      real(DP), intent(inout) :: py(:, :)
      real(DP), intent(inout) :: pz(:, :)
      real(DP), intent(in) :: m(:, :)
      real(DP) :: ran(natom * 3 * nwalk)
      integer :: iat, iw, pom

      langham = langham + ekin_p(px, py, pz, m, natom, nwalk)

      pom = 1
      call gautrg(ran, natom * 3 * nwalk)

      ! This is a local version of the thermostat (PILE-L)
      ! TODO: implement global version according to equations 51-54
      do iw = 1, nwalk
         do iat = 1, natom
            px(iat, iw) = c1(iw) * px(iat, iw) + c2(iw) * ran(pom) * dsqrt(m(iat, iw))
            py(iat, iw) = c1(iw) * py(iat, iw) + c2(iw) * ran(pom + 1) * dsqrt(m(iat, iw))
            pz(iat, iw) = c1(iw) * pz(iat, iw) + c2(iw) * ran(pom + 2) * dsqrt(m(iat, iw))
            pom = pom + 3
         end do
      end do

      langham = langham - ekin_p(px, py, pz, m, natom, nwalk)
   end subroutine pile_step

   subroutine finalize_pile()
      if (allocated(c1)) deallocate (c1, c2)
   end subroutine finalize_pile

   subroutine print_gle_header(inose, ipimd, inormalmodes)
      integer, intent(in) :: inose
      integer, intent(in) :: ipimd
      integer, intent(in) :: inormalmodes

      write (stdout, *) ""
      write (stdout, *) "Initialization of GLE thermostat."
      write (stdout, *) "Please cite the following review paper about the GLE methodology"
      write (stdout, *) "M. Ceriotti, G. Bussi and M. Parrinello"
      write (stdout, *) "J. Chem. Theory Comput. 6, 1170 (2010)"
      write (stdout, *) "Colored-noise thermostats a la carte"
      write (stdout, *) "http://dx.doi.org/10.1021/ct900563s"
      write (stdout, *) ""
      write (stdout, *) "The citations for specific GLE use cases is at"
      write (stdout, *) "http://gle4md.org/index.html?page=refs"
      write (stdout, *) ""

      if (inose == 2 .and. ipimd == 0) then
         write (stdout, *) "Using Quantum Thermostat, please cite"
         write (stdout, *) "M. Ceriotti, G. Bussi and M. Parrinello"
         write (stdout, *) "Phy. Rev. Lett. 103, 030603 (2009)"
         write (stdout, *) "Nuclear quantum effects in solids using a colored-noise thermostat"
         write (stdout, *) "http://dx.doi.org/10.1103/PhysRevLett.103.030603"
         write (stdout, *) ""
      else if (inose == 2 .and. ipimd == 1 .and. inormalmodes == 0) then
         write (stdout, *) "Using PI+GLE thermostat, please cite:"
         write (stdout, *) "M. Ceriotti, D. E. Manolopoulos, and M. Parrinello"
         write (stdout, *) "J. Chem. Phys. 134, 084104 (2011)"
         write (stdout, *) "Accelerating the convergence of path integral dynamics"
         write (stdout, *) "with a generalized Langevin equation"
         write (stdout, *) "http://dx.doi.org/10.1063/1.3556661"
         write (stdout, *) ""
      else if (inose == 2 .and. ipimd == 1 .and. inormalmodes == 1) then
         write (stdout, *) "Using PIGLET thermostat, please cite"
         write (stdout, *) "M. Ceriotti, D. E. Manolopoulos, Rev. Lett. 109, 100604 (2012)"
         write (stdout, *) "Efficient first-principles calculation of the"
         write (stdout, *) "quantum kinetic energy and momentum distribution of nuclei"
         write (stdout, *) "http://dx.doi.org/10.1103/PhysRevLett.109.100604"
         write (stdout, *) ""
      end if
   end subroutine print_gle_header

   integer function read_ns(fname)
      use mod_error, only: fatal_error
      character(len=*), intent(in) :: fname
      integer :: u, i
      open (newunit=u, file=fname, action="read", status="old")
      read (u, *, iostat=i) read_ns
      if (i /= 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Could not read ns from the first line of file '//fname)
      end if
      close (u)
   end function read_ns

   subroutine read_gle_matrix(ns, fname, M)
      use mod_error, only: fatal_error
      integer, intent(in) :: ns
      character(len=*), intent(in) :: fname
      real(DP), intent(out) :: M(:, :)
      integer :: ns_read
      integer :: u, i

      open (newunit=u, file=fname, action='read', form='formatted', access='sequential')

      read (u, *, iostat=i) ns_read
      if (i /= 0) then
         close (u)
         call fatal_error(__FILE__, __LINE__, &
            & 'Could not read ns from the first line of file '//fname)
         return
      end if
      if (ns_read /= ns) then
         close (u)
         call fatal_error(__FILE__, __LINE__, &
            & 'ns read from the first line of file '//fname//' does not match')
         return
      end if

      do i = 1, ns + 1
         read (u, *) M(i, :)
      end do
      close (u)
   end subroutine read_gle_matrix

   subroutine gle_init(dt)
      use mod_const, only: AUtoEV
      use mod_utils, only: append_rank
      use mod_error, only: fatal_error
      use mod_mpi, only: get_mpi_rank
      use mod_general, only: natom, nwalk, ipimd, inormalmodes
      use mod_nhc, only: temp, inose
      real(DP), intent(in) :: dt
      real(DP), allocatable :: gA(:, :), gC(:, :)
      character(len=100) :: glea, glec, glea_centroid, glec_centroid
      integer :: i, iw

      call print_gle_header(inose, ipimd, inormalmodes)

      langham = 0.D0 ! sets to zero accumulator for langevin 'conserved' quantity

      glea = append_rank('GLE-A')
      glec = append_rank('GLE-C')
      glea_centroid = append_rank('GLE-A.centroid')
      glec_centroid = append_rank('GLE-C.centroid')

      ns = read_ns(glea)

      allocate (gA(ns + 1, ns + 1))
      allocate (gC(ns + 1, ns + 1))

      write (stdout, *) 'Reading A-matrix. Expecting a.u. units!'
      call read_gle_matrix(ns, glea, gA)

      ! read C (in eV!), or init to kT
      if (inose == 4) then
         write (stdout, *) "Using canonical-sampling, Cp=kT"
         gC = 0.0D0
         do i = 1, ns + 1
            gC(i, i) = temp
         end do
      else
         write (stdout, *) 'Reading specialized Cp matrix. Expecting eV units!'
         call read_gle_matrix(ns, glec, gC)
         gC = gC / AUtoEV
      end if

      ! Propagator matrices
      allocate (gS(ns + 1, ns + 1))
      allocate (gT(ns + 1, ns + 1))
      call compute_propagator(gA, gC, ns, dt, gT, gS)

      ! This is the main array holding auxiliary GLE momemta
      allocate (ps(natom * 3, ns, nwalk))
      ! Temporary array
      allocate (gp(natom * 3, ns + 1))
      ! Initialize the auxiliary vectors.
      ! We keep general - as we might be using non-diagonal C
      ! to break detailed balance - and we use cholesky decomposition of C
      ! since one would like to initialize correctly the velocities in
      ! gA used as a temporary array here
      call cholesky(gC, gA, ns + 1)
      do iw = 1, nwalk
         call initialize_momenta(gA, iw, natom, ps)
      end do

      ! For PIGLET, we repeat the procedure with separate
      ! GLE matrices for the Path Integral centroid.
      if (inormalmodes == 1) then
         write (stdout, *) 'Reading A matrix for centroid. Expecting a.u. units!'
         call read_gle_matrix(ns, glea_centroid, gA)

         write (stdout, *) 'Reading specialized Cp matrix for centroid. Expecting eV units!'
         call read_gle_matrix(ns, glec_centroid, gC)
         gC = gC / AUtoEV

         allocate (gS_centroid(ns + 1, ns + 1))
         allocate (gT_centroid(ns + 1, ns + 1))
         call compute_propagator(gA, gC, ns, dt, gT_centroid, gS_centroid)

         call cholesky(gC, gA, ns + 1)
         call initialize_momenta(gA, 1, natom, ps)
      end if

      deallocate (gA)
      deallocate (gC)

      ! To ensure numerical stability of End-to-end GLE tests,
      ! we read precomputed propagators.
      if (gle_test) then
         call read_propagator(dt, ns, 'GLE-S.bin', gS)
         call read_propagator(dt, ns, 'GLE-T.bin', gT)
         if (inormalmodes == 1) then
            call read_propagator(dt, ns, 'GLE-S.centroid.bin', gS_centroid)
            call read_propagator(dt, ns, 'GLE-T.centroid.bin', gT_centroid)
         end if
      end if
   end subroutine gle_init

   subroutine initialize_momenta(C, iw, natom, ps)
      !use mod_arrays,  only: px, py, pz, vx, vy, vz, amt
      use mod_error, only: fatal_error
      use mod_random, only: gautrg
      real(DP), intent(in) :: C(:, :)
      integer, intent(in) :: iw, natom
      real(DP), intent(inout) :: ps(:, :, :)
      real(DP), allocatable :: gr(:)
      integer :: i, j

      allocate (gr(ns + 1))

      do j = 1, natom * 3
         call gautrg(gr, ns + 1)
         gp(j, :) = matmul(C, gr)
      end do

      ! TODO, actually pass this to initialize momenta
      ! case of generic C, we use an extra slot for gp for the physical momentum
      !do i = 1, natom
      !   px(i,iw) = gp(i,1)
      !   py(i,iw) = gp(i + natom,1)
      !   pz(i,iw) = gp(i + 2*natom,1)
      !end do
      !vx = px / amt
      !vy = py / amt
      !vz = pz / amt

      do j = 1, natom * 3
         do i = 1, ns
            ps(j, i, iw) = gp(j, i + 1)
         end do
      end do

      deallocate (gr)
   end subroutine initialize_momenta

   subroutine pile_restin(u)
      ! Restart file unit
      integer, intent(in) :: u
      read (u, *) langham
   end subroutine pile_restin

   subroutine pile_restout(u)
      integer, intent(in) :: u
      write (u, '(ES24.16E3)') langham
   end subroutine pile_restout

   subroutine gle_restin(u)
      use mod_general, only: natom, nwalk
      ! Restart file unit
      integer, intent(in) :: u
      integer :: iat, iw, is

      do iw = 1, nwalk
         do iat = 1, natom * 3
            do is = 1, ns - 1
               read (u, '(ES25.16E3)', advance="no") ps(iat, is, iw)
            end do
            read (u, '(ES25.16E3)') ps(iat, ns, iw)
         end do
      end do
      read (u, *) langham
   end subroutine gle_restin

   subroutine gle_restout(u)
      use mod_general, only: natom, nwalk
      ! Restart file unit
      integer, intent(in) :: u
      character(len=50) :: chformat
      integer :: iat, iw, is

      write (chformat, '(A1,I0,A)') '(', ns, 'ES25.16E3)'
      do iw = 1, nwalk
         do iat = 1, natom * 3
            write (u, fmt=trim(chformat)) (ps(iat, is, iw), is=1, ns)
         end do
      end do
      write (u, *) langham
   end subroutine gle_restout

   subroutine finalize_gle()
      if (allocated(gS)) then
         deallocate (gS, gT, ps, gp)
      end if
      if (allocated(gS_centroid)) then
         deallocate (gS_centroid, gT_centroid)
      end if
   end subroutine finalize_gle

   subroutine write_propagator(M, dt, ns, fname)
      real(DP), intent(in) :: M(:, :)
      real(DP), intent(in) :: dt
      integer, intent(in) :: ns
      character(len=*), intent(in) :: fname
      integer :: u

      open (newunit=u, file=fname, action="write", access="sequential", form="unformatted")
      write (u) dt, ns
      write (u) M
      close (u)
   end subroutine write_propagator

   subroutine read_propagator(dt, ns, fname, M)
      use mod_error, only: fatal_error
      real(DP), intent(in) :: dt
      integer, intent(in) :: ns
      character(len=*), intent(in) :: fname
      real(DP), intent(out) :: M(:, :)
      real(DP) :: dt_read
      integer :: ns_read
      integer :: u

      open (newunit=u, file=fname, action="read", status="old", access="sequential", form="unformatted")
      read (u) dt_read, ns_read
      if (dt /= dt_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "dt read from file "//fname//" does not match")
         return
      end if
      if (ns /= ns_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "ns read from file "//fname//" does not match")
         return
      end if
      read (u) M
      close (u)
   end subroutine read_propagator

   subroutine compute_propagator(A, C, ns, dt, T, S)
      real(DP), intent(in) :: A(ns + 1, ns + 1)
      real(DP), intent(in) :: C(ns + 1, ns + 1)
      integer, intent(in) :: ns
      real(DP), intent(in) :: dt
      real(DP), intent(out) :: T(ns + 1, ns + 1)
      real(DP), intent(out) :: S(ns + 1, ns + 1)
      real(DP) :: TMP(ns + 1, ns + 1)

      ! the deterministic part of the propagator
      call matrix_exp(-dt * A, ns + 1, 15, 15, T)

      ! the stochastic part
      TMP = C - matmul(T, matmul(C, transpose(T)))
      call cholesky(TMP, S, ns + 1)
   end subroutine compute_propagator

   ! The GLE propagator.
   ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
   ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
   ! gp(j,2:ns+1) the current values of additional momenta.
   ! the matrix multiplies are performed on the whole array at once, and the new momentum
   ! is passed back to the caller, while the new s's are kept stored in ps.
   subroutine gle_step(px, py, pz, m)
      use mod_general, only: natom, nwalk, inormalmodes
      use mod_utils, only: ekin_p
      real(DP), intent(inout) :: px(:, :)
      real(DP), intent(inout) :: py(:, :)
      real(DP), intent(inout) :: pz(:, :)
      real(DP), intent(in) :: m(:, :)
      integer :: i, iat, iw

      langham = langham + ekin_p(px, py, pz, m, natom, nwalk)

      do iw = 1, nwalk

         do i = 1, ns
            do iat = 1, natom * 3
               gp(iat, i + 1) = ps(iat, i, iw)
            end do
         end do

         do iat = 1, natom
            gp(iat, 1) = px(iat, iw)
            gp(iat + natom, 1) = py(iat, iw)
            gp(iat + natom * 2, 1) = pz(iat, iw)
         end do

         if (inormalmodes == 1 .and. iw == 1) then
            ! PIGLET centroid propagation
            call gle_propagate(gp, gT_centroid, gS_centroid, m, iw)
         else
            call gle_propagate(gp, gT, gS, m, iw)
         end if

         ! Pass the results back
         do iat = 1, natom
            px(iat, iw) = gp(iat, 1)
            py(iat, iw) = gp(iat + natom, 1)
            pz(iat, iw) = gp(iat + 2 * natom, 1)
         end do

         do i = 1, ns
            do iat = 1, natom * 3
               ps(iat, i, iw) = gp(iat, i + 1)
            end do
         end do

      end do

      langham = langham - ekin_p(px, py, pz, m, natom, nwalk)
   end subroutine gle_step

   subroutine gle_propagate(p, T, S, mass, iw)
      use mod_general, only: natom
      use mod_random, only: gautrg
      real(DP), intent(inout) :: p(:, :)
      real(DP), intent(in) :: T(:, :), S(:, :), mass(:, :)
      integer, intent(in) :: iw
      real(DP) :: ngp(natom * 3, ns + 1)
      real(DP) :: ran(natom * 3)
      real(DP) :: sqm
      integer :: i, j

      ngp = transpose(matmul(T, transpose(p)))

      ! now, must compute random part.
      ! first, fill up p of random n
      do i = 1, ns + 1
         call gautrg(ran, natom * 3)
         do j = 1, natom
            sqm = dsqrt(mass(j, iw))
            p(j, i) = ran(j) * sqm
            p(j + natom, i) = ran(j + natom) * sqm
            p(j + natom * 2, i) = ran(j + natom * 2) * sqm
         end do
      end do

      p = ngp + transpose(matmul(S, transpose(p)))
   end subroutine gle_propagate

   ! matrix exponential by scale & square.
   ! one can diagonalize with lapack, but it's not worth it, as
   ! we call this routine only once!
   subroutine matrix_exp(M, n, j, k, EM)
      integer, intent(in) :: n, j, k
      real(DP), intent(in) :: M(n, n)
      real(DP), intent(out) :: EM(n, n)

      real(DP) :: tc(j + 1), SM(n, n)
      integer :: p, i
      tc(1) = 1
      do i = 1, j
         tc(i + 1) = tc(i) / dble(i)
      end do

      !scale
      SM = M * (1.D0 / 2.D0**k)
      EM = 0.D0
      do i = 1, n
         EM(i, i) = tc(j + 1)
      end do

      !taylor exp of scaled matrix
      do p = j, 1, -1
         EM = matmul(SM, EM); 
         do i = 1, n
            EM(i, i) = EM(i, i) + tc(p)
         end do
      end do

      !square
      do p = 1, k
         EM = matmul(EM, EM)
      end do
   end subroutine matrix_exp

   ! TODO: replace by more stable procedure from i-Py.
   ! Brute-force "stabilized" cholesky decomposition.
   ! in practice, we compute LDL^T decomposition, and force
   ! to zero negative eigenvalues.
   subroutine cholesky(SST, S, n)
      use mod_files, only: stderr
      integer, intent(in) :: n
      real(DP), intent(in) :: SST(n, n)
      real(DP), intent(out) :: S(n, n)
      real(DP) :: L(n, n), D(n, n)
      integer :: i, j, k
      S = 0.D0
      L = 0.D0
      D = 0.D0
      do i = 1, n
         L(i, i) = 1.0D0
         do j = 1, i - 1
            L(i, j) = SST(i, j); 
            do k = 1, j - 1
               L(i, j) = L(i, j) - L(i, k) * L(j, k) * D(k, k)
            end do
            if (D(j, j) > 1.0D-10) then
               L(i, j) = L(i, j) / D(j, j)
            else
               write (stderr, '(A)') "Warning: zero eigenvalue in LDL^T decomposition."
               L(i, j) = 0.D0
            end if
         end do
         D(i, i) = SST(i, i)
         do k = 1, i - 1
            D(i, i) = D(i, i) - L(i, k)**2 * D(k, k)
         end do
      end do
      do i = 1, n
         if (D(i, i) >= 0.0D0) then
            D(i, i) = dsqrt(D(i, i))
         else
            write (stderr, *) "WARNING: negative eigenvalue (", D(i, i), ")in LDL^T decomposition."
            D(i, i) = 0.0D0
         end if
      end do
      S = matmul(L, D)
   end subroutine cholesky

end module mod_gle_private

! Public interface
module mod_gle
   use mod_gle_private
   implicit none
   private
   ! Input keywords
   public :: tau0_langevin
   public :: gle_test
   ! Public subroutines
   public :: pile_init, gle_init
   public :: finalize_gle, finalize_pile
   public :: gle_step, pile_step
   public :: gle_restout, gle_restin
   public :: pile_restout, pile_restin
   ! Public functions
   public :: get_langham
end module mod_gle
