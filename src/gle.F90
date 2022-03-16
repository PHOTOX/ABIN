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

! The original code was modified to allow for PI+GLE and PIGLET simulation
! in ABIN by Daniel Hollas, danekhollas at gmail dot com

module mod_gle
   use mod_const, only: DP
   use mod_random, only: gautrg
   implicit none
   private
   public :: readQT, ns, ps, langham, tau0_langevin
   public :: gle_test
   public :: gle_step, pile_step, pile_init, gle_init, finalize_gle, finalize_pile
   ! Public for unit tests
   public :: read_propagator, write_propagator
   real(DP), allocatable :: gS(:, :), gT(:, :)
   real(DP), allocatable :: ps(:, :, :)
   ! For PIGLET
   real(DP), allocatable :: gS_centroid(:, :), gT_centroid(:, :)
   ! the following are only temporary helper arrays
   real(DP), allocatable :: gp(:, :), ngp(:, :), ran(:)
   ! Relaxation time of the white noise Langevin thermostat (PILE)
   ! In ABIN input, it should be set in picoseconds
   real(DP) :: tau0_langevin = -1.0D0
   real(DP) :: langham = 0.0D0
   integer :: ns, ns_centroid, readQT = 1
   logical :: gle_test
   real(DP), allocatable :: c1(:), c2(:)
   save

contains
   ! Initialize white-noise PILE thermostat,
   ! which can be used both for classical MD and PIMD.
   subroutine pile_init(dt, tau0)
      use mod_const, only: PI, AUtoFS
      use mod_general, only: natom, nwalk, ipimd, inormalmodes
      use mod_nhc, only: temp
      use mod_utils, only: abinerror
      real(DP), intent(in) :: dt, tau0
      real(DP) :: gam, omega, tau
      integer :: iw

      ! Sanity check
      if (ipimd == 1 .and. inormalmodes <= 0) then
         write (*, *) 'ERROR: You must use normal mode coordinates with PILE thermostat.'
         call abinerror('pile_init')
      end if

      if (tau0 <= 0.0D0) then
         write (*, *) 'ERROR: tau0_langevin for PILE thermostat was not set or was negative.'
         write (*, *) 'Set "tau0_langevin" in picoseconds in the ABIN input in section "nhcopt".'
         call abinerror('pile_init')
      end if

      langham = 0.D0 ! sets to zero accumulator for langevin 'conserved' quantity
      allocate (ran(natom * 3 * nwalk)) !allocating arrays for random numbers produced by gautrg

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

      write (*, *) 'C1', (c1(iw), iw=1, nwalk)
      write (*, *) 'C2', (c2(iw), iw=1, nwalk)

   end subroutine

   ! white-noise propagator. time-step has been set in wn_init
   subroutine pile_step(px, py, pz, m)
      use mod_general, only: natom, nwalk
      real*8, intent(inout) :: px(:, :)
      real*8, intent(inout) :: py(:, :)
      real*8, intent(inout) :: pz(:, :)
      real*8, intent(in) :: m(:, :)
      integer :: iat, iw, pom

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
   end subroutine

   subroutine finalize_pile()
      if (allocated(ran)) deallocate (ran, c1, c2)
   end subroutine finalize_pile

   subroutine print_gle_header(inose, ipimd, inormalmodes)
      integer, intent(in) :: inose
      integer, intent(in) :: ipimd
      integer, intent(in) :: inormalmodes
      print*, ""
      print*, "Initialization of GLE thermostat."
      print*, "Please cite the following review paper about the GLE methodology"
      print*, "M. Ceriotti, G. Bussi and M. Parrinello"
      print*, "J. Chem. Theory Comput. 6, 1170 (2010)"
      print*, "Colored-noise thermostats a la carte"
      print*, "http://dx.doi.org/10.1021/ct900563s"
      print*, ""
      print*, "The citations for specific GLE use cases is at"
      print*, "http://gle4md.org/index.html?page=refs"
      print*, ""

      if (inose == 2 .and. ipimd == 0) then
         print *, "Using Quantum Thermostat, please cite"
         print*, "M. Ceriotti, G. Bussi and M. Parrinello"
         print*, "Phy. Rev. Lett. 103, 030603 (2009)"
         print*, "Nuclear quantum effects in solids using a colored-noise thermostat"
         print*, "http://dx.doi.org/10.1103/PhysRevLett.103.030603"
         print*, ""
      else if (inose == 2 .and. ipimd == 1 .and. inormalmodes == 0) then
         print *, "Using PI+GLE thermostat, please cite:"
         print*, "M. Ceriotti, D. E. Manolopoulos, and M. Parrinello"
         print*, "J. Chem. Phys. 134, 084104 (2011)"
         print*, "Accelerating the convergence of path integral dynamics"
         print*, "with a generalized Langevin equation"
         print*, "http://dx.doi.org/10.1063/1.3556661"
         print*, ""
      else if (inose == 2 .and. ipimd == 1 .and. inormalmodes == 1) then
         print *, "Using PIGLET thermostat, please cite"
         print*, "M. Ceriotti, D. E. Manolopoulos, Rev. Lett. 109, 100604 (2012)"
         print*, "Efficient first-principles calculation of the"
         print*, "quantum kinetic energy and momentum distribution of nuclei"
         print*, "http://dx.doi.org/10.1103/PhysRevLett.109.100604"
         print*, ""
      end if
   end subroutine print_gle_header

   subroutine gle_init(dt)
      use mod_const, only: AUtoEV
      use mod_general, only: natom, nwalk, ipimd, inormalmodes, my_rank, iremd
#if __GNUC__ == 0
      use mod_general, only: irest
#endif
      use mod_utils, only: abinerror
      use mod_nhc, only: temp, inose
      implicit none
      real(DP), intent(in) :: dt
      real(DP), allocatable :: gA(:, :), gC(:, :), gr(:)
      real(DP), allocatable :: gA_centroid(:, :), gC_centroid(:, :)
      integer :: i, cns, ios, iw
      character(len=10) :: glea, glec
      character(len=2) :: char_my_rank

      if (my_rank == 0) then
         call print_gle_header(inose, ipimd, inormalmodes)
      end if

      ! reads in matrices
      ! reads A (in a.u. units)
      if (iremd == 1) then
         if (my_rank == 0) then
            write (*, *) "iremd=1 - Expecting matrices in form:"//&
                        &" GLE-A.id_of_replica, GLE-C.id_of_replica (ie ./GLE-A.00)"
         end if
         write (char_my_rank, '(I0.2)') my_rank
         glea = 'GLE-A.'//char_my_rank
         glec = 'GLE-C.'//char_my_rank
      else
         glea = 'GLE-A'
         glec = 'GLE-C'
      end if

      open (121, file=glea, status='OLD', iostat=ios, action='read')
      if (ios /= 0) then
         write (0, *) "Error: could not read GLE-A file!"
         write (0, *) "Exiting..."
         call abinerror('gle_init')
      end if

      ! try to open GLE-C file
      open (122, file=glec, status='OLD', action='read', iostat=ios)
      if (ios /= 0 .and. inose == 2) then
         write (0, *) "Error: could not read GLE-C file!"
         call abinerror("gle_init")
      end if

      if (inormalmodes == 1) then
         read (121, *) ns
         allocate (gA_centroid(ns + 1, ns + 1))
         allocate (gC_centroid(ns + 1, ns + 1))
         allocate (gS_centroid(ns + 1, ns + 1))
         allocate (gT_centroid(ns + 1, ns + 1))
         if (my_rank == 0) write (6, *) '# Reading A-matrix for centroid. Expecting a.u. units!'
         do i = 1, ns + 1
            read (121, *) gA_centroid(i, :)
         end do
         ns_centroid = ns

         ! Read C matrix for centroid
         read (122, *) ns
         if (ns /= ns_centroid) then
            write (*, *) 'ERROR: Inconsistent size of A and C matrices for centroid.'
            call abinerror("gle_init")
         end if
         if (my_rank == 0) write (6, *) '# Reading C-matrix for centroid. Expecting eV units!'
         do i = 1, ns + 1
            read (122, *) gC_centroid(i, :)
            gC_centroid(i, :) = gC_centroid(i, :) / AUtoEV
         end do

         call compute_propagator(gA_centroid, gC_centroid, gT_centroid, gS_centroid, dt)

         ! For initialization of momenta, see below
         gA_centroid = gC_centroid
         call cholesky(gA_centroid, gC_centroid, ns + 1)
      end if

      read (121, *) ns

      if (inormalmodes == 1 .and. ns /= ns_centroid) then
         write (*, *) 'ERROR: Size of A matrix for centroid and other normal modes does not match!'
         write (*, *) 'Please, double check file GLE-A!'
         call abinerror('gle_init')
      end if

      !allocate everything we need
      if (natom * 3 > ns + 1) then
         allocate (ran(natom * 3)) !allocating arrays for random numbers produced by gautrg
      else
         allocate (ran(ns + 1))
      end if

      allocate (gA(ns + 1, ns + 1))
      allocate (gC(ns + 1, ns + 1))
      allocate (gS(ns + 1, ns + 1))
      allocate (gT(ns + 1, ns + 1))
      allocate (gp(natom * 3, ns + 1))
      allocate (ngp(natom * 3, ns + 1))
      allocate (gr(ns + 1))
      allocate (ps(natom * 3, ns, nwalk)) !each bead has to have its additional momenta

      if (my_rank == 0) print*, 'Reading A-matrix. Expecting a.u. units!'
      do i = 1, ns + 1
         read (121, *) gA(i, :)
      end do

      close (121)

      ! reads C (in eV!), or init to kT
      if (inose == 4) then
         if (my_rank == 0) print*, "Using canonical-sampling, Cp=kT"
         gC = 0.0D0
         do i = 1, ns + 1
            gC(i, i) = temp
         end do
      else
         if (my_rank == 0) print*, 'Reading specialized Cp matrix. Expecting eV units!'
         read (122, *) cns
         if (cns /= ns) then
            write (*, *) "Error: size mismatch matrices in GLE-A and GLE-C!"
            call abinerror('gle_init')
         end if
         do i = 1, ns + 1
            read (122, *) gC(i, :)
            gC(i, :) = gC(i, :) / AUtoEV
         end do
      end if
      close (122)

      ! WARNING: gA is rewritten here
      ! TODO: do not overwrite gA
      if (gle_test) then
         call read_propagator(gT, gS, dt, ns)
      else
         call compute_propagator(gA, gC, gT, gS, dt)
      end if

      ! Initialize the auxiliary vectors.
      ! we keep general - as we might be using non-diagonal C
      ! to break detailed balance - and we use cholesky decomposition of C
      ! since one would like to initialize correctly the velocities in

      ! DH: ps rewritten in init.f90 if irest.eq.1

      ! TODO: Do not overwrite gA
      gA = gC
      call cholesky(gA, gC, ns + 1)

      do iw = 1, nwalk
         ! Centroid was already initialized
         if (inormalmodes == 1 .and. iw == 1) then
            call initialize_momenta(gC_centroid, 1)
         else
            call initialize_momenta(gC, iw)
         end if
      end do

      langham = 0.D0 ! sets to zero accumulator for langevin 'conserved' quantity

      deallocate (gA)
      deallocate (gC)
      if (inormalmodes == 1) then
         deallocate (gA_centroid)
         deallocate (gC_centroid)
      end if
#if __GNUC__ == 0
      if (irest /= 0 .and. ns > 6) then
         write (*, *) 'ERROR: Restarting GLE thermostat with ns>6&
          & is not supported with non-GNU compilers.'
         call abinerror('gle_init')
      end if
#endif
   end subroutine gle_init

   subroutine initialize_momenta(C, iw)
      !use mod_arrays,  only: px, py, pz, vx, vy, vz, amt
      use mod_error, only: fatal_error
      use mod_general, only: natom
      use mod_utils, only: abinerror
      real(DP), intent(in) :: C(:, :)
      integer, intent(in) :: iw
      real(DP), allocatable :: gr(:)
      integer :: i, j
      ! WARNING: this routine must be called after arrays are allocated!
      if (.not. allocated(ps)) then
         call fatal_error(__FILE__, __LINE__, "oh no, programming error in gle.F90!")
      end if

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

   subroutine finalize_gle()
      if (allocated(gS)) then
         deallocate (gS, gT, ps, gp, ngp, ran)
      end if
      if (allocated(gS_centroid)) then
         deallocate (gS_centroid, gT_centroid)
      end if

      if (allocated(ran)) deallocate (ran)
   end subroutine finalize_gle

   subroutine write_propagator(T, S, dt, ns)
      real(DP), intent(out) :: T(:, :), S(:, :)
      real(DP), intent(in) :: dt
      integer, intent(in) :: ns
      integer :: u

      open (newunit=u, file='GLE-T', action="write", access="sequential", form="unformatted")
      write (u) dt, ns
      write (u) T
      close (u)

      open (newunit=u, file='GLE-S', action="write", access="sequential", form="unformatted")
      write (u) dt, ns
      write (u) S
      close (u)
   end subroutine write_propagator

   subroutine read_propagator(T, S, dt, ns)
      use mod_error, only: fatal_error
      real(DP), intent(out) :: T(:, :), S(:, :)
      real(DP), intent(in) :: dt
      integer, intent(in) :: ns
      real(DP) :: dt_read
      integer :: ns_read
      integer :: u

      open (newunit=u, file='GLE-T', action="read", status="old", access="sequential", form="unformatted")
      read (u) dt_read, ns_read
      if (dt /= dt_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "dt read from GLE-T does not match")
         return
      end if
      if (ns /= ns_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "ns read from GLE-T does not match")
         return
      end if
      read (u) T
      close (u)

      open (newunit=u, file='GLE-S', action="read", status="old", access="sequential", form="unformatted")
      read (u) dt_read, ns_read
      if (dt /= dt_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "dt read from GLE-S does not match")
         return
      end if
      if (ns /= ns_read) then
         close (u)
         call fatal_error(__FILE__, __LINE__, "ns read from GLE-S does not match")
         return
      end if
      read (u) S
      close (u)
   end subroutine read_propagator

   ! Matrix A is rewritten on output
   subroutine compute_propagator(A, C, T, S, dt)
      real(DP), intent(inout) :: A(:, :)
      real(DP), intent(in) :: C(:, :)
      real(DP), intent(out) :: T(:, :), S(:, :)
      real(DP), intent(in) :: dt

      ! the deterministic part of the propagator
      call matrix_exp(-dt * A, ns + 1, 15, 15, T)

      ! the stochastic part, we use A as a temporary array
      ! TODO: Do not overwrite A, makes things confusing
      A = C - matmul(T, matmul(C, transpose(T)))
      call cholesky(A, S, ns + 1)

   end subroutine compute_propagator

   ! The GLE propagator.
   ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
   ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
   ! gp(j,2:ns+1) the current values of additional momenta.
   ! the matrix multiplies are performed on the whole array at once, and the new momentum
   ! is passed back to the caller, while the new s's are kept stored in ps.
   subroutine gle_step(px, py, pz, m)
      use mod_general, only: natom, nwalk, inormalmodes
      real(DP), intent(inout) :: px(:, :)
      real(DP), intent(inout) :: py(:, :)
      real(DP), intent(inout) :: pz(:, :)
      real(DP), intent(in) :: m(:, :)
      integer :: i, iat, iw

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

         !iw enddo
      end do

   end subroutine gle_step

   subroutine gle_propagate(p, T, S, mass, iw)
      use mod_general, only: natom
      real(DP), intent(inout) :: p(:, :)
      real(DP), intent(in) :: T(:, :), S(:, :), mass(:, :)
      real(DP) :: sqm
      integer, intent(in) :: iw ! this one is stupid
      integer :: i, j
      ! ran and ngp arrays are allocated in gle_init,
      ! maybe we should move the allocation here

      ngp = transpose(matmul(T, transpose(p)))

      ! now, must compute random part.
      ! first, fill up p of random n
      do i = 1, ns + 1
         call gautrg(ran, natom * 3)
         do j = 1, natom
            sqm = dsqrt(mass(j, iw))
            !<-- if m!= 1, alternatively one could perform the scaling here (check also init!)
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
      real*8, intent(in) :: M(n, n)
      real*8, intent(out) :: EM(n, n)

      real*8 :: tc(j + 1), SM(n, n)
      integer p, i
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
      integer, intent(in) :: n
      real*8, intent(in) :: SST(n, n)
      real*8, intent(out) :: S(n, n)
      real*8 :: L(n, n), D(n, n)
      integer i, j, k
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
               write (*, *) "Warning: zero eigenvalue in LDL^T decomposition."
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
            write (*, *) "Warning: negative eigenvalue (", D(i, i), ")in LDL^T decomposition."
            D(i, i) = 0.0D0
         end if
      end do
      S = matmul(L, D)
   end subroutine cholesky

end module mod_gle
