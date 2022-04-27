!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!  The following code was kindly provided by prof. Marius Lewerenz.
!  Version of this code is included in the program MCMC2.
!  See the following article for details.
!  Bonhommeau, D. A.; Lewerenz, M. and Gaigeot, M. 
!  MCMC2 (version 1.1): A {Monte Carlo} code for multiply-charged clusters
!  Comp. Phys. Comm. , 185, 2014, pp. 1188-1191
!  doi: 10.1016/j.cpc.2013.09.026
!
!  set of random number routines extracted from toolbox
!                                          m. lewerenz may/2001
!-----------------------------------------------------------------------
module mod_random
!---------------------------- ranlfg.in! -------------------------------
!---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
!
!      parameters for lagged fibonacci generators and  generator state
!
!-----------------------------------------------------------------------
!
!      possible (np,nq) values, (np,np-nq) is also valid:
!              (17,5), (250,103), (521,158), (1279,418),
!              (2281,715), (4423,1393), (1279,1063)
!      ref.: Bhanot et al., phys. rev b 33, 7841 (1986);
!            Zierler, inf. control 15, 67 (1961)
!
!      mersenne prime primitive trinomials:
!      (heringa et al. int.j.mod.phys.c 3, 561 (1992))
!
!      (89,38)
!      (127,1), (127,7), (127,15), (127,30), (127,63)
!      (521,32), (521,48), (521,158), (521,168)
!      (607,105), (607,147), (607, 273)
!      (1279,216), (1279,418)
!      (2281,715), (2281,915), (2281,1029)
!      (3217,67), (3217,576)
!      (4423,271), (4423,369), (4423,370), (4423,649), (4424,1393),
!                  (4423,1419), (4423,2098)
!      (9689,84), (9689,471), (9689,1836), (9689,2444), (9689,4187)
!      (19937,881), (19937,7083), (19937,9842)
!      (23209,1530), (23209,6619), (23209,9739)
!      (44497,8575), (44497,21034)
!      (110503,25230), (110503,53719)
!      (132049,7000), (132049,33912), (132049,41469), (132049,52549),
!                     (132049,54454)
!
!      another pair from brent92 who recommends q=0.618p : (258,175) 
!      brent's ranu4 uses (132049,79500)

!-----------------------------------------------------------------------
!     parameter (np=250,nq=103)
!     parameter (np=2281,nq=715)
!     parameter (np=4423,nq=1393)
      use mod_const, only: DP
      use mod_error, only: fatal_error
      use mod_files, only: stdout, stderr
      private
      public :: gautrg, vranf
      public :: write_prng_state, read_prng_state
      integer,parameter :: np=1279, nq=418
      real(DP)  :: x(np)
      integer   :: last, init
      ! gautrg variables determining its state
      real(DP)  :: gsave
      integer   :: isave = -1
      ! For restart file
      character(len=*), parameter :: chprng='PRNG STATE (OPTIONAL)'
      save
      contains


      subroutine gautrg(gran, nran, iseed)
!     DH WARNING: initialiazation from gautrg and vranf are different
!     gautrg gran is shifted by 1
!     Apparently, nran must be > 1
!
!      vectorized portable unit standard deviation gaussian random
!      number generator using the box-mueller transformation method.
!      this method is faster than the central limit method, when the
!      uniform random numbers are comparatively expensive.
!      version working on a square with sine and cosine functions. the
!      same sequence is produced independent of the sequence of calls
!      if no calls to other vranf based prn generators are made.
!
!      gran  : vector of length nran returning the random numbers
!      nran  : number of desired random numbers, with nran=0 and iseed
!              not 0 only generator initialization is done. no action
!              when both are zero.
!      iseed : (OPTIONAL) a positive non-zero integer to start generator.
!
!      times for the generation of 10**6 prn`s:
!      rs6000/350 41 mhz    2.31 s 
!      r8000      75 mhz    1.43 s
!      axp21064  200 mhz    1.20 s
!      axp21164  250 mhz    0.572 s
!      t3e-900   450 mhz    0.52 s
!      cray-t90  450 mhz    0.0728 s   1.28*10**7 prn/s
!
!      subroutines called: vranf, r1mach
!      m. lewerenz 6/may/90, modified 17/jun/91, mar/95
 
      implicit none
      integer, intent(in) :: nran
      integer, intent(in), optional :: iseed
      real(DP), intent(inout) :: gran(nran)

      real(DP),parameter  :: two=2.d0,twom=-two,one=1.d0
      integer,parameter :: npoly=11
      real(DP) :: scf(npoly),xran(1),ccf(npoly)
      real(DP) :: cx, cxi, sx, sxi, y
      real(DP) :: fac, trig
      real(DP) :: tiny, twopi, pi4
      integer :: i, newran
      save tiny,twopi,pi4,scf,ccf

!      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
!      FOR COS(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14

       DATA ccf/ 0.9999999999999986D+00, 0.6612476846390664D-13, &
               -0.4999999999989523D+00,-0.5434658088910759D-10, &
                0.4166666737609693D-01,-0.4648977428396692D-08, &
               -0.1388871052129944D-02,-0.4228394738587799D-07, &
                0.2486361945804866D-04,-0.5317743184071916D-07, &
               -0.2539224676809412D-06/

!      POLYNOMIAL FROM CHEBYSHEV APPROXIMATION ON [ 0.000, 0.790]
!      FOR SIN(X) WITH ABSOLUTE ERROR LESS THAN 0.2220E-14

      DATA scf/-0.9443414574112906D-15, 0.1000000000000244D+01, &
               -0.1224236196202217D-10,-0.1666666664242968D+00, &
               -0.2471495821870120D-08, 0.8333348067492644D-02, &
               -0.5482536616811601D-07,-0.1982815612039858D-03, &
               -0.2017619095413939D-06, 0.2948964053761139D-05, &
               -0.1051448397925916D-06/

      if(isave.lt.0) then
        isave=0
        tiny = r1mach()**2
        pi4=atan(one)
        twopi=pi4*8
      end if

      if (present(iseed)) then
         if (iseed < 0) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Random number seed must be a positive integer!')
            return
         end if
         call vranf(gran, 0, iseed)
      end if

      if(nran.gt.0) then
        newran=nran-isave
        if(isave.eq.1) gran(1)=gsave
        call vranf(gran(isave+1), newran)
        do 100 i=1,newran-1,2
        fac = sqrt(twom*log(gran(isave+i)+tiny))
        !DHmod: renamed x to y to avoid conflict with x(np) from module
        y=pi4*gran(isave+i+1)
        cx=(((((((((ccf(11)*y+ccf(10))*y+ccf(9))*y+ccf(8))*y &
                            +ccf(7))*y+ccf(6))*y+ccf(5))*y  &
                            +ccf(4))*y+ccf(3))*y+ccf(2))*y+ccf(1)
        sx=(((((((((scf(11)*y+scf(10))*y+scf(9))*y+scf(8))*y &
                            +scf(7))*y+scf(6))*y+scf(5))*y  &
                            +scf(4))*y+scf(3))*y+scf(2))*y+scf(1)

        sxi=(two*cx)*sx
        cxi=(two*cx)*cx-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        sxi=(two*cxi)*sxi
        cxi=(two*cxi)*cxi-one
        gran(isave+i)=fac*sxi
        gran(isave+i+1)=fac*cxi
  100   continue

        if(mod(newran,2).eq.1) then
          call vranf(xran, 1)
          fac=sqrt(twom*log(gran(nran)+tiny))
          trig=twopi*xran(1)
          gran(nran)=fac*sin(trig)
          gsave=fac*cos(trig)
          isave=1
        else
          isave=0
        end if
      end if
      return
      end subroutine gautrg
!-----------------------------------------------------------------------
!-------------------- ranf/vranf uniform random package ----------------
!-----------------------------------------------------------------------

      subroutine vranf(ranv, nran, iseed)

!      machine independent portable uniform random number generator for
!      the interval [0,1) based on a floating point subtractive lagged
!      fibonacci method similar to the feedback shift register method
!      proposed by kirkpatrick/stoll. subtractive (or additive) version
!      gives much better prn quality than the original xor operation.
!      an additive variant is given in v r a n f a.

!      v r a n f  and r a n f use the same method and can be used to
!      work together on the same random sequence. either of them can be
!      used for initialization. the state of the generator can be saved
!      or retrieved for restart with routine r s a v e f.
! ----------------------------------------------------------------------

!      ranv  : vector of length nran for random numbers; output
!      nran  : number of desired random numbers. nran=0 and iseed.ne.0
!              -> only generator initialization. no action when both 
!              are zero.; input
!      iseed : (OPTIONAL) positive non-zero integer to start generator
!
! ----------------------------------------------------------------------
!      method:
!      x(k+np)=x(k)-x(k+np-nq), initial x array generated by xuinit.
!      ieee standard requires real(DP) to have at least 48 mantissa bits.
!      with nbit=48 this generator is entirely machine independent
!      and will always give the same random sequence. you can change
!      the period of the generator by setting a different nbit value or
!      changing np and nq appropriately.

! ----------------------------------------------------------------------
!      this is a floating implementation of a generator described in:
!      M.H. KALOS, P.A. WHITLOCK, MONTE CARLO METHODS, APPENDIX,
!                  WILEY 1986
!      D.W. HEERMANN, COMPUTER SIMULATION METHODS, 2ND ED.,SPRINGER 199
!                  APPENDIX A1
!      d. stauffer, f.w. hehl, v. winkelmann, j.g. zabolitzky,
!                  computer simulation & computer algebra, section 2.2

!      original references:
!      s. kirkpatrick, e.p. stoll, j. comput. phys. 40, 517 (1981)
!      r.c. tausworthe, random numbers generated by linear recurrence
!                  modulo 2, math. comp. 19, 201 (1965)
!      t.g. lewis, w.h. payne, generalized feedback shift register
!                  pseudorandom number algorithm, j. acm 20, 456 (1973)

! ----------------------------------------------------------------------
!      other (np,nq) values: (17,5), (250,103), (521,158), (1279,418),
!                            (2281,715), (4423,1393), (1279,1063)
!          ref.: Bhanot et al., phys. rev b 33, 7841 (1986);
!                Zierler, inf. control 15, 67 (1961)
! ----------------------------------------------------------------------
!      alternative additive formulation bypassing if statements:
!      temp=x(k)+x(k+np-nq)
!      x(k)=temp-float(int(temp))

!      alternative subtractive formulation bypassing if statements:
!      temp=x(k)-x(k+np-nq)
!      x(k)=(temp+one)-float(int(temp+one))
!----------------------------------------------------------------------
!      subroutines called: xuinit            m. lewerenz may/91 & nov/93

      implicit none
      real(DP), intent(out) :: ranv(nran)
      integer, intent(in) :: nran
      integer, intent(in), optional :: iseed
      integer, parameter :: nratio=np/nq, nexec=4
      real(DP), parameter :: zero = 0.0D0, one = 1.0D0
      real(DP) :: x1
      integer :: i, j, k, left, loop, limit

      if (present(iseed)) then
         if (iseed <= 0) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Random number seed must be a positive integer!')
            return
         end if
         ! table initialization by xuinit
         call xuinit(x, np, nq, 0, nexec, iseed, init, last)
      end if

      if (nran > 0 .and. init == 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Incorrect initialization in vranf')
         return
      end if

      !  fibonacci generator updates elements of x in a cyclic fashion
      !  and copies them into ranv in blocks of max. length np.
      !  loop split into chunks of max. length nq to avoid recurrence.

      if (nran > 0) then
         j = 0
         left = nran
10       continue
         loop = min(nq, left + last) - last

         do 500 i = last + 1, last + loop
            x1 = x(i) - x(i + np - nq)
            if (x1 < zero) x1 = x1 + one
            x(i) = x1
            j = j + 1
            ranv(j) = x1
  500    continue

         if (last < nratio * nq) then
            do 650 k = 1, nratio - 1
               limit = min((k + 1) * nq, left + last)
               loop = limit - max(k * nq, last)

               do 600 i = max(k * nq, last) + 1, max(k * nq, last) + loop
                  x1 = x(i) - x(i - nq)
                  if (x1 < zero) x1 = x1 + one
                  x(i) = x1
                  j = j + 1
                  ranv(j) = x1
  600      continue
  650      continue
         end if

         limit = min(np, left + last)
         loop = limit - max(nratio * nq, last)

         do 700 i = max(nratio * nq, last) + 1, max(nratio * nq, last) + loop
            x1 = x(i) - x(i - nq)
            if (x1 < zero) x1 = x1 + one
            x(i) = x1
            j = j + 1
            ranv(j) = x1
  700    continue

         last = mod(limit, np)
         left = nran - j
         if (left > 0) goto 10
      end if
      end subroutine vranf

!-----------------------------------------------------------------------

      subroutine xuinit(y, np, nq, mode, nexec, iseed, init, last)

!      initializes a (np,nq) lagged fibonacci generator table
!      with random bits generated by a congruential generator using
!      l'ecuyers decomposition. ref.: bratley et al. p. 214

!      y     : vector of length np for initial random number table;
!              output
!      np,nq : parameters p and q of feedback shift register generator;
!              input
!      mode  : operation for lfg:
!              mode=<0 -> subtractive generator, mode=1 additive
!      nexec : number of warm up cycles for the table. nexec*nbit*np
!              random numbers are generated and discarded; input
!      iseed : integer seed for congruential generator generating
!              the bits of the initial table entries; input
!      init  : returns updated seed of congruential generator.
!              0 if table was not initialized, > 0 if ok; output
!      last  : pointer to the last used number in the table; output
!      subroutines called : none             m. lewerenz mar/93, mar/98

      implicit none
      real(DP), intent(inout) :: y(np)
      integer, intent(in) :: np, nq, mode, nexec, iseed
      integer, intent(inout) :: init, last
      integer, parameter :: ia=40692, ib=52774, ic=3791, ip=2147483399
      integer, parameter :: iphalf = ip / 2, nbit = 48
      real(DP), parameter :: half = 0.5d0, zero = 0.d0
      real(DP) :: add
      integer :: ix, i, k1, j
      logical :: high

      if (nq >= np .or. iseed == 0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Illegal seed parameter(s) in xuinit')
         return
      end if

      ! set table to zero and exercise the bit generator a little

      ix = iabs(iseed)

      do i = 1, np
         x(i) = zero
         k1 = ix / ib
         ix = ia * (ix - k1 * ib) - k1 * ic
         if (ix < 0) ix = ix + ip
      end do

      !  assemble np numbers with mantissa length nbit from random bits
      ! 'high' toggle compensates for bias from odd ip

      high = .true.
      do i = 1, np
         add = half
         do j = 1, nbit
            k1 = ix / ib
            ix = ia * (ix - k1 * ib) - k1 * ic
            if (ix < 0) ix = ix + ip
            if (high) then
               if (ix >= iphalf) y(i) = y(i) + add
               high = .false.
            else
               if (ix > iphalf) y(i) = y(i) + add
               high = .true.
            end if
            add = add * half
         end do
      end do

      if (nexec > 0) then
         call xuwarm(y, np, nq, mode, nbit * nexec)
      end if

      init = ix
      last = 0
      end subroutine xuinit

!-----------------------------------------------------------------------

      subroutine xuwarm(y, np, nq, mode, nexec)

!      warms up (p,q) lagged fibonacci generators (lfg) by nexec rounds

!      y     : vector of length np for initial random number table;
!              output
!      np,nq : parameters p and q of feedback shift register generator;
!              input
!      mode  : operation for lfg:
!              mode=<0 -> subtractive generator, mode=1 additive
!      nexec : number of warm up cycles for the table. nexec*nbit*np
!              random numbers are generated and discarded; input
!      subroutines called: x                         m. lewerenz mar/98

      implicit none
      ! DHmod, renamed x to y to prevent collision with module
      real(DP), intent(inout) :: y(np)
      real(DP), parameter :: zero = 0.d0, one = 1.d0
      integer, intent(in) :: mode, nexec, np, nq
      integer :: i, k

      if (nq.ge.np.or.np.eq.0.or.nq.eq.0) then
         call fatal_error(__FILE__, __LINE__, &
            & 'Illegal table parameter(s) in xuwarm')
         return
      else

!      exercise the generator for nexec rounds of np prn's
!      separate sections for subtractive or additive version

          if(mode.le.0) then
            do k=1,nexec
              do i=1,nq
                y(i)=y(i)-y(i+np-nq)
                if(y(i).lt.zero) y(i)=y(i)+one
              end do
              do i=nq+1,np
                y(i)=y(i)-y(i-nq)
                if(y(i).lt.zero) y(i)=y(i)+one
              end do
            end do
          else
            do k=1,nexec
              do i=1,nq
                y(i)=y(i)+y(i+np-nq)
                if(y(i).ge.one) y(i)=y(i)-one
              end do
              do i=nq+1,np
                y(i)=y(i)+y(i-nq)
                if(y(i).ge.one) y(i)=y(i)-one
              end do
            end do
          end if
      end if
      return
      end subroutine xuwarm


      function R1MACH()
      real(DP), parameter :: ONE=1.D0, TWO=2.D0, HALF=0.5D0
      real(DP) :: R1MACH
      real(DP) :: eps
      real(DP) :: U, COMP
      integer :: icall
      SAVE ICALL,EPS
      DATA ICALL,EPS/0,ONE/
! ---------------------------------------------------------------------
! THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE.
! THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
! U SUCH THAT 1.0 + U .NE. 1.0E0
! ---------------------------------------------------------------------
      IF(ICALL.EQ.0) THEN
        ICALL=1
        U = ONE
  10    U = U*HALF
        COMP = ONE + U
        IF(COMP .NE. ONE) GO TO 10
        EPS = U*TWO
      END IF
      R1MACH = EPS
      RETURN
      END FUNCTION R1MACH

      ! Write PRNG state into opened restart file
      subroutine write_prng_state(uout)
         integer, intent(in) :: uout
         integer :: i

         write (uout,'(A)') chprng
         write (uout, '(I0,X,I0)') init, last
         write (uout, '(I0,X,ES24.16E3)') isave, gsave
         do i = 1, np
            write (uout, '(ES24.16E3)') x(i)
         end do
      end subroutine write_prng_state

      ! Read PRNG state from opened restart file
      subroutine read_prng_state(uin)
         integer, intent(in) :: uin
         character(len=len(chprng) + 20) :: readstring
         integer :: i, iost

         read (uin, '(A)', iostat=iost) readstring
         ! The assumption here is that the PRNG state is at the
         ! end of the restart file. It is okay if it is missing,
         ! we simply use the seed from input file.
         if (iost /= 0) then
            write (stderr, *) 'WARNING: PRNG state not present in restart.xyz'
            write (stderr, *) 'Random seed from input file will be used.'
         else if (chprng /= trim(adjustl(readstring))) then
            call fatal_error(__FILE__, __LINE__, &
               & 'Unexpected line in restart file when trying to read prng state.'//&
               & new_line('a')//'Expected: '//chprng//&
               & new_line('a')//'Got: '//trim(adjustl(readstring)))
         else
            read (uin, *) init, last
            read (uin, *) isave, gsave
            do i = 1, np
               read (uin, *) x(i)
            end do
         end if
      end subroutine read_prng_state

end module mod_random

module mod_prng_init
   use, intrinsic :: iso_fortran_env, only: INT64
   use mod_const, only: DP
   use mod_files, only: stdout, stderr
   ! We initialize our PRNG via call to gautrg
   ! with optional seed parameter.
   use mod_random, only: gautrg
   implicit none
   private
   public :: initialize_prng
   ! Used in utils/abin-randomint.f90
   public :: get_random_seed, initialize_fortran_prng, random_ints

contains

   ! PRNG INITIALIZATION
   ! We want the simulations to be repeatable so we ask the user
   ! to provide a single random seed. Single random seed is all we need
   ! to initialize our PRNG implemented in the vranf routine in mod_random.
   !
   ! We also allow the user to not provide a seed at all (or provide negative seed),
   ! in which case we determine it from /dev/urandom.
   ! We print it to the output so that the simulation can be exactly repeated if needed.
   ! (in this case repeatability is not possible for REMD).
   !
   ! Things get complicated for REMD, where we need to initialize
   ! nreplica independent prngs, but we still want only a single seed
   ! from the user, and we want repeatability. We solve this by generating
   ! additional seeds by a another much more simple PRNG in function lcg(),
   ! taken from the Gfortran documentation.
   !
   ! We also implement additional high-quality integer PRNG, random_ints(),
   ! which uses the random_number() subroutine, as defined in Fortran standard.
   ! Unfortunately, we cannot use this routine for deterministic generation of REMD seeds,
   ! because random_number() itself needs to be seeded by a call to random_seed().
   ! Here's the kicker: random_seed() in general requires an array of seeds :-(
   ! https://gcc.gnu.org/onlinedocs/gcc-4.9.4/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
   ! Also, random_number() implementation depends on the compiler,
   ! so it is not a good fit for out test suite.
   subroutine initialize_prng(seed, mpi_rank)
      integer, intent(inout) :: seed
      integer, intent(in) :: mpi_rank
      integer(INT64) :: new_seed
      integer :: i
      real(DP) :: drans(1)

      write (stdout, *) 'Initializing pseudo-random number generator'

      ! If the user does not provide the seed, or provides negative one,
      ! we will determine it automatically. Note that the seed is
      ! passed back so that it can be printed in init.F90
      if (seed <= 0) then
         seed = get_random_seed()
         if (mpi_rank > 0) then
            write (*, '(A,I0,A,I0)') 'MPI rank = ', mpi_rank, ' Seed = ', seed
         else
            write (stdout, '(A,I0,A,I0)') 'MPI rank = ', mpi_rank, ' Seed = ', seed
         end if
      end if

      ! Generate different random number seeds for different MPI processes.
      if (mpi_rank > 0) then
         new_seed = int(seed, kind(new_seed))
         do i = 1, mpi_rank
            new_seed = lcg(new_seed)
         end do
         seed = int(new_seed, kind(seed))
      end if

      ! Initializes PRNG for integers, random_ints(),
      ! currently not in use in ABIN, but used in utils/abin-randomint.f90
      call initialize_fortran_prng(seed)

      ! Initialize our custom pseudo-random number generator.
      ! This call has to happen before we read the restart file,
      ! to allocate internal arrays.
      ! If we are restarting, the PRNG state initialized here is overwritten.
      call gautrg(drans, 0, seed)
   end subroutine initialize_prng

   ! If the user doesn't provide initial random seed via
   ! 'irandom' variable in the input file, we will take it
   ! from unix kernel via /dev/urandom.
   ! https://linux.die.net/man/4/urandom
   ! As fallback, we use current date/time and PID.
   !
   ! Code taken from:
   ! https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
   !
   ! For REMD, this will make the seed unique for each replica.
   ! Unfortunatelly, we don't currently provide a way for a user
   ! to set custom seed for each replica separately so the approach
   ! used here makes the REMD simulations non-repeatable.
   ! To have repeatable REMD simulations, user should provide irandom
   ! in the input file such that this function is not called.
   integer function get_random_seed() result(seed)
      integer :: un, istat
      integer :: getpid, pid
      integer, dimension(8) :: dt
      integer(INT64) :: rate
      integer(INT64) :: t

      open (newunit=un, file="/dev/urandom", access="stream", &
         & form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         write (stdout, *) 'Getting random seed from /dev/urandom'
         read (un) seed
         close (un)
         seed = iabs(seed)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         write (stdout, *) 'Could not open /dev/urandom'//new_line('A')//&
            & 'Using date and time to get random seed'
         ! https://gcc.gnu.org/onlinedocs/gfortran/SYSTEM_005fCLOCK.html
         call system_clock(count=t, count_rate=rate)
         if (rate == 0) then
            call date_and_time(values=dt)
            ! Seconds of Unix time
            t = (dt(1) - 1970) * 365_INT64 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_INT64 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24_INT64 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
         end if
         pid = getpid()
         seed = int(ieor(t, int(pid, kind(t))))
         seed = abs(seed)
      end if
   end function get_random_seed

   ! Initializing PRNG subroutine random_number() defined in Fortran standard
   ! https://stackoverflow.com/questions/51893720/correctly-setting-random-seeds-for-repeatability
   subroutine initialize_fortran_prng(user_seed)
      integer, intent(in) :: user_seed
      integer, allocatable :: seeds(:)
      integer :: i, seed_size
      integer(INT64) :: s
      real(DP) :: drans(100)

      call random_seed(size=seed_size)
      allocate (seeds(seed_size))

      ! We're use a simple PRNG defined below for the initial seed state
      s = int(user_seed, kind(s))
      do i = 1, seed_size
         s = lcg(s)
         seeds(i) = int(s)
      end do

      call random_seed(put=seeds)
      ! Prime the prng by discarding first 100 values
      call random_number(drans)
   end subroutine initialize_fortran_prng

   ! PRNG for integers, based on Fortran standard subroutine random_number()
   ! https://stackoverflow.com/questions/23057213/how-to-generate-integer-random-number-in-fortran-90-in-the-range-0-5
   subroutine random_ints(iran, n)
      integer, intent(out) :: iran(:)
      integer, intent(in) :: n
      real(DP) :: dran(n)

      call random_number(dran)
      iran = floor(dran * huge(n))
   end subroutine random_ints

   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
   ! Taken from:
   ! https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
   integer function lcg(s)
      integer(INT64) :: s

      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_INT64)
      end if
      s = mod(s * 279470273_INT64, 4294967291_INT64)
      lcg = int(mod(s, int(huge(0), INT64)), kind(0))
   end function lcg

end module mod_prng_init
