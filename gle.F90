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
! in ABIN by Daniel Hollas

module mod_gle
   use mod_const, only: DP
   use mod_random, only: gautrg
   implicit none
   private
   public :: readQT, ns, ps, langham
   public :: gle_step, wn_step, wn_init, gle_init, finalize_gle
   real(DP), allocatable :: gS(:,:), gT(:,:)
   real(DP), allocatable :: ps(:,:,:)
   ! For PIGLET
   real(DP), allocatable :: gS_centroid(:,:), gT_centroid(:,:)
   ! the following are only temporary helper arrays
   real(DP), allocatable :: gp(:,:), ngp(:,:), ran(:)
   real(DP)         :: wnt, wns, langham
   integer          :: ns, ns_centroid, readQT=1
   save

contains

   ! initialize white-noise thermostat. 
   ! the init and the propagator here are written with the same philosophy
   ! used for the full-fledged colored-noise stuff.
   subroutine wn_init(dt,wopt)
   use mod_general,only:natom
   use mod_nhc,only:temp
   real(DP), intent(in) :: dt,wopt
   real(DP) g
   g = 2.d0 * wopt
   wnt = exp(-dt*g)
   wns = sqrt(temp*(1.d0-wnt*wnt))
   langham = 0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
   allocate( ran(natom*3) )  !allocating arrays for random numbers produced by gautrg
   end subroutine

   ! white-noise propagator. time-step has been set in wn_init
   subroutine wn_step(px,py,pz,m)
   use mod_general, only: natom
   real*8, intent(inout)  :: px(:,:)
   real*8, intent(inout)  :: py(:,:)
   real*8, intent(inout)  :: pz(:,:)
   real*8, intent(in)     :: m(:,:)
   integer :: iat,iw,pom
   iw=1
   pom=1
   call gautrg(ran,natom*3,0,6)
   do iat=1,natom
      px(iat,iw) = wnt*px(iat,iw)+wns*ran(pom)*sqrt(m(iat,iw))  
      py(iat,iw) = wnt*py(iat,iw)+wns*ran(pom+1)*sqrt(m(iat,iw))
      pz(iat,iw) = wnt*pz(iat,iw)+wns*ran(pom+2)*sqrt(m(iat,iw))
      pom = pom + 3
   end do
   end subroutine
  
   subroutine gle_init(dt)
   use mod_const, only: AUtoEV
   use mod_general,only: natom, nwalk, inormalmodes
   use mod_utils, only: abinerror
   use mod_nhc,only: temp,inose
   implicit none
   real(DP), intent(in)  :: dt
   real(DP), allocatable :: gA(:,:), gC(:,:), gr(:)
   real(DP), allocatable :: gA_centroid(:,:), gC_centroid(:,:)
   integer :: i, j, cns, ios, iw
   

   write(6,*) "# Initialization of GLE thermostat.                           "
   write(6,*) "# Please cite the relevant works among:                       "
   write(6,*) "#                                                             "
   write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
   write(6,*) "# Phy. Rev. Lett. 102, 020601 (2009)                          "
   write(6,*) "#                                                             "
   write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
   write(6,*) "# Phy. Rev. Lett. 103, 030603 (2009)                          "
   write(6,*) "#                                                             "

   !reads in matrices
   !reads A (in a.u. units)
   open(121,file='GLE-A',status='OLD',iostat=ios,action='read')
   if (ios.ne.0)then
      write(0,*) "Error: could not read GLE-A file!"
      write(0,*) "Exiting..."
      call abinerror('gle_init')
   end if

   ! try to open GLE-C file
   open(122,file='GLE-C',status='OLD',action='read',iostat=ios)
   if (ios.ne.0.and.inose.eq.2)then
      write(0,*) "Error: could not read GLE-C file!"
      call abinerror("gle_init")
   end if

   if(inormalmodes.eq.1)then
      read(121,*) ns
      allocate(gA_centroid(ns+1,ns+1))
      allocate(gC_centroid(ns+1,ns+1))
      allocate(gS_centroid(ns+1,ns+1))
      allocate(gT_centroid(ns+1,ns+1))
      write(6,*)'# Reading A-matrix for centroid. Expecting a.u. units!!!!'
      do i=1,ns+1
         read(121,*) gA_centroid(i,:)
      enddo
      ns_centroid = ns

      ! Read C matrix for centroid
      read(122,*) ns
      if(ns.ne.ns_centroid)then
         write(*,*)'ERROR: Inconsistent size of A and C matrices for centroid.'
         call abinerror("gle_init")
      end if
      write(6,*)'# Reading C-matrix for centroid. Expecting a.u. units!!!!'
      do i=1, ns+1
         read(122,*) gC_centroid(i,:)
         gC_centroid(i,:) = gC_centroid(i,:) / AUtoEV
      enddo

      call compute_propagator(gA_centroid, gC_centroid, gT_centroid, gS_centroid, dt)
      
      ! For initialization of momenta, see below
      gA_centroid = gC_centroid   
      call cholesky(gA_centroid, gC_centroid, ns+1)
   end if

   read(121,*) ns

   if(inormalmodes.eq.1.and.ns.ne.ns_centroid)then
      write(*,*)'ERROR: Size of A matrix for centroid and other normal modes does not match!'
      write(*,*)'Double check file GLE-A!!'
      call abinerror('gle_init')
   end if

   !allocate everything we need
   if(natom*3.gt.ns+1)then 
      allocate( ran(natom*3) )  !allocating arrays for random numbers produced by gautrg
   else
      allocate( ran(ns+1) )  
   endif

   allocate(gA(ns+1,ns+1))
   allocate(gC(ns+1,ns+1))
   allocate(gS(ns+1,ns+1))
   allocate(gT(ns+1,ns+1))
   allocate(gp(natom*3,ns+1))   
   allocate(ngp(natom*3,ns+1))   
   allocate(gr(ns+1))
   allocate(ps(natom*3, ns, nwalk)) !each bead has to have its additional momenta

   
   write(6,*)'# Reading A-matrix. Expecting a.u. units!!!!'
   do i=1,ns+1
      read(121,*) gA(i,:)
   enddo


   close(121)

   !gamma for a WN langevin will be 1/tau, which would make it optimal for w=1/(2tau) angular frequency. 
   !DHmod: we DONT scale gA (which is expected to be fitted by http://gle4md.berlios.de/ ) 
   !A matrix is expected to be in atomic units of time 
   !gA=gA

   ! reads C (in eV!), or init to kT
   if(inose.eq.4)then
      write(6,*) "# Using canonical-sampling, Cp=kT"
      gC=0.0d0
      do i=1,ns+1
         gC(i,i)=temp
      enddo
   else    
      write(6,*)"# Reading specialized Cp matrix."
      write(6,*)'# Expecting eV units!!!'
      read(122,*) cns
      if (cns.ne.ns)then
          write(0,*) " Error: size mismatch matrices in GLE-A and GLE-C!"
          call abinerror('gle_init')
      endif
      do i=1, ns+1
         read(122,*) gC(i,:)
         gC(i,:) = gC(i,:) / AUtoEV
      enddo
   end if
   close(122)

   call compute_propagator(gA, gC, gT, gS, dt)


   ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
   ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
   ! of C. again, since one would like to initialize correctly the velocities in 
   ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
   ! could then use it to initialize the momentum in the calling code

   ! DH: ps rewritten in init.f90 if irest.eq.1
   ! always initialize, maybe rewritten in restart

   gA = gC   
   call cholesky(gA, gC, ns+1)

   do iw = 1, nwalk
      ! Centroid was already initialized
      if(inormalmodes.eq.1.and.iw.eq.1)then
         call initialize_momenta(gC_centroid, 1)
      else
         call initialize_momenta(gC, iw)
      end if
   end do

   langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
   
   deallocate(gA)
   deallocate(gC)
   if(inormalmodes.eq.1)then
      deallocate(gA_centroid)
      deallocate(gC_centroid)
   end if
   end subroutine gle_init

   subroutine initialize_momenta(C, iw)
   use mod_arrays,  only: px, py, pz
   use mod_general, only: natom
   use mod_utils,   only: abinerror
   real(DP), intent(in) :: C(:,:)
   integer, intent(in)  :: iw
   real(DP),allocatable :: gr(:)
   integer  :: i, j
   ! WARNING: this routine must be called after arrays are allocated!
   if(.not.allocated(ps))then
      write(*,*)"WHOOOPS: Fatal programming error int gle.F90!!"
      call abinerror("initialize_momenta")
   end if

   allocate(gr(ns+1))

   do j=1,natom*3
      call gautrg(gr,ns+1,0,6)
      ! TODO, actually pass this to initialize momenta
      gp(j,:) = matmul(C, gr)
   end do

   do j=1,natom*3
      do i=1,ns
         ps(j,i,iw) = gp(j,i+1)
      enddo
   enddo

   deallocate(gr)
   end subroutine initialize_momenta


   subroutine finalize_gle()
   if(allocated(gS))then
      deallocate(gS, gT, ps, gp, ngp, ran)
   end if
   if(allocated(gS_centroid))then
      deallocate(gS_centroid, gT_centroid)
   end if

   end subroutine finalize_gle


   subroutine compute_propagator(A, C, T, S, dt)
      ! Matrix A is rewritten on output
   real(DP), intent(inout) :: A(:,:) 
   real(DP), intent(in)    :: C(:,:)
   real(DP), intent(out)   :: T(:,:), S(:,:)
   real(DP), intent(in)    :: dt

   ! the deterministic part of the propagator is obtained in a second
   call matrix_exp(-dt*A, ns+1, 15, 15, T)
    
   ! the stochastic part is just as easy. we use gA as a temporary array
   A = C - matmul(T, matmul(C, transpose(T)) )
   call cholesky(A, S, ns+1)

   end subroutine compute_propagator


   ! The GLE propagator. 
   ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
   ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
   ! gp(j,2:ns+1) the current values of additional momenta. 
   ! the matrix multiplies are performed on the whole array at once, and the new momentum
   ! is passed back to the caller, while the new s's are kept stored in ps.
   subroutine gle_step(px, py, pz, m)
   use mod_general,only:natom, nwalk, inormalmodes
   real(DP), intent(inout)  :: px(:,:)
   real(DP), intent(inout)  :: py(:,:)
   real(DP), intent(inout)  :: pz(:,:)
   real(DP), intent(in)     :: m(:,:)
   integer                :: i, j, iat, iw

!   call printf(px,py,pz)

   do iw=1,nwalk

      do i=1,ns
         do iat=1,natom*3
            gp(iat,i+1) = ps(iat,i,iw)
         enddo
      enddo
      
      do iat=1,natom
         gp(iat,1) = px(iat,iw)         
         gp(iat+natom,1) = py(iat,iw)   
         gp(iat+natom*2,1) = pz(iat,iw) 
      enddo

      if(inormalmodes.eq.2.and.iw.eq.1)then
         ! PIGLET centroid propagation
         call gle_propagate(gp, gT_centroid, gS_centroid, m, iw)
      else
         call gle_propagate(gp, gT, gS, m, iw)
      end if

      ! Pass the results back
      do iat=1,natom
         px(iat,iw) = gp(iat,1)           
         py(iat,iw) = gp(iat+natom,1)     
         pz(iat,iw) = gp(iat+2*natom,1)   
      end do

      do i=1,ns
         do iat=1,natom*3
            ps(iat,i,iw) = gp(iat,i+1)
         enddo
      enddo

   !iw enddo
   enddo

   end subroutine gle_step


   subroutine gle_propagate(p, T, S, m, iw)
   use mod_general, only: natom, inormalmodes
   real(DP), intent(inout) :: p(:,:)
   real(DP), intent(in) :: T(:,:), S(:,:), m(:,:)
   integer, intent(in)  :: iw  ! this one is stupid
   integer :: i, j
   ! ran and ngp arrays are allocated in gle_init,
   ! maybe we should move the allocation here

#ifdef USELIBS
   call dgemm('n','t',natom*3,ns+1,ns+1,1.0d0,p,natom*3,T,ns+1,0.0d0,ngp,natom*3)
#else
   ngp = transpose(matmul(gT,transpose(p)))
#endif

   ! now, must compute random part. 
   ! first, fill up gp of random n
   do i=1,ns+1
      call gautrg(ran,natom*3,0,6)
      do j=1,natom
         !<-- if m!= 1, alternatively one could perform the scaling here (check also init!)
         p(j,i)=ran(j)*sqrt(m(j,iw))
         p(j+natom,i)=ran(j+natom)*sqrt(m(j,iw))
         p(j+natom*2,i)=ran(j+natom*2)*sqrt(m(j,iw))
      end do
   end do

#ifdef USELIBS    
   call dgemm('n','t',natom*3,ns+1,ns+1,1.0d0,p,natom*3,S,ns+1,1.0d0,gp,natom*3)
   p = ngp
#else
   p = ngp + transpose(matmul(S,transpose(p)))
#endif

   end subroutine gle_propagate


  ! matrix exponential by scale & square.
  ! one can diagonalize with lapack, but it's not worth it, as 
  ! we call this routine only once!      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real*8, intent(in)   :: M(n,n)
    real*8, intent(out)  :: EM(n,n)
    
    real *8 :: tc(j+1), SM(n,n)
    integer p, i
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1.d0/2.d0**k)
    EM=0.d0
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo
    
    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo
    
    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo
  end subroutine matrix_exp
  
  ! TODO: replace by more stable procedure from i-Py???
  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
  subroutine cholesky(SST, S, n)
    integer, intent(in)  :: n
    real*8, intent(in)   :: SST(n,n)
    real*8, intent(out)   :: S(n,n)
    real *8 :: L(n,n), D(n,n) 
    integer i,j,k
    S=0.d0
    L=0.d0
    D=0.d0
    do i=1,n
       L(i,i)=1.0d0
       do j=1,i-1
          L(i,j)=SST(i,j);
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k,k)
          enddo
          if (D(j,j).gt. 1.0d-10) then
            L(i,j)=L(i,j)/D(j,j) 
          else
            write(0,*) "Warning: zero eigenvalue in LDL^T decomposition."
            L(i,j)=0.d0
          end if
       enddo
       D(i,i)=SST(i,i)
       do k=1,i-1
          D(i,i)=D(i,i)-L(i,k)**2*D(k,k)
       end do
    enddo
    do i=1,n
       if ( D(i,i).ge. 0.0d0 ) then
         D(i,i)=sqrt(D(i,i))
       else
         write(0,*) "Warning: negative eigenvalue (",D(i,i),")in LDL^T decomposition."
         D(i,i)=0.0d0
       end if
    end do
    S=matmul(L,D)
end subroutine cholesky

real*8 function ran2(idum)  !we don't use this anymore
  implicit none
  real *8 x
  integer, intent(inout) :: idum
  integer iseed,i
  integer, allocatable :: seed(:)
    
  if(idum.le.0) then 
     idum=-idum
     call random_seed(size=iseed) 
     allocate(seed(iseed))
     do i=1,iseed  !ugly. once again, this is just a stub. you should use a GOOD prng!
       seed(i)=idum+i
     end do
     call random_seed(put=seed)
  endif
  call random_number(harvest=x)
  ran2=x
  return
end function ran2

real*8 function rang(idum) !we don't use this anymore
  implicit none
  integer, intent(inout) :: idum
  integer iset
  real(8) gset,v1,v2,fac,rsq
  save iset,gset
  data iset/0/
  if (iset.eq.0) then
1    v1=2.d0*ran2(idum)-1.d0
     v2=2.d0*ran2(idum)-1.d0
     rsq=v1**2+v2**2  
     if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1 
     fac=sqrt(-2.d0*log(rsq)/rsq) 
     gset=v1 * fac
     rang=v2 * fac
     iset=1
  else 
     rang=gset
     iset=0 
  endif
  return
end function rang

end module mod_gle
