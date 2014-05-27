! *********************************************************
! * here we define the functions to be called for         *
! * colored-noise thermostatting. incidentally, we also   *
! * write here a couple of functions for white-noise.     *
! *                                                       *
! * code is licensed under GPLv3 [www.gnu.org]            *
! * please consider citing the relevant papers (listed    *
! * below) if you use GLE in your simulations.            *
! *                                                       *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************

module mod_gle
  implicit none
  real*8, allocatable, save ::  gS(:,:), gT(:,:), gp(:,:), ngp(:,:)
  real*8, allocatable :: ran(:)
  real*8, allocatable :: ps(:,:,:)
  real*8 wnt, wns, langham
  integer, save :: ns,irand
contains

  ! initialize white-noise thermostat. 
  ! the init and the propagator here are written with the same phylosophy 
  ! used for the full-fledged colored-noise stuff.
  subroutine wn_init(dt,wopt,irandom)
    use mod_general,only:natom
    use mod_nhc,only:temp
    implicit none
    real*8, intent(in) :: dt,wopt
    real*8 g
    integer irandom
    g=2.*wopt
    wnt=exp(-dt*g)
    wns=sqrt(temp*(1.-wnt*wnt))
    langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
    allocate( ran(natom*3) )  !allocating arrays for random numbers produced by gautrg
!    call gautrg(ran,0,irandom,6)  !inicializace prng
  end subroutine

  ! white-noise propagator. time-step has been set in wn_init
  subroutine wn_step(px,py,pz,m)
    use mod_array_size
    use mod_general
    use mod_nhc, only:temp
    implicit none
    real*8, intent(inout)  :: px(npartmax,nwalkmax)
    real*8, intent(inout)  :: py(npartmax,nwalkmax)
    real*8, intent(inout)  :: pz(npartmax,nwalkmax)
    real*8, intent(in)  :: m(npartmax,nwalkmax)
    integer :: iat,iw,pom
    !DH: nemusime tady prenasobovat hmotnostma?
    iw=1
    pom=1
    call gautrg(ran,natom*3,0,6)
    do iat=1,natom
      px(iat,iw)=wnt*px(iat,iw)+wns*ran(pom)*sqrt(m(iat,iw))  !zmenit generator...jak to
      py(iat,iw)=wnt*py(iat,iw)+wns*ran(pom+1)*sqrt(m(iat,iw))  !zmenit generator...jak to
      pz(iat,iw)=wnt*pz(iat,iw)+wns*ran(pom+2)*sqrt(m(iat,iw))  !zmenit generator...jak to
      pom=pom+3
    end do
  end subroutine
  
  ! initialize gle_init
  subroutine gle_init(dt,irandom)
    use mod_array_size
    use mod_general,only: natom,nwalk
    use mod_nhc,only: temp
    implicit none
    real*8, intent(in)  :: dt
    real *8, allocatable :: gA(:,:), gC(:,:), gr(:)
    integer i, j, k, h, cns, ios, irandom,iw
    

    write(6,*) "# Initialization of GLE thermostat.                           "
    write(6,*) "# Please cite the relevant works among:                       "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phy. Rev. Lett. 102, 020601 (2009)                          "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phy. Rev. Lett. 102, 020601 (2009)                          "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phy. Rev. Lett. 103, 030603 (2009)                          "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phy. Rev. ???????  in press (2009)                          "
    write(6,*) "#                                                             "

    !reads in matrices
    !reads A (in units of the "optimal" frequency of the fitting)
    open(121,file='GLE-A',status='OLD',iostat=ios)
    read(121,*) ns

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
    if(nwalk.gt.1) allocate(ps(natom*3,ns,nwalk))!each bead has to have its additional momenta

    !WARNING: TODO: zkontrolovat, ze to robime spravne
!    call gautrg(ran,0,IRandom,6)  !inicializace prng
    
    if (ios.ne.0) write(0,*) "Error: could not read GLE-A file!"
    write(6,*)'Reading A-matrix. Expecting a.u. units!!!!'
    do i=1,ns+1
       read(121,*) gA(i,:)
    enddo
    close(121)

    ! gamma for a WN langevin will be 1/tau, which would make it optimal for w=1/(2tau) angular frequency. 
    !DHmod: we DONT scale gA (which is expected to be fitted by http://gle4md.berlios.de/ ) 
    !A metrix is expected to be in atomic units of time 
    !gA=gA

    ! reads C (in au!), or init to kT
    ! C matrix from web must be scaled accordingly!!!(ie from eV to au)
    open(121,file='GLE-C',status='OLD',iostat=ios)
    if (ios.ne.0) then            
       write(6,*) "# Using canonical-sampling, Cp=kT"
       gC=0.
       do i=1,ns+1
          gC(i,i)=temp
       enddo
    else    
       write(6,*) "# Reading specialized Cp matrix"
       write(6,*)'#Expecting eV units!!!!'
       read(121,*) cns
       if (cns.ne.ns)then
           write(0,*) " Error: size mismatch between given GLE-A and GLE-C!"
           stop 1
       endif
       do i=1,ns+1
          read(121,*) gC(i,:)
          gC(i,:)=gC(i,:)/autoev
       enddo
    end if
    close(121)

    ! the deterministic part of the propagator is obtained in a second
    call matrix_exp(-dt*gA, ns+1,15,15,gT)
    
    ! the stochastic part is just as easy. we use gA as a temporary array
    gA=gC-matmul(gT,matmul(gC,transpose(gT)))
    call cholesky(gA, gS, ns+1)


    ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
    ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
    ! of C. again, since one would like to initialize correctly the velocities in 
    ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
    ! could then use it to initialize the momentum in the calling code

    !DH: ps or gp rewritten in init.F90 if irest.eq.1
    gA=gC   
    call cholesky(gA, gC, ns+1)
    
    do iw=1,nwalk

    do j=1,natom*3
      call gautrg(ran,ns+1,0,6)
      do i=1,ns+1
        gr(i)=ran(i)
      enddo
      gp(j,:)=matmul(gC,gr)
    end do

    if(nwalk.gt.1)then
     do j=1,natom*3
      do i=1,ns
       ps(j,i,iw)=gp(j,i+1)
      enddo
     enddo
    endif
!nwalk ENDDO
    enddo
    
    deallocate(gA)
    deallocate(gC)
    deallocate(gr)
    langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
  end subroutine gle_init

  ! the GLE propagator. 
  ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
  ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
  ! gp(j,2:ns+1) the current values of additional momenta. 
  ! the matrix multiplies are performed on the whole array at once, and the new momentum
  ! is passed back to the caller, while the new s's are kept stored in gp.
  ! please note that one can avoid the double conversion between mass-scaled and actual
  ! momentum/velocity (as used by the calling code) by scaling with the mass the white-noise
  ! random numbers. 
  subroutine gle_step(px,py,pz,m)
    use mod_array_size
    use mod_general,only:natom,nwalk
    implicit none
    real*8, intent(inout)  :: px(npartmax,nwalkmax)
    real*8, intent(inout)  :: py(npartmax,nwalkmax)
    real*8, intent(inout)  :: pz(npartmax,nwalkmax)
    real*8, intent(in)     :: m(npartmax,nwalkmax)
    integer i, j, iat, iw
    real*8 mfac, totm

    do iw=1,nwalk
!for GLE+PIMD, we store additional momenta in ps 3d matrices
    if(nwalk.gt.1)then
     do i=1,ns
      do iat=1,natom*3
       gp(iat,i+1)=ps(iat,i,iw)
      enddo
     enddo
    endif
      

    do iat=1,natom
     gp(iat,1)=px(iat,iw)         !<-- if m!= 1, here a proper scaling must be performed
     gp(iat+natom,1)=py(iat,iw)   !<-- if m!= 1, here a proper scaling must be performed
     gp(iat+natom*2,1)=pz(iat,iw) !<-- if m!= 1, here a proper scaling must be performed
    enddo

#ifdef USELIBS
    call dgemm('n','t',ndim,ns+1,ns+1,1.0d0,gp,natom*3,gT,ns+1,0.0d0,ngp,ndim)
#else
    ngp=transpose(matmul(gT,transpose(gp)))
#endif

    !now, must compute random part. 
    !first, fill up gp of random n
    do i=1,ns+1
     call gautrg(ran,natom*3,0,6)
     do j=1,natom
        gp(j,i)=ran(j)*sqrt(m(j,iw))     !<-- if m!= 1, alternatively one could perform the scaling here (check also init!)
        gp(j+natom,i)=ran(j+natom)*sqrt(m(j,iw))     
        gp(j+natom*2,i)=ran(j+natom*2)*sqrt(m(j,iw))     
      end do
    end do

#ifdef USELIBS    
    call dgemm('n','t',natom*3,ns+1,ns+1,1.0d0,gp,natom*3,gS,ns+1,1.0d0,ngp,natom*3)
    gp=ngp
#else
    gp=ngp+transpose(matmul(gS,transpose(gp)))
#endif

    do iat=1,natom
     px(iat,iw)=gp(iat,1)           !<-- if m!= 1, here a proper inverse scaling must be performed
     py(iat,iw)=gp(iat+natom,1)     !<-- if m!= 1, here a proper inverse scaling must be performed
     pz(iat,iw)=gp(iat+2*natom,1)   !<-- if m!= 1, here a proper inverse scaling must be performed
    end do

    if(nwalk.gt.1)then
     do i=1,ns
      do iat=1,natom*3
       ps(iat,i,iw)=gp(iat,i+1)
      enddo
     enddo
    endif

    !iw enddo
    enddo
  end subroutine gle_step

  ! matrix exponential by scale & square.
  ! one can diagonalize with lapack, but it's not worth it, as 
  ! we call this routine only once!      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real*8, intent(in)   :: M(n,n)
    real*8, intent(out)   :: EM(n,n)
    
    real *8 :: tc(j+1), tmp(n,n), SM(n,n)
    integer p, i
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1./2.**k)
    EM=0.
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
  
  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
  subroutine cholesky(SST, S, n)
    integer, intent(in)  :: n
    real*8, intent(in)   :: SST(n,n)
    real*8, intent(out)   :: S(n,n)
    real *8 :: L(n,n), D(n,n) 
    integer i,j,k
    S=0.
    L=0.
    D=0.
    do i=1,n
       L(i,i)=1.0
       do j=1,i-1
          L(i,j)=SST(i,j);
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k,k)
          enddo
          if (D(j,j).ne. 0.0) then
            L(i,j)=L(i,j)/D(j,j) 
          else
            write(0,*) "Warning: zero eigenvalue in LDL^T decomposition."
            L(i,j)=0.
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
         D(i,i)=0.0
       end if
    end do
    S=matmul(L,D)
  end subroutine cholesky

real*8 function ran2(idum)
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

real*8 function rang(idum)
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
end module
