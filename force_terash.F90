module mod_terampi_sh
! ----------------------------------------------------------------
! Interface for TeraChem based Surface Hopping.
! Based on the FMS interface from FMS90
! (files TerachemModule.f90 and )
!
! Original Author: Unknown 
!
! Date: December 2015
! ----------------------------------------------------------------
   use mod_const, only: DP
#ifdef MPI
   use mod_terampi, only: newcomm, qmcharges, qmcoords, dxyz_all

   implicit none
   private
   public :: force_terash, init_terash, send_terash
   real(DP), allocatable :: CIvecs(:,:), MO(:,:), CMO(:,:), NAC(:)
   real(DP), allocatable :: CIvecs_old(:,:), MO_old(:,:), CMO_old(:,:)
   integer :: civec, nbf 
   save
contains


subroutine force_terash(x, y, z, fx, fy, fz, eclas)
   use mod_const, only: DP
   include 'mpif.h'
   real(DP),intent(in)      ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   ::  eclas

   call send_terash(x, y, z, fx, fy, fz)

   call receive_terash(fx, fy, fz, eclas)


end subroutine force_terash

subroutine receive_terash(fx, fy, fz, eclas)
   use mod_const, only: DP, ANG
   use mod_array_size, only: NSTMAX
   use mod_general, only: idebug, natom
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror
   use mod_sh, only: check_CIVector, en_array, nstate, istate, nacx, nacy, nacz, tocalc
   include 'mpif.h'
   real(DP),intent(inout) :: fx(:,:), fy(:,:), fz(:,:)
   real(DP),intent(inout) :: eclas
   real(DP) :: dip(NSTMAX*3), tdip((NSTMAX-1)*3) ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iat,iw, ist1, ist2, itrj, ipom, i
   integer  :: bufints(20)

   itrj = 1
   iw = 1

!  Receive energies from TC
   if (idebug>0) write(*, '(a)') 'Receiving energies from TC.'
!  DH WARNING this will only work if itrj = 1
   call MPI_Recv( en_array, nstate, MPI_DOUBLE_PRECISION, &
           MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

   eclas = en_array(istate(itrj), itrj)

   if (idebug>0) write(*, '(a)') 'Receiving transition dipoles from TC.'
   call MPI_Recv( TDip, (nstate-1)*3,  &
        MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
!   do i=1, nstate-1
!      T_FMS%ElecStruc%TransDipole(i+1,:)=TDip(3*(i-1)+1:3*(i-1)+3)
!   end do

!    Receive dipole moment from TC
   if (idebug>0) write(*, '(a)') 'Receiving dipole moments from TC.'
   call MPI_Recv( Dip,nstate*3, &
          MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
!  do i=1, nstate
!     T_FMS%ElecStruc%Dipole(i,1:3)=Dip(3*(i-1)+1:3*(i-1)+3)
!  end do

!   Receive partial charges from TC
   if (idebug>0) write(*, '(a)') 'Receiving charges from TC.'
   call MPI_Recv( qmcharges, natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

!  Receive gradients from TC
   if (idebug>0) write(*, '(a)') 'Receiving gradients from TC.'
   call MPI_Recv( dxyz_all, 3*natom, MPI_DOUBLE_PRECISION, &
        MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
   if (idebug>0) write(*, *)dxyz_all(1,1)
  

   do iat=1,natom
      fx(iat,iw)=-dxyz_all(1,iat)
      fy(iat,iw)=-dxyz_all(2,iat)
      fz(iat,iw)=-dxyz_all(3,iat)
   end do

!  Receive nbf, CI len from TC
!  assume here that they are the same as initially   
   if (idebug>0) write(*, '(a)') 'Receiving nbf and civec from TC.'
   call MPI_Recv( bufints, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
     MPI_ANY_TAG, newcomm, status, ierr)

!  Receive MOs from TC
   if (idebug>0) write(*, '(a)') 'Receiving MOs from TC.'
   call MPI_Recv( MO, nbf*nbf, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
     MPI_ANY_TAG, newcomm, status, ierr)

!   T_FMS%ElecStruc%OldOrbitals=MO

   if (idebug>0) write(*, '(a)') 'Receiving canonical MOs from TC.'
!      if (tc_nml%cis == 'yes') then
!        call MPI_Recv( CIvecs_tmp, tc_nml%cisnumstates*(ci_length), 
!     $      MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm,
!     $      status, ierr)
!      else if (tc_nml%casscf == 'yes') then
   call MPI_Recv( CMO, nbf*nbf, MPI_DOUBLE_PRECISION, &
          MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

   if (idebug>0) write(*, '(a)') 'Receiving CI vectors from TC.'
   call MPI_Recv( CIvecs, nstate*civec,  &
           MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
!      end if
        
!      if (tc_nml%cis == 'yes') then
!        CIvecs(1:ci_length, 2:tc_nml%cisnumstates+1) = CIvecs_tmp
!        CIvecs(ci_length+1, 1) = 1.D0
!      end if

   i = Check_CIVector(CIvecs, CIvecs_old, civec, nstate)

   CIVecs_old = Civecs

   do ist1=1, nstate-1
      do ist2=ist1+1, nstate
      ! At this point, TeraChem sends all couplings
!         if (tocalc(ist1, ist2).eq.0) cycle
         if (idebug>0) write(*, '(a,i3,i3)') 'Receiving NA couplings between states.'&
         ,ist1, ist2

         call MPI_Recv( NAC, 3*natom, MPI_DOUBLE_PRECISION, &
              MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

         if (idebug>0) write(*, *)(NAC(i),i=1,3*natom)

         ipom = 1
         do iat=1,natom
            nacx(iat, itrj, ist1, ist2) = NAC(ipom)
            nacy(iat, itrj, ist1, ist2) = NAC(ipom+1)
            nacz(iat, itrj, ist1, ist2) = NAC(ipom+2)
            nacx(iat,itrj,ist2,ist1) = -nacx(iat,itrj,ist1,ist2)
            nacy(iat,itrj,ist2,ist1) = -nacy(iat,itrj,ist1,ist2)
            nacz(iat,itrj,ist2,ist1) = -nacz(iat,itrj,ist1,ist2)
            ipom=ipom+3
         end do
      end do
   end do

end subroutine receive_terash

subroutine send_terash(x, y, z, vx, vy, vz)
   use mod_array_size, only: NSTMAX
   use mod_const, only: DP, ANG, AUTOFS
   use mod_general, only: natom, idebug, it
   use mod_system, only: names
   use mod_qmmm,  only: natqm
   use mod_utils, only: abinerror
   use mod_sh,    only: istate, tocalc, nstate 
   include 'mpif.h'
   real(DP),intent(in)      ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   ::  vx(:,:),vy(:,:),vz(:,:)
   real(DP)                 ::  bufdoubles(100), vels(3, 1000) !DH initial hack
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iw, iat, itrj, FMSInit,i, ist1, ist2
   integer  :: bufints(NSTMAX*(NSTMAX-1)/2)

   itrj = 1
   iw = 1
   do iat=1, natqm
      qmcoords(1,iat) = x(iat,iw) / ANG
      qmcoords(2,iat) = y(iat,iw) / ANG
      qmcoords(3,iat) = z(iat,iw) / ANG
      vels(1,iat) = vx(iat,iw) 
      vels(2,iat) = vy(iat,iw) 
      vels(3,iat) = vz(iat,iw) 
   end do

   ! Send ESinit
   FMSinit=0
   bufints(1)=FMSinit
   bufints(2)=natom
   ! Maybe TODO: send doCoup only if tocalc contains at least one non-zero element
   bufints(3)=1      ! doCoup
   bufints(4)=itrj   ! T_FMS%TrajID
   bufints(5)=0 ! T_FMS%CentID(1)
   bufints(6)=0 ! T_FMS%CentID(2)
   bufints(7)=istate(itrj)-1  ! T_FMS%StateID ! currently not used in fms.cpp
   bufints(8)=civec  ! ci_length
   bufints(9)=istate(itrj)-1  ! iCalcState-1 ! TC Target State
   bufints(10)=istate(itrj)-1 ! jCalcState-1
   bufints(11)=0 ! first_call
   bufints(12)=0 ! glirestTC

   call MPI_Send(bufints, 12, MPI_INTEGER, 0, 2, newcomm, ierr )

   ! temporary hack
   bufdoubles(1) = it * AUtoFS
   ! Send Time 
   call MPI_Send(bufdoubles, 1, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )

!  Send atom types
   call MPI_Send(names, 2*natom, MPI_CHARACTER, 0, 2, newcomm, ierr )

!  Send coordinates
   call MPI_Send(qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
   write(*, '(a)') 'Sent initial coordinates to TeraChem.'

!  Send previous MOs
   if(idebug.gt.0) write(*,*)'Sending previous orbitals.', nbf*nbf
   call MPI_Send(MO, nbf*nbf, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

!  Send previous CI vecs and MOs
!   if (tc_nml%cis == 'yes') then
!        call MPI_Send( CIvecs_tmp, ci_length*tc_nml%cisnumstates, 
!     $      MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
!      else if (tc_nml%casscf == 'yes') then
   if(idebug.gt.0) write(*,*)'Sending canonical orbitals.', nbf*nbf
   call MPI_SSend(CMO, nbf*nbf, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
   if(idebug.gt.0) write(*,*)'Sending CI vector of size ', civec*nstate
   call MPI_SSend(CIvecs, civec*nstate, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)
!  end if

   if(idebug.gt.0) write(*,*)'Sending velocities'
!  Send velocity vectors
!  Not really needed for CASSCF, so send 0 instead
   vels = 0.0d0
   call MPI_Send(vels, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
! Not sure why second set of velocities is needed. (maybe for centroids in FMS)
   call MPI_SSend(vels , 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )

   i=1
   do ist1=1,nstate-1
      do ist2=ist1+1,nstate
         bufints(i)=tocalc(ist1, ist2)
         i=i+1
      end do
   end do

   if(idebug.gt.0)then
      write(*,*)'Sending tocalc array.'
      write(*,*)(bufints(i),i=1,nstate*(nstate-1)/2)
   end if
   call MPI_Send(bufints, nstate*(nstate-1)/2, MPI_INTEGER, 0, 2, newcomm, ierr )

   if(idebug.gt.0) write(*,*)'Succesfully sent all data to TeraChem-FMS'

end subroutine send_terash

subroutine init_terash(x, y, z)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, nwalk, DP, natom
   use mod_system, only: names
   use mod_qmmm, only: natqm
   use mod_sh, only: nstate
   include 'mpif.h'
   real(DP),intent(in)  ::  x(:,:), y(:,:), z(:,:)
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, FMSinit, iat, iw
   integer  :: bufints(20)

   iw = 1
   do iat=1, natqm
      qmcoords(1,iat) = x(iat,iw) / ANG
      qmcoords(2,iat) = y(iat,iw) / ANG
      qmcoords(3,iat) = z(iat,iw) / ANG
   end do

!  Send 0 ints for ESinit
   bufints(1:12) = 0

   FMSinit = 1
   bufints(1) = FMSinit
   bufints(2) = natom
   call MPI_Send( bufints, 12, MPI_INTEGER, 0, 2, newcomm, ierr )

!  send initial time step (so far 0, no restart yet)
   bufints(1) = 0.0d0 ! it
   call MPI_Send( bufints, 1, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )

!  Send atom types
   call MPI_Send( names, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )

!  Send coordinates
   call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
   write(*, '(a)') 'Sent initial coordinates to TeraChem.'

!  Receive nbf and CI length
   call MPI_Recv( bufints, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
     MPI_ANY_TAG, newcomm, status, ierr)

!JWS: Begin Edit
!   if (tc_nml%cis == 'yes') then 
!     esLCiVec = bufints(1)+1        ! Include GS CSF
!   else if (tc_nml%casscf == 'yes') then
!   esLCiVec = bufints(1)
!   end if
   civec = bufints(1)

   nbf = bufints(2)
   write(*,*)'CI vector and nuber of AOs:', CiVec, nbf

   allocate(MO(nbf, nbf))
   allocate(MO_old(nbf, nbf))
   allocate(CMO(nbf, nbf))
   allocate(CMO_old(nbf, nbf))
   allocate(CiVecs(civec,nstate))
   allocate(CiVecs_old(civec,nstate))
   allocate(NAC(natom*3))

end subroutine init_terash

#endif

end module mod_terampi_sh

