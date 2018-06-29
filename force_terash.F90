module mod_terampi_sh
!----------------------------------------------------------------
! Interface for TeraChem based Surface Hopping.
! Based on the FMS interface from FMS90 (TerachemModule.f90)
!
! Original Authors: Basile Curchod, J. Snyder and Ed Hohenstein
!----------------------------------------------------------------
   use mod_const, only: DP
   use mod_terampi, only: chsys_sleep
   implicit none
   private
#ifdef MPI
   public :: init_terash, send_terash
#endif
   public :: force_terash, finalize_terash
   public :: write_wfn, read_wfn, move_new2old_terash, move_old2new_terash
   real(DP), allocatable :: CIvecs(:,:), MO(:,:), blob(:), NAC(:)
   real(DP), allocatable :: CIvecs_old(:,:), MO_old(:,:), blob_old(:)
   real(DP), allocatable :: SMatrix(:)
   integer :: civec, nbf, blobsize, oldWfn = 0 
   save

CONTAINS

subroutine force_terash(x, y, z, fx, fy, fz, eclas)
   use mod_const, only: DP
   use mod_terampi, only: newcomms
   real(DP),intent(in)      ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)   ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(inout)   ::  eclas

   ! for SH, we use only one TC server...
   ! might be changes if we ever implement more elaborate SH schemes
#ifdef MPI
   call send_terash(x, y, z, fx, fy, fz, newcomms(1))

   call receive_terash(fx, fy, fz, eclas, newcomms(1))
#else
   write(*,*) "FATAL ERROR: ABIN not compiled with MPI, cannot connect to TeraChem"
   stop 1
#endif

end subroutine force_terash


#ifdef MPI

subroutine receive_terash(fx, fy, fz, eclas, newcomm)
   use mod_const, only: DP, ANG
   use mod_array_size, only: NSTMAX
   use mod_general, only: idebug, natom, en_restraint
   use mod_qmmm, only: natqm
   use mod_utils, only: abinerror
   use mod_io, only: print_charges, print_dipoles, print_transdipoles
   use mod_sh_integ, only: nstate
   use mod_sh, only: check_CIVector, en_array, istate, nacx, nacy, nacz, tocalc
   include 'mpif.h'
   real(DP),intent(inout) :: fx(:,:), fy(:,:), fz(:,:)
   real(DP),intent(inout) :: eclas
   integer, intent(in)    :: newcomm
   real(DP) :: dip(NSTMAX*3), tdip((NSTMAX-1)*3) ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
   real(DP) :: qmcharges( size(fx,1) )
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iat,iw, ist1, ist2, itrj, ipom, i
   integer  :: bufints(20)
   logical  :: ltest

   itrj = 1
   iw = 1

!  Receive energies from TC
   if (idebug>0) write(*, '(a)') 'Receiving energies from TC.'

   ! DH reduce cpu usage comming from MPI_Recv() via system call to 'sleep'.
   ! Not elegant, but MPICH apparently does not currently provide better solution.
   ! Based according to an answer here:
   ! http://stackoverflow.com/questions/14560714/probe-seems-to-consume-the-cpu

   ltest = .false.
   do while(.not.ltest)
      call MPI_IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG,newcomm,ltest, status, ierr)
      call system(chsys_sleep)
   end do

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
   ! TODO: these things should be printed in analysis.F90
   ! TODO: move charges and dipoles to array module and make them universal
   ! TODO: move TDIP to surface hopping module
   ! allow reading this stuff from other programs as well
   call print_transdipoles(TDip, istate(itrj), nstate-1 )

!  Receive dipole moment from TC
   if (idebug>0) write(*, '(a)') 'Receiving dipole moments from TC.'
   call MPI_Recv( Dip,nstate*3, &
          MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
!  do i=1, nstate
!     T_FMS%ElecStruc%Dipole(i,1:3)=Dip(3*(i-1)+1:3*(i-1)+3)
!  end do
   call print_dipoles(Dip, iw, nstate )

!  Receive partial charges from TC
   if (idebug>0) write(*, '(a)') 'Receiving atomic charges from TC.'
   call MPI_Recv( qmcharges, natqm, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

   call print_charges(qmcharges, istate(itrj) )

!  Receive MOs from TC
   if (idebug>0) write(*, '(a)') 'Receiving MOs from TC.'
   call MPI_Recv( MO, nbf*nbf, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
     MPI_ANY_TAG, newcomm, status, ierr)

!   T_FMS%ElecStruc%OldOrbitals=MO

   if (idebug>0) write(*, '(a)') 'Receiving CI vectors from TC.'
   call MPI_Recv( CIvecs, nstate*civec,  &
           MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)
        
   if (idebug>0) write(*,*) "Receiving wavefunction overlap via MPI."
   call MPI_Recv(SMatrix, nstate*nstate, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr);

   ! Should change the following according to what is done in TeraChem
   i = Check_CIVector(CIvecs, CIvecs_old, civec, nstate)

   CIVecs_old = Civecs

   if (idebug>0) write(*, '(a)') 'Receiving blob from TC.'
   call MPI_Recv( blob, blobsize,  &
           MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

   if (idebug>0) write(*, '(a)') 'Receiving gradients and NACME from TC.'
   do ist1=1, nstate
      do ist2=ist1, nstate

         if (idebug>0) write(*, '(a,i3,i3)') 'Receiving derivatives between states.'&
         ,ist1, ist2


         call MPI_Recv( NAC, 3*natom, MPI_DOUBLE_PRECISION, &
              MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr)

         if (idebug>0) write(*, *)(NAC(i),i=1,3*natom)


         ipom = 1
         if (ist1.eq.ist2.and.istate(itrj).eq.ist1)then
         ! GRADIENTS
            do iat=1,natom
               fx(iat,iw)=-NAC(ipom)
               fy(iat,iw)=-NAC(ipom+1)
               fz(iat,iw)=-NAC(ipom+2)
               ipom = ipom + 3
            end do
         else if (ist1.eq.ist2)then
            ! DH2Jirka: here we will read excited state forces..
            ! perhaps we can use the iw index for the excited state force e.g.
            ! (this assumes, that the initial state is ground state)
            if (en_restraint.ge.1)then
               if(ist1.gt.2)then
                  write(*,*)'ERROR: Energy restraint not implemented for more than 2 states!'
                  call abinerror('receive_terash')
               end if
               do iat=1,natom
                  fx(iat,2)=-NAC(ipom)
                  fy(iat,2)=-NAC(ipom+1)
                  fz(iat,2)=-NAC(ipom+2)
                  ipom = ipom + 3
               end do
            else
               cycle
            end if
         else
         ! NACME
            do iat=1,natom
               nacx(iat, itrj, ist1, ist2) = NAC(ipom)
               nacy(iat, itrj, ist1, ist2) = NAC(ipom+1)
               nacz(iat, itrj, ist1, ist2) = NAC(ipom+2)
               nacx(iat,itrj,ist2,ist1) = -nacx(iat,itrj,ist1,ist2)
               nacy(iat,itrj,ist2,ist1) = -nacy(iat,itrj,ist1,ist2)
               nacz(iat,itrj,ist2,ist1) = -nacz(iat,itrj,ist1,ist2)
               ipom = ipom + 3
            end do
         end if

      end do
   end do

   oldWfn = 1

end subroutine receive_terash


subroutine send_terash(x, y, z, vx, vy, vz, newcomm)
   use mod_array_size, only: NSTMAX
   use mod_const, only: DP, ANG, AUTOFS
   use mod_general, only: natom, idebug, sim_time, en_restraint
   use mod_system, only: names
   use mod_qmmm,  only: natqm
   use mod_utils, only: abinerror
   use mod_sh_integ, only: nstate
   use mod_sh,    only: istate, tocalc, ignore_state 
   include 'mpif.h'
   real(DP),intent(in)     :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout)  :: vx(:,:),vy(:,:),vz(:,:)
   integer, intent(in)     :: newcomm
   real(DP)                ::  bufdoubles(100)
   real(DP) :: qmcoords(3, size(x,1)), vels(3,size(vx,1) )
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iw, iat, itrj, i, ist1, ist2
   integer  :: bufints(NSTMAX*(NSTMAX-1)/2+NSTMAX)
   integer, parameter :: FMSInit = 0

   itrj = 1
   iw = 1
   do iat=1, natqm
      qmcoords(1,iat) = x(iat,iw)
      qmcoords(2,iat) = y(iat,iw)
      qmcoords(3,iat) = z(iat,iw)
      vels(1,iat) = vx(iat,iw)
      vels(2,iat) = vy(iat,iw)
      vels(3,iat) = vz(iat,iw)
   end do

   ! Send ESinit
   bufints(1)=FMSinit
   bufints(2)=natom
   bufints(3)=1   ! doCoup
   bufints(4)=0   ! TrajID=0 for SH
   bufints(5)=0   ! T_FMS%CentID(1)
   bufints(6)=0   ! T_FMS%CentID(2)
   bufints(7)=istate(itrj)-1  ! T_FMS%StateID ! currently not used in fms.cpp
   bufints(8)=oldWfn  ! does ABIN have info about WF?
   bufints(9)=istate(itrj)-1  ! iCalcState-1 ! TC Target State
   bufints(10)=istate(itrj)-1 ! jCalcState-1
   bufints(11)=0 ! first_call, not used
   bufints(12)=0 ! FMSRestart, not used

   call MPI_Send(bufints, 12, MPI_INTEGER, 0, 2, newcomm, ierr )

   ! The following bit is not in FMS code
   ! let ABIN decide which derivatives should TC compute
   i=1
   if (ignore_state.gt.0)then
      do ist1=1, nstate
         tocalc(ist1,ignore_state)  = 0
         tocalc(ignore_state, ist1) = 0
      end do
   end if
   do ist1=1,nstate
      do ist2=ist1,nstate
         if(ist1.eq.ist2.and.ist1.eq.istate(itrj))then
            bufints(i) = 1
         else if (ist1.eq.ist2)then
            ! DH hack for jirka
            ! this will work only if we compute only S0 and S1 states
            if(en_restraint.ge.1)then
               bufints(i) = 1
            else
               bufints(i) = 0
            end if
         else
            bufints(i) = tocalc(ist1, ist2)
         end if
         i=i+1
      end do
   end do

   if(idebug.gt.0)then
      write(*,*)'Sending derivative matrix logic.'
      write(*,*)(bufints(i),i=1,nstate*(nstate-1)/2+nstate)
   end if
   call MPI_SSend(bufints, nstate*(nstate-1)/2+nstate, MPI_INTEGER, 0, 2, newcomm, ierr )

   ! temporary hack
   bufdoubles(1) = sim_time ! * AUtoFS !* dt
   ! Send Time 
   call MPI_Send(bufdoubles, 1, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )

!  Send coordinates
   call MPI_Send(qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
   if(idebug.gt.0) write(*, '(a)') 'Sent coordinates to TeraChem.'

!  Send previous diabatic MOs
   if(idebug.gt.0) write(*,*)'Sending previous orbitals.', nbf*nbf
   call MPI_Send(MO, nbf*nbf, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

!  Send previous CI vecs
   if(idebug.gt.0) write(*,*)'Sending CI vector of size ', civec*nstate
   call MPI_Send(CIvecs, civec*nstate, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

   if(idebug.gt.0) write(*,*)'Sending blob.'
   call MPI_Send(blob, blobsize, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr)

   if(idebug.gt.0) write(*,*)'Sending velocities'
!  Only needed for numerical NACME, so send 0 instead for now
   vels = 0.0d0
   call MPI_Send(vels, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
!  Imaginary velocities for FMS, not needed here, sending zeros...
   call MPI_SSend(vels , 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )


   if(idebug.gt.0) write(*,*)'Succesfully sent all data to TeraChem-FMS'

end subroutine send_terash


subroutine init_terash(x, y, z)
   use mod_const, only: DP, ANG
   use mod_general, only: idebug, nwalk, DP, natom
   use mod_system, only: names
   use mod_qmmm, only: natqm
   use mod_sh_integ, only: nstate
   use mod_terampi, only: newcomms, natmm_tera
   include 'mpif.h'
   real(DP),intent(in)  ::  x(:,:), y(:,:), z(:,:)
   real(DP) :: qmcoords(3, size(x,1))
   integer  :: status(MPI_STATUS_SIZE)
   integer  :: ierr, iat, iw, newcomm
   integer, parameter :: FMSinit = 1
   integer  :: bufints(3)

   ! use only one TC server !
   newcomm = newcomms(1)

   iw = 1
   do iat=1, natqm
      qmcoords(1,iat) = x(iat,iw)
      qmcoords(2,iat) = y(iat,iw)
      qmcoords(3,iat) = z(iat,iw)
   end do

   bufints(1) = FMSinit
   bufints(2) = natom-natmm_tera
   bufints(3) = natmm_tera
   call MPI_SSend( bufints, 3, MPI_INTEGER, 0, 2, newcomm, ierr )
   if (idebug.gt.0) write(*, '(a)') 'Sent initial FMSinit.'

!  Send atom types
   call MPI_Send( names, 2*natqm, MPI_CHARACTER, 0, 2, newcomm, ierr )
   if (idebug.gt.0) write(*, '(a)') 'Sent initial atom types.'

!  Send coordinates
   call MPI_Send( qmcoords, 3*natom, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr )
   if (idebug.gt.0) write(*, '(a)') 'Sent initial coordinates to TeraChem.'

!-- START RECEIVING INFO FROM TeraChem ------!

!  Receive nbf,CI length and blob size
   call MPI_Recv( bufints, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
     MPI_ANY_TAG, newcomm, status, ierr)

   civec = bufints(1)
   nbf = bufints(2)
   blobsize = bufints(3)

   write(*,*)'size of CI vector, number of AOs, blob size:', CiVec, nbf, blobsize

   allocate(MO(nbf, nbf))
   allocate(MO_old(nbf, nbf))
   allocate(CiVecs(civec,nstate))
   allocate(CiVecs_old(civec,nstate))
   allocate(NAC(natom*3))
   allocate(blob(blobsize))
   allocate(blob_old(blobsize))
   allocate(SMatrix(nstate*nstate))
   blob = 0.0d0
   blob_old = 0.0d0

end subroutine init_terash

#endif

subroutine finalize_terash()
   if( allocated(MO) )then
      deallocate(MO, MO_old)
      deallocate(blob, blob_old)
      deallocate(CiVecs, CiVecs_old)
      deallocate(NAC, SMatrix)
   end if
end subroutine finalize_terash

subroutine write_wfn()
   use mod_files, only: UWFN
   use mod_general, only: it, sim_time, iremd, my_rank, narchive
   use mod_sh_integ, only: nstate
   use mod_utils, only: archive_file
   character(len=200)    :: chout, chsystem
   logical  :: file_exists

   if(iremd.eq.1)then
      write(chout, '(A,I2.2)')'wfn.bin.', my_rank
   else
      chout='wfn.bin'
   end if

   INQUIRE(FILE=chout, EXIST=file_exists)
   chsystem='mv '//trim(chout)//'  '//trim(chout)//'.old'
   if(file_exists) call system(chsystem)

   open(UWFN, file=chout, action='WRITE',status="NEW",access="Sequential",form="UNFORMATTED")

   write(UWFN)it, sim_time
   write(UWFN)nbf
   write(UWFN)MO
   write(UWFN)civec, nstate
   write(UWFN)Civecs
   write(UWFN)blobsize
   write(UWFN)blob

   close(UWFN)

   if(modulo(it,narchive).eq.0)  call archive_file('wfn.bin',it)

end subroutine write_wfn


subroutine read_wfn()
   use mod_files, only: UWFN
   use mod_general, only: iremd, my_rank, iknow, it
   use mod_chars, only: chknow
   use mod_utils, only: abinerror, archive_file
   use mod_sh_integ,    only: nstate
   character(len=200)    :: chout, chsystem
   logical  :: file_exists
   integer  :: temp, temp2, time_step
   real(DP) :: stime

   if(iremd.eq.1)then
      write(chout, '(A,I2.2)')'wfn.bin.', my_rank
   else
      chout='wfn.bin'
   end if

   INQUIRE(FILE=chout, EXIST=file_exists)
   if(.not.file_exists)then
      write(*,*)'ERROR: wavefunction restart file does not exist! ', chout
      write(*,*)chknow
      if(iknow.ne.1) call abinerror('read_wfn')
      RETURN
   end if

   open(UWFN, file=chout, action='READ',status="OLD",access="Sequential",form="UNFORMATTED")

   read(UWFN)time_step, stime
   read(UWFN)temp
   if(temp.ne.nbf)then
      write(*,*)'ERROR: Number of MOs in restart file is inconsistent!'
      GO TO 10
   end if
   read(UWFN)MO
   read(UWFN)temp, temp2
   if(temp.ne.civec.or.temp2.ne.nstate)then
      write(*,*)'ERROR: Number and/or size of the CI vectors in restart file is inconsistent!'
      GO TO 10
   end if
   read(UWFN)CIVecs
   read(UWFN)temp
   if(temp.ne.blobsize)then
      write(*,*)'ERROR: Size of blob in restart file is inconsistent!'
      GO TO 10
   end if
   read(UWFN)blob

   close(UWFN)

   oldWFN = 1
   call archive_file('wfn.bin',it)

   RETURN

10 close(UWFN)
   write(*,*)'If you want to proceed, delete file "wfn.bin" and then...'
   write(*,*)chknow
   call abinerror('read_wfn')

end subroutine read_wfn


subroutine move_new2old_terash
   MO_old = MO
   CIVecs_old = CIVecs
   blob_old = blob
end subroutine move_new2old_terash


subroutine move_old2new_terash
   MO = MO_old
   CIVecs = CIVecs_old
   blob = blob_old
end subroutine move_old2new_terash


end module mod_terampi_sh

