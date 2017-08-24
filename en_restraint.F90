module mod_en_restraint
  use mod_const,  only: DP
  use mod_nhc,    only: temp
  use mod_general,  only: natom, nwalk
  implicit none
  public
  real(DP):: restr
  real(DP), allocatable :: fxr(:,:), fyr(:,:), fzr(:,:)
  character(len=200) :: restrain_pot

CONTAINS

! ENERGY RESTRAINT subroutine
! Uses method of Lagrange multipliers to fix excitation energy of molecule - adapts forces based on knowledge
! of ground and excited state gradients

! V1 - using energy gradient as approximation for next step energy change

subroutine en_rest_init(dt)
   use mod_general,  only: natom
   implicit none
   real(DP),intent(in) :: dt
   integer             :: iat

   allocate(fxr(natom,2)) ! two states, not beads.
   allocate(fyr(natom,2))
   allocate(fzr(natom,2))

end subroutine en_rest_init

subroutine energy_restraint(x, y, z, fx, fy, fz, eclas)
   use mod_general,  only: natom, it, dt0 ! dt0 is the time step
   use mod_utils,    only: abinerror
   use mod_const,    only: AMU
   use mod_system,   only: am
   use mod_sh,       only: en_array
   use mod_terampi_sh, only: force_terash
   implicit none
   real(DP),intent(inout)  :: x(:,:),y(:,:),z(:,:),eclas
   real(DP),intent(inout)  :: fx(:,:),fy(:,:),fz(:,:)
   ! DH: I am surprised this works, since natom is not a constant
   real(DP),dimension(natom):: fxgs,fygs,fzgs,fxes,fyes,fzes
   real(DP)                :: eclasexc, eclasground, Egrad, deltaE, lambda, lsum
   integer                 :: ios,iat,iat2,iw
   character(len=30)       :: chforce_ground, chforce_exc

do iw=1,nwalk

   ! DH call the extra forces and potentials here
   if(restrain_pot.eq.'_tera_')then
      call force_terash(x, y, z, fxr, fyr, fzr, eclasground)
      fxgs(iat) = -fxr(iat,1)
      fygs(iat) = -fyr(iat,1)
      fzgs(iat) = -fzr(iat,1)
      fxes(iat) = -fxr(iat,2)
      fyes(iat) = -fyr(iat,2)
      fzes(iat) = -fzr(iat,2)
      eclasground = en_array(1, 1)
      eclasexc = en_array(2, 1)
      ! TODO: need to figure out how to get eclasexc

   else

      ! Should be done as a separate call eventually..
!      call force_abin(x, y, z, fxr, fyr, fzr, energy)

    write(chforce_ground,'(A,I3.3)')'engrad.ground.dat.',iw
    write(chforce_exc,'(A,I3.3)')'engrad.exc.dat.',iw

!-----READING energy of groud state (engrad.ground.dat)
    open(901,file=chforce_ground,status='OLD',iostat=ios,action='read')
    if (ios.ne.0)then
       write(*,*)'Error: could not read engrad.ground.dat.001! Check if you use proper &
       external script for restrained energy dynamics.'
       call abinerror('energy_restraint')
    end if
    read(901,*) eclasground

!----READING gradient of ground state
!----This 
      do iat=1,natom
         read(901,*,IOSTAT=ios)fxgs(iat), fygs(iat), fzgs(iat)
         if(ios.ne.0)then
            write(*,*)'Fatal problem with reading gradients from file engrad.ground.dat '
            call abinerror('energy_restraint')
         endif
       !---Conversion to forces
         fxgs(iat) = -fxgs(iat)
         fygs(iat) = -fygs(iat)
         fzgs(iat) = -fzgs(iat)
      end do
      close(901)

!-----READING energy of excited state (engrad.exc.dat)
    open(901,file=chforce_exc,status='OLD',iostat=ios,action='read')
    if (ios.ne.0)then
       write(*,*)'Error: could not read engrad.exc.dat.001! Check if you use proper external script for restrained energy dynamics.'
       call abinerror('energy_restraint')
    end if
    read(901,*) eclasexc
    !Saving energy difference
    Egrad = eclasexc-eclasground
    write(*,*)'Energy difference ES-GS (',iw,') is',Egrad*2625.5697

!----READING gradient of excited state
      do iat=1,natom
         read(901,*,IOSTAT=ios)fxes(iat), fyes(iat), fzes(iat)
         if(ios.ne.0)then
            write(*,*)'Fatal problem with reading gradients from file engrad.exc.dat '
            call abinerror('energy_restraint')
         endif
       !---Conversion to forces
         fxes(iat) = -fxes(iat)
         fyes(iat) = -fyes(iat)
         fzes(iat) = -fzes(iat)
      end do
      close(901)

   ! restraint_pot endif
   end if

!---Computing lagrange multiplier lambda

!     DH: better use more descriptive variable names..restr does not say much
      deltaE=(Egrad-restr)

!write(*,*)'fx old',fx

      lsum=0        
      do iat2=1, natom
         lsum = lsum + 1/(am(iat2)) * (fxes(iat2) - fxgs(iat2))**2
         lsum = lsum + 1/(am(iat2)) * (fyes(iat2) - fygs(iat2))**2
         lsum = lsum + 1/(am(iat2)) * (fzes(iat2) - fzgs(iat2))**2
      end do
!      lambda =  deltaE / (lsum * dt0 * dt0)
! DH: the following line was uncommented, assuming this was just a debug
!      lambda = 0
      !---Applying new forces
      do iat=1,natom
         fx(iat,iw) = fx(iat,iw) + lambda * (fxes(iat) - fxgs(iat))
         fy(iat,iw) = fy(iat,iw) + lambda * (fyes(iat) - fygs(iat))
         fz(iat,iw) = fz(iat,iw) + lambda * (fzes(iat) - fzgs(iat))
      end do
        

!---
!write(*,*)'fx new',fx
!write(*,*)'lsum is',lsum
write(*,*)'deltaE',deltaE
write(*,*)'Lambda multiplier:',lambda
end do
write(*,*)'--using 3 state (MD,GS,ES) version--'
write(*,*)''


end subroutine energy_restraint

end module

