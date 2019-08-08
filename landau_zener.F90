! Driver routines for Landau-Zener excited state dynamics
!    Based on PHYSICAL REVIEW A 84, 014701 (2011)
!    Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
!    Andrey K. Belyaev1,2 and Oleg V. Lebedev1
! J. Suchan, J. Chalabala

! No TeraChem functionality YET
! No proper restart functionality (energy array)
! No energy checks

module mod_lz
   use mod_const, only: DP
   use mod_utils, only: abinerror
   use mod_array_size, only: NSTMAX, NTRAJMAX
   implicit none
   !private 
   public :: lz_init, lz_hop, lz_finalize !Routines
   public :: initstate_lz, nstate_lz, deltaE_lz !User defined variables
   public :: en_array_lz, tocalc_lz, istate_lz !Routine variables

   real(DP) :: deltaE_lz
   ! Initial electronic state
   integer :: initstate_lz = 1
   integer :: tocalc_lz(NSTMAX)
   integer :: nstate_lz

   !Module variables
   integer  :: istate_lz
   real(DP),allocatable :: en_array_lz(:,:)
   real(DP),allocatable :: fx_old(:), fy_old(:), fz_old(:)
   save

   CONTAINS

   !Initialization
   subroutine lz_init(x, y, z, vx, vy, vz, dt)
   use mod_general,  ONLY: natom
   real(DP),intent(in)    :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
   real(DP),intent(in)    :: dt
   integer :: ist1

   !Initial state
   istate_lz = initstate_lz

   !Which gradients to compute
   do ist1=1, nstate_lz
       if(ist1.eq.istate_lz)then
          tocalc_lz(ist1) = 1
       else
          tocalc_lz(ist1) = 0
       end if
   end do

   !Allocate energy arrays
   allocate(en_array_lz(nstate_lz, 3)) !last 3 energies (1: current, 2: n-1, 3: n-3)
   allocate(fx_old(natom),fy_old(natom),fz_old(natom))
   en_array_lz=0.0d0

   end subroutine lz_init


   !LZ singlets hop
   subroutine lz_hop(x, y, z, vx, vy, vz, fxc, fyc, fzc, amt, dt, eclas)
   use mod_const,    ONLY: ANG, AUTOFS, PI, AUTOEV
   use mod_files,    ONLY: UPOP,UPROB,UPES
   use mod_general,  ONLY: natom, pot, nwrite, idebug, it, sim_time
   use mod_random,   ONLY: vranf
   use mod_kinetic,  ONLY: ekin_v
   use mod_interfaces,ONLY: force_clas
   real(DP),intent(inout)    :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) :: vx(:,:),vy(:,:),vz(:,:)
   real(DP),intent(inout) :: fxc(:,:), fyc(:,:), fzc(:,:)
   real(DP),intent(in)    :: amt(:,:)
   real(DP),intent(inout) :: eclas
   real(DP),intent(in)    :: dt

   real(DP) :: prob(NSTMAX)
   real(DP) :: ran(10)
   real(DP) :: en_diff(3), second_der
   integer  :: ihop, ihist, ist1, ist2, iat
   integer  :: ist                     ! =istate_lz
   real(DP) :: Ekin, Epot, hop_rdnum, vel_rescale, stepfs
   character(len=20) :: formt

   !TODO: energy drift
   !call check_energy(vx_old, vy_old, vz_old, vx, vy, vz, itrj)
   !call check_energydrift(vx, vy, vz, itrj)

   !Current state
   ist = istate_lz
   !Compute probabilities (all states except current)
   prob=0.0d0
   
   do ist1=1, nstate_lz
     if (ist1.eq.ist) cycle
     do ihist=1, 3
        en_diff(ihist) =  abs(en_array_lz(ist, ihist) - en_array_lz(ist1, ihist))
     end do
     ! Three point minima of adiabatic splitting Zjk
     if ((en_diff(1).gt.en_diff(2)).and.(en_diff(2).lt.en_diff(3)).and.(it.gt.2))then
        second_der = ((en_diff(3) - 2 * en_diff(2) + en_diff(1)) / dt**2)
        prob(ist1) = exp(-PI/2*(sqrt(en_diff(2)**3 / second_der)))
        !write(*,*)"Hop? (",ist,"->",ist1 ,")dE/a.u.", en_diff(2), "Probability:", prob(ist1)
        if(prob(ist1).gt.1)then
           call abinerror('landau_zener_prob')
        end if

     end if
   end do 
   !write(*,*) "diff1",en_diff(1),"diff2",en_diff(2),"diff3",en_diff(3)

   !Hop?
   ihop=0
   call vranf(ran,1,0,6)
   hop_rdnum = ran(1)

   ! Determine, whether we hopped or not
   do ist1=1, nstate_lz
      if (ist1.eq.ist) cycle
      if(hop_rdnum.lt.prob(ist1))then
         ihop = ist1
         exit
      end if
   end do

   if(ihop.ne.0)then
    Ekin = ekin_v(vx,vy,vz)
    Epot = en_diff(1)
    ! Energy conservation criteria
    if ((Epot.lt.Ekin).and.(abs(Epot * AUTOEV).lt.deltaE_lz))then
        !HOP
        write(*,*)"Hop! (",istate_lz,"->",ihop ,")dE/a.u.", en_diff(2), "Probability:", prob(ihop), "Random n:",hop_rdnum
        istate_lz = ihop
        !Get new forces
        do iat=1, natom
            fx_old(iat) = fxc(iat,1)
            fy_old(iat) = fyc(iat,1)
            fz_old(iat) = fzc(iat,1)
        end do                     
        call force_clas(fxc, fyc, fzc, x, y, z, eclas, pot)
        !Adjust velocities from previous step to new forces
        do iat=1, natom
           vx(iat,1) = vx(iat,1) + (dt/2.0d0) * (-fx_old(iat)+fxc(iat,1)) / amt(iat,1)
           vy(iat,1) = vy(iat,1) + (dt/2.0d0) * (-fy_old(iat)+fyc(iat,1)) / amt(iat,1)
           vz(iat,1) = vz(iat,1) + (dt/2.0d0) * (-fz_old(iat)+fzc(iat,1)) / amt(iat,1)
        end do 

        !Simple scaling velocities
        !vel_rescale = sqrt(1-(Epot/Ekin))

        !do iat=1, natom
        !    vx(iat,1) = vx(iat,1) * vel_rescale
        !    vy(iat,1) = vy(iat,1) * vel_rescale
        !    vz(iat,1) = vz(iat,1) * vel_rescale
        !end do

        !Gradient to compute
         do ist1=1, nstate_lz
           if(ist1.eq.istate_lz)then
              tocalc_lz(ist1) = 1
           else                 
              tocalc_lz(ist1) = 0
           end if
         end do                        
     else
        write(*,*)"NO Hop (",istate_lz,"->",ihop ,")dE/a.u.", en_diff(2), "Probability:", prob(ihop), "Random n:",hop_rdnum
     end if

   end if                  

   !Write
   if(modulo(it,nwrite).eq.0)then
     stepfs = sim_time * AUtoFS

     write(formt,'(A7,I3,A7)')'(F15.2,',nstate_lz,'E20.10)'
     write(UPES,fmt=formt)stepfs,(en_array_lz(ist1,1),ist1=1,nstate_lz)
   end if

   end subroutine lz_hop

   subroutine lz_finalize()
   ! Deallocate arrays
       deallocate(en_array_lz, fx_old, fy_old, fz_old)
   end subroutine lz_finalize

end module mod_lz
