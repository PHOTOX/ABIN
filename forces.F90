
! A wrapper routine for getting forces and energies
! Ab-initio programs are called from force_abin routine
subroutine force_clas(fx,fy,fz,x,y,z,energy,chpot)
   use mod_const,    only: DP
   use mod_general,  only: natom, nwalk, istage, inormalmodes, iqmmm, it, &
                           pot, pot_ref, imini, idebug
   use mod_qmmm,     only: force_LJCoul
   use mod_nab,      only: ipbc,wrap,nsnb,force_nab
   use mod_sbc,      only: force_sbc, isbc !,ibag
   use mod_system,   only: conatom
   use mod_nhc,      only: inose
   use mod_transform
   use mod_interfaces, only: force_wrapper
   use mod_plumed,   only: iplumed, plumedfile, force_plumed
   implicit none
   real(DP),intent(inout) :: x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(out)   :: energy
   character(len=*),intent(in) :: chpot
   real(DP)  :: transx(size(x,1),size(x,2))
   real(DP)  :: transy(size(y,1),size(y,2))
   real(DP)  :: transz(size(z,1),size(z,2))
   real(DP)  :: fxab(size(fx,1),size(fx,2))
   real(DP)  :: fyab(size(fy,1),size(fy,2))
   real(DP)  :: fzab(size(fz,1),size(fz,2))
   integer   :: iat,iw
   real (DP) :: eclas

! Initialization
   do iw=1,nwalk
      do iat=1,natom
         fx(iat,iw)=0.0d0
         fy(iat,iw)=0.0d0
         fz(iat,iw)=0.0d0
         fxab(iat,iw)=0.0d0
         fyab(iat,iw)=0.0d0
         fzab(iat,iw)=0.0d0
      enddo
   enddo

   eclas = 0.0d0

! Back stage transformation, Cartesian coordinates are kept in trans
! matrices (even if staging is OFF!)
   if(istage.eq.1)then
      call QtoX(x,y,z,transx,transy,transz)
   else if(inormalmodes.gt.0)then
      if(idebug.gt.0) write(*,*)'Transforming coordinates back to cartesian'
      call UtoX(x,y,z,transx,transy,transz)
   else
      do iat=1,natom
         do iw=1,nwalk
            transx(iat,iw) = x(iat,iw)
            transy(iat,iw) = y(iat,iw)
            transz(iat,iw) = z(iat,iw)
         enddo
      enddo 
   endif

!  wraping molecules back to the box 
   if (chpot.eq.'nab'.and.ipbc.eq.1.and.modulo(it,nsnb).eq.0) call wrap(transx,transy,transz)


   ! LET'S GET FORCES! Ab initio interface is still deeper in force_abin
   call force_wrapper(transx, transy, transz, fxab, fyab, fzab, eclas, chpot, nwalk)

!  Spherical harmonic potential
   if (isbc.eq.1) call force_sbc(transx,transy,transz,fxab,fyab,fzab)
!  if (ibag.eq.1) call force_bag(transx,transy,transz,fxab,fyab,fzab)

!---------QMMM SECTION-----------------
!  ONIOM method (iqmmm=1) is called in force_abin
!  The following are not working at the moment
   if(iqmmm.eq.2) call force_nab(transx, transy, transz, fxab, fyab, fzab, eclas, nwalk)
   if(iqmmm.eq.3) call force_LJCoul(transx, transy, transz, fxab, fyab, fzab, eclas)

!--------PLUMED SECTION---------------
   ! not sure there should be imini, might produce weird results in US
   if(iplumed.eq.1.and.it.gt.imini) call force_plumed(transx,transy,transz,fxab,fyab,fzab,eclas)
!------------------------------------


!  For reference potential and ring-polymer contraction
   if(pot_ref.ne.'none'.and.chpot.eq.pot)then
      ! fxab now holds the full potential,
      ! but we need the difference force on the output
      fx = fxab; fy = fyab; fz = fzab
      fxab = 0.0d0; fyab=0.0d0; fzab=0.0d0
      energy = eclas
      eclas = 0.0d0

      call force_wrapper(transx, transy, transz, fxab, fyab, fzab, eclas, pot_ref, nwalk)
      if(isbc.eq.1)  call force_sbc(transx,transy,transz,fxab,fyab,fzab)
      if(iqmmm.eq.2) call force_nab(transx, transy, transz, fxab, fyab, fzab, eclas, nwalk)
      if(iqmmm.eq.3) call force_LJCoul(transx, transy, transz, fxab, fyab, fzab, eclas)
      if(iplumed.eq.1.and.it.gt.imini) call force_plumed(transx,transy,transz,fxab,fyab,fzab,eclas)

      fxab = fx - fxab
      fyab = fy - fyab
      fzab = fz - fzab
      ! fxab now holds the difference force
      ! we return the difference forces, but full energy
!      eclas = energy - eclas
      eclas = energy
   end if


!--TRANSFORMING FORCES FROM CARTESIAN TO STAGING or NORMAL MODE COORDS--!
!--forces are divided by nwalk inside the FXtoFQ routine!---------------!
   if(istage.eq.1)then

      call FXtoFQ(fxab,fyab,fzab,fx,fy,fz)

   else if(inormalmodes.gt.0)then
      ! for PIGLET

      call XtoU(fxab, fyab, fzab, fx, fy, fz)

   else if (inose.eq.2)then
      ! for PI+GLE 
      do iw=1,nwalk 
         do iat=1,natom
            fx(iat,iw) = fxab(iat,iw)
            fy(iat,iw) = fyab(iat,iw) 
            fz(iat,iw) = fzab(iat,iw) 
         enddo 
      enddo

   else

      ! no transformation, classical MD
      do iw=1,nwalk 
         do iat=1,natom
            fx(iat,iw) = fxab(iat,iw) / nwalk
            fy(iat,iw) = fyab(iat,iw) / nwalk 
            fz(iat,iw) = fzab(iat,iw) / nwalk 
         enddo 
      enddo

   endif
    
   ! Constraining atoms
   ! Warning, this kills energy conservation!
   if(conatom.gt.0)then
      do iw=1,nwalk
         do iat=1,conatom
            fx(iat,iw)=0.0d0
            fy(iat,iw)=0.0d0
            fz(iat,iw)=0.0d0
         enddo
       enddo
   endif

   energy = eclas

!--for PBC we do wrapping of molecules back to the box    
!  We should probably be doing this somewhere else, 
!  this is confusing
   if(chpot.eq.'nab'.and.ipbc.eq.1)then

! Stage transformation,
      if(istage.eq.1)then
         call XtoQ(transx,transy,transz,x,y,z)
      else if(inormalmodes.gt.0)then
         call XtoU(transx,transy,transz,x,y,z)
      else
         x = transx
         y = transy
         z = transz
      endif

   endif

end subroutine force_clas


subroutine force_wrapper(x, y, z, fx, fy, fz,  e_pot, chpot, walkmax)
   use mod_const,    only: DP
   use mod_interfaces, only: force_abin
   use mod_general,  only: natom, ipimd
   use mod_water,    only: watpot
   use mod_qmmm,     only: force_LJCoul
   use mod_nab,      only: force_nab
   use mod_harmon,   only: force_harmon,force_2dho,force_morse,force_doublewell
   use mod_guillot,  only: force_guillot
   use mod_cp2k,     only: force_cp2k
#ifdef MPI
   use mod_terampi,     only: force_tera
   use mod_terampi_sh,  only: force_terash
#endif
   implicit none
   real(DP),intent(in)    ::  x(:,:),y(:,:),z(:,:)
   real(DP),intent(inout) ::  fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(out)   ::  e_pot
   character(len=*),intent(in) :: chpot
   integer, intent(in)    :: walkmax
   real(DP)               :: eclas

   eclas = 0.0d0
! Here we decide which forces we want.
! By default we call external program in force_abin routine
   SELECT CASE (chpot)
     case ("mm")
       call force_LJCoul(x, y, z, fx, fy, fz, eclas)
     case ("mmwater")
       call force_water(x, y, z, fx, fy, fz, eclas, natom, walkmax, watpot)
     case ("harm")
       call force_harmon(x, y, z, fx, fy, fz, eclas)
     case ("2dho")
       call force_2dho(x, y, z, fx, fy, fz, eclas)
     case ("morse")
       call force_morse(x, y, z, fx, fy, fz, eclas)
     case ("guillot")
       call force_guillot(x, y, z, fx, fy, fz, eclas)
     case ("doublewell")
       call force_doublewell(x, y, z, fx, fy, fz, eclas)
     case ("nab")
       call force_nab(x, y, z, fx, fy, fz, eclas, walkmax)
     case ("_cp2k_")
       call force_cp2k(x, y, z, fx, fy, fz, eclas, walkmax)
#ifdef MPI
      case ("_tera_")
         if(ipimd.eq.2.or.ipimd.eq.4)then
            call force_terash(x, y, z, fx, fy, fz, eclas)
         else
            call force_tera(x, y, z, fx, fy, fz, eclas, walkmax)
         end if
#endif
      case DEFAULT
        call force_abin(x, y, z, fx, fy, fz, eclas, chpot, walkmax)
        eclas = eclas / walkmax
   END SELECT

   e_pot = eclas
end subroutine force_wrapper


subroutine force_quantum(fx, fy, fz, x, y, z, amg, energy)
   use mod_const,      only: DP
   use mod_array_size
   use mod_general, only: nwalk, inormalmodes, istage, natom, idebug
   use mod_nhc,   only: temp, inose
   use mod_utils, only: printf
   implicit none
   real(DP), intent(in) :: x(:,:),y(:,:),z(:,:)
   real(DP), intent(in) :: amg(:,:)
   real(DP), intent(inout) :: fx(:,:),fy(:,:),fz(:,:)
   real(DP), intent(out) :: energy
   real(DP) :: ak(size(x,1),size(x,2))
   real(DP) :: equant
   integer  :: iat,iw,i,j,kplus,kminus

   fx = 0.0d0
   fy = 0.0d0
   fz = 0.0d0


!  TODO: ak parametry se nemuseji pocitat pokazde znova
!  Setting the quantum force constants
!  ak is defined is m*P/(beta^2*hbar^2)
   do iw=1,nwalk
      do i=1,natom
         ak(i,iw) = nwalk * amg(i,iw) * TEMP**2
      enddo
   enddo

!  for PI+GLE we have different hamiltonian
   if(inose.eq.2)then
     ak = nwalk * ak
   endif

!   if(inormalmodes.eq.2)then
!      ak = ak / dsqrt(nwalk*1.0d0)
!   end if
   ! Tuckerman normal modes Hamiltonian
   if (inormalmodes.eq.2) then
     ak = NWALK * TEMP**2 * amg / sqrt(nwalk*1.0d0)
   end if


   equant=0.0d0


!  If the staging transformation is not used
   if(istage.eq.0.and.inormalmodes.eq.0)then
      do j=1,natom 
         do i=1,nwalk 
            kplus=i+1
            kminus=i-1
            if(i.eq.1)then
               kminus=nwalk
            endif
            if(i.eq.nwalk)then
               kplus=1
            endif
            fx(j,i)=(x(j,i)-x(j,kplus))
            fx(j,i)=fx(j,i)+(x(j,i)-x(j,kminus))
            fx(j,i)=-fx(j,i)*ak(j,i)
            fy(j,i)=(y(j,i)-y(j,kplus))
            fy(j,i)=fy(j,i)+(y(j,i)-y(j,kminus))
            fy(j,i)=-fy(j,i)*ak(j,i)
            fz(j,i)=(z(j,i)-z(j,kplus))
            fz(j,i)=fz(j,i)+(z(j,i)-z(j,kminus))
            fz(j,i)=-fz(j,i)*ak(j,i)
            equant = equant + 0.5*ak(j,i)*(x(j,i)-x(j,kplus))**2
            equant = equant + 0.5*ak(j,i)*(y(j,i)-y(j,kplus))**2
            equant = equant + 0.5*ak(j,i)*(z(j,i)-z(j,kplus))**2
         enddo
      enddo
   endif

!  If the staging or normal mode transformation is used
   if(istage.eq.1.or.inormalmodes.gt.0)then
      do j=1,natom 
         do i=1,nwalk 
            fx(j,i) = -ak(j,i) * x(j,i)   
            fy(j,i) = -ak(j,i) * y(j,i)       
            fz(j,i) = -ak(j,i) * z(j,i) 
            equant = equant + 0.5d0*ak(j,i)*x(j,i)**2
            equant = equant + 0.5d0*ak(j,i)*y(j,i)**2
            equant = equant + 0.5d0*ak(j,i)*z(j,i)**2
         enddo
      enddo
   endif 

   energy = equant

   if(idebug.eq.1)then
      write(*,*)'EQUANT', equant
   end if
      
end subroutine force_quantum

                                        
! Based on:
! Efficient stochastic thermostatting of path integral molecular dynamics
! J. Chem. Phys. 133, 124104 ?2010?
subroutine propagate_nm(x, y, z, px, py, pz, m, dt)
   use mod_const, only: DP, PI
   use mod_array_size
   use mod_general,   only: nwalk, natom
   use mod_nhc,       only: temp, inose
   implicit none
   real(DP), intent(inout) :: x(:,:), y(:,:), z(:,:)
   real(DP), intent(inout) :: px(:,:), py(:,:), pz(:,:)
   real(DP), intent(in)    :: m(:,:), dt
   real(DP) :: omega(size(x,2)), omega_n, om, omt, c, s
   real(DP) :: x_old(size(x,1),size(x,2))
   real(DP) :: y_old(size(x,1),size(x,2))
   real(DP) :: z_old(size(x,1),size(x,2))
   real(DP) :: px_old(size(px,1),size(px,2))
   real(DP) :: py_old(size(px,1),size(px,2))
   real(DP) :: pz_old(size(px,1),size(px,2))
   integer  :: iat, iw, k

   ! expecting m = am = physical masses

   x_old = x
   y_old = y
   z_old = z
   px_old = px
   py_old = py
   pz_old = pz

   omega_n = TEMP * NWALK

   do iw=2, nwalk
      k = iw -1
      omega(iw) = 2 * omega_n * sin(k * PI / nwalk) 
   enddo

   ! First, propagate centroid
   iw = 1
   do iat = 1, natom
      X(iat,iw) = X(iat,iw) + dt * PX(iat,iw) / M(iat,iw)
      Y(iat,iw) = Y(iat,iw) + dt * PY(iat,iw) / M(iat,iw)
      Z(iat,iw) = Z(iat,iw) + dt * PZ(iat,iw) / M(iat,iw)
   end do

   ! eq 23 from J. Chem. Phys. 133, 124104 2010
   ! exact propagation of a free ring polymer in normal mode coordinates
   do iw = 2, nwalk
      om = omega(iw)
      omt = omega(iw) * dt
      c = cos(omt)
      s = sin(omt)
      do iat = 1, natom
         ! Propagate positions
         x(iat, iw) = x_old(iat, iw) * c &
            + px_old(iat, iw) * s / m(iat,iw) / om
         y(iat, iw) = y_old(iat, iw) * c &
            + py_old(iat, iw) * s / m(iat,iw) / om
         z(iat, iw) = z_old(iat, iw) * c &
            + pz_old(iat, iw) * s / m(iat,iw) / om

         ! propagate momenta
         px(iat, iw) = px_old(iat, iw) * c &
            - x_old(iat, iw) * s * m(iat,iw) * om
         py(iat, iw) = py_old(iat, iw) * c &
            - y_old(iat, iw) * s * m(iat,iw) * om
         pz(iat, iw) = pz_old(iat, iw) * c &
            - z_old(iat, iw) * s * m(iat,iw) * om

      end do

   end do

end subroutine propagate_nm


