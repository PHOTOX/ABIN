module mod_plumed
  use mod_const, only: DP
  use mod_nhc,  only:temp
  implicit none
  public
  integer :: plumed
  character(len=40) :: plumedfile


CONTAINS

#ifdef PLUM
subroutine plumed_init(dt)
 use mod_general, only: natom
 use mod_const,   only: ANG,AUTOFS
 implicit none
 real, parameter:: kB=8.3144621d-3
 real(DP) :: dt, plumed_kbt, plumed_energyUnits, plumed_lengthUnits, plumed_timeUnits
 character(256) :: plumedoutfile

 !variables computation
 plumed_kbt=kB*temp
 plumed_energyUnits=2625.5697  !conversion Ha -> kJ/mol
 plumed_lengthUnits=0.1/ANG    !conversion Angstrom*ANG -> nm
 plumed_timeUnits=AUTOFS*0.001 ! au time -> ps
 plumedoutfile = 'plumed.out'

 !plumed init
 call plumed_f_gcreate()
 call plumed_f_gcmd("setRealPrecision"//char(0),8)
 call plumed_f_gcmd("setMDEnergyUnits"//char(0),plumed_energyUnits)
 call plumed_f_gcmd("setMDLengthUnits"//char(0),plumed_lengthUnits)
 call plumed_f_gcmd("setMDTimeUnits"//char(0),plumed_timeUnits)
 call plumed_f_gcmd("setPlumedDat"//char(0),trim(adjustl(plumedfile))//char(0))
 call plumed_f_gcmd("setNatoms"//char(0),natom)
 call plumed_f_gcmd("setMDEngine"//char(0),"abin")
 call plumed_f_gcmd("setTimestep"//char(0),dt)
 call plumed_f_gcmd("setKbT"//char(0),plumed_kbt)
 call plumed_f_gcmd("setLogFile"//char(0),trim(adjustl(plumedoutfile))//char(0)) 
 
 call plumed_f_gcmd("init"//char(0),0);
end subroutine plumed_init
#endif

#ifdef PLUM
subroutine plumed_stepinfo(it,amt)
 use mod_general, only: natom,nwalk
 use mod_const, only:AMU
integer,intent(in)  ::it
real(DP),intent(in) ::amt(:,:)
real(DP),allocatable::amm(:)
integer             ::iat,iw

allocate( amm(natom) )
      do iw=1,nwalk
       do iat=1,natom
        amm(iat)=amt(iat,iw)/AMU
       enddo
      enddo

call plumed_f_gcmd("setStep"//char(0),it)
call plumed_f_gcmd("setMasses"//char(0),amm)

end subroutine plumed_stepinfo


subroutine force_plumed(x,y,z,fx,fy,fz,eclas)
 use mod_general, only: natom,it
 use mod_const, only:AMU
 use mod_system, only:am
implicit none
real(DP),intent(inout)  :: x(:,:), y(:,:), z(:,:)
real(DP),intent(inout)  :: fx(:,:),fy(:,:),fz(:,:)
real(DP),intent(in)     :: eclas
real(DP),allocatable    :: xx(:),yy(:),zz(:),fxx(:),fyy(:),fzz(:),amm(:)
real(DP)                :: pbox(3,3),pcharges(natom),plumvirial(3,3),fconv,gpp(natom*3)
integer                 :: iat,a

allocate(xx(natom),yy(natom),zz(natom),fxx(natom),fyy(natom),fzz(natom),amm(natom))

      do a=1,natom
       amm(a)=am(a)/AMU
      enddo

      fconv=1.88972599                       ! force from Eh/bohr to Eh/ang
!     fconv=2625.5*18.89725                  ! force from Eh/bohr to kJ/mol/nm
      do iat=1,natom
        xx(iat)=x(iat,1)
        yy(iat)=y(iat,1)
        zz(iat)=z(iat,1)
        fxx(iat)=fx(iat,1)*fconv
        fyy(iat)=fy(iat,1)*fconv
        fzz(iat)=fz(iat,1)*fconv
        pcharges(iat)=0
      enddo

write(*,*)'Input forces x are:',fx
write(*,*)'Forces x for PLUMED are:',fxx

      pbox=0.0
      plumvirial=0.0

call plumed_f_gcmd("setStep"//char(0),it)
call plumed_f_gcmd("setPositionsX"//char(0),xx)
call plumed_f_gcmd("setPositionsY"//char(0),yy)
call plumed_f_gcmd("setPositionsZ"//char(0),zz)
call plumed_f_gcmd("setMasses"//char(0),amm)
call plumed_f_gcmd("setCharges"//char(0),pcharges)  ! setting zeros, abin is not capable of this
call plumed_f_gcmd("setVirial"//char(0),plumvirial) ! setting zeros, could cause problem in NVP
call plumed_f_gcmd("setBox"//char(0),pbox)
call plumed_f_gcmd("setEnergy"//char(0),eclas)
call plumed_f_gcmd("setForcesX"//char(0),fxx)
call plumed_f_gcmd("setForcesY"//char(0),fyy)
call plumed_f_gcmd("setForcesZ"//char(0),fzz)

call plumed_f_gcmd("calc"//char(0),0); 

write(*,*)'Obtained forces x are:',fxx

      do iat=1,natom
!        x(iat,1)=xx(iat)
        fx(iat,1)=fxx(iat)/fconv
        fy(iat,1)=fyy(iat)/fconv
        fz(iat,1)=fzz(iat)/fconv
      enddo
write(*,*)'Output forces to MD x are:',fx

end subroutine force_plumed
#endif

end module

