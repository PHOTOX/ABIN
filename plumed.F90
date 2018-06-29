module mod_plumed
  use mod_const, only: DP
  use mod_nhc,   only: temp
  implicit none
  public
  integer :: iplumed
  character(len=40) :: plumedfile=''


CONTAINS

subroutine plumed_init(dt)
   use mod_general, only: natom, nwalk, irest
   use mod_const,   only: ANG, AUTOFS
   implicit none
   real, parameter :: kB=8.3144621d-3
   real(DP) :: dt, plumed_kbt, plumed_energyUnits, plumed_lengthUnits, plumed_timeUnits
   character(256) :: plumedoutfile
   integer  :: api_version

   !variables computation
   plumed_kbt = temp * 8.3144598d0     ! in kJ/mol
   plumed_energyUnits = 2625.5697d0    !conversion Ha -> kJ/mol
   plumed_lengthUnits = 0.1d0 / ANG    !conversion Angstrom*ANG -> nm
   plumed_timeUnits = AUTOFS * 0.001d0   ! au time -> ps
   plumedoutfile = 'plumed.out'

#ifdef PLUM
   ! plumed init
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

    ! TODO: 
     
    ! Require at least version 2
    ! plumed_f_gcmd("getApiVersion"//char(0),api_version);   ! Pass the api version that plumed is using

    ! This is valid only if API VERSION > 2
    ! plumed_cmd(plumedmain,"setRestart",&res);                      // Pointer to
    ! an integer saying if we are restarting (zero means no, one means yes)

 
   call plumed_f_gcmd("init"//char(0),0);
#endif
end subroutine plumed_init


subroutine force_plumed(x, y, z, fx, fy, fz, eclas)
   use mod_general,  only: natom, it, nwrite
   use mod_const,    only: AMU
   use mod_system,   only: am
   implicit none
   real(DP),intent(in)     :: x(:,:), y(:,:), z(:,:)
   real(DP),intent(inout)  :: fx(:,:),fy(:,:),fz(:,:)
   real(DP),intent(in)     :: eclas
   real(DP)    :: xx(size(x,1)), yy(size(y,1)),zz(size(z,1))
   real(DP)    :: fxx(size(fx,1)),fyy(size(fy,1)),fzz(size(fz,1))
   real(DP)    :: amm(natom)
   real(DP)    :: pbox(3,3),pcharges(natom),plumvirial(3,3),gpp(natom*3)
   integer     :: iat, iw

   ! Apply plumed only to centroid in PIMD
   ! Normal modes needed
   iw = 1

   do iat = 1, natom
      amm(iat) = am(iat) / AMU
   enddo

   do iat=1,natom
      xx(iat) = x(iat,iw)
      yy(iat) = y(iat,iw)
      zz(iat) = z(iat,iw)
      fxx(iat) = fx(iat,iw)
      fyy(iat) = fy(iat,iw)
      fzz(iat) = fz(iat,iw)
      pcharges(iat) = 0
   enddo

   pbox = 0.0d0
   plumvirial = 0.0d0

#ifdef PLUM

   call plumed_f_gcmd("setStep"//char(0),it)
   call plumed_f_gcmd("setPositionsX"//char(0),xx)
   call plumed_f_gcmd("setPositionsY"//char(0),yy)
   call plumed_f_gcmd("setPositionsZ"//char(0),zz)
   call plumed_f_gcmd("setMasses"//char(0),amm)
   call plumed_f_gcmd("setCharges"//char(0),pcharges)  ! setting zeros, abin is not capable of this
   call plumed_f_gcmd("setVirial"//char(0),plumvirial) ! setting zeros, could cause problem in NpT
   call plumed_f_gcmd("setBox"//char(0),pbox)
   call plumed_f_gcmd("setEnergy"//char(0),eclas)
   call plumed_f_gcmd("setForcesX"//char(0),fxx)
   call plumed_f_gcmd("setForcesY"//char(0),fyy)
   call plumed_f_gcmd("setForcesZ"//char(0),fzz)
   call plumed_f_gcmd("calc"//char(0),0); 
! http://plumed.github.io/doc-v2.0/developer-doc/html/class_p_l_m_d_1_1_plumed_main.html#a4da81a378216f0739fe50bf620d72655
   ! TODO: call plumed_cmd(plumedmain,"fflush"//char(0)); every nwrite

#endif


   do iat=1,natom
      fx(iat,iw) = fxx(iat)
      fy(iat,iw) = fyy(iat)
      fz(iat,iw) = fzz(iat)
   enddo

end subroutine force_plumed

end module

