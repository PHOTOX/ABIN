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
 plumed_energyUnits=2600       !conversion Ha -> kJ/mol
 plumed_lengthUnits=0.1/ANG    !conversion Angstrom*ANG -> nm
 plumed_timeUnits=AUTOFS*0.001 ! au time -> ps
 plumedoutfile = 'plumed.out'

 !plumed init
 call plumed_f_gcreate()
 call plumed_f_gcmd("setRealPrecision"//char(0),8)
 call plumed_f_gcmd("setKbT"//char(0),plumed_kbt)
 call plumed_f_gcmd("setMDEnergyUnits"//char(0),plumed_energyUnits)
 call plumed_f_gcmd("setMDLengthUnits"//char(0),plumed_lengthUnits)
 call plumed_f_gcmd("setMDTimeUnits"//char(0),plumed_timeUnits)
 call plumed_f_gcmd("setPlumedDat"//char(0),trim(adjustl(plumedfile))//char(0))
 call plumed_f_gcmd("setNatoms"//char(0),natom)
 call plumed_f_gcmd("setMDEngine"//char(0),"abin")
 call plumed_f_gcmd("setTimestep"//char(0),dt)
 call plumed_f_gcmd("setLogFile"//char(0),trim(adjustl(plumedoutfile))//char(0)) 
 
 call plumed_f_gcmd("init"//char(0),0);
end subroutine plumed_init
#endif

!subroutine force_plumed(x,y,z,fx,fy,fz,eclas)

!end subroutine force_plumed





end module

