This is a sample input file for Surface-Hopping simulation in ABIN 

&general
pot='molpro'		! where do we obtain forces? options are:g09,orca,tera,turbo,molpro,nab,harm,morse,guillot,2dho 
natom=8,		!number of atoms
irest=0,		! should we restart from restart.xyz? (ignoring mini.dat)

ipimd=2,		! classical simulation 0, quantum simulation 1, surface-hopping 2, steepest descent 3 
nstep=1000,		! number of steps
dt=20.,			! timestep [au]
irandom=85631316,  	! random seed

nwrite=1,		! how often some output should be printed (estimators etc.)
nwritex=1,		! how often should we print coordinates?
nrest=1,		! how often we print restart files?
nwritev=0,		! how often we print velocities?
/

&nhcopt
inose=0,		! Thermostating: Nose-Hoover 1, microcanonical 0,GLE 2, Langevin 3
!temp=0.00,		! Usually, you would take initial velocities from WIGNER or set them to zero
/

&sh
istate_init=2,		! initial electronic state
nstate=3,		! number of electronic states
deltaE=2.0,		! maximum energy difference [eV], for which we calculate NA coupling
PopThr=0.001,           ! minimum population of either state, for which we compute NA coupling
EnergyDifThr=1.00,      ! maximum energy difference between two consecutive steps
EnergyDriftThr=1.00,    ! maximum energy drift from initial total energy
substep=100,		! number of substeps for solving ESCH
/

----------------------END OF INPUT---------------------------------------------
