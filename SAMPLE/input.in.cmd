This is a sample input file for ABIN for CMD simulation

&general
pot='g09'		! where do we obtain forces? options are:g09,orca,tera,turbo,molpro,orca,qchem 
natom=27,		!number of atoms
ipimd=0,		! classical simulation 0, quantum simulation 1, surface-hopping 2, steepest descent 3 
imini=2000,		! equilibration period (should be at least 2000), or number of steps for minimization
nstep=50000,		! number of steps
dt=20.,			! timestep [au]
irandom=1651563,  	! random seed
irest=0,		! should we restart from restart.xyz? (ignoring mini.dat)

nwrite=1,		! how often some output should be printed (estimators etc.)
nwritex=5,		! how often should we print coordinates?
nrest=5,		! how often we print restart files?
nwritev=0,		! how often we print velocities?
/


&nhcopt
temp=298.15,		! temperature [K] for Maxwell-Boltzmann sampling and thermostat
inose=1,		! Thermostating: Nose-Hoover 1, microcanonical 0,GLE 2
nchain=4,		! number of nose-hoover chains
tau0=0.001,		! relaxation time of NHC thermostat
nrespnose=3,		! number of inner respa steps for NH thermostat
nyosh=7,		! number of steps in suzuki-yoshida scheme,can be 1,3 or 7
/


----------------------END OF INPUT---------------------------------------------
!!!!!EVERYTHING BELOW IS IGNORED!!!!
