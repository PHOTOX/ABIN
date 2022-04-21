This is a sample input file for ABIN for Classical MD simulation

&general
pot='g09'		! where do we obtain forces? options are:g09,orca,tera,turbo,molpro,orca,qchem 
mdtype='md',		! md = classical MD
nstep=50000,		! number of steps
dt=20.,			! timestep [au]

irest=0,		! should we restart from restart.xyz? (ignoring input coordinates and velocities)

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
