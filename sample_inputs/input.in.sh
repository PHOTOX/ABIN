This is a sample input file for Surface-Hopping simulation in ABIN 

&general
pot='molpro-sh'		! where do we obtain forces?
irest=0,		! should we restart from restart.xyz?

mdtype='sh',		! sh = surface-hopping MD
nstep=1000,		! number of steps
dt=20.,			! timestep [au]

nwrite=1,		! how often some output should be printed (energies, surface hopping probabilities etc.)
nwritex=1,		! how often to print coordinates?
nrest=1,		! how often to print restart files?
nwritev=0,		! how often to print velocities?
/

&nhcopt
inose=0,		! NVE ensemble (thermostat turned OFF)
!temp=0.00,		! Usually, you would read initial velocities from previous ground state sampling (-v option)
/

&sh
istate_init=2,		! initial electronic state (1 is ground state)
nstate=3,		! number of electronic states
deltaE=2.0,		! maximum energy difference (eV) between states for which we calculate NA coupling
PopThr=0.001,           ! minimum population of either state, for which we compute NA coupling
EnergyDifThr=0.50,      ! maximum energy difference between two consecutive steps
EnergyDriftThr=0.50,    ! maximum energy drift from initial total energy
/
