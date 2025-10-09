This is a sample input file for ABIN for Classical MD simulation --
-- plain old Newton's equations of motion.

&general
pot='g09'		! where do we obtain forces and energies?
mdtype='md',		! md = classical MD
nstep=50000,		! number of steps
dt=20.,			! timestep (a.u.)

irest=0,		! should we restart from restart.xyz? (ignoring input coordinates and velocities)

nwrite=1,		! how often some output should be printed (energies etc.)
                        ! 1 = every time step
                        ! 10 = each 10th time step
nwritex=5,		! how often should we print coordinates?
nrest=5,		! how often we print restart files?
nwritev=0,		! how often we print velocities? (by default ABIN does not print velocities)
/


&thermostat
temp=298.15,		! temperature [K] for Maxwell-Boltzmann sampling and thermostat
therm='nhc',		! Thermostat options: 
                        ! 'none'= NVE( microcanonical)
                        ! 'nhc' = Nose-Hoover Chains
                        ! 'gle' = Generalized Langevin Equation (GLE)
tau0=0.001,		! relaxation time of NHC thermostat in picoseconds
/
