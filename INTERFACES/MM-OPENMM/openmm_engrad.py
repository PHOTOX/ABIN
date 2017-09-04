# --------------------------------------------
# Single point MM energy and forces using OpenMM. 
# To be used with ABIN

# Before running the actual MD simulation, 
# you should run this script with your template PDB file to verify that everything works
# ---------------------------------------------

from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Starting configuration to use.
# Need to make this an input variable
pdb = PDBFile('input001.pdb')

# Force field to use. OPENMM has AMBER force fields by default
forcefield = ForceField('parm99.xml','tip3p.xml')
# Note that for the qTIP4P/F model we need to solve the issue of the virtual site
# Probably, we need to set the virtual site correctly in ABIN at each step
# (could do it externally in r.mm, but will need it eventually in ABIN anyway
# for proper QMMM interface with TeraChem

# See this for the parameters of qTIP4P/Fw
# https://web.stanford.edu/group/markland/software.html

# Initialize the parameters of the simulation. 
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,nonbondedCutoff=2.0*nanometer,constraints=None,rigidWater=False)

# But we will not do actuall simulation
integrator = VerletIntegrator(0.0002*picoseconds)

# By default, let's just use CPU's
platform = Platform.getPlatformByName('CPU')
properties = {}
# Use CUDA Platform in mixed precision.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'mixed'}

# Setup simulation
simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

state = simulation.context.getState(getEnergy=True, getForces=True)
energy = state.getPotentialEnergy().value_in_unit(kilojoules/mole)
forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)

AUtoKCAL = 627.50946943
KJMOLtoAU = 1 / 4.184 / AUtoKCAL
NMtoBOHR = 10 * 1.889726132873

# Energy in atomic units
print(energy * KJMOLtoAU)

# Forces in atomic units
for f in forces:
   print('%g %g %g' % (-f[0]*KJMOLtoAU/NMtoBOHR, -f[1]*KJMOLtoAU/NMtoBOHR, -f[2]*KJMOLtoAU/NMtoBOHR))


