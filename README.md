[![DOI](https://zenodo.org/badge/28882168.svg)](https://zenodo.org/badge/latestdoi/28882168)

----------------
1. What is ABIN?
----------------

ABIN is a program for performing ab initio molecular dynamics. It was designed specifically to deal with quantum nuclear effects. It can do path integral simulations and also utilizes recently developed quantum thermostat based on General Langevin equation. It can also simulate non-adiabatic events using Surface-hoping algorithm.

The basic philosophy of this program is simple. While the program itself handles the propagation of the system according to the equations of motion, the forces and energies are taken from the external electronic structure program such as ORCA or TeraChem. The call to the chosen external program is handled via a simple shell script. Therefore, writing a new interface is rather straightforward and can be done without any changes to the ABIN code or ab initio code.

The code is provided under the GNU General Public License.
A full text of the license can be found in the file LICENCE.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

A full documentation can be found in the folder DOC.

---------------
2. INSTALLATION
---------------
To compile the code, you need a gfortran and gcc compiler with a version at least 4.3.
However, version 4.6 and higher is recommended.
If your distribution does not have the required version of the compilers, it can be downloaded from this website:
https://gcc.gnu.org/wiki/GFortranBinaries
The site offers compiled binaries and required libraries so compilation is typically not needed.

For some features, you will also need to install the FFTW library. It is usually provided with your Linux distribution.
It can also be downloaded here: http://www.fftw.org/

The compilation can be as easy as: 
$ cd src; make

You can test the installation by:
$ make test

If you modify the source code and want to recompile, you should always clean up by:
$ make  clean

-------------------------------
3. STRUCTURE OF THE SOURCE CODE
-------------------------------

Here you find the information about the key source files.

SAMPLE/     - Contains sample input files.

INTERFACES/ - Contains bash interfaces to common ab initio codes.

UTIL/       - Contains some handy scripts that might be useful in conjuction with the MD code.


abin.f90       - Main routine.

modules.f90    - Modules containing all basic variables.

init.f03       - Big ugly initialization routine. It reads the input parameters, restart files and checks for errors in input.

mdstep.f03     - Routines for propagation of equations of motion (velocity Verlet and RESPA).

forces.f03     - Main routine for getting forces (force_clas) and quantum forces in PIMD (force_quant).

force_abin.f90 - Routine that calls ab initio BASH interfaces.

nosehoover.f03 - Nosé-Hoover thermostat.

surfacehop.f90 - Surface Hopping dynamics

sh_integ.f90   - Propagation of electronic WF in SH dynamics 

gle.F90        - GLE thermostat

shake.f03      - Contraints using SHAKE algorithm.

analysis.f90   - Driver routine for printing output.

NAB/           - Contains force field routines taken from AmberTools package.

EWALD/         - Code for Ewald summation taken from Macsimus package. Not tested.


other files:

minimizer.f90  - Naive steepest-descent minimizer.

random.f90     - Random number generator.

arrays.f90     - Module with all dynamically allocated arrays.

vinit.f90      - Initializes initial velocities.

potentials.f90 - Contains analytical potentials and hessians for harmonic and Morse oscillators.

ekin.f90       - Evaluates current kinetic energy and temperature.

density.f90    - Evaluates distance, angle and dihedrals distributions during dynamics.

estimators.f90 - Energy and heat capacity estimators for PIMD.

transform.F90  - Coordinate transformations for PIMD.

utils.f90      - Some general purpose routines, e.g. Get_distance() and abinerror()    

