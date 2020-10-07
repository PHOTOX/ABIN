[![DOI](https://zenodo.org/badge/28882168.svg)](https://zenodo.org/badge/latestdoi/28882168)
[![CI](https://github.com/PHOTOX/ABIN/workflows/GFortran%20CI/badge.svg?branch=master&event=push)](https://github.com/PHOTOX/ABIN/actions?query=workflow%3A%22GFortran+CI%22)
[![codecov](https://codecov.io/gh/PHOTOX/ABIN/branch/master/graph/badge.svg)](https://codecov.io/gh/PHOTOX/ABIN)

----------------
1. What is ABIN?
----------------

ABIN is a program for performing ab initio molecular dynamics.
It was designed specifically to deal with quantum nuclear effects.
It can do path integral simulations and also utilizes quantum thermostat based on General Langevin equation.
It can also simulate non-adiabatic events using Surface-hoping algorithm.

The basic philosophy of this program is simple.
While the program itself handles the propagation of the system according to the equations of motion,
the forces and energies are taken from the external electronic structure program such as ORCA or TeraChem.
The call to the chosen external program is handled via a simple shell script.
Therefore, writing a new interface is rather straightforward
and can be done without any changes to ABIN or the ab initio code.

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
To compile the code, you need a gfortran and g++ compiler version >=4.3.
Version 4.6 and higher is recommended.

For some features, you will also need to install the FFTW library.
It is usually provided with your Linux distribution.
It can also be downloaded here: http://www.fftw.org/

The compilation can be as easy as:

`$ ./configure && make`

You can test the installation by:

`$ make test`

If you modify the source code and want to recompile, you should always clean up by:

`$ make clean`

-------------------------------
3. STRUCTURE OF THE SOURCE CODE
-------------------------------

Short descriptions of key source files.

| Path     | Description |
|----------|-------------|
| SAMPLE/         | Sample input files. |
| INTERFACES/     | BASH interfaces to common ab initio codes. |
| UTIL/           | Handy scripts that might be useful in conjuction with the MD code. |
| abin.F90        | Main routine. |
| modules.F90     | Modules containing all basic variables. |
| init.F90        | Big ugly init routine. Reads input parameters, restart files and checks for errors in input. |
| mdstep.f90      | Routines for propagation of equations of motion (velocity Verlet and RESPA). |
| forces.f90      | Main routine for getting forces (force_clas) and quantum forces in PIMD (force_quant). |
| force_abin.f90  | Routine that calls ab initio BASH interfaces. |
| nosehoover.f90  | Nosé-Hoover thermostat. |
| surfacehop.f90  | Surface Hopping dynamics |
| sh_integ.f90    | Propagation of electronic wavefunction in SH dynamics | 
| landau_zener.f90| Landau-Zener dynamics |
| gle.F90         | Generalized Langevin Equation thermostat |
| shake.f90       | Constraints using SHAKE algorithm. |
| analysis.f90    | Driver routine for printing output. |
| arrays.f90      | Module containing all dynamically allocated arrays.| 
| random.f90      | Random number generator. |
| vinit.f90       | Initial velocities. |
| potentials.f90  | Analytical potentials and hessians for harmonic and Morse oscillators and others. |
| density.f90     | Evaluates distance, angle and dihedrals distributions during dynamics. | 
| estimators.f90  | Energy and heat capacity estimators for PIMD. | 
| transform.F90   | Coordinate transformations for PIMD. |
