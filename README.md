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
To compile the code, you'll need a GFortran and g++ compiler version >=7.0.
Earlier versions might work, but always run test to verify.
Intel compiler might work, but there is a know issue with restarting
GLE based simulations, so this functionality is automatically disabled
when using ifort.

For some features, you will also need to install the FFTW library.
It is usually provided with your Linux distribution.
It can also be downloaded here: http://www.fftw.org/

The compilation can be as easy as:

`$ ./configure && make`

You can test the installation by:

`$ make test`

If you modify the source code and want to recompile,
you should always clean up before the recompilation:

`$ make clean && make`

-------------------------------
3. STRUCTURE OF THE SOURCE CODE
-------------------------------

Short descriptions of key source files.

| Path     | Description |
|----------|-------------|
| sample_inputs   | Sample input files. |
| interfaces/     | BASH interfaces to common ab initio codes. |
| utils/          | Handy scripts that might be useful in conjuction with the MD code. |
| unit_tests/     | Unit tests, run by `make unittest`
| tests/          | End-to-End tests, run by `make e2etest`
| src/            | ABIN source code
| dev_scripts/    | Setup for ABIN devs + install scripts for optional libraries needed by ABIN
