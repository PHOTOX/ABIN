[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1228462.svg)](https://zenodo.org/badge/latestdoi/28882168)
[![CI](https://github.com/PHOTOX/ABIN/workflows/GFortran%20CI/badge.svg?branch=master&event=push)](https://github.com/PHOTOX/ABIN/actions?query=workflow%3A%22GFortran+CI%22)
[![codecov](https://codecov.io/gh/PHOTOX/ABIN/branch/master/graph/badge.svg)](https://codecov.io/gh/PHOTOX/ABIN)

## What is ABIN?

ABIN is a program for performing ab initio molecular dynamics.
It was designed specifically to deal with quantum nuclear effects.
It can do path integral simulations and also utilizes quantum thermostat based on General Langevin equation.
It can also simulate non-adiabatic events using Surface-hoping algorithm.

The basic philosophy of ABIN program is simple.
While the program itself handles the propagation of the system according to the equations of motion,
the forces and energies are taken from an external electronic structure program such as ORCA or TeraChem.
The call to the chosen external program is handled via a simple shell script interface.
Therefore, writing a new interface is rather straightforward
and can be done without any changes to ABIN or the ab initio code.

The code is provided under the GNU General Public License.
A full text of the license can be found in the file LICENCE.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

A full documentation can be found in the folder DOC.

## Installation

To compile the code, you'll need a Fortran and C++ compiler.
GNU compilers are tested the most, GFortran and g++ compiler versions >=7.0.
Versions >=5.4 likely work as well, but always run the test suite to verify.
Intel compiler (ifort) is supported since version >=2018,
including the newly open-sourced versions (Intel OneAPI).

The compilation can be as easy as:

`$ ./configure && make`

Always test the installation by running the test suite:

`$ make test`

If you modify the source code and want to recompile,
you should always clean up before the recompilation:

`$ make clean && make`

## Optional dependencies

Some functionality relies on external libraries. These are optional,
and the code automatically recognizes which feature is supported for a given build.
Run `configure -h` to see all the options and how to configure them.

To install the libraries, you can use the install scripts in `dev_scripts/`.
We use these in our Continuous Integration testing suite on Github using the Ubuntu 18.04 image.

The optional libraries are:
 - [MPICH](https://www.mpich.org/): An MPI implementation used for Replica Exchange MD and MPI interface with TeraChem.
      - If you just need REMD you can also use other MPI libraries such as OpenMPI or IntelMPI.
 - [FFTW](http://www.fftw.org/): Fast Fourier Transform library used for normal mode transformation in Path Integral MD.
 - [PLUMED](https://www.plumed.org/): A collection of very useful tools for free energy calculations (MetaDynamics, Umbrella Sampling etc).

For some features, you will also need to install the FFTW library.
It is usually provided with your Linux distribution,
but can also be downloaded from http://www.fftw.org/


## Structure of the repository

| Path             | Description                                  |
|------------------|----------------------------------------------|
| src/             | ABIN source code
| sample\_inputs   | Sample input files.
| interfaces/      | BASH interfaces to common _ab initio_ codes.
| utils/           | Handy scripts that might be useful in conjuction with the MD code.
| unit\_tests/     | Unit tests; run by `make unittest` (needs pFUnit library installed)
| tests/           | End-to-End tests; run by `make e2etest`
| dev\_scripts/    | Setup for ABIN devs and install scripts for optional libraries.
