[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1228462.svg)](https://zenodo.org/badge/latestdoi/28882168)
[![CI](https://github.com/PHOTOX/ABIN/workflows/GFortran%20CI/badge.svg?branch=master&event=push)](https://github.com/PHOTOX/ABIN/actions?query=workflow%3A%22GFortran+CI%22)
[![codecov](https://codecov.io/gh/PHOTOX/ABIN/branch/master/graph/badge.svg)](https://codecov.io/gh/PHOTOX/ABIN)

## What is ABIN?

ABIN is a program for performing ab initio molecular dynamics.
It is a general purpose program that was initially designed to model nuclear quantum effects (NQE).
NQE can be most rigirously captured with path integral MD (PIMD), but also within the Quantum Thermostat based on General Langevin Equation framework developed by Michele Cerriotti.
ABIN can also simulate non-adiabatic events using Surface-hoping algorithm, using either the classical fewest-switches algorithm (FSSH) or simpler Landau-Zener approach which does not require non-adiabatic couplings. The LZ approach can also capture singlet-triplet transitions.

The basic philosophy of ABIN program is simple — 
while the program itself handles the propagation of the system according to the equations of motion,
the forces and energies are taken from an external electronic structure program such as ORCA or TeraChem.
The call to the chosen external program is handled via a simple shell script interface.
Therefore, writing a new interface is straightforward
and can be done without any changes to ABIN or the ab initio code.

The code is provided under the GNU General Public License.
A full text of the license can be found in the file LICENCE.

The documentation (work-in-progress) can be found in the folder `docs/`.

## Installation

To compile the code, you'll need a Fortran and C++ compiler.
GNU compilers are tested the most, GFortran and g++ compiler versions >=7.0.
Intel compiler (ifort) is supported since version >=2018,
including the newly open-sourced versions (Intel OneAPI).

The compilation can be as easy as:

```console
$ ./configure && make
```

Always test the installation by running the test suite:

```console
$ make test
```
If you modify the source code and want to recompile,
you should always clean up before the recompilation:

```console
$ make clean && make
```

## Running the code

After compilation, the executable binary should be available in `bin/abin`.
You can copy the binary somewhere in your PATH, or add bin/ to your PATH.
To execute an MD simulation using an input file `input.in` and initial XYZ coordinates in file `geom.xyz`, run:

```console
bin/abin -i input.in -x geom.xyz
```

Run `bin/abin --help` to see other options.
Example input files for various types of simulations can be found in `sample_inputs/`.

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
 - [TCPB-CPP](https://github.com/mtzgroup/tcpb-cpp): [EXPERIMENTAL] TCPB interface to TeraChem


## Structure of the repository

| Path             | Description                                  |
|------------------|----------------------------------------------|
| src/             | ABIN source code
| sample\_inputs   | Sample input files.
| interfaces/      | BASH interfaces to common quantum chemistry codes.
| utils/           | Handy scripts that might be useful in conjuction with the MD code.
| unit\_tests/     | Unit tests; run by `make unittest` (needs pFUnit library installed)
| tests/           | End-to-End tests; run by `make e2etest`
| dev\_scripts/    | Setup for ABIN devs and install scripts for optional libraries.

## For developers

Contributions are very much welcome! Please see our [contribution guide](CONTRIBUTING.md).
