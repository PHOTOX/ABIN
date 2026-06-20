# MACE MPI Interface for ABIN

This directory contains the Python-based server for the MACE (Machine Learning Atomic Cluster Expansion) potential.

## Requirements

- Python >= 3.8
- PyTorch >= 1.12
- mace-torch
- ase
- mpi4py
- numpy

You can install all dependencies using:
```console
pip install mace-torch ase mpi4py numpy
```

## Usage

The MACE server is typically launched alongside ABIN using MPI. ABIN communicates with the server to obtain energies and forces.

### Automatic Launch

It is recommended to use the provided launch script in `utils/`:
```console
./utils/run.mace_mpi_abin.sh
```
