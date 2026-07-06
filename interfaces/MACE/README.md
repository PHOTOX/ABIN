# MACE MPI Interface for ABIN

This directory contains the Python-based server for the MACE (Machine Learning Atomic Cluster Expansion) potential.

## Requirements

- Python >= 3.8
- PyTorch >= 1.12
- mace-torch
- ase
- mpi4py
- numpy

The exact dependencies are specified as inline metadata in `mace_server.py`.
We highly recommend installing the dependencies in a fresh virtual environment
using the [uv package manager](https://github.com/astral-sh/uv):

```console
# Install uv first, https://github.com/astral-sh/uv#installation
uv venv  # Creates a new virtual environment in .venv folder
uv pip install -r mace_server.py --torch-backend=auto
source .venv/bin/activate   # Activates the environment
```

> [!IMPORTANT]
> The tricky part is to install the correct binary version of PyTorch, as it needs to match your CUDA version.
> The `--torch-backend=auto` should autodetect your environment and install the correct version.
> In HPE clusters, make sure you load your CUDA environment before running the installation.
> Please read the [uv PyTorch documentation](https://docs.astral.sh/uv/guides/integration/pytorch/#automatic-backend-selection) for more information.
> If your machine doesn't have a GPU accelerator, the installation will automatically pick up the CPU-only build.

ABIN itself must be compiled using the MPICH compiler, see top-level [README.md](../../README.md) for instructions.

After installation, run the MACE tests to make sure the basic communication works.
```console
$ make test TEST="MACE MACE_ERROR MACE_ERROR2"
Running tests in directories:
MACE MACE_ERROR MACE_ERROR2
MACE	PASSED
=======================
MACE_ERROR	PASSED
=======================
MACE_ERROR2	PASSED
=======================
 
3 tests PASSED.
```

These tests use a fake potential, not the actual MACE implementation.
Make sure you test the energy conservation using a short NVE dynamics
before running production calculations!

## Usage

The MACE server must be launched alongside ABIN using mpirun.
It is recommended to use the provided launch script `utils/run.mace_mpi_abin.sh`.

Rough steps:

1. Specify `pot=_mace_` in the ABIN input file, everything else is the same.
2. Copy `interfaces/MACE/mace_server.py` and `utils/run.mace_mpi_abin.sh` into your working directory.
3. In the launch script, specify the ABIN input files, and paths to MPICH and Python installation.
4. In the launch script, specify path to your MACE model file, and whether to evaluate the model on CPU or CUDA GPU.
5. Profit!
