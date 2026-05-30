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
Or use the provided script in the root directory:
```console
./dev_scripts/install_mace.sh
```

## Usage

The MACE server is typically launched alongside ABIN using MPI. ABIN communicates with the server to obtain energies and forces.

### Manual Launch (OpenMPI)

When using OpenMPI, you need an `ompi-server` running for the connection handshake:

```console
# 1. Start ompi-server
ompi-server --no-daemonize -r ompi_uri.txt &

# 2. Start MACE server
mpirun --ompi-server file:ompi_uri.txt -n 1 python3 mace_server.py &

# 3. Start ABIN (in another terminal or same script)
mpirun --ompi-server file:ompi_uri.txt -n 1 abin -i input.in -x geom.xyz
```

### Automatic Launch

It is recommended to use the provided launch script in `utils/`:
```console
./utils/run.mace_mpi_abin.sh
```

## Configuration

MACE settings are controlled via the `&mace` namelist in the ABIN input file (`input.in`). Below is the full list of available parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mace_model` | String | `'MACE-OFF23_medium.model'` | Path to a local `.model` file (TorchScript or PyTorch checkpoint) or a foundation model name (e.g., `MACE-OFF23_medium.model`, `medium`). Foundation models are automatically downloaded. |
| `mace_device` | String | `'cpu'` | Execution device: `'cpu'` or `'cuda'`. |
| `mace_default_dtype` | String | `'float64'` | Default tensor precision: `'float32'` or `'float64'`. |
| `mace_batch_size` | Integer | `64` | Batch size for evaluations. |
| `mace_compute_stress` | Logical | `.false.` | Whether to compute the stress tensor. |
| `mace_return_contributions` | Logical | `.false.` | Whether to return energy contributions per body order. |
| `mace_info_prefix` | String | `'MACE_'` | Prefix for MACE-related metadata in the ASE output. |
| `mace_head` | String | `''` | Model head identifier (for multi-head models). |
| `mace_max_mpi_wait_time` | Real | `60.0` | Maximum time (seconds) to wait for the MACE port file to appear. |
| `mace_mpi_milisleep` | Integer | `50` | Sleep interval (milliseconds) when polling for the MACE port or results. |

See `sample_inputs/input.in.mace` for a template.

## MACE License
The provided MACE model (MACE-OFF23_medium.model) is downloaded from the [MACE-OFF GitHub repository](https://github.com/ACEsuit/mace-off) and is licensed under the [Academic Software License Agreement (ASLA)](https://github.com/ACEsuit/mace-off/blob/main/LICENSE.md).

Note that this is not our model and we are not responsible for its performance or any issues related to its use. The purpose of its inclusion here is to provide a working example of how to use MACE with ABIN.

For full terms and conditions, please refer to the [MACE-OFF GitHub repository](https://github.com/ACEsuit/mace-off).