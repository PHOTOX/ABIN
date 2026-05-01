"""
Mock MACE MPI server for E2E testing.

This server mimics the MPI protocol of the real mace_server.py
but returns hardcoded energy and forces instead of running the
actual MACE ML model. This avoids the need for mace-torch/PyTorch.

The returned values use a simple harmonic potential for water:
energy = 0.5 * sum(forces^2) (not actually, just fixed values)
"""

import sys
import numpy as np
from mpi4py import MPI

MACE_TAG_EXIT = 666
MACE_TAG_DATA = 2

# Atomic unit conversions
BOHR_TO_ANG = 0.529177249
EV_TO_HARTREE = 1.0 / 27.211399


def main():
    comm = MPI.COMM_WORLD

    # Open MPI port and write to file
    port_name = MPI.Open_port()
    print(f"[MockMACE]: Port opened: {port_name}", flush=True)

    port_file = "mace_port.txt.1"
    with open(port_file, "w") as f:
        f.write(port_name)

    # Accept connection from ABIN
    print("[MockMACE]: Waiting for ABIN...", flush=True)
    abin_comm = comm.Accept(port_name)
    print("[MockMACE]: Connected!", flush=True)

    # Receive number of atoms
    natom_buf = np.empty(1, dtype=np.intc)
    abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
    natom = int(natom_buf[0])
    print(f"[MockMACE]: natom = {natom}", flush=True)

    # Receive atom types
    atom_type_buf = bytearray(natom * 2)
    abin_comm.Recv([atom_type_buf, MPI.CHAR], source=0, tag=MACE_TAG_DATA)
    atom_types_str = atom_type_buf.decode('ascii')
    print(f"[MockMACE]: atom types = {atom_types_str}", flush=True)

    # Receive config string
    config_status = MPI.Status()
    abin_comm.Probe(source=0, tag=MACE_TAG_DATA, status=config_status)
    config_len = config_status.Get_count(MPI.CHAR)
    config_buf = bytearray(config_len)
    abin_comm.Recv([config_buf, MPI.CHAR], source=0, tag=MACE_TAG_DATA)
    config_str = config_buf.decode('ascii')
    print(f"[MockMACE]: config = {config_str}", flush=True)

    # Main loop
    step = 0
    max_steps = 100
    while step < max_steps:
        status = MPI.Status()
        abin_comm.Probe(source=0, tag=MPI.ANY_TAG, status=status)

        if status.Get_tag() == MACE_TAG_EXIT:
            print("[MockMACE]: Exit signal received", flush=True)
            abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_EXIT)
            break

        # Receive natom per step
        abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)

        # Receive coordinates (3*natom doubles in Bohr)
        coords = np.empty((3, natom), dtype=np.float64)
        abin_comm.Recv([coords, MPI.DOUBLE], source=0, tag=MACE_TAG_DATA)

        step += 1

        # Compute mock energy and forces using simple harmonic potential
        # E = 0.5 * k * sum(r^2) where r is displacement from origin
        k = 0.01  # force constant in Hartree/Bohr^2
        coords_t = coords.T  # (natom, 3)

        # Simple harmonic energy around the center of mass
        com = np.mean(coords_t, axis=0)
        displ = coords_t - com
        energy = 0.5 * k * np.sum(displ ** 2)

        # Forces = -gradient = -k * displacement
        forces = -k * displ  # (natom, 3)

        # Send energy
        energy_buf = np.array([energy], dtype=np.float64)
        abin_comm.Send([energy_buf, MPI.DOUBLE], dest=0, tag=0)

        # Send forces (3, natom) for Fortran column-major
        forces_send = forces.T.copy()
        abin_comm.Send([forces_send, MPI.DOUBLE], dest=0, tag=0)

        print(f"[MockMACE]: Step {step}, energy = {energy:.10f} Hartree", flush=True)

    # Cleanup
    abin_comm.Disconnect()
    MPI.Close_port(port_name)
    print("[MockMACE]: Server stopped.", flush=True)


if __name__ == "__main__":
    main()
