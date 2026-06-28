#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "ase>=3.18.0",
#     "mace-torch>=0.3.10",
#     "mpi4py>=4.1.2",
#     "numpy>=1.26.0",
# ]
# ///
"""
MACE MPI SERVER

This server communicates with ABIN via MPI (using mpi4py).

Usage:
  mpirun -n 1 python mace_server.py

  The server writes its MPI port to 'mace_port.txt.1' for ABIN to read.
"""

import argparse
import functools
import sys
import time
from traceback import print_tb
from pathlib import Path

LOG_NAME = "MaceMPIServer"

# MPI Tags (must match Fortran module mod_mace_mpi)
MACE_TAG_EXIT = 666
MACE_TAG_DATA = 2


def log(message, should_print=True):
    msg_formatted = f"[{LOG_NAME}]: {message!s}"
    if should_print:
        print(msg_formatted, flush=True)
    return msg_formatted


def parse_cmd():
    desc = "MACE MPI server for ground state MD with ABIN"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "--model-path",
        type=str,
        required=True,
        help="Path to MACE model file",
    )
    parser.add_argument(
        "--device",
        type=str,
        choices=("cpu", "cuda"),
        required=True,
        help="Device for model inference",
    )
    config = parser.parse_args()
    if not Path(config.model_path).is_file():
        sys.exit(f"ERROR: file '{config.model_path}' not found")
    return config


class MaceModel:
    """
    Manages the MACE ML model for evaluating atomic configurations.
    Configuration is received from ABIN via MPI.
    """

    def __init__(self, config):
        from mace.calculators import MACECalculator

        log("initializing MACE model")

        # Set ASE calculator
        self.calculator = MACECalculator(
            model_paths=config.model_path,
            device=config.device,
        )

    def evaluate(self, atom_types, coords_bohr):
        """
        Evaluate energy and forces for a single configuration.

        Parameters:
            atom_types: list of atomic symbols (e.g. ['H', 'O', 'H'])
            coords_bohr: numpy array of shape (natom, 3) in Bohr

        Returns:
            energy_hartree: energy in Hartree
            forces_hartree_bohr: forces in Hartree/Bohr, shape (natom, 3)
        """
        import ase

        # Unit conversions
        bohr_to_ang = 0.529177249
        ev_to_hartree = 1.0 / 27.211399
        ev_per_ang_to_hartree_per_bohr = ev_to_hartree * bohr_to_ang

        # Convert coordinates from Bohr to Angstrom
        coords_ang = coords_bohr * bohr_to_ang

        # Create ASE atoms object
        pbc = (False, False, False)
        # cell_size = 100.0  # Angstroms
        # cell = ((cell_size, 0, 0), (0, cell_size, 0), (0, 0, cell_size))
        atoms = ase.Atoms(symbols=atom_types, positions=coords_ang, pbc=pbc)
        atoms.set_calculator(self.calculator)

        if self.head is not None:
            atoms.info["head"] = self.head

        energy_hartree = atoms.get_potential_energy() * ev_to_hartree
        forces_hartree_bohr = atoms.get_forces() * ev_per_ang_to_hartree_per_bohr

        return energy_hartree, forces_hartree_bohr

def connect_to_abin():
    """Establish initial connection to ABIN"""
    from mpi4py import MPI

    # Open MPI port and write to file for ABIN to read
    port_name = MPI.Open_port()
    log(f"MPI port opened: {port_name}")

    port_file = "mace_port.txt.1"
    with open(port_file, "w") as f:
        f.write(port_name)
    log(f"Port written to {port_file}")

    # Accept connection from ABIN
    log("Waiting for ABIN to connect...")
    abin_comm = MPI.COMM_WORLD.Accept(port_name)
    log("Connection from ABIN accepted!")
    return port_name, abin_comm


def exception_handler(shutdown_callback, exception_type, exception, traceback):
    """Try to gracefully shutdown communication with ABIN upon uncaught exceptions"""
    print(f"Unexpected {exception_type.__name__}: {exception}")
    print_tb(traceback)
    # Restore original exception handling to prevent endless loop
    # in case of uncaught excpetion during shutdown
    sys.excepthook = None
    shutdown_callback()
    sys.exit(1)


def main(config):
    import numpy as np
    from mpi4py import MPI

    port_name, abin_comm = connect_to_abin()

    def shutdown_communication():
        """Gracefully shutdown communication with ABIN"""
        log("Shutting down communication with ABIN...")
        try:
            abin_comm.Disconnect()
        except BaseException as e:  # noqa: E722
            log(e)
        else:
            log("Disconnected")

        try:
            MPI.Close_port(port_name)
        except BaseException as e:  # noqa: E722
            log(e)
        else:
            log("Port {port_name} close")


    def error_shutdown():
        log("Sending ERROR tag to ABIN")
        error_energy = np.array([0.0], dtype=np.float64)
        # This is best effort only, since ABIN might be dead already, ignore any errors here
        try:
            abin_comm.Send([error_energy, MPI.DOUBLE], dest=0, tag=1)
        except Exception as e:
            log(e)
            pass

        shutdown_communication()

        sys.exit(1)

    # Call error_shutdown upon any unhandled exception
    sys.excepthook = functools.partial(exception_handler, error_shutdown)

    # Receive number of atoms
    natom_buf = np.empty(1, dtype=np.intc)
    abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
    natom = int(natom_buf[0])
    log(f"Received number of atoms: {natom}")

    # Receive atom types
    atom_type_buf = bytearray(natom * 2)
    abin_comm.Recv([atom_type_buf, MPI.CHAR], source=0, tag=MACE_TAG_DATA)
    atom_types_str = atom_type_buf.decode('ascii')
    atom_types = [atom_types_str[i:i + 2].strip() for i in range(0, len(atom_types_str), 2)]
    log(f"Received atom types: {atom_types}")

    # Load MACE model
    mace_model = MaceModel(config)
    log("MACE model ready. Entering main loop.")

    # Main loop: receive coordinates, compute, send results
    eval_count = 0
    while True:
        # Receive natom (sent each step for protocol consistency)
        status = MPI.Status()
        abin_comm.Probe(source=0, tag=MPI.ANY_TAG, status=status)

        if status.Get_tag() == MACE_TAG_EXIT:
            log("Received exit signal from ABIN")
            # Consume the message
            try:
                abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_EXIT)
            except Exception as e:
                log(e)
            break

        abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
        natom_step = int(natom_buf[0])

        if natom_step != natom:
            log(f"ERROR: Received natom={natom_step}, expected {natom}")
            error_shutdown()

        # Receive coordinates (3*natom doubles, in Bohr)
        coords = np.empty((natom, 3), dtype=np.float64)
        abin_comm.Recv([coords, MPI.DOUBLE], source=0, tag=MACE_TAG_DATA)

        coords_bohr = coords.copy()

        eval_count += 1
        log(f"Evaluation {eval_count}")

        energy, forces = mace_model.evaluate(atom_types, coords_bohr)
        log(f"Evaluation {eval_count}: energy = {energy:.15f} Hartree")

        # Send energy (1 double, in Hartree)
        energy_buf = np.array([energy], dtype=np.float64)
        abin_comm.Send([energy_buf, MPI.DOUBLE], dest=0, tag=0)

        # Send forces (3*natom doubles, in Hartree/Bohr)
        # Transpose back to (3, natom) to match Fortran column-major layout
        if forces.dtype != np.float64:
            forces_send = forces.T.astype(np.float64)
        else:
            forces_send = forces.T.copy()
        abin_comm.Send([forces_send, MPI.DOUBLE], dest=0, tag=0)

    shutdown_communication()


if __name__ == "__main__":
    config = parse_cmd()
    main(config)
