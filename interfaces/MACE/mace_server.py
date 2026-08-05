#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
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

  The server writes its MPI port to 'mace_port.txt' for ABIN to read.
"""
# ruff: file-ignore[blind-except]

import argparse
import functools
import logging
import sys
import warnings
from pathlib import Path
from time import perf_counter

def setup_logger(debug=True):
    """Configure standard library logging to output to sys.stdout."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
    logging.captureWarnings(True)
    warnings.formatwarning = (
        lambda msg, cat, fname, lineno, line=None: f"{fname}:{lineno}: {cat.__name__}: {str(msg).strip()}"
    )


# Configure logger at module import time (debug=True by default)
setup_logger(debug=True)


def log_environment_info(config):
    """Log debug information about Python, PyTorch, CUDA, dependencies, and environment."""
    import platform
    import ase
    import mpi4py
    import numpy as np
    import torch

    logging.debug("=== MACE Server Environment ===")
    logging.debug(f"Python Executable : {sys.executable}")
    logging.debug(f"Python Version    : {sys.version.replace('\n', ' ')}")
    logging.debug(f"Platform / OS     : {platform.platform()}")
    logging.debug(f"Working Directory : {Path.cwd()}")
    logging.debug(f"Model Path        : {config.model_path}")
    logging.debug(f"Configured Device : {config.device}")

    # Dependency versions
    logging.debug(f"PyTorch Version   : {torch.__version__}")
    logging.debug(f"ASE Version       : {ase.__version__}")
    logging.debug(f"NumPy Version     : {np.__version__}")
    logging.debug(f"mpi4py Version    : {mpi4py.__version__}")

    try:
        import mace

        logging.debug(f"MACE Version      : {getattr(mace, '__version__', 'unknown')}")
    except ImportError:
        logging.debug("MACE Version      : Not installed")

    # CUDA & Hardware details
    cuda_avail = torch.cuda.is_available()
    logging.debug(f"CUDA Available    : {cuda_avail}")
    if cuda_avail:
        logging.debug(f"CUDA PyTorch Build: {torch.version.cuda}")
        logging.debug(f"GPU Device Count  : {torch.cuda.device_count()}")
        logging.debug(f"GPU Device Name   : {torch.cuda.get_device_name(0)}")
        mem_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        logging.debug(f"GPU Total Memory  : {mem_gb:.2f} GB")
    logging.debug("======================================================")


# MPI Tags (must match Fortran module mod_mace_mpi)
MACE_TAG_EXIT = 666
MACE_TAG_DATA = 2
MACE_TAG_ERROR = 13

MACE_PORT_FILE = "mace_port.txt"


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
    parser.add_argument(
        "--debug",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable or disable debug logging output",
    )
    config = parser.parse_args()
    if not config.debug:
        setup_logger(debug=False)
    model = config.model_path
    if not model.startswith("__MOCK_") and not Path(model).is_file():
        logging.error(f"File '{config.model_path}' not found")
        sys.exit(1)
    return config


class MaceModel:
    """
    Manages the MACE ML model for evaluating atomic configurations.
    Configuration is received from ABIN via MPI.
    """

    def __init__(self, config):
        from mace.calculators import MACECalculator

        logging.info("initializing MACE model")

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
        coords_ang = coords_bohr.copy() * bohr_to_ang

        # Create ASE atoms object
        pbc = (False, False, False)
        # cell_size = 100.0  # Angstroms
        # cell = ((cell_size, 0, 0), (0, cell_size, 0), (0, 0, cell_size))
        atoms = ase.Atoms(symbols=atom_types, positions=coords_ang, pbc=pbc)
        atoms.calc = self.calculator

        energy_hartree = atoms.get_potential_energy() * ev_to_hartree
        forces_hartree_bohr = atoms.get_forces() * ev_per_ang_to_hartree_per_bohr

        return energy_hartree, forces_hartree_bohr


class HarmonicModel:
    """
    This is just a mock for testing, using harmonic potential.
    """

    k = 0.01  # force constant in Hartree/Bohr^2

    def __init__(self, config):
        self.model_path = config.model_path
        logging.info("Using Harmonic Mock Model")

    def evaluate(self, _atom_types, coords_bohr):
        """
        Compute mock energy and forces using simple harmonic potential
        E = 0.5 * k * sum(r^2) where r is displacement from origin
        """
        import numpy as np

        if config.model_path == "__MOCK_ERROR__":
            raise RuntimeError("Simulating error condition")

        coords_t = coords_bohr.T  # (natom, 3)

        # Simple harmonic energy around the center of mass
        com = np.mean(coords_t, axis=0)
        displ = coords_t - com
        energy = 0.5 * self.k * np.sum(displ**2)

        # Forces = -gradient = -k * displacement
        forces = -self.k * displ  # (natom, 3)
        return energy, forces


def connect_to_abin():
    """Establish initial connection to ABIN"""
    from mpi4py import MPI

    # Open MPI port and write to file for ABIN to read
    port_name = MPI.Open_port()
    logging.info(f"MPI port opened: {port_name}")

    with open(MACE_PORT_FILE, "w", encoding="utf-8") as f:
        f.write(port_name)
    logging.info(f"Port written to {MACE_PORT_FILE}")

    # Accept connection from ABIN
    logging.info("Waiting for ABIN to connect...")
    abin_comm = MPI.COMM_WORLD.Accept(port_name)
    logging.info("Connection from ABIN accepted!")
    return port_name, abin_comm


# https://docs.python.org/3/library/sys.html#sys.excepthook
def exception_handler(shutdown_callback, exception_type, exception, traceback):
    """Try to gracefully shutdown communication with ABIN upon uncaught exceptions"""
    logging.error(
        f"Unexpected {exception_type.__name__}: {exception}",
        exc_info=(exception_type, exception, traceback),
    )
    # Restore original exception handling to prevent endless loop
    # in case of uncaught excpetion during shutdown
    sys.excepthook = sys.__excepthook__
    shutdown_callback()
    sys.exit(1)


def main(config):
    import numpy as np
    from mpi4py import MPI

    log_environment_info(config)
    port_name, abin_comm = connect_to_abin()

    def shutdown_communication():
        """Gracefully shutdown communication with ABIN"""
        logging.info("Shutting down communication with ABIN...")
        try:
            abin_comm.Disconnect()
        except Exception as e:
            logging.error(f"Error disconnecting ABIN communicator: {e}")
        else:
            logging.info("ABIN communicator disconnected")

        try:
            MPI.Close_port(port_name)
        except Exception as e:
            logging.error(f"Error closing port {port_name}: {e}")
        else:
            logging.info(f"Port {port_name} closed")

    def error_shutdown():
        logging.warning("Sending ERROR tag to ABIN")
        # This is best effort only, since ABIN might be dead already
        try:
            abin_comm.Send([MPI.BOTTOM, MPI.INT], dest=0, tag=MACE_TAG_ERROR)
        except Exception as e:
            logging.error(f"Error sending ERROR tag to ABIN: {e}")

        shutdown_communication()

        sys.exit(1)

    def check_incoming_msg():
        """Check whether ABIN is sending ERROR or EXIT message"""
        status = MPI.Status()
        abin_comm.Probe(source=0, tag=MPI.ANY_TAG, status=status)

        if (tag := status.Get_tag()) in (MACE_TAG_EXIT, MACE_TAG_ERROR):
            if tag == MACE_TAG_EXIT:
                logging.info("Received graceful exit signal from ABIN")
                exit_code = 0
            else:
                logging.warning("Received ERROR signal from ABIN. Stopping server")
                exit_code = 1

            try:
                abin_comm.Recv([MPI.BOTTOM, MPI.INT], source=0, tag=tag)
            except Exception as e:
                logging.error(f"Error receiving tag payload: {e}")

            shutdown_communication()
            sys.exit(exit_code)

    # Call error_shutdown upon any unhandled exception
    sys.excepthook = functools.partial(exception_handler, error_shutdown)

    # Receive number of atoms
    check_incoming_msg()
    natom_buf = np.empty(1, dtype=np.intc)
    abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
    natom = int(natom_buf[0])
    logging.info(f"Received number of atoms: {natom}")

    # Receive atom types
    check_incoming_msg()
    byte_buf = bytearray(natom * 2)
    mpi_status = MPI.Status()
    abin_comm.Recv([byte_buf, MPI.CHAR], source=0, tag=MACE_TAG_DATA, status=mpi_status)
    assert mpi_status.Get_count() == natom * 2
    # TODO: Check that size of received data!
    atom_types_str = byte_buf.decode("ascii")
    atom_types = [
        atom_types_str[i : i + 2].strip() for i in range(0, len(atom_types_str), 2)
    ]
    logging.info(f"Received atom types: {atom_types}")
    assert len(atom_types) == natom

    # Load MACE model
    if config.model_path.startswith("__MOCK_"):
        # This is for testing purposes only
        mace_model = HarmonicModel(config)
    else:
        mace_model = MaceModel(config)

    logging.info("MACE model ready. Entering main loop.")

    # Main loop: receive coordinates, compute, send results
    eval_count = 0
    while True:
        start_loop = perf_counter()
        check_incoming_msg()

        # Receive coordinates (3*natom doubles, in Bohr)
        coords = np.empty((natom, 3), dtype=np.float64)
        abin_comm.Recv(
            [coords, MPI.DOUBLE], source=0, tag=MACE_TAG_DATA, status=mpi_status
        )
        assert mpi_status.Get_elements(MPI.DOUBLE) == natom * 3

        start = perf_counter()

        energy, forces = mace_model.evaluate(atom_types, coords)

        end = perf_counter()
        time_ms = (end - start) * 1000
        logging.info(f"Step {eval_count} done in {time_ms:.3f} miliseconds")
        logging.info(f"Energy = {energy:.15f} Hartree")

        # Send energy (1 double, in Hartree)
        energy_buf = np.array([energy], dtype=np.float64)
        abin_comm.Send([energy_buf, MPI.DOUBLE], dest=0, tag=MACE_TAG_DATA)

        # Send forces (3*natom doubles, in Hartree/Bohr)
        # Transpose back to (3, natom) to match Fortran column-major layout
        if forces.dtype != np.float64:
            forces_send = forces.T.astype(np.float64)
        else:
            forces_send = forces.T.copy()
        abin_comm.Send([forces_send, MPI.DOUBLE], dest=0, tag=MACE_TAG_DATA)

        end_loop = perf_counter()
        loop_ms = (end_loop - start_loop) * 1000
        # Note: Communication overhead includes the ABIN propagation time,
        # which should however be negligible.
        logging.info(f"Communication overhead = {loop_ms - time_ms:.3f} ms")

        eval_count += 1


if __name__ == "__main__":
    config = parse_cmd()
    main(config)
