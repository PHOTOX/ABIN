######################
# MACE MPI SERVER
######################
# Version: 6.0.0
#
# This server communicates with ABIN via MPI (using mpi4py).
# MACE configuration (model path, device, etc.) is received from ABIN's
# &mace namelist section at initialization time.
#
# Dependencies:
#   pip install mpi4py mace-torch torch ase numpy
#
# Usage:
#   mpirun -n 1 python mace_server.py
#
# The server writes its MPI port to 'mace_port.txt.1' for ABIN to read.

import os
import time

LOG_NAME = "MaceMPIServer"

# MPI Tags (must match Fortran module mod_mace_mpi)
MACE_TAG_EXIT = 666
MACE_TAG_DATA = 2


def log(message, should_print=True):
    msg_formatted = "[{}]: {}".format(LOG_NAME, str(message))
    if should_print:
        print(msg_formatted, flush=True)
    return msg_formatted


def parse_config_string(config_str):
    """Parse key=value config string sent from ABIN's &mace namelist."""
    config = {}
    for pair in config_str.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            config[key.strip()] = value.strip()
    return config


class MaceModel:
    """
    Manages the MACE ML model for evaluating atomic configurations.
    Configuration is received from ABIN via MPI.
    """

    def __init__(self, config):
        from mace.calculators import MACECalculator

        log("initializing MACE model")
        model_path = config.get('model', '')
        if not model_path:
            raise RuntimeError("MACE model path not specified")

        device = config.get('device', 'cpu')
        head_val = config.get('head', '')
        self.head = head_val if head_val else None

        if not os.path.isfile(model_path):
            raise RuntimeError(f"file '{model_path}' not found")

        # Set ASE calculator
        self.calculator = MACECalculator(
            model_paths=model_path,
            device=device,
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


def main():
    import numpy as np
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

    # Receive configuration string
    config_status = MPI.Status()
    abin_comm.Probe(source=0, tag=MACE_TAG_DATA, status=config_status)
    config_len = config_status.Get_count(MPI.CHAR)
    config_buf = bytearray(config_len)
    abin_comm.Recv([config_buf, MPI.CHAR], source=0, tag=MACE_TAG_DATA)
    config_str = config_buf.decode('ascii').strip()
    log(f"Received config: {config_str}")

    config = parse_config_string(config_str)
    default_dtype = config.get('default_dtype', 'float64')

    # Load MACE model
    mace_model = MaceModel(config)
    log("MACE model ready. Entering main loop.")

    def shutdown_server():
        log("Shutting down...")
        abin_comm.Disconnect()
        MPI.Close_port(port_name)
        log("Server stopped.")


    # Main loop: receive coordinates, compute, send results
    eval_count = 0
    while True:
        # Receive natom (sent each step for protocol consistency)
        status = MPI.Status()
        abin_comm.Probe(source=0, tag=MPI.ANY_TAG, status=status)

        if status.Get_tag() == MACE_TAG_EXIT:
            log("Received exit signal from ABIN")
            # Consume the message
            abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_EXIT)
            break

        abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
        natom_step = int(natom_buf[0])

        if natom_step != natom:
            log(f"ERROR: Received natom={natom_step}, expected {natom}")
            break

        # Receive coordinates (3*natom doubles, in Bohr)
        coords = np.empty((natom, 3), dtype=np.float64)
        abin_comm.Recv([coords, MPI.DOUBLE], source=0, tag=MACE_TAG_DATA)

        coords_bohr = coords.copy()

        eval_count += 1
        log(f"Evaluation {eval_count}")

        try:
            energy, forces = mace_model.evaluate(atom_types, coords_bohr)
        except Exception as e:
            log(f"ERROR during evaluation: {e}")
            # Send error tag
            error_energy = np.array([0.0], dtype=np.float64)
            abin_comm.Send([error_energy, MPI.DOUBLE], dest=0, tag=1)
            shutdown_server()
            raise e
        else:
            log(f"Evaluation {eval_count}: energy = {energy:.15f} Hartree")
            # log(f"Evaluation {eval_count}: forces = \n{forces}")

        try:
            # Send energy (1 double, in Hartree)
            energy_buf = np.array([energy], dtype=np.float64)
            abin_comm.Send([energy_buf, MPI.DOUBLE], dest=0, tag=0)

            # Send forces (3*natom doubles, in Hartree/Bohr)
            # Transpose back to (3, natom) to match Fortran column-major layout
            if default_dtype != 'float64':
                forces_send = forces.T.astype(np.float64)
            else:
                forces_send = forces.T.copy()
            #log(f"Energy sent to ABIN ({eval_count}): {energy_buf}")
            #log(f"Forces sent to ABIN ({eval_count}): {forces_send!r}")
            abin_comm.Send([forces_send, MPI.DOUBLE], dest=0, tag=0)
        except Exception as e:
            log(f"ERROR when sending energy and forces to ABIN: {e}")
            error_energy = np.array([0.0], dtype=np.float64)
            abin_comm.Send([error_energy, MPI.DOUBLE], dest=0, tag=1)
            shutdown_server()
            raise e

    shutdown_server()


if __name__ == "__main__":
    main()
