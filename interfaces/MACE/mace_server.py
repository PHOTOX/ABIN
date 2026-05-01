######################
# MACE MPI SERVER
######################
# Version: 6.0.0
# Evaluation logic based on MACE/eval_config [mace.cli.eval_config]
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

import sys
import time
from io import StringIO

import numpy as np
from mpi4py import MPI

# These imports are deferred until after config is received
# so the server can start quickly and report port
import ase.data

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
        pair = pair.strip()
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
        import torch
        from mace.tools import torch_tools, utils
        from mace import data as mace_data

        self.mace_data = mace_data
        self.torch_tools = torch_tools
        self.utils = utils
        self.torch = torch

        log("initializing MaceModel...")
        start_time = time.time()

        model_path = config.get('model', '')
        if not model_path:
            raise RuntimeError("MACE model path not specified")

        device = config.get('device', 'cpu')
        default_dtype = config.get('default_dtype', 'float64')
        self.batch_size = int(config.get('batch_size', '64'))
        self.compute_stress = config.get('compute_stress', 'false').lower() == 'true'
        self.return_contributions = config.get('return_contributions', 'false').lower() == 'true'
        self.info_prefix = config.get('info_prefix', 'MACE_')
        head_val = config.get('head', '')
        self.head = head_val if head_val else None

        torch_tools.set_default_dtype(default_dtype)
        self.device = str(torch_tools.init_device(device))

        import os
        if os.path.isfile(model_path):
            # Load from local file
            try:
                self.model = torch.jit.load(f=model_path, map_location=device).to(device)
            except Exception:
                log("Failed to load as TorchScript, trying as regular PyTorch model...")
                self.model = torch.load(f=model_path, map_location=device, weights_only=False).to(device)
        else:
            # Try loading as a MACE foundation model (auto-downloads)
            self.model = self._load_foundation_model(model_path, device)

        for param in self.model.parameters():
            param.requires_grad = False

        log("model loaded in %.2f s" % (time.time() - start_time))

    def _load_foundation_model(self, model_name, device):
        """
        Load a MACE foundation model by name.
        Supported names include:
          - MACE-OFF23_medium.model, MACE-OFF23_small.model, MACE-OFF23_large.model
          - medium, small, large (shorthand for MACE-OFF23)
          - MACE-MP-0_medium.model, etc.
        """
        from mace.calculators.foundations_models import mace_off, mace_mp

        name = model_name.replace('.model', '').strip()
        log(f"Loading foundation model: {name}")

        # Parse model family and size
        name_lower = name.lower()
        if name_lower in ('small', 'medium', 'large'):
            size = name_lower
            calc = mace_off(model=size, device=device, return_raw_model=True)
        elif 'mace-off23' in name_lower or 'mace_off23' in name_lower:
            size = name.split('_')[-1].lower()
            if size not in ('small', 'medium', 'large'):
                size = 'medium'
            calc = mace_off(model=size, device=device, return_raw_model=True)
        elif 'mace-mp' in name_lower or 'mace_mp' in name_lower:
            size = name.split('_')[-1].lower()
            if size not in ('small', 'medium', 'large'):
                size = 'medium'
            calc = mace_mp(model=size, device=device, return_raw_model=True)
        else:
            raise RuntimeError(
                f"Unknown model '{model_name}'. Provide a path to a local .model file, "
                f"or use a foundation model name like 'MACE-OFF23_medium.model' or 'medium'."
            )

        return calc

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
        from mace.tools import torch_geometric

        # Unit conversions
        bohr_to_ang = 0.529177249
        ev_to_hartree = 1.0 / 27.211399
        ev_per_ang_to_hartree_per_bohr = ev_to_hartree * bohr_to_ang

        # Convert coordinates from Bohr to Angstrom
        coords_ang = coords_bohr * bohr_to_ang

        # Create ASE atoms object
        atoms = ase.Atoms(symbols=atom_types, positions=coords_ang)

        if self.head is not None:
            atoms.info["head"] = self.head

        configs = [self.mace_data.config_from_atoms(atoms)]

        # Prepare dataset
        z_table = self.utils.AtomicNumberTable([
            int(z) for z in self.model.atomic_numbers
        ])

        try:
            heads = self.model.heads
        except AttributeError:
            heads = None

        data_loader = torch_geometric.dataloader.DataLoader(dataset=[
            self.mace_data.AtomicData.from_config(
                config,
                z_table=z_table,
                cutoff=float(self.model.r_max),
                heads=heads
            )
            for config in configs
        ], batch_size=self.batch_size, shuffle=False, drop_last=False)

        # Evaluate
        for batch in data_loader:
            batch = batch.to(self.device)
            output = self.model(
                batch.to_dict(),
                compute_stress=self.compute_stress
            )

            energy_ev = self.torch_tools.to_numpy(output["energy"])[0]
            forces_ev_ang = np.split(
                self.torch_tools.to_numpy(output["forces"]),
                indices_or_sections=batch.ptr[1:],
                axis=0,
            )[0]

        # Convert to atomic units
        energy_hartree = energy_ev * ev_to_hartree
        forces_hartree_bohr = forces_ev_ang * ev_per_ang_to_hartree_per_bohr

        return energy_hartree, forces_hartree_bohr


def main():
    comm = MPI.COMM_WORLD

    # Open MPI port and write to file for ABIN to read
    port_name = MPI.Open_port()
    log(f"MPI port opened: {port_name}")

    port_file = "mace_port.txt.1"
    with open(port_file, "w") as f:
        f.write(port_name)
    log(f"Port written to {port_file}")

    # Accept connection from ABIN
    log("Waiting for ABIN to connect...")
    abin_comm = comm.Accept(port_name)
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
            abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_EXIT)
            break

        abin_comm.Recv([natom_buf, MPI.INT], source=0, tag=MACE_TAG_DATA)
        natom_step = int(natom_buf[0])

        if natom_step != natom:
            log(f"ERROR: Received natom={natom_step}, expected {natom}")
            break

        # Receive coordinates (3*natom doubles, in Bohr)
        coords = np.empty((3, natom), dtype=np.float64)
        abin_comm.Recv([coords, MPI.DOUBLE], source=0, tag=MACE_TAG_DATA)
        # Transpose to (natom, 3) for ASE
        coords_bohr = coords.T.copy()

        eval_count += 1
        log(f"Evaluation {eval_count}: evaluating...")

        try:
            energy, forces = mace_model.evaluate(atom_types, coords_bohr)

            # Send energy (1 double, in Hartree)
            energy_buf = np.array([energy], dtype=np.float64)
            abin_comm.Send([energy_buf, MPI.DOUBLE], dest=0, tag=0)

            # Send forces (3*natom doubles, in Hartree/Bohr)
            # Transpose back to (3, natom) to match Fortran column-major layout
            forces_send = forces.T.copy()
            abin_comm.Send([forces_send, MPI.DOUBLE], dest=0, tag=0)

            log(f"Evaluation {eval_count}: energy = {energy:.10f} Hartree")

        except Exception as e:
            log(f"ERROR during evaluation: {e}")
            # Send error tag
            error_energy = np.array([0.0], dtype=np.float64)
            abin_comm.Send([error_energy, MPI.DOUBLE], dest=0, tag=1)
            break

    # Cleanup
    log("Shutting down...")
    abin_comm.Disconnect()
    MPI.Close_port(port_name)
    log("Server stopped.")


if __name__ == "__main__":
    main()
