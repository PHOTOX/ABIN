#!/usr/bin/env python3
"""Check NVE drift against classical MD and electronic rotation against an analytic solution.

The NaI crossing also checks population transfer and bounds total-energy drift.
Coverage is limited to the chosen initial conditions and NVE timesteps.
"""
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / 'bin/abin'

# Two degenerate, flat surfaces with constant derivative coupling. Nuclear
# NVE motion must match classical MD, while the electronic rotation is
# determined exactly by the actual displacement.
GRADIENT_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
import sys
# Ehrenfest requests every state. Classical MD needs one energy and gradient.
request = Path('eh_state.dat')
n = int(request.read_text().splitlines()[0]) if request.exists() else 1
Path('engrad.dat.' + sys.argv[2]).write_text('0\\n' * n + '0 0 0\\n' * (2 * n))
'''
NAC_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
n = int(Path('eh_state.dat').read_text().splitlines()[0])
# Write two opposite atomic vectors for each upper-triangle state pair.
Path('nacm.dat').write_text(''.join(f'{0.4 * k} 0 0\\n{-0.4 * k} 0 0\\n'
                                 for k in range(1, n * (n - 1) // 2 + 1)))
'''


def table(work, name):
    """Read a nonempty numeric output table and require all values to be finite."""
    rows = [[float(x) for x in line.split()]
            for line in (work / name).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert rows and all(math.isfinite(x) for row in rows for x in row), name
    return rows


def run_case(mdtype, pot='mock', nstep=12, crossing=False, nstate=2, nwrite=1):
    """Run one NVE trajectory and return its final positions, velocities, and diagnostics."""
    with tempfile.TemporaryDirectory(prefix='abin-eh-nve-') as folder:
        work = Path(folder)
        # Install the mock adapter at the path ABIN uses for this potential.
        adapter = work / 'MOCK'
        adapter.mkdir()
        for name, script in [('r.mock', GRADIENT_SCRIPT), ('r.mock.nacm', NAC_SCRIPT)]:
            path = adapter / name
            path.write_text(script)
            path.chmod(0o755)
        separation = 4.0 if pot == 'mock' else 13.1
        if crossing:
            separation = 6.0 * 1.889726132873  # Start inside the NaI avoided crossing (bohr).
        (work / 'mini.xyz').write_text(f'2\nBohr\nNa 0 0 0\nI {separation} 0 0\n')
        # Give the flat surfaces steady drift and send NaI outward through the crossing.
        velocities = 'Na 0.0002 0.0001 0\nI -0.0001 0 0.0001\n'
        if crossing:
            velocities = 'Na -0.0017 0 0\nI 0.0003 0 0\n'
        (work / 'velocities.in').write_text('2\nAtomic units\n' + velocities)
        # Disable the thermostat and velocity adjustments to preserve the supplied velocities.
        (work / 'input.in').write_text(f'''&general
mdtype='{mdtype}', pot='{pot}', irest=0, irandom=123,
nstep={nstep}, dt=5.0, nwrite={nwrite}, nwritex=1, nwritev=1,
xyz_units='bohr',
/
&thermostat
inose=0,
temp=300.0, rem_comvel=.false., rem_comrot=.false., scaleveloc=0,
/
&eh
nstate={nstate}, istate_init={2 if crossing else 1}, substep=100,
/
&system
/''')
        result = subprocess.run([str(BINARY), '-x', 'mini.xyz', '-v', 'velocities.in'],
                                cwd=work, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr
        # Read the final nuclear state separately from the sampled diagnostics.
        restart = (work / 'restart.xyz').read_text().splitlines()
        coords = [[float(x) for x in row.split()] for row in restart[2:4]]
        velocities = [[float(x) for x in row.split()] for row in restart[5:7]]
        temperatures = table(work, 'temper.dat')
        energies = table(work, 'energies.dat')
        assert len(temperatures) == len(energies) == nstep // nwrite
        if mdtype == 'ehrenfest':
            populations = table(work, 'pop.dat')
            dotproducts = table(work, 'dotprod.dat')
            pes = table(work, 'PES.dat')
            # Each write must include synchronized electronic diagnostics and a unit norm.
            assert len(populations) == len(dotproducts) == len(pes) == nstep // nwrite
            assert all(abs(row[-1] - 1) < 1e-9 for row in populations)
            for row, coupling, energy in zip(populations, dotproducts, pes):
                assert len(row) == nstate + 3
                # The state label must identify the largest reported population.
                assert row[1] == 1 + max(range(nstate), key=lambda i: row[2 + i])
                assert row[0] == coupling[0] == energy[0]
                assert len(coupling) == 1 + nstate * (nstate - 1) // 2
                if pot == 'mock':
                    # Constant drift velocity and distinct couplings check signs,
                    # atomic units, and upper-triangle pair ordering.
                    assert all(abs(value - 0.4 * k * 0.0003) < 1e-12
                               for k, value in enumerate(coupling[1:], 1))
            if pot == 'mock' and nstate == 2:
                # Integrate the constant coupling over relative displacement
                # to obtain the electronic rotation angle.
                angle = 0.4 * (coords[0][0] - (coords[1][0] - separation))
                assert abs(populations[-1][3] - math.sin(angle)**2) < 1e-9, (
                    populations[-1], angle)
                assert populations[-1][3] > 1e-6  # Exercise electronic evolution.
            if crossing:
                # A successful run and conserved norm alone also pass with a
                # frozen wavefunction. Require transfer through the NaI crossing.
                final_distance = math.dist(coords[0], coords[1])
                assert final_distance > 7.5 * 1.889726132873, final_distance
                assert max(row[2] for row in populations) > 0.15, populations[-1]
                assert populations[-1][2] > 0.05, populations[-1]
                assert max(row[3] for row in energies) - min(row[3] for row in energies) < 1e-7
        return coords + velocities, temperatures, energies


# Compare the final nuclear state and sampled NVE diagnostics on the flat surfaces.
classical, temperatures, energies = run_case('md')
eh, eh_temperatures, eh_energies = run_case('ehrenfest')
assert max(abs(x - y) for erow, crow in zip(eh, classical)
           for x, y in zip(erow, crow)) < 1e-12
assert eh_temperatures == temperatures
assert [row[:4] for row in eh_energies] == [row[:4] for row in energies]
assert all(len(row) == 3 for row in temperatures)

# Run the physical two-state NaI model in the supported NVE ensemble too.
run_case('ehrenfest', pot='_nai_', nstep=20)
run_case('ehrenfest', pot='_nai_', nstep=400, crossing=True)
# Exercise all three state pairs with a write interval longer than one timestep.
run_case('ehrenfest', nstate=3, nstep=12, nwrite=3)

print('Ehrenfest NVE trajectories, electronic drift velocities and NaI population transfer passed.')
