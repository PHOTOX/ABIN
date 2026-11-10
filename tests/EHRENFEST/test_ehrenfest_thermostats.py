#!/usr/bin/env python3
"""Check shared nuclear thermostats and the Ehrenfest drift-velocity coupling."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / 'bin/abin'

# Two degenerate, flat surfaces with constant derivative coupling. Nuclear
# motion must match classical MD, while the electronic rotation is determined
# exactly by the actual displacement, including all thermostat impulses.
GRADIENT_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
import sys
request = Path('eh_state.dat')
n = int(request.read_text().splitlines()[0]) if request.exists() else 1
Path('engrad.dat.' + sys.argv[2]).write_text('0\\n' * n + '0 0 0\\n' * (2 * n))
'''
NAC_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
Path('nacm.dat').write_text('0.4 0 0\\n-0.4 0 0\\n')
'''


def table(work, name):
    rows = [[float(x) for x in line.split()]
            for line in (work / name).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert rows and all(math.isfinite(x) for row in rows for x in row), name
    return rows


def run_case(mdtype, thermostat, section='thermostat', pot='mock', nstep=12):
    with tempfile.TemporaryDirectory(prefix='abin-eh-thermostat-') as folder:
        work = Path(folder)
        adapter = work / 'MOCK'
        adapter.mkdir()
        for name, script in [('r.mock', GRADIENT_SCRIPT), ('r.mock.nacm', NAC_SCRIPT)]:
            path = adapter / name
            path.write_text(script)
            path.chmod(0o755)
        # Nonzero auxiliary coupling exercises the colored-noise GLE path.
        # A is in atomic units and C is in eV, as required by gle_init.
        (work / 'GLE-A').write_text('1\n0.01 0.003\n-0.003 0.02\n')
        (work / 'GLE-C').write_text('1\n0.02585 0\n0 0.02585\n')
        separation = 4.0 if pot == 'mock' else 13.1
        (work / 'mini.xyz').write_text(f'2\nBohr\nNa 0 0 0\nI {separation} 0 0\n')
        (work / 'velocities.in').write_text('2\nAtomic units\nNa 0.0002 0.0001 0\nI -0.0001 0 0.0001\n')
        (work / 'input.in').write_text(f'''&general
mdtype='{mdtype}', pot='{pot}', irest=0, irandom=123,
nstep={nstep}, dt=5.0, nwrite=1, nwritex=1, nwritev=1,
xyz_units='bohr',
/
&{section}
{thermostat}
temp=300.0, rem_comvel=.false., rem_comrot=.false., scaleveloc=0,
/
&eh
nstate=2, istate_init=1, substep=100,
/
&system
/
''')
        result = subprocess.run([str(BINARY), '-x', 'mini.xyz', '-v', 'velocities.in'],
                                cwd=work, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr
        restart = (work / 'restart.xyz').read_text().splitlines()
        coords = [[float(x) for x in row.split()] for row in restart[2:4]]
        velocities = [[float(x) for x in row.split()] for row in restart[5:7]]
        temperatures = table(work, 'temper.dat')
        energies = table(work, 'energies.dat')
        assert len(temperatures) == len(energies) == nstep
        if mdtype == 'ehrenfest':
            populations = table(work, 'pop.dat')
            assert len(populations) == nstep
            assert all(abs(row[-1] - 1) < 1e-9 for row in populations)
            if pot == 'mock':
                angle = 0.4 * (coords[0][0] - (coords[1][0] - separation))
                assert abs(populations[-1][2] - math.sin(angle)**2) < 1e-9, (
                    thermostat, populations[-1], angle)
                assert populations[-1][2] > 1e-6  # Exercise electronic evolution.
        return coords + velocities, temperatures, energies


cases = [
    ('NVE', 'inose=0,', 'thermostat'),
    ('global NHC', "therm='nhc', imasst=0, tau0=0.002,", 'thermostat'),
    ('massive NHC', 'inose=1, imasst=1, tau0=0.002,', 'nhcopt'),
    ('GLE', "therm='gle', gle_test=.false.,", 'thermostat'),
    ('canonical GLE', 'inose=4, gle_test=.false.,', 'nhcopt'),
    ('Langevin', 'inose=3, tau0_langevin=0.002,', 'thermostat'),
]

baseline = None
for name, thermostat, section in cases:
    eh, temperatures, energies = run_case('ehrenfest', thermostat, section)
    classical, classical_temperatures, classical_energies = run_case('md', thermostat, section)
    assert max(abs(x - y) for erow, crow in zip(eh, classical)
               for x, y in zip(erow, crow)) < 1e-12, name
    assert temperatures == classical_temperatures, name
    assert [row[:4] for row in energies] == [row[:4] for row in classical_energies], name
    if baseline is None:
        baseline = eh
        assert all(len(row) == 3 for row in temperatures)
    else:
        assert all(len(row) == 4 for row in temperatures), name
        assert max(abs(row[-1] - temperatures[0][-1]) for row in temperatures) < 1e-8, name
        assert max(abs(row[2] - energies[0][2]) for row in energies) > 1e-7, name
        assert max(abs(x - y) for erow, brow in zip(eh, baseline)
                   for x, y in zip(erow, brow)) > 1e-7, name

# Also run the physical two-state NaI model with nonzero energy gaps and forces.
for name, thermostat, section in cases[1:5]:
    run_case('ehrenfest', thermostat, section, pot='_nai_', nstep=20)

print('Ehrenfest thermostats passed: global/massive NHC, GLE, canonical GLE, Langevin, '
      'classical trajectories, heat accounting, electronic drift velocities and NaI.')
