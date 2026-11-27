#!/usr/bin/env python3
"""Exercise Ehrenfest input validation and the shell contract with synthetic exports."""
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / 'bin/abin'

GRADIENT_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
import sys
n = int(input())
# The gradient request must select every state and use the Ehrenfest state file.
assert list(map(int, input().split())) == [1] * n
assert Path('eh_state.dat').exists() and not Path('state.dat').exists()
with open('calls', 'a') as f:
    print(f'gradients {n}', file=f)
if Path('stale_gradients').exists() and sys.argv[1] == '1':
    sys.exit(0)
# Emit all energies first, then two Cartesian gradient rows per state.
rows = [str((i + 1) * 0.01) for i in range(n)]
for i in range(n):
    rows += [f'{(i + 1) * 0.0001} 0 0', f'{-(i + 1) * 0.0001} 0 0']
if Path('truncate_gradients').exists():
    rows = rows[:-1]
Path('engrad.dat.' + sys.argv[2]).write_text('\\n'.join(rows) + '\\n')
'''
NAC_SCRIPT = '''#!/usr/bin/env python3
from pathlib import Path
import sys
n = int(input())
# Expect one request flag per upper-triangle state pair and record the requested accuracy.
assert list(map(int, input().split())) == [1] * (n * (n - 1) // 2)
step, _, accuracy = sys.argv[1:]
with open('calls', 'a') as f:
    print(f'nac {accuracy}', file=f)
zero_nac = '0 0 0\\n' * (n * (n - 1))
if Path('stale_gradients').exists() and step == '0':
    # The initial gradient has been consumed. Leave a valid file for the next step.
    energies = ''.join(f'{(i + 1) * 0.01}\\n' for i in range(n))
    Path('engrad.dat.001').write_text(energies + '0 0 0\\n' * (2 * n))
if Path('stale_nac').exists():
    if accuracy == '7':
        Path('nacm.dat').write_text(zero_nac)
        sys.exit(1)
    # Report success on retry while producing no replacement.
    sys.exit(0)
if Path('coupled').exists() or Path('flip_nac').exists():
    assert n == 2
    sign = -1 if Path('flip_nac').exists() and step == '1' else 1
    with open('nac_signs', 'a') as f:
        print(sign, file=f)
    Path('nacm.dat').write_text(f'{sign * 0.4} {sign * 0.2} {-sign * 0.3}\\n'
                                f'{-sign * 0.4} {-sign * 0.2} {sign * 0.3}\\n')
    sys.exit(0)
# A one-row result forces a retry at accuracy 5 or fails both attempts.
if Path('fail_nac').exists() or (Path('retry_nac').exists() and accuracy == '7'):
    Path('nacm.dat').write_text('0 0 0\\n')
else:
    Path('nacm.dat').write_text(zero_nac)
'''


def run_case(nstates=3, initial_state=2, *, mode='normal', general='', system='',
             eh_input=None, fallback=5, diagnostic=None, launches=(), restart_at=None):
    """Run one shell-backend trajectory and check its outputs or expected failure."""
    with tempfile.TemporaryDirectory(prefix='abin-eh-backend-') as folder:
        work = Path(folder)
        # Give each run its own mock scripts, call logs, and trajectory outputs.
        (work / 'MOCK-EH').mkdir()
        for name, script in [('r.mock-eh', GRADIENT_SCRIPT), ('r.mock-eh.nacm', NAC_SCRIPT)]:
            executable = work / 'MOCK-EH' / name
            executable.write_text(script)
            executable.chmod(0o755)
        (work / 'mini.xyz').write_text('2\nTest\nI 0 0 0\nNa 4 0 0\n')
        coupled = mode in ('coupled', 'flip_nac')
        nstep = 3 if coupled else 1
        velocities = 'I 0.008 0.004 0.002\nNa -0.004 -0.003 0.001\n' if coupled else 'I 0 0 0\nNa 0 0 0\n'
        (work / 'velocities.in').write_text('2\nTest\n' + velocities)
        if eh_input is None:
            eh_input = f'''&eh
nstate={nstates}, istate_init={initial_state}, nac_accu1=7, nac_accu2={fallback},
/
'''
        input_text = f'''&general
mdtype='ehrenfest', pot='mock-eh', irandom=123,
irest={{irest}}, nstep={{nstep}}, dt=0.1, nwrite=1, nwritef=1, nwritex=1, nwritev=0,
{general}
/
&thermostat
inose=0, temp=0.0, rem_comvel=.false., rem_comrot=.false., scaleveloc=0,
/
{eh_input}&system
{system}
/
'''
        input_path = work / 'input.in'
        first_steps = nstep if restart_at is None else restart_at
        input_path.write_text(input_text.format(irest=0, nstep=first_steps))
        if mode != 'normal':
            (work / mode).touch()
        result = subprocess.run([str(BINARY), '-x', 'mini.xyz', '-v', 'velocities.in'],
                                cwd=work, capture_output=True, text=True, timeout=30)
        if restart_at is not None:
            assert result.returncode == 0, result.stdout + result.stderr
            input_path.write_text(input_text.format(irest=1, nstep=nstep))
            result = subprocess.run([str(BINARY), '-x', 'mini.xyz'], cwd=work,
                                    capture_output=True, text=True, timeout=30)
        calls = work / 'calls'
        actual_calls = calls.read_text().splitlines() if calls.exists() else []
        details = f'{mode}, {general}\n{result.stdout[-1500:]}{result.stderr}'
        if diagnostic is not None:
            assert result.returncode != 0, details
            assert (work / 'ERROR').is_file(), details
            assert (work / 'ERROR').read_text().strip() == diagnostic, details
            assert actual_calls == list(launches), (diagnostic, actual_calls)
            return
        assert result.returncode == 0, details
        # Initialization and each step request gradients followed by NACs.
        # The retry case must fall back from accuracy 7 to 5 at each geometry.
        per_geometry = [f'gradients {nstates}', 'nac 7']
        if mode == 'retry_nac':
            per_geometry.append(f'nac {fallback}')
        assert actual_calls == per_geometry * (nstep + 1), actual_calls
        if coupled:
            signs = (work / 'nac_signs').read_text().splitlines()
            assert signs == (['1', '-1', '1', '1'] if mode == 'flip_nac' else ['1'] * 4), signs
            # Gradients are small, so the prescribed drift gives the NAC projection to this precision.
            projection = 0.4 * 0.012 + 0.2 * 0.007 - 0.3 * 0.001
            products = [row for row in (work / 'dotprod.dat').read_text().splitlines()
                        if row.strip() and not row.lstrip().startswith('#')]
            assert len(products) == nstep, products
            assert all(abs(float(row.split()[1]) - projection) < 1e-6 for row in products), (
                'Unexpected NAC projection', mode, products)
            final_populations = (work / 'pop.dat').read_text().splitlines()[-1].split()
            assert float(final_populations[3]) > 1e-8, final_populations
            return {name: (work / name).read_text()
                    for name in ('dotprod.dat', 'pop.dat', 'forces.xyz', 'PES.dat', 'energies.dat', 'movie.xyz')}
        # Zero NACs keep the chosen pure state populated throughout the trajectory.
        population = [float(v) for v in (work / 'pop.dat').read_text().splitlines()[-1].split()]
        expected = [0.0] * nstates
        expected[initial_state - 1] = 1.0
        assert population[1] == initial_state
        assert population[2:-1] == expected and abs(population[-1] - 1) < 1e-12
        dotproducts = [float(v) for v in (work / 'dotprod.dat').read_text().splitlines()[-1].split()]
        assert dotproducts == [population[0]] + [0.0] * (nstates * (nstates - 1) // 2)
        pes = [float(v) for v in (work / 'PES.dat').read_text().splitlines()[-1].split()]
        assert all(abs(e - (i + 1) * .01) < 1e-12 for i, e in enumerate(pes[1:]))
        # The two forces must be the negatives of the selected state's gradients.
        forces = (work / 'forces.xyz').read_text().splitlines()
        fx1 = float(forces[2].split()[1])
        fx2 = float(forces[3].split()[1])
        assert abs(fx1 + initial_state * .0001) < 1e-12
        assert abs(fx2 - initial_state * .0001) < 1e-12


# Select every gradient block in both the two-state and three-state layouts.
for state in (1, 2, 3):
    run_case(3, state)
for state in (1, 2):
    run_case(2, state)

# Each invalid input must fail at the Ehrenfest check before calling chemistry.
eh_error = 'ERROR in ehrenfest.F90: '
for options, message in (
    (dict(eh_input=''), 'Missing or invalid &eh namelist'),
    (dict(eh_input='&eh\nnstate=invalid,\n/\n'), 'Missing or invalid &eh namelist'),
    (dict(general='nwalk=2,'), 'nwalk must be 1 for Ehrenfest dynamics'),
    (dict(general='iqmmm=1,'), 'iqmmm must be 0 for Ehrenfest dynamics'),
    (dict(system='nshake=1,'), 'nshake must be 0 for Ehrenfest dynamics'),
    (dict(general='en_restraint=1,'), 'en_restraint must be 0.0 for Ehrenfest dynamics'),
    (dict(general="pot_ref='_morse_',"), 'pot_ref must be _none_ for Ehrenfest dynamics'),
    (dict(general='conatom=1,'), 'conatom must be 0 for Ehrenfest dynamics'),
    (dict(general='iplumed=1,'), 'iplumed must be 0 for Ehrenfest dynamics'),
    (dict(general='isbc=1,'), 'isbc must be 0 for Ehrenfest dynamics'),
    (dict(nstates=1, initial_state=1), 'nstate must be >= 2 for Ehrenfest dynamics since for'
                                      'nstate=1the Ehrenfest dynamics reduces to classical dynamics'),
):
    run_case(**options, diagnostic=eh_error + message)

# Check both retry outcomes and the deletion of files left during the trajectory.
backend_error = 'ERROR in eh_backend.F90: '
nac_error = backend_error + 'Invalid or missing Ehrenfest NAC vectors'
run_case(mode='retry_nac')
run_case(mode='truncate_gradients', diagnostic=backend_error + 'Invalid Ehrenfest energies or state gradients',
         launches=('gradients 3',))
run_case(mode='fail_nac', diagnostic=nac_error, launches=('gradients 3', 'nac 7', 'nac 5'))
run_case(mode='fail_nac', fallback=7, diagnostic=nac_error, launches=('gradients 3', 'nac 7'))
run_case(mode='stale_gradients', diagnostic=backend_error + 'Missing Ehrenfest energies and gradients',
         launches=('gradients 3', 'nac 7', 'gradients 3'))
run_case(mode='stale_nac', diagnostic=nac_error, launches=('gradients 3', 'nac 7', 'nac 5'))

baseline = run_case(2, 1, mode='coupled')
reversed_nac = run_case(2, 1, mode='flip_nac')
restarted = run_case(2, 1, mode='flip_nac', restart_at=1)
for name in baseline:
    assert baseline[name] == reversed_nac[name], f'NAC sign reversal changed {name}'
    assert reversed_nac[name] == restarted[name], f'Ehrenfest restart changed {name}'
print('Ehrenfest shell backend, input validation, NAC retries, sign continuity, and restart passed.')
