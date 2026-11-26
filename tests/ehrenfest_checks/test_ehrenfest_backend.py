#!/usr/bin/env python3
"""Check the all-state shell contract with synthetic energies and gradients.

Zero NACs isolate the selected state's gradient. Separate cases check stale
or truncated output, retry and retry exhaustion, and single-state rejection.
Coverage here is limited to ABIN's shell protocol and pure-state force selection.
The real-adapter matrix covers program-specific input and export handling.
"""
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
with open('gradient_calls', 'a') as f:
    f.write(str(n) + '\\n')
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
with open('nac_calls', 'a') as f:
    f.write(sys.argv[3] + '\\n')
# A one-row result forces a retry at accuracy 5 or fails both attempts.
if Path('fail_nac').exists() or (Path('retry_nac').exists() and sys.argv[3] == '7'):
    Path('nacm.dat').write_text('0 0 0\\n')
else:
    Path('nacm.dat').write_text('0 0 0\\n' * (n * (n - 1)))
'''


def run_case(nstates, initial_state, mode='normal'):
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
        (work / 'velocities.in').write_text('2\nTest\nI 0 0 0\nNa 0 0 0\n')
        (work / 'input.in').write_text(f'''&general
mdtype='ehrenfest', pot='mock-eh', irest=0, irandom=123,
nstep=1, dt=0.1, nwrite=1, nwritef=1, nwritex=1, nwritev=0,
/
&nhcopt
inose=0, temp=0.0,
/
&eh
nstate={nstates}, istate_init={initial_state}, nac_accu1=7, nac_accu2=5,
/
&system
/
''')
        if mode != 'normal':
            (work / mode).touch()
        # A stale complete result must never mask truncated new output.
        if mode == 'truncate_gradients':
            (work / 'engrad.dat.001').write_text('0\n' * 100)
        result = subprocess.run([str(BINARY), '-x', 'mini.xyz', '-v', 'velocities.in'],
                                cwd=work, capture_output=True, text=True, timeout=30)
        # Single-state input must fail before either adapter is called.
        if nstates == 1:
            assert result.returncode != 0, 'Single-state Ehrenfest input was accepted'
            assert 'nstate must be >= 2' in result.stderr, result.stderr
            assert not (work / 'gradient_calls').exists()
            assert not (work / 'nac_calls').exists()
            return
        if mode in ('truncate_gradients', 'fail_nac'):
            assert result.returncode != 0, 'Malformed backend output was accepted'
            assert 'Invalid Ehrenfest' in result.stderr or 'Invalid or missing Ehrenfest NAC' in result.stderr
            return
        assert result.returncode == 0, result.stdout[-2000:] + result.stderr
        # Initialization and the first step each request gradients and NACs.
        # The retry case must fall back from accuracy 7 to 5 at both geometries.
        calls = (work / 'gradient_calls').read_text().splitlines()
        assert calls == [str(nstates)] * 2, calls
        accuracies = (work / 'nac_calls').read_text().splitlines()
        assert accuracies == (['7', '5'] * 2 if mode == 'retry_nac' else ['7'] * 2), accuracies
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
# Exercise input rejection, recovery after a failed NAC, and malformed output.
run_case(1, 1)
run_case(3, 2, 'retry_nac')
run_case(3, 2, 'truncate_gradients')
run_case(3, 2, 'fail_nac')
print('Independent Ehrenfest shell backend: all states, gradient ordering, NAC retries, '
      'invalid-output checks and single-state rejection passed.')
