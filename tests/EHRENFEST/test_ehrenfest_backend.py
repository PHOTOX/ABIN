#!/usr/bin/env python3
"""Exercise the independent all-state shell contract without a QM installation."""
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
assert list(map(int, input().split())) == [1] * n
assert Path('eh_state.dat').exists() and not Path('state.dat').exists()
with open('gradient_calls', 'a') as f:
    f.write(str(n) + '\\n')
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
assert list(map(int, input().split())) == [1] * (n * (n - 1) // 2)
with open('nac_calls', 'a') as f:
    f.write(sys.argv[3] + '\\n')
if Path('fail_nac').exists() or (Path('retry_nac').exists() and sys.argv[3] == '7'):
    Path('nacm.dat').write_text('0 0 0\\n')
else:
    Path('nacm.dat').write_text('0 0 0\\n' * (n * (n - 1)))
'''


def run_case(nstates, initial_state, mode='normal'):
    with tempfile.TemporaryDirectory(prefix='abin-eh-backend-') as folder:
        work = Path(folder)
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
        if mode in ('truncate_gradients', 'fail_nac'):
            assert result.returncode != 0, 'Malformed backend output was accepted'
            assert 'Invalid Ehrenfest' in result.stderr or 'Invalid or missing Ehrenfest NAC' in result.stderr
            return
        assert result.returncode == 0, result.stdout[-2000:] + result.stderr
        calls = (work / 'gradient_calls').read_text().splitlines()
        assert calls == [str(nstates)] * 2, calls
        if nstates == 1:
            assert not (work / 'nac_calls').exists()
        else:
            accuracies = (work / 'nac_calls').read_text().splitlines()
            assert accuracies == (['7', '5'] * 2 if mode == 'retry_nac' else ['7'] * 2), accuracies
        population = [float(v) for v in (work / 'pop.dat').read_text().splitlines()[-1].split()]
        expected = [0.0] * nstates
        expected[initial_state - 1] = 1.0
        assert population[1:-1] == expected and abs(population[-1] - 1) < 1e-12
        pes = [float(v) for v in (work / 'PES.dat').read_text().splitlines()[-1].split()]
        assert all(abs(e - (i + 1) * .01) < 1e-12 for i, e in enumerate(pes[1:]))
        forces = (work / 'forces.xyz').read_text().splitlines()
        fx1 = float(forces[2].split()[1])
        fx2 = float(forces[3].split()[1])
        assert abs(fx1 + initial_state * .0001) < 1e-12
        assert abs(fx2 - initial_state * .0001) < 1e-12


for state in (1, 2, 3):
    run_case(3, state)
run_case(1, 1)
run_case(3, 2, 'retry_nac')
run_case(3, 2, 'truncate_gradients')
run_case(3, 2, 'fail_nac')
print('Independent Ehrenfest shell backend: all states, gradient ordering, NAC retries and invalid-output checks passed.')
