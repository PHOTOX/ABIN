#!/usr/bin/env python3
"""Test the real BAGEL-EH scripts using deterministic BAGEL-format exports."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / 'bin/abin'

FAKE_BAGEL = '''import json
import os
from pathlib import Path
import sys

document = json.loads(Path(sys.argv[1]).read_text())
blocks = document['bagel']
forces = next(block for block in blocks if block['title'] == 'forces')
assert forces['export'] is True and forces['export_single'] is False
assert blocks[0]['angstrom'] is True
method = forces['method'][0]
nstates = method['nstate']
natoms = len(blocks[0]['geometry'])
gradients = [grad for grad in forces['grads'] if grad['title'] == 'force']
nacs = [grad for grad in forces['grads'] if grad['title'] == 'nacme']
assert [grad['target'] for grad in gradients] == list(range(nstates))
assert [(grad['target'], grad['target2']) for grad in nacs] == [
    (i, j) for i in range(nstates - 1) for j in range(i + 1, nstates)]
assert all(grad['nacmtype'] == 'full' for grad in nacs)
with Path('requests.jsonl').open('a') as out:
    out.write(json.dumps(document) + '\\n')
error = os.environ.get('BAGEL_TEST_ERROR', '')
if error == 'exit':
    sys.exit(7)
if error == 'exception':
    print('EXCEPTION RAISED: SMITH convergence not reached')
    sys.exit(0)
suffix = os.environ.get('BAGEL_TEST_SUFFIX', '.out')
energies = [f'{-75 + i * .01:.12f}' for i in range(nstates)]
if error == 'energy': energies[-1] = 'NaN'
Path('ENERGY' + suffix).write_text('\\n'.join(energies) + '\\n')
Path('orbitals.archive').write_text('mock orbitals')
for grad in forces['grads']:
    i = grad['target']
    if grad['title'] == 'force':
        stem = f'FORCE_{i}'
        scale = (i + 1) * .0001
        if error == 'missing' and i == 1: continue
    else:
        j = grad['target2']
        stem = f'NACME_{i}_{j}'
        if error == 'missing_nac' and (i, j) == (1, 2): continue
        scale = ((i + 1) * 10 + j + 1) * .001
    rows = ['                    0                   1                   2']
    for atom in range(natoms):
        value = scale * (1 if atom == 0 else -1)
        rows.append(f'{atom} {value:.10e} {value * 2:.10e} {value * 3:.10e}')
    if error == 'truncated': rows.pop()
    if error == 'reordered': rows[1:] = reversed(rows[1:])
    # No trailing blank line: the final atom must never be discarded.
    Path(stem + suffix).write_text('\\n'.join(rows) + '\\n')
'''


@contextmanager
def setup(method='sa_casscf', suffix='.out', error=''):
    with tempfile.TemporaryDirectory(prefix='abin bagel eh ') as folder:
        work = Path(folder)
        adapter = work / 'BAGEL-EH'
        shutil.copytree(ROOT / 'interfaces/BAGEL-EH', adapter)
        with (adapter / 'bagel.inp').open('a') as config:
            config.write(f'\nmethod={method}\nnact=2\nnclosed=4\n')
        fake = work / 'mock bagel.py'
        fake.write_text(FAKE_BAGEL)
        env = dict(os.environ, BAGELEXE=shlex.join([sys.executable, str(fake)]),
                   BAGEL_TEST_SUFFIX=suffix, BAGEL_TEST_ERROR=error)
        (work / 'geom.dat.001').write_text('H 0 0 0\nF 0 0 1.0\n')
        yield work, adapter, env


def invoke(work, adapter, env, mode, accuracy=7):
    flags = '1 1 1'  # 3 gradients or the 3 upper-triangular pairs.
    script = adapter / ('r.bagel-eh' if mode == 'gradients' else 'r.bagel-eh.nacm')
    return subprocess.run([str(script), '0', '001', str(accuracy)], input=f'3\n{flags}\n',
                          cwd=work, env=env, capture_output=True, text=True, timeout=30)


def require_success(result):
    assert result.returncode == 0, result.stdout + result.stderr


for method in ('sa_casscf', 'ms_caspt2', 'xms_caspt2'):
    for suffix in ('.out', '.txt'):
        with setup(method, suffix) as (work, adapter, env):
            require_success(invoke(work, adapter, env, 'gradients'))
            rows = [[float(x) for x in line.split()] for line in (work / 'engrad.dat.001').read_text().splitlines()]
            assert len(rows) == 9
            assert all(abs(rows[i][0] - (-75 + .01 * i)) < 1e-12 for i in range(3))
            for state in range(3):
                assert abs(rows[3 + 2 * state][0] - (state + 1) * .0001) < 1e-14
                assert rows[4 + 2 * state] == [-x for x in rows[3 + 2 * state]]
            require_success(invoke(work, adapter, env, 'nac', accuracy=5))
            nacs = [[float(x) for x in line.split()] for line in (work / 'nacm.dat').read_text().splitlines()]
            assert len(nacs) == 6
            assert [row[0] for row in nacs[::2]] == [.012, .013, .023]
            for i in range(0, 6, 2):
                assert nacs[i + 1] == [-x for x in nacs[i]]
            documents = [json.loads(line) for line in (adapter / 'requests.jsonl').read_text().splitlines()]
            assert documents[0]['bagel'][1]['title'] == 'rohf'
            assert len(documents) == 1  # NAC entry point must not launch BAGEL.
            grad_method = documents[0]['bagel'][2]['method'][0]
            assert abs(grad_method['thresh'] - 1e-8) < 1e-20
            assert ('smith' in grad_method) == (method != 'sa_casscf')
            if method != 'sa_casscf':
                assert grad_method['smith']['ms'] is True
                assert grad_method['smith']['xms'] == (method == 'xms_caspt2')
            # A second gradient call must reuse the saved reference.
            require_success(invoke(work, adapter, env, 'gradients'))
            latest = json.loads((adapter / 'input001.json').read_text())
            assert latest['bagel'][1]['title'] == 'load_ref'

for error in ('exit', 'exception', 'missing', 'truncated', 'reordered', 'energy', 'missing_nac'):
    with setup(error=error) as (work, adapter, env):
        # Complete old exports must never mask missing or malformed new data.
        (adapter / 'FORCE_1.out').write_text('0 1 2\n0 1 2 3\n1 4 5 6\n')
        (work / 'engrad.dat.001').write_text('stale output')
        result = invoke(work, adapter, env, 'gradients')
        assert result.returncode != 0, error
        assert not (work / 'engrad.dat.001').exists(), error
        assert 'ERROR: BAGEL-EH:' in result.stderr, result.stderr

with setup() as (work, adapter, env):
    result = invoke(work, adapter, env, 'nac')
    assert result.returncode != 0 and 'combined gradient' in result.stderr
    require_success(invoke(work, adapter, env, 'gradients'))
    # ABIN deletes public outputs before invoking the scripts.
    (work / 'engrad.dat.001').unlink()
    (work / 'nacm.dat').write_text('stale output')
    require_success(invoke(work, adapter, env, 'nac'))
    result = invoke(work, adapter, env, 'nac', accuracy=8)
    assert result.returncode != 0 and 'accuracy' in result.stderr
    assert not (work / 'nacm.dat').exists()
    assert len((adapter / 'requests.jsonl').read_text().splitlines()) == 1

# A changed request must not return cached NACs, even at the same timestep.
for change in ('geometry', 'settings', 'timestep', 'states', 'bead'):
    with setup() as (work, adapter, env):
        require_success(invoke(work, adapter, env, 'gradients'))
        if change == 'geometry':
            (work / 'geom.dat.001').write_text('H 0 0 0\nF 0 0 1.1\n')
        elif change == 'settings':
            with (adapter / 'bagel.inp').open('a') as config:
                config.write('shift=0.3\n')
        timestep = '1' if change == 'timestep' else '0'
        bead = '002' if change == 'bead' else '001'
        request = '2\n1\n' if change == 'states' else '3\n1 1 1\n'
        result = subprocess.run([str(adapter / 'r.bagel-eh.nacm'), timestep, bead, '5'],
                                input=request, cwd=work, env=env, capture_output=True, text=True)
        assert result.returncode != 0, change
        assert not (work / 'nacm.dat').exists(), change

# A failed new calculation invalidates the previous successful NAC cache.
with setup() as (work, adapter, env):
    require_success(invoke(work, adapter, env, 'gradients'))
    env['BAGEL_TEST_ERROR'] = 'missing_nac'
    assert invoke(work, adapter, env, 'gradients').returncode != 0
    assert not (work / 'engrad.dat.001').exists()
    assert invoke(work, adapter, env, 'nac').returncode != 0
    assert not (work / 'nacm.dat').exists()

# Single-state gradients have no NAC pairs.
with setup() as (work, adapter, env):
    result = subprocess.run([str(adapter / 'r.bagel-eh'), '0', '001', '7'],
                            input='1\n1\n', cwd=work, env=env, capture_output=True, text=True)
    require_success(result)
    assert len((work / 'engrad.dat.001').read_text().splitlines()) == 3

# Run the adapter with only Bash and its Unix tools on PATH. The mock BAGEL
# also uses shell builtins, so neither gradient nor NAC jobs can invoke Python.
with setup() as (work, adapter, env):
    tools = work / 'shell tools'
    tools.mkdir()
    for name in ('bash', 'awk', 'dirname', 'rm', 'cp', 'grep', 'cat', 'mktemp', 'mv', 'cmp'):
        (tools / name).symlink_to(shutil.which(name))
    fake = work / 'mock bagel.sh'
    fake.write_text('''#!/bin/bash
printf '%s\\n' '-7.5D+01' '-7.499D+01' '-7.498D+01' > ENERGY.out
printf '%s\\n' 'mock orbitals' > orbitals.archive
for stem in FORCE_0 FORCE_1 FORCE_2 NACME_0_1 NACME_0_2 NACME_1_2; do
    printf '%s\\n' '0 1 2' '0 1D-3 -2d-3 .3' '1 -1D-3 2d-3 -.3' > "$stem.out"
done
''')
    fake.chmod(0o755)
    env.update(PATH=str(tools), BAGELEXE=shlex.quote(str(fake)))
    assert shutil.which('python3', path=env['PATH']) is None
    for mode, output, count in (('gradients', 'engrad.dat.001', 9), ('nac', 'nacm.dat', 6)):
        require_success(invoke(work, adapter, env, mode))
        rows = (work / output).read_text().splitlines()
        assert len(rows) == count
        assert [float(x) for x in rows[-2].split()] == [.001, -.002, .3]

# Exercise the adapter through ABIN's independent Ehrenfest backend too.
for nac_accuracy in (5, 9):
    with setup() as (work, adapter, env):
        (work / 'mini.xyz').write_text('2\nTest\nH 0 0 0\nF 0 0 1.0\n')
        (work / 'velocities.in').write_text('2\nTest\nH 0 0 0\nF 0 0 0\n')
        (work / 'input.in').write_text('''&general
mdtype='ehrenfest', pot='bagel-eh', irest=0, irandom=123,
nstep=1, dt=0.01, nwrite=1, nwritex=1, nwritev=0,
/
&nhcopt
inose=0, temp=0.0,
/
&eh
nstate=3, istate_init=2, nac_accu1=NAC_ACCURACY,
/
&system
/
'''.replace('NAC_ACCURACY', str(nac_accuracy)))
        require_success(subprocess.run([str(BINARY), '-x', 'mini.xyz', '-v', 'velocities.in'],
                                       cwd=work, env=env, capture_output=True, text=True, timeout=30))
        requests = (adapter / 'requests.jsonl').read_text().splitlines()
        assert len(requests) == 2  # One combined calculation per endpoint.
        threshold = json.loads(requests[0])['bagel'][2]['method'][0]['thresh']
        assert abs(threshold - 10 ** -(max(7, nac_accuracy) + 1)) < 1e-20
        population = [float(x) for x in (work / 'pop.dat').read_text().splitlines()[-1].split()]
        assert abs(population[-1] - 1) < 1e-10

print('BAGEL-EH interface passed: combined derivatives, cache validation, accuracy, export formats, orbital reuse and ABIN integration.')
