#!/usr/bin/env python3
"""BAGEL's all-state gradient/NAC interface for Ehrenfest dynamics.

Input: https://nubakery.org/grad/force.html
Exports: https://nubakery.org/grad/dynamics.html
"""
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys


def number(value):
    result = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(result):
        raise ValueError('BAGEL returned a non-finite value')
    return result


def export_path(stem):
    for suffix in ('.out', '.txt'):
        path = Path(stem + suffix)
        if path.is_file():
            return path
    raise ValueError(f'BAGEL did not export {stem}.out or {stem}.txt')


def read_vector(stem, natoms):
    path = export_path(stem)
    vectors = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or line.lstrip().startswith('#'):
            continue
        if not vectors and fields == ['0', '1', '2']:
            continue  # BAGEL's Cartesian column header.
        if len(fields) != 4 or int(fields[0]) != len(vectors):
            raise ValueError(f'Invalid or unordered atom row in {path}: {line}')
        vectors.append([number(value) for value in fields[1:]])
    if len(vectors) != natoms:
        raise ValueError(f'{path} contains {len(vectors)} atom rows; expected {natoms}')
    return vectors


def method_input(nstates, accuracy):
    method = os.environ['method']
    if method not in ('sa_casscf', 'ms_caspt2', 'xms_caspt2'):
        raise ValueError(f'Unsupported BAGEL-EH method: {method}')
    # BAGEL has no separate NAC accuracy keyword; follow the existing BAGEL
    # interface convention by controlling the CASSCF threshold with n+1.
    cas = {
        'title': 'casscf' if method == 'sa_casscf' else 'caspt2',
        'charge': int(os.environ['charge']),
        'nspin': int(os.environ['nspin']),
        'nstate': nstates,
        'nact': int(os.environ['nact']),
        'nclosed': int(os.environ['nclosed']),
        'thresh': 10.0 ** (-accuracy - 1),
        'maxiter': int(os.environ['maxiter_CASSCF']),
    }
    if method != 'sa_casscf':
        cas['smith'] = {
            'method': 'caspt2', 'ms': True, 'xms': method == 'xms_caspt2',
            'thresh': number(os.environ['thresh_CASPT2']),
            'maxiter': int(os.environ['maxiter_CASPT2']),
            'shift': number(os.environ['shift']), 'shift_imag': True,
        }
    return cas


def build_input(mode, nstates, accuracy, geometry):
    blocks = [{
        'title': 'molecule', 'basis': os.environ['basis'],
        'df_basis': os.environ['df_basis'], 'angstrom': True, 'geometry': geometry,
    }]
    if Path('orbitals.archive').is_file():
        blocks.append({'title': 'load_ref', 'file': 'orbitals', 'continue_geom': False})
    elif mode == 'nac':
        raise ValueError('NAC calculation requires orbitals.archive from the gradient calculation')
    else:
        blocks.append({
            'title': os.environ['hf_variant'], 'charge': int(os.environ['charge']),
            'nspin': int(os.environ['nspin']),
        })
    maxziter = int(os.environ['maxziter'])
    if mode == 'gradients':
        grads = [{'title': 'force', 'target': state, 'maxziter': maxziter} for state in range(nstates)]
    else:
        # The Ehrenfest force uses derivative couplings, not unweighted
        # interstate couplings; preserve the i<j orientation in the export.
        grads = [{'title': 'nacme', 'target': i, 'target2': j,
                  'nacmtype': 'full', 'maxziter': maxziter}
                 for i in range(nstates - 1) for j in range(i + 1, nstates)]
    blocks.append({'title': 'forces', 'export': True, 'export_single': False,
                   'grads': grads, 'method': [method_input(nstates, accuracy)]})
    if mode == 'gradients':
        blocks.append({'title': 'save_ref', 'file': 'orbitals'})
    return {'bagel': blocks}


def main():
    if len(sys.argv) != 5 or sys.argv[1] not in ('gradients', 'nac'):
        raise ValueError('Expected mode timestep bead accuracy')
    mode, timestep, bead, accuracy = sys.argv[1:]
    int(timestep)
    if not bead.isdigit():
        raise ValueError('Invalid bead index')
    accuracy = int(accuracy)
    if accuracy < 1:
        raise ValueError('Accuracy must be positive')
    output_path = Path(f'../engrad.dat.{bead}' if mode == 'gradients' else '../nacm.dat')
    output_path.unlink(missing_ok=True)
    nstates = int(sys.stdin.readline())
    flags = [int(value) for value in sys.stdin.readline().split()]
    expected = nstates if mode == 'gradients' else nstates * (nstates - 1) // 2
    if nstates < 1 or flags != [1] * expected or expected < 1:
        raise ValueError('Ehrenfest requires all states or all upper-triangular NAC pairs')
    geometry = []
    for line in Path(f'../geom.dat.{bead}').read_text().splitlines():
        atom, x, y, z = line.split()
        geometry.append({'atom': atom, 'xyz': [number(x), number(y), number(z)]})
    if not geometry:
        raise ValueError('Empty BAGEL geometry')
    document = build_input(mode, nstates, accuracy, geometry)
    prefix = f'input{bead}' + ('.nacm' if mode == 'nac' else '')
    input_file, log_file = Path(prefix + '.json'), Path(prefix + '.out')
    input_file.write_text(json.dumps(document, indent=2) + '\n')
    # Export files have no timestep stamp. Remove both supported suffixes
    # before every invocation so a failed calculation cannot reuse old data.
    for pattern in ('ENERGY.out', 'ENERGY.txt', 'FORCE_*.out', 'FORCE_*.txt', 'NACME_*.out', 'NACME_*.txt'):
        for path in Path('.').glob(pattern):
            path.unlink()
    command = shlex.split(os.environ['BAGELEXE'])
    if not command:
        raise ValueError('BAGELEXE is empty')
    with log_file.open('w') as log:
        log.write(input_file.read_text() + '\n')
        log.flush()
        result = subprocess.run(command + [str(input_file)], stdout=log, stderr=subprocess.STDOUT)
    with Path(prefix + '.out.all').open('a') as history:
        history.write(f'\nTIMESTEP = {timestep}\n' + log_file.read_text())
    if result.returncode != 0 or 'EXCEPTION RAISED' in log_file.read_text().upper():
        raise ValueError(f'BAGEL failed; inspect {log_file} (exit status {result.returncode})')
    rows = []
    if mode == 'gradients':
        energies = [number(value) for value in export_path('ENERGY').read_text().split()]
        if len(energies) != nstates:
            raise ValueError(f'Expected {nstates} BAGEL state energies, got {len(energies)}')
        rows.extend([[energy] for energy in energies])
        for state in range(nstates):
            rows.extend(read_vector(f'FORCE_{state}', len(geometry)))
    else:
        for i in range(nstates - 1):
            for j in range(i + 1, nstates):
                rows.extend(read_vector(f'NACME_{i}_{j}', len(geometry)))
    text = ''.join(' '.join(f'{value:.17e}' for value in row) + '\n' for row in rows)
    temporary = output_path.with_name(output_path.name + '.tmp')
    temporary.write_text(text)
    temporary.replace(output_path)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print(f'ERROR: BAGEL-EH: {error}', file=sys.stderr)
        sys.exit(2)
