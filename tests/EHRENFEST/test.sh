#!/bin/bash
set -euo pipefail

if [[ $1 == clean ]]; then
    exit 0
fi

export OMP_NUM_THREADS=1
exec > abin.out 2>&1

# Run the standalone checks before creating the outputs used for comparison.
python3 ../ehrenfest_checks/test_ehrenfest_backend.py "$1"
python3 ../ehrenfest_checks/test_ehrenfest_nve.py "$1"
"$1" -x mini.xyz -v velocities.in
python3 interface_matrix.py "$1"
