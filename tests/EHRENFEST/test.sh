#!/bin/bash
set -euo pipefail

if [[ $1 == clean ]]; then
    exit 0
fi

# Keep the four-state shell-contract regression alongside the real-adapter matrix.
# The standard runner compares both sets of outputs with the .ref files here.
"$1" -x mini.xyz -v velocities.in > abin.out 2>&1
python3 interface_matrix.py "$1"
