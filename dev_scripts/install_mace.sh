#!/bin/bash
set -euo pipefail

REPO_URL="https://github.com/ACEsuit/mace.git"
REPO_DIR="$HOME/mace"
if [[ "$#" -eq 1 && ! -z $1 ]];then
   REPO_DIR=$1
fi

if [[ -e $REPO_DIR ]];then
  echo "ERROR: $REPO_DIR already exists."
  exit 1
fi

git clone "${REPO_URL}" "${REPO_DIR}" && cd "$REPO_DIR"

pip install --upgrade pip
pip install .
pip install mpi4py

echo "
Successfully installed MACE and mpi4py.

To use with ABIN, set MACE_PYTHON in utils/run.mace_mpi_abin.sh
to point to the Python interpreter that has these packages, e.g.:

export MACE_PYTHON=$(which python3)
"
