#!/bin/bash
# Launch script for ABIN + MACE MPI interface.
#
# This script launches both the Python MACE server and ABIN,
# connecting them via MPI.
#
# Prerequisites:
#   - ABIN compiled with MPI=TRUE
#   - Python environment with: mpi4py, mace-torch, torch, ase, numpy
#   - MPICH (not OpenMPI)
#
# Usage:
#   ./run.mace_mpi_abin.sh

set -uo pipefail

# ABIN SETUP
ABIN_OUT=abin.out
ABIN_IN=input.in
GEOM_IN=mini.xyz
VELOC_IN=

MODEL_PATH=/path/to/mace.model
MACE_DEVICE=cpu   # 'cpu' or 'cuda'

# Path to Python with MACE dependencies (mpi4py, mace-torch, torch, ase, numpy).
# Uncomment and set to your Python interpreter:
# export MACE_PYTHON=/path/to/my/conda/envs/mace/bin/python
MACE_PYTHON="${MACE_PYTHON:-python3}"

# Path to MACE server script
MACE_SERVER=MACE/mace_server.py
# Path to ABIN binary
ABINEXE=./abin
if [[ -z ${MPI_PATH-} ]];then
  MPIRUN=mpirun
else
  MPIRUN=$MPI_PATH/bin/mpirun
fi


##### END OF INPUT #####

function files_exist() {
   local error=""
   for file in "$@"; do
      if [[ ! -f $file ]];then
         echo "ERROR: Cannot find file $file" >&2
         error=1
      fi
   done
   if [[ -n ${error-} ]];then
      exit 1
   fi
}

function validate_inputs() {
  files_exist $ABIN_IN $GEOM_IN $MACE_SERVER $ABINEXE

  # Check pot='_mace_' in ABIN input
  test=$(grep -E -o -e "^[^!]*pot[[:space:]]*=[[:space:]]*['\"]_mace_[\"']" $ABIN_IN || true)
  if [[ -z $test ]];then
    echo "ERROR: You did not specify pot='_mace_' in $ABIN_IN." >&2
    exit 1
  fi
}

# Cleanup function to stop the background processes
function cleanup {
  if [[ -n ${mace_pid-} ]] && kill -0 $mace_pid >& /dev/null; then
    echo "ERROR: MACE server $mace_pid is still running!"
    kill ${mace_pid} &> /dev/null
  fi

  if [[ -n ${abin_pid-} ]] && kill -0 $abin_pid >& /dev/null; then
    echo "ERROR: ABIN process $abin_pid is still running!"
    kill ${abin_pid} &> /dev/null
  fi
}

function wait_for_portfile {
  # Wait 10s for the MACE server to write the port file
  MAX_WAIT=20
  i=0
  while [[ ! -f mace_port.txt ]]; do
    if [[ $i -gt $MAX_WAIT ]]; then
      echo "ERROR: MACE server did not write port file"
      exit 1
    fi
    sleep 0.5
    (( i++ ))
  done
}

# Validate input files exist
validate_inputs

# Automatically call the cleanup function when the script
# exits or is interrupted by a signal
trap cleanup INT ABRT TERM EXIT

# LAUNCH MACE SERVER
$MPIRUN $MACE_PYTHON $MACE_SERVER --device $MACE_DEVICE --model-path $MODEL_PATH > mace_server.out 2>&1 &
mace_pid=$!
echo "Launched MACE server (PID: ${mace_pid})"

wait_for_portfile

# LAUNCH ABIN
ABIN_CMD="$ABINEXE -i $ABIN_IN -x $GEOM_IN"
if [[ -n $VELOC_IN ]];then
   ABIN_CMD="$ABIN_CMD -v $VELOC_IN"
fi
$MPIRUN $ABIN_CMD &> $ABIN_OUT &
abin_pid=$!
echo "Launched ABIN (PID: ${abin_pid})"
echo "(Monitor abin.out and mace_server.out for progress)"

# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $abin_pid >& /dev/null) && (kill -0 $mace_pid >& /dev/null) ); do sleep 1; done

echo "Simulation finished."
