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
ABIN_IN=input.in
GEOM_IN=mini.xyz
VELOC_IN=
ABIN_OUTPUT=abin.out
ABINEXE="./abin"  # Path to ABIN binary
# Path to MPICH installation. Leave commented if you're using system-wide installation
# or if you're loading it via e.g. `load module mpich`
# This must match the installation that you used for ABIN compilation!
# MPI_PATH=/home/${USER}/software/mpich/4.3.2/install

# MACE SETUP
MODEL_PATH=/path/to/mace.model
MACE_DEVICE=cpu   # 'cpu' or 'cuda'

# If you're on a system with multiple GPUs, select which GPU to use
# with CUDA_VISIBLE_DEVICES
#export CUDA_VISIBLE_DEVICES=0

# Path to Python executable with MACE dependencies (mpi4py, mace-torch, torch, ase, numpy).
# Typically, you'll provide an absolute path to your conda or virtual environment
MACE_PYTHON="../.venv/bin/python"

# Path to MACE server script (copy from interfaces/MACE)
MACE_SERVER=./mace_server.py
MACE_OUTPUT=mace_server.out


##### END OF USER INPUT #####
#

# Convert to absolute path if not provided as such
MACE_PYTHON=$(command -v ${MACE_PYTHON:-python3})
if [[ -z ${MACE_PYTHON} ]]; then
  echo "ERROR: Could not find python executable!"
  exit 1
fi
if [[ -z ${MPI_PATH-} ]]; then
  MPIRUN=$(command -v mpirun)
  if [[ -z $MPIRUN ]]; then
    echo "Command mpirun not found. Please specify MPI_PATH variable"
    exit 1
  fi
else
  MPIRUN="$MPI_PATH/bin/mpirun"
  export LD_LIBRARY_PATH="${MPI_PATH}/lib:${LD_LIBRARY_PATH-}"
fi

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
  files_exist "$ABIN_IN" "$GEOM_IN" "$MACE_SERVER" "$ABINEXE" "$MACE_PYTHON" "$MPIRUN"

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
    echo "WARNING: MACE server $mace_pid is still running!"
    kill ${mace_pid}
  fi

  if [[ -n ${abin_pid-} ]] && kill -0 $abin_pid >& /dev/null; then
    echo "WARNING: ABIN process $abin_pid is still running!"
    kill ${abin_pid}
  fi
}

function wait_for_portfile {
  # Wait 20s for the MACE server to write the port file
  MAX_WAIT=20
  i=0
  while [[ ! -f mace_port.txt ]]; do
    if [[ $i -gt $MAX_WAIT ]] || ! kill -0 $mace_pid >& /dev/null; then
      echo "ERROR: MACE server did not write port file 'mace_port.txt'"
      echo "cat $MACE_OUTPUT"
      cat $MACE_OUTPUT
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
$MPIRUN $MACE_PYTHON $MACE_SERVER --device "$MACE_DEVICE" --model-path "$MODEL_PATH" > "$MACE_OUTPUT" 2>&1 &
mace_pid=$!
echo "Launched MACE server (PID: ${mace_pid})"

wait_for_portfile

# LAUNCH ABIN
ABIN_CMD="$ABINEXE -i $ABIN_IN -x $GEOM_IN"
if [[ -n $VELOC_IN ]];then
   ABIN_CMD="$ABIN_CMD -v $VELOC_IN"
fi
$MPIRUN $ABIN_CMD &> $ABIN_OUTPUT &
abin_pid=$!
echo "Launched ABIN (PID: ${abin_pid})"
echo "(Monitor $ABIN_OUTPUT and $MACE_OUTPUT for progress)"

# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $abin_pid >& /dev/null) && (kill -0 $mace_pid >& /dev/null) ); do sleep 1; done
sleep 2 # grace time for program termination

echo "Simulation finished."
