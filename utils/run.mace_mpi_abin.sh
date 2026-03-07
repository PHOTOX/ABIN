#!/bin/bash
# Launch script for ABIN + MACE MPI interface.
#
# This script launches both the Python MACE server and ABIN,
# connecting them via MPI.
#
# Prerequisites:
#   - ABIN compiled with MPI=TRUE
#   - Python environment with: mpi4py, mace-torch, torch, ase, numpy
#     Install with: pip install mpi4py mace-torch torch ase numpy
#   - MPICH (not OpenMPI)
#
# Usage:
#   bash run.mace_mpi_abin.sh

set -euo pipefail

# ABIN SETUP
ABIN_OUT=abin.out
ABIN_IN=input.in
GEOM_IN=mini.xyz
VELOC_IN=

# Path to Python with MACE dependencies (mpi4py, mace-torch, torch, ase, numpy).
# Uncomment and set to your Python interpreter:
# export MACE_PYTHON=/path/to/my/conda/envs/mace/bin/python
MACE_PYTHON="${MACE_PYTHON:-python3}"

# Path to MACE server script
MACE_SERVER=MACE/mace_server.py

################

LAUNCH_DIR=$PWD

function files_exist() {
   local error=""
   for file in $* ;
   do
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
  files_exist $ABIN_IN $GEOM_IN $MACE_SERVER

  # Check pot='_mace_' in ABIN input
  test=$(egrep -o -e "^[^!]*pot[[:space:]]*=[[:space:]]*['\"]_mace_[\"']" $ABIN_IN || true)
  if [[ -z $test ]];then
    echo "ERROR: You did not specify pot='_mace_' in $ABIN_IN." >&2
    exit 1
  fi
}

# Validate input files exist
validate_inputs

# Determine MPIRUN
if [[ -z ${MPI_PATH-} ]];then
  MPIRUN=mpirun
else
  MPIRUN=$MPI_PATH/bin/mpirun
fi

MPIRUN_ABIN="$MPIRUN -n 1"
MPIRUN_MACE="$MPIRUN -n 1"

# Determine ABIN executable location
if [[ -z ${ABINEXE-} ]];then
  ABINEXE=bin/abin
fi

echo "Starting MACE MPI simulation"
echo "=============================="

declare -A job_pids

# LAUNCH MACE SERVER
$MPIRUN_MACE $MACE_PYTHON $MACE_SERVER > mace_server.out 2>&1 &
job_pids[mace]=$!
echo "Launched MACE server (PID: ${job_pids[mace]})"

# Wait a moment for the port file to be created
sleep 1

# LAUNCH ABIN
ABIN_CMD="$ABINEXE -i $ABIN_IN -x $GEOM_IN"
if [[ -n $VELOC_IN ]];then
   ABIN_CMD=$ABIN_CMD" -v $VELOC_IN"
fi
$MPIRUN_ABIN $ABIN_CMD > $ABIN_OUT 2>&1 &
job_pids[abin]=$!
echo "Launched ABIN (PID: ${job_pids[abin]})"

# PERIODICALLY CHECK WHETHER ABIN AND MACE ARE RUNNING
function join_by { local IFS="$1"; shift; echo "$*"; }
regex=`join_by \| ${job_pids[@]}`
while true;do
   njobs=$(ps -eo pid | egrep "$regex" | wc -l)
   if [[ $njobs -eq 0 ]];then
      echo "Both ABIN and MACE server stopped"
      break
   elif [[ $njobs -lt ${#job_pids[@]} ]];then
      echo "One of ABIN or MACE server died. Killing the rest." >&2
      kill -9 ${job_pids[@]} 2>/dev/null || true
      break
   fi
   sleep 3
done

echo "Simulation finished."
