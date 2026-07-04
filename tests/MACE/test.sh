#!/bin/bash
set -uo pipefail

ABINEXE=$1

ABINOUT=abin.out
ABININ=input.in
ABINGEOM=mini.xyz

MACE_SERVER=../../interfaces/MACE/mace_server.py
MACE_OUT=mace_server.out

ABIN_CMD="$ABINEXE -i $ABININ -x $ABINGEOM"
MACE_CMD="python3 $MACE_SERVER --device cpu --model-path __MOCK_HARMONIC__"

# If $1 = "clean"; exit early.
if [[ "${1-}" = "clean" ]]; then
  rm -f $MACE_OUT $ABINOUT ./*.dat ./*.diff
  rm -f restart.xyz velocities.xyz forces.xyz movie.xyz restart.xyz.old
  rm -f mace_port.txt ERROR ompi_uri.txt
  exit 0
fi

# Determine MPI paths
if [[ -z ${MPI_PATH-} ]]; then
  export MPIRUN=mpirun
else
  export MPIRUN=$MPI_PATH/bin/mpirun
fi

# Detect OpenMPI vs MPICH
IS_OPENMPI=false
if $MPIRUN --version 2>&1 | grep -q "Open MPI"; then
  IS_OPENMPI=true
fi

ompi_server_pid=""
MPIRUN_EXTRA_ARGS=""

if [[ "$IS_OPENMPI" = "true" ]]; then
  # OpenMPI requires ompi-server for MPI_Comm_connect/accept
  OMPI_SERVER=${MPI_PATH-}/bin/ompi-server
  if [[ -z ${MPI_PATH-} ]]; then
    OMPI_SERVER=ompi-server
  fi

  if ! which $OMPI_SERVER &> /dev/null; then
    echo "Skipping MACE test: ompi-server not found (required for OpenMPI)"
    exit 0
  fi

  OMPI_URI_FILE="$PWD/ompi_uri.txt"
  $OMPI_SERVER --no-daemonize -r "$OMPI_URI_FILE" &
  ompi_server_pid=$!
  sleep 1

  if [[ ! -f "$OMPI_URI_FILE" ]]; then
    echo "ERROR: ompi-server did not create URI file" >&2
    kill $ompi_server_pid 2>/dev/null
    exit 1
  fi

  MPIRUN_EXTRA_ARGS="--ompi-server file:$OMPI_URI_FILE"
fi

MPIRUN_CMD="$MPIRUN -n 1 $MPIRUN_EXTRA_ARGS"

# Cleanup function to stop the background processes
function cleanup {
  if [[ -n ${macepid-} ]] && kill -0 $macepid >& /dev/null; then
    echo "ERROR: MACE server $macepid is still running!" >> ERROR
    kill ${macepid-} &> /dev/null || true
  fi

  if [[ -n ${abinpid-} ]] && kill -0 $abinpid >& /dev/null; then
    echo "ERROR: ABIN process $abinpid is still running!" >> ERROR
    kill ${abinpid-} &> /dev/null || true
  fi

  if [[ -n "${ompi_server_pid-}" ]]; then
    kill $ompi_server_pid &> /dev/null || true
  fi
}

# Automatically call the cleanup function when the script
# exits or is interrupted by a signal
trap cleanup INT ABRT TERM EXIT

# Launch mock MACE server
$MPIRUN_CMD $MACE_CMD > $MACE_OUT 2>&1 &
macepid=$!

# Wait for the server to write the port file
MAX_WAIT=15
i=0
while [[ ! -f mace_port.txt && $i -lt $MAX_WAIT ]]; do
  sleep 0.5
  let ++i
done

if [[ ! -f mace_port.txt ]]; then
  echo "ERROR: MACE server did not write port file" >> ERROR
  cat $MACE_OUT 2>/dev/null || true
  exit 1
fi

# Launch ABIN
$MPIRUN_CMD $ABIN_CMD > $ABINOUT 2>&1 &
abinpid=$!

# Give both processes 10 seconds to finish
MAX_ITER=20
iter=0
# Note about 'kill -0' https://unix.stackexchange.com/questions/169898/what-does-kill-0-do
while ( (kill -0 $abinpid >& /dev/null) || (kill -0 $macepid >& /dev/null) ); do
  if [[ $iter -gt $MAX_ITER ]]; then
    echo "Test did not finish in time" >> ERROR
    break
  fi
  sleep 0.5
  let iter++
done

# Any errors in this script should be echoed to file ERROR
# so that the overall test fails when comparing to empty ERROR.ref file
# Here we create an empty one in case no errors actually occured, as is expected.
touch ERROR
